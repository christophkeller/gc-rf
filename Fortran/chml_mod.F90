!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chml_mod.F90
!
! !DESCRIPTION: Module Chml\_Mod contines arrays and routines for the
!  chemical integrator based on a random forest regressor. 
!\\
!\\
! !INTERFACE: 
!
MODULE Chml_Mod
!
! !USES:
!
  USE Precision_Mod            ! For GEOS-Chem Precision (fp)
  USE GIGC_ErrCode_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_Chml
  PUBLIC  :: Run_Chml
  PUBLIC  :: Cleanup_Chml
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Eval_Chml
  PRIVATE :: ReadTxtArrI
  PRIVATE :: ReadTxtArrF
  PRIVATE :: CreateTrees
!    
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Generic type for tree
  TYPE :: TREE
     INTEGER                  :: n_nodes                     ! Number of nodes
     INTEGER                  :: n_leafs                     ! Number of leafs
     INTEGER,  POINTER        :: node_feature(:)   => NULL() ! Node features
     REAL(fp), POINTER        :: node_threshold(:) => NULL() ! Node threshold
     INTEGER,  POINTER        :: node_left   (:)   => NULL() ! Node left ID 
     INTEGER,  POINTER        :: node_right  (:)   => NULL() ! Node right ID
     REAL(fp), POINTER        :: leaf_value(:)     => NULL() ! Leaf value
  END TYPE

  ! Generic type for random forest of one species
  TYPE :: RFtype
     LOGICAL                  :: hasRF               ! Is NN available?
     CHARACTER(LEN=1023)      :: srcFile             ! RF source file
     CHARACTER(LEN=255)       :: KPP_SpcName         ! KPP species name
     INTEGER                  :: Kpp_SpcID           ! KPP species ID 
     CHARACTER(LEN=255)       :: GC_SpcName          ! GEOS-Chem species name
     INTEGER                  :: GC_SpcID            ! GEOS-Chem species ID 
     INTEGER                  :: IdxInCin            ! Species index in input vector 
     INTEGER                  :: n_trees             ! # of trees 
     TYPE(TREE), POINTER      :: forest(:) => NULL() ! list of trees
     INTEGER                  :: out_type            ! output type
     REAL(fp)                 :: corr                ! correction factor 
  END TYPE

  ! Generic type for random forest
  TYPE :: RF
     INTEGER                     :: n_features
     INTEGER,  POINTER           :: TRCid(:)
     INTEGER,  POINTER           :: JVALid(:)
     INTEGER                     :: NUMDEN_ID
     INTEGER                     :: TEMP_ID 
     INTEGER                     :: PRESS_ID 
     INTEGER                     :: SUNCOS_ID 
     INTEGER                     :: RH_ID
     INTEGER                     :: QLIQ_ID
     INTEGER                     :: QICE_ID
     CHARACTER(LEN=255), POINTER :: feature_names(:)
     INTEGER                     :: NSPEC = 0
     TYPE(RFtype), POINTER       :: SpcRF(:) => NULL()
  END TYPE

  ! Random forest structure
  TYPE(RF), POINTER              :: MyRF => NULL()

  ! Location of source files. Hardcoded for now 
  CHARACTER(LEN=511)             :: srcFileRoot = '/discover/nobackup/cakelle2/ML/chml/sandbox_v4/chml_models/'

  ! Output types
  INTEGER, PARAMETER             :: OutTypeTend = 1  ! Tendency
  INTEGER, PARAMETER             :: OutTypeConc = 2  ! Concentration
  INTEGER, PARAMETER             :: OutTypeNorm = 3  ! Normalized (Conc2/Conc1)
  INTEGER, PARAMETER             :: OutTypeLogc = 4  ! Log concentration 
  INTEGER, PARAMETER             :: OutTypeDrop = 5  ! Log concentration 


  ! Maximum number of trees in forest
  INTEGER,   PARAMETER           :: max_trees = 30

  ! Time step used when training forest
  REAL(fp), PARAMETER            :: DT_RF = 15.0*60.0_fp

  ! Adjust for NOx?
  LOGICAL,  PARAMETER            :: DoNOx = .TRUE.
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: run_chml
!
! !DESCRIPTION: Subroutine Run\_Chml is the driver subroutine for
!  the random forest forward integrator. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_Chml ( am_I_Root, nin, CIN, JVAL, NUMDEN, TEMP, PRESS, &
                        SUNCOS, RH, QLIQ, QICE, DT, COUT, TOTT, TREET, RC, verb )
!
! !USES:
!
    USE CMN_FJX_MOD,      ONLY : JVN_
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root  ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: nin 
    REAL(fp),       INTENT(IN)    :: CIN(nin)   
    REAL(fp),       INTENT(IN)    :: JVAL(JVN_)
    REAL(fp),       INTENT(IN)    :: NUMDEN
    REAL(fp),       INTENT(IN)    :: TEMP 
    REAL(fp),       INTENT(IN)    :: PRESS
    REAL(fp),       INTENT(IN)    :: SUNCOS
    REAL(fp),       INTENT(IN)    :: RH  
    REAL(fp),       INTENT(IN)    :: QLIQ
    REAL(fp),       INTENT(IN)    :: QICE
    REAL(fp),       INTENT(IN)    :: DT
    LOGICAL,        INTENT(IN)    :: verb
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(INOUT) :: COUT(nin)
    REAL(fp),       INTENT(OUT)   :: TOTT
    REAL(fp),       INTENT(OUT)   :: TREET 
    INTEGER,        INTENT(INOUT) :: RC         ! Success or failure
! 
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(RFtype), POINTER :: ThisRF
    INTEGER               :: I, N, nft, nout_
    INTEGER               :: IDX
    INTEGER               :: idNO, idNO2 
    LOGICAL               :: NOxTend, NOxSet
    REAL(fp), ALLOCATABLE :: cin_(:), cout_(:)
    REAL(fp), ALLOCATABLE :: cin__(:)
    REAL(fp)              :: t0, t1, treet_, treetsum, ncalls
    REAL(fp)              :: NOx, NOrat, NOxin
    CHARACTER(LEN=255)    :: ErrMsg,   ThisLoc

    !=======================================================================
    ! Run_Chml begins here!
    !=======================================================================

    ! Initialize
    RC        =  GIGC_SUCCESS
    ErrMsg    =  ''
    ThisLoc   =  ' -> at Run_Chml (in module GeosCore/chml_mod.F)'
    ThisRF    => NULL()

    ! Set inputs/output
    nft   = MyRF%n_features 
    nout_ = 1
    ALLOCATE(cin_(nft),cin__(nft),cout_(nout_))
    cin_  =    0.0_fp
    cin__ =    0.0_fp
    cout_ =    0.0_fp
    NOx   = -999.9_fp
    NOrat = -999.9_fp

    ! Map inputs 
    DO I=1,nft
       IF ( MyRF%TRCid(I) > 0 ) THEN
          cin_(I) = CIN(MyRF%TRCid(I)) / NUMDEN * 1.0e15
       ELSEIF ( MyRF%JVALid(I) > 0 ) THEN
          cin_(I) = JVAL(MyRF%JVALid(I))
       ELSEIF ( MyRF%TEMP_ID == I ) THEN
          cin_(I) = TEMP
       ELSEIF ( MyRF%NUMDEN_ID == I ) THEN
          cin_(I) = NUMDEN 
       ELSEIF ( MyRF%SUNCOS_ID == I ) THEN
          cin_(I) = SUNCOS
       ELSEIF ( MyRF%RH_ID == I ) THEN
          cin_(I) = RH
       ELSEIF ( MyRF%PRESS_ID == I ) THEN
          cin_(I) = PRESS
       ELSEIF ( MyRF%QLIQ_ID == I ) THEN
          cin_(I) = QLIQ
       ELSEIF ( MyRF%QICE_ID == I ) THEN
          cin_(I) = QICE
       ELSE
          IF (am_I_Root) THEN
             write(*,*) 'Cannot match feature ',I,nft
          ENDIF
          RC = GIGC_FAILURE
          EXIT
          cin_(I) = 1e-20
       ENDIF
    ENDDO

    ! Error trap
    IF ( RC==GIGC_FAILURE ) RETURN

    ! Run time
    CALL CPU_TIME(t0)
    treetsum = 0.0_fp 
    ncalls   = 0.0_fp

    ! Memorize NO and NO2 index
    idNO  = -1
    idNO2 = -1
    NOxTend = .FALSE.
    NOxSet  = .FALSE.

    ! Run instances
    DO I=1,MyRF%NSPEC
       ThisRF => MyRF%SpcRF(I)

       ! Is NO or NO2?
       IF ( TRIM(ThisRF%GC_SpcName) == 'NO'  ) idNO  = I
       IF ( TRIM(ThisRF%GC_SpcName) == 'NO2' ) idNO2 = I

       ! Compute forests 
       IF ( ThisRF%hasRF ) THEN

          ! Index in CIN/COUT
          IDX = ThisRF%IdxInCin

          ! Input features
          cin__(:) = cin_(:)
          IF ( ThisRF%out_type == OutTypeDrop ) THEN
             DO N=1,nft
                IF ( MyRF%TRCid(N) == I ) THEN
                   cin__(N) = 0.0e0
                   IF ( verb ) write(*,*) 'set input conc to zero: ',N,TRIM(MyRF%feature_names(N)),TRIM(ThisRF%GC_SpcName)
                   EXIT
                ENDIF
             ENDDO
          ENDIF

          ! Run random forest
          cout_ = 0.0_fp
          CALL Eval_Chml ( am_I_Root, ThisRF, nft, nout_, cin_, cout_, treet_, RC )
          IF ( RC /= GIGC_SUCCESS ) EXIT

          ! Accumulate statistics for timer
          treetsum = treetsum + treet_
          ncalls   = ncalls + 1.0

          IF ( verb ) THEN
           write(*,*) 'CHML raw out: ',I,trim(ThisRF%GC_SpcName),cout_(1)
          ENDIF

          ! Normalize output
          IF ( TRIM(ThisRF%GC_SpcName) == 'NOratio' ) THEN 
             NOrat = cout_(1) / 1.0e3
             NOrat = MIN(MAX(NOrat,0.0),1.0)

          ELSEIF ( TRIM(ThisRF%GC_SpcName) == 'NOx' ) THEN 
             NOx   = cout_(1) / 1.0e15 * NUMDEN
             IF ( ThisRF%out_type == OutTypeTend ) THEN
                NOxTend = .TRUE.
             ELSE
                NOx   = MAX(0.0_fp,NOx)
             ENDIF
             NOxSet = .TRUE.
          ELSE
             ! Normalize back
             cout_(1) = cout_(1) / 1.0e15 * NUMDEN
             ! Apply tendency to species array, use provided time step DT
             IF ( ThisRF%out_type == OutTypeTend ) THEN
                COUT(IDX) = CIN(IDX) + ( cout_(1) * DT )
                IF ( verb ) THEN
                 write(*,*) 'CHML applied tend: ',I,IDX,trim(ThisRF%GC_SpcName),COUT(IDX),CIN(IDX),cout_(1)
                ENDIF
             ELSEIF ( ThisRF%out_type == OutTypeConc .OR. &
                      ThisRF%out_type == OutTypeDrop ) THEN
                !COUT(I) = 10**cout_(1) 
                COUT(IDX) = cout_(1) 
                IF ( verb ) THEN
                 write(*,*) 'CHML applied conc: ',I,IDX,trim(ThisRF%GC_SpcName),COUT(IDX),CIN(IDX),cout_(1)
                ENDIF
             ELSEIF ( ThisRF%out_type == OutTypeNorm ) THEN
                COUT(IDX) = CIN(IDX) * 10**cout_(1)
             ELSEIF ( ThisRF%out_type == OutTypeLogc ) THEN
                COUT(IDX) = 10**cout_(1)
             ENDIF

             ! Make sure value is not negative
             COUT(IDX) = MAX(0.0_fp,COUT(IDX))
          ENDIF

       ! If not found keep output concentration as is. 
       !ELSE
       !   IF ( TRIM(ThisRF%GC_SpcName) /= 'NOratio' .AND. &
       !        TRIM(ThisRF%GC_SpcName) /= 'NOx'            ) THEN 
       !      COUT(I) = COUT(I)
       !   ENDIF
       ENDIF
    ENDDO

    ! Adjust for NOx if available 
    ! NOx is in kg N/m3 (after normalization). NOratio is the ratio of NO / NOx 
    IF ( DoNOx .and. NOxSet ) THEN
       IF ( verb ) THEN
          write(*,*) 'Before applying NOx ratios: '
          write(*,*) 'NO : ',idNO,CIN(idNO),COUT(idNO)
          write(*,*) 'NO2: ',idNO2,CIN(idNO2),COUT(idNO2)
       ENDIF
       IF ( NOxTend ) THEN
          NOxin       = CIN(idNO)*14.0_fp/30.0_fp + CIN(idNO2)*14.0_fp/46.0_fp
          NOx         = NOxin + ( NOx * DT )
          NOx         = MAX(0.0_fp,NOx)
          COUT(idNO ) = NOrat * NOx * 30.0_fp / 14.0_fp
          COUT(idNO2) = ( 1.0 - NOrat ) * NOx * 46.0_fp / 14.0_fp
       ELSE
          COUT(idNO ) = NOrat * NOx * 30.0_fp / 14.0_fp
          COUT(idNO2) = ( 1.0 - NOrat ) * NOx * 46.0_fp / 14.0_fp
          IF ( verb ) THEN
             write(*,*) 'Applied NOx ratios: ',NOrat,NOx
             write(*,*) 'NO : ',idNO,CIN(idNO),COUT(idNO)
             write(*,*) 'NO2: ',idNO2,CIN(idNO2),COUT(idNO2)
          ENDIF
       ENDIF
    ENDIF

    ! Run time
    CALL CPU_TIME(t1)
    TOTT  = ( t1 - t0 ) * 1000.0
    TREET = treetsum / ncalls

    ! Verbose
    IF ( verb ) THEN
     write(*,*) 'CHML debug:'
     write(*,*) 'raw input conc:'
     do I=1,nin
      write(*,*) i,cin(i)
     enddo
     write(*,*) 'converted input features:'
     do I=1,nft
      write(*,*) i,cin_(i)
     enddo
     write(*,*) 'output conc: ',size(cout,1)
     do I=1,nin
      write(*,*) i,cout(i)
     enddo

    ENDIF

    ! Cleanup 
    DEALLOCATE(cin_,cin__,cout_)

  END SUBROUTINE Run_Chml 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: eval_chml
!
! !DESCRIPTION: Subroutine Eval\_Chml runs an instance of the random forest.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Eval_Chml ( am_I_Root, ThisRF, nft, nout, cin, cout, treet, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root  ! Is this the root CPU?
    TYPE(RFtype),   INTENT(IN)    :: ThisRF     ! Random forest instance for species
    INTEGER,        INTENT(IN)    :: nft        ! # of input features
    INTEGER,        INTENT(IN)    :: nout       ! # of outputs
    REAL(fp),       INTENT(IN)    :: cin(nft)   ! input(s)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT)   :: cout(nout) ! output(s)
    REAL(fp),       INTENT(OUT)   :: treet      ! average tree execution time 
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
! 
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: itree, ntr, nnd, nlf, depth
    INTEGER               :: inode, ileft, iright, iftr
    INTEGER               :: leafindex
    REAL(fp)              :: ithreshold
    REAL(fp)              :: t0, t1 
    REAL(fp), ALLOCATABLE :: guesses(:)
    REAL(fp)              :: mean_guess
    TYPE(TREE), POINTER   :: ThisTree => NULL()
    CHARACTER(LEN=255)    :: ErrMsg,   ThisLoc

    !=======================================================================
    ! Eval_Chml begins here!
    !=======================================================================

    ! Initialize
    RC        =  GIGC_SUCCESS
    ErrMsg    =  ''
    ThisLoc   =  ' -> at Run_Instance (in module GeosCore/chml_mod.F)'

    ! Number of trees & max. number of nodes
    ntr = ThisRF%n_trees 

    ! Time
    CALL CPU_TIME(t0)

    ! Get estimate for every tree
    ALLOCATE(guesses(ntr))
    guesses = 0.0_fp
    DO itree = 1, ntr
       ! Get three
       ThisTree => ThisRF%forest(itree)

       ! Get number of nodes and leafs on this tree
       nnd = ThisTree%n_nodes
       nlf = ThisTree%n_leafs

       ! start at first node
       inode = 1
       depth = 1
       ! repeat the following until we reach leaf. Leaf nodes
       ! have identical left and right children node IDs.
       DO
          ! If this is a leaf the node index is larger than the
          ! # of nodes. Get leaf value and go to next tree 
          IF ( inode > nnd ) THEN
             leafindex = inode - nnd
             EXIT
          ENDIF

          ! if we get here we are on a node:

          ! get children nodes of this node 
          !ileft  = ThisTree%node_left (inode)
          !iright = ThisTree%node_right(inode)

          ! get index of feature, plus feature threshold   
          iftr       = ThisTree%node_feature(inode) 
          !ithreshold = ThisTree%node_threshold(inode)

          ! Next node: left if feature is true, right otherwise
          !IF ( cin(iftr) <= ithreshold ) THEN
          IF ( cin(iftr) <= ThisTree%node_threshold(inode) ) THEN
             inode = ThisTree%node_left (inode)
          ELSE
             inode = ThisTree%node_right(inode)
          ENDIF

          ! count tree depth: this is primarily for debugging
          !depth = depth + 1
          !IF ( depth > 200 ) THEN
          !   write(*,*) 'Your tree is too deep!'
          !   EXIT
          !ENDIF
       ENDDO

       ! Assign leaf value to guess
       guesses(itree) = ThisTree%leaf_value(leafindex)
    ENDDO

    ! Average across all trees
    cout(1) = sum(guesses) / ntr 

    ! Apply correction factor
    cout(1) = cout(1) * ThisRF%corr

    ! Time
    CALL CPU_TIME(t1)
    treet = ( t1 - t0 ) * 1000.0 / ntr 

    ! Cleanup & return
    DEALLOCATE(guesses)
    RC = GIGC_SUCCESS

  END SUBROUTINE Eval_Chml
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_chml
!
! !DESCRIPTION: Subroutine Init\_Chml initializes the random forest. 
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE Init_Chml ( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE inquireMod,           ONLY : findFreeLUN
    USE GIGC_Input_Opt_Mod,   ONLY : OptInput
    USE GIGC_State_Chm_Mod,   ONLY : ChmState
!
! !INPUT PARAMETERS:
!    
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!    
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!    
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, N, GCID
    INTEGER                :: nft, ntr, ntruse, nnd, tmp 
    INTEGER                :: lun, ios, idx
    REAL(fp)               :: corr
    LOGICAL                :: exists
    CHARACTER(LEN=1023)    :: srcFile
    CHARACTER(LEN=8191)    :: line
    CHARACTER(LEN= 31)     :: tpe 
    CHARACTER(LEN=255)     :: ivar, spc, ErrMsg, ThisLoc

    !=======================================================================
    ! Init_Chml begins here!
    !=======================================================================

    ! Assume success
    RC       = GIGC_SUCCESS

    !=======================================================================
    ! Initialize variables
    !=======================================================================
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_Chml (in module GeosCore/chml_mod.F90)'

    ! This should not happen 
    IF ( ASSOCIATED(MyRF) ) THEN
       IF ( am_I_Root ) THEN
          WRITE(*,*) 'Warning: initialize Chml multiple times ' // &
                     '- cleanup previously defined random forest!!'
       ENDIF
       CALL Cleanup_Chml( am_I_Root, RC )
    ENDIF

    !=======================================================================
    ! Read header and define generic features of random forest 
    !=======================================================================
    ALLOCATE(MyRF)
    MyRF%NSPEC = Input_Opt%N_TRACERS + 2 ! last 2 are NOx and NOratio 

    ! -----------------------
    ! Read header information 
    srcFile = TRIM(srcFileRoot)//'forest_header.txt' 
    INQUIRE(file=trim(srcFile),exist=exists)
    IF ( .NOT. exists ) THEN
       IF(am_I_Root) WRITE(*,*) 'Error: ',TRIM(srcFile),' not found'
       RC = GIGC_FAILURE
       RETURN 
    ENDIF
    ! Open file
    lun = FindFreeLun()
    OPEN(unit=lun,file=TRIM(srcFile))
    ! Number of inputs 
    READ(lun,'(a)',IOSTAT=IOS) line 
    IF ( IOS < 0 ) THEN
       IF(am_I_Root) WRITE(*,*) 'Error reading ',TRIM(srcFile)
       RC = GIGC_FAILURE
       RETURN 
    ENDIF
    ! Get number of features
    idx = INDEX(trim(line),'Features:')
    IF ( idx <= 0 ) THEN
       IF ( am_I_Root ) &
          WRITE(*,*) '`Features:` not found - error reading ',TRIM(srcFile),'::',TRIM(line)
       RC = GIGC_FAILURE
       RETURN 
    ENDIF
    READ(line(idx+9:len(line)),*) nft 
    MyRF%n_features = nft

    ! Allocate index vectors
    ALLOCATE(MyRF%TRCId(nft))
    ALLOCATE(MyRF%JVALid(nft))
    ALLOCATE(MyRF%feature_names(nft))
    MyRF%TRCid     = -1
    MyRF%JVALid    = -1
    MyRF%NUMDEN_ID = -1
    MyRF%TEMP_ID   = -1
    MyRF%PRESS_ID  = -1
    MyRF%SUNCOS_ID = -1
    MyRF%RH_ID     = -1
    MyRF%QLIQ_ID   = -1
    MyRF%QICE_ID   = -1
    MyRF%feature_names = ''

    ! Loop over all input features
    DO I=1,nft
       READ(lun,'(a)',IOSTAT=IOS) ivar
       IF ( IOS < 0 ) THEN
          IF(am_I_Root) WRITE(*,*) 'Error reading ',TRIM(srcFile), ' - line ',I
          RC = GIGC_FAILURE
          EXIT 
       ENDIF
       ! Get index
       IF ( TRIM(ivar)=='KPP_AIRDEN' ) THEN
          MyRF%NUMDEN_ID = I
       ELSEIF ( TRIM(ivar)=='KPP_TEMP' ) THEN
          MyRF%TEMP_ID = I
       ELSEIF ( TRIM(ivar)=='KPP_PRESS' ) THEN
          MyRF%PRESS_ID = I
       ELSEIF ( TRIM(ivar)=='KPP_SUNCOS' ) THEN
          MyRF%SUNCOS_ID = I
       ELSEIF ( TRIM(ivar)=='KPP_RH' ) THEN
          MyRF%RH_ID = I
       ELSEIF ( TRIM(ivar)=='KPP_QLIQ' ) THEN
          MyRF%QLIQ_ID = I
       ELSEIF ( TRIM(ivar)=='KPP_QICE' ) THEN
          MyRF%QICE_ID = I
       ELSEIF ( INDEX(ivar,'KPP_BEFORE_CHEM') > 0 ) THEN
          spc  = ivar(17:len(ivar))
          MyRF%TRCid(I) = 0
          DO N=1,Input_Opt%N_TRACERS
             IF ( TRIM(spc)==TRIM(Input_Opt%TRACER_NAME(N)) ) THEN
                MyRF%TRCid(I) = N
                EXIT
             ENDIF 
          ENDDO
          IF ( MyRF%TRCid(I) == 0 ) THEN
             IF(am_I_Root) WRITE(*,*) 'Warning: species not in KPP: ',TRIM(spc)
          ENDIF 
       ELSEIF ( INDEX(ivar,'KPP_JVAL') > 0 ) THEN
          READ(ivar(10:12),*) MyRF%JVALid(I)
       ELSE
          IF(am_I_Root) WRITE(*,*) 'Error: do not understand ',TRIM(ivar)
          RC = GIGC_FAILURE
          EXIT 
       ENDIF
       ! Add feature name for reference
       MyRF%feature_names(I) = trim(ivar)
    ENDDO
    IF ( RC==GIGC_FAILURE ) RETURN

    ! Close file
    CLOSE(lun)

    ! Verbose
    IF ( am_I_Root ) THEN
       write(*,*) 'Random forest: # of input features: ',nft
       DO I=1,nft
          write(*,*) 'feature ',I,trim(MyRF%feature_names(I)),MyRF%TRCid(I),MyRF%JVALid(I)
       ENDDO 
       write(*,*) 'Will use NOx output? ',DoNOx
    ENDIF

    !=======================================================================
    ! Define random forest for every species
    !=======================================================================
    ALLOCATE(MyRF%SpcRF(MyRF%NSPEC))
    DO N = 1, MyRF%NSPEC
       ! Initialize variables & arrays to empty values
       MyRF%SpcRF(N)%hasRF       = .FALSE.
       MyRF%SpcRF(N)%srcFile     = '' 
       MyRF%SpcRF(N)%KPP_SpcID   = -1 
       MyRF%SpcRF(N)%KPP_SpcName = '' 
       MyRF%SpcRF(N)%GC_SpcID    = -1 
       MyRF%SpcRF(N)%GC_SpcName  = '' 
       MyRF%SpcRF(N)%out_type    = -1
       MyRF%SpcRF(N)%n_trees     = -1
       MyRF%SpcRF(N)%corr        =  1.0
       MyRF%SpcRF(N)%forest      => NULL()
    
       ! Extract species name
       IF ( N <= Input_Opt%N_TRACERS ) THEN 
          spc = TRIM(Input_Opt%TRACER_NAME(N))
       ELSEIF ( N == Input_Opt%N_TRACERS + 1 ) THEN
          spc = 'NOx' 
       ELSEIF ( N == Input_Opt%N_TRACERS + 2 ) THEN
          spc = 'NOratio'
       ENDIF

       ! Always fill GC species name
       MyRF%SpcRF(N)%GC_SpcName  = TRIM(spc)

       ! -----------------------
       ! Read header 
       srcFile = TRIM(srcFileRoot)//TRIM(spc)//'/forest_'//TRIM(spc)//'_header.txt' 
       INQUIRE(file=trim(srcFile),exist=exists)
       IF ( .NOT. exists ) THEN
          IF ( am_I_Root ) &
             WRITE(*,*) 'Skip ',TRIM(spc),' - header file not found: ',TRIM(srcFile)
          CYCLE
       ENDIF
       ! Open file
       IF ( am_I_Root ) write(*,*) 'Reading ',TRIM(srcFile)
       lun = FindFreeLun()
       OPEN(unit=lun,file=TRIM(srcFile))
       ! Double-check # of features
       READ(lun,'(a)',IOSTAT=IOS) line
       IF ( IOS < 0 ) THEN
          CLOSE(lun)
          CYCLE 
       ENDIF
       READ(line(10:len(line)),*) tmp 
       IF ( tmp /= myRF%n_features ) THEN
          IF ( am_I_Root ) &
             WRITE(*,*) 'Features /= n_features: ',TRIM(spc),'::',tmp,myRF%n_features
          CLOSE(lun)
          CYCLE
       ENDIF
       ! Double-check species (predictor) 
       READ(lun,'(a)',IOSTAT=IOS) line
       IF ( IOS < 0 ) THEN
          CLOSE(lun)
          CYCLE 
       ENDIF
       IF ( TRIM(ADJUSTL(line(12:len(line)))) /= TRIM(spc) ) THEN
          IF ( am_I_Root ) &
             WRITE(*,*) 'Species mismatch ',TRIM(spc),TRIM(ADJUSTL(line(12:len(line))))
          CLOSE(lun)
          CYCLE
       ENDIF
       ! Get output operation (output type) 
       READ(lun,'(a)',IOSTAT=IOS) line
       IF ( IOS < 0 ) THEN
          CLOSE(lun)
          CYCLE 
       ENDIF
       tpe = TRIM(ADJUSTL(line(8:len(line))))
       IF ( trim(tpe)/="T" .AND. trim(tpe)/="L" .AND. trim(tpe)/="C" .AND. trim(tpe)/="D" ) THEN
          WRITE(*,*) 'Wrong output type - must be T, L, or C: ',TRIM(tpe),'::',TRIM(line)
          RETURN 
       ENDIF

       ! Get number of trees
       READ(lun,'(a)',IOSTAT=IOS) line
       IF ( IOS < 0 ) THEN
          CLOSE(lun)
          CYCLE 
       ENDIF
       idx = INDEX(trim(line),'Trees:')
       IF ( idx <= 0 ) THEN
          IF ( am_I_Root ) &
             WRITE(*,*) '`Trees:` not found - error reading ',TRIM(srcFile),'::',TRIM(line)
          RC = GIGC_FAILURE
          RETURN 
       ENDIF
       READ(line(idx+6:len(line)),*) ntr

       ! Get max. number of nodes 
       READ(lun,'(a)',IOSTAT=IOS) line
       IF ( IOS < 0 ) THEN
          CLOSE(lun)
          CYCLE 
       ENDIF
       idx = INDEX(trim(line),'Nodes:')
       IF ( idx <= 0 ) THEN
          IF ( am_I_Root ) &
             WRITE(*,*) '`Nodes:` not found - error reading ',TRIM(srcFile),'::',TRIM(line)
          RC = GIGC_FAILURE
          RETURN 
       ENDIF
       READ(line(idx+6:len(line)),*) nnd 

       ! Get correction factor 
       corr = 1.0
       READ(lun,'(a)',IOSTAT=IOS) line
       IF ( IOS < 0 ) THEN
          CLOSE(lun)
          CYCLE 
       ENDIF
       idx = INDEX(trim(line),'corr_factor:')
       IF ( idx <= 0 ) THEN
          IF ( am_I_Root ) THEN 
             WRITE(*,*) '`corr_factor:` not found - will use correction factor of 1.0'
          ENDIF
       ELSE 
           READ(line(idx+13:len(line)),*) corr 
       ENDIF

       ! Close file
       CLOSE(lun)

       ! ----------------------------
       ! Fill type
       IF ( max_trees > 0 ) THEN
          ntruse = min(max_trees,ntr)
       ELSE
          ntruse = ntr
       ENDIF
       MyRF%SpcRF(N)%n_trees     = ntruse
       MyRF%SpcRF(N)%srcFile     = TRIM(srcFileRoot)
       MyRF%SpcRF(N)%KPP_SpcID   = N
       MyRF%SpcRF(N)%KPP_SpcName = TRIM(spc)
       GCID                      = N
       MyRF%SpcRF(N)%GC_SpcID    = GCID 
       MyRF%SpcRF(N)%GC_SpcName  = TRIM(spc)
       SELECT CASE ( TRIM(tpe) ) 
          CASE ( "T" )
             MyRF%SpcRF(N)%out_type = OutTypeTend
          CASE ( "L" )
             MyRF%SpcRF(N)%out_type = OutTypeLogc
          CASE ( "C" )
             MyRF%SpcRF(N)%out_type = OutTypeConc
          CASE ( "D" )
             MyRF%SpcRF(N)%out_type = OutTypeDrop
       END SELECT
       MyRF%SpcRF(N)%corr        = corr   
       ! Find index in input vector
       MyRF%SpcRF(N)%IdxInCin = -1
       DO I=1,MyRF%n_features
          IF ( MyRF%TRCId(I) == N ) THEN
             MyRF%SpcRF(N)%IdxInCin = N
             EXIT
          ENDIF
       ENDDO
       IF ( MyRF%SpcRF(N)%IdxInCin < 0 ) THEN
          IF ( am_I_Root ) THEN
             write(*,*) 'Warning: species is not element of input vector: ',N,GCID,TRIM(spc)
          ENDIF
       ENDIF 

       ! ----------------------------
       ! Create all trees 
       CALL CreateTrees( am_I_Root, N, spc, nnd, ntr, RC )
       IF ( RC /= GIGC_SUCCESS ) EXIT

       ! ----------------------------
       ! Verbose
       IF ( am_I_Root ) THEN
          write(*,*) 'Defined random forest for species ',TRIM(spc),': '
          write(*,*) '# of features    : ',nft
          write(*,*) 'Max. # of nodes  : ',nnd
          write(*,*) '# of trees       : ',ntruse
          write(*,*) 'index in input   : ',MyRF%SpcRF(N)%IdxInCin
          write(*,*) 'output type      : ',TRIM(tpe)
          write(*,*) 'correction factor: ',MyRF%SpcRF(N)%corr
          write(*,*) 'forests read from: ',TRIM(srcFileRoot) 
       ENDIF

       ! Mark this forest as being active
       MyRF%SpcRF(N)%hasRF = .TRUE.
    ENDDO

    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE Init_Chml
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: CreateTrees
!
! !DESCRIPTION: Helper routine to create trees of a forest
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE CreateTrees( am_I_Root, SpcID, spc, nnd, ntr, RC )
!
! !USES:
!
    USE inquireMod,       ONLY : findFreeLUN
!
! !INPUT PARAMETERS:
!    
    LOGICAL,          INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
    INTEGER,          INTENT(IN)  :: SpcID
    CHARACTER(LEN=*), INTENT(IN)  :: spc
    INTEGER,          INTENT(IN)  :: nnd 
    INTEGER,          INTENT(IN)  :: ntr 
!
! !INPUT/OUTPUT PARAMETERS:
!    
    INTEGER,          INTENT(OUT) :: RC          ! Success or failure?
!    
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Local variables
    CHARACTER(LEN=1023)    :: srcFile
    INTEGER                :: itree
    INTEGER                :: I, J
    INTEGER                :: nid, lid
    INTEGER                :: ntruse
    INTEGER                :: oleft, nleft, oright, nright, inode
    INTEGER,  ALLOCATABLE  :: newid(:), isnode(:)
    INTEGER,  ALLOCATABLE  :: left(:,:)
    INTEGER,  ALLOCATABLE  :: right(:,:)
    INTEGER,  ALLOCATABLE  :: feature(:,:)
    REAL(fp), ALLOCATABLE  :: threshold(:,:)
    REAL(fp), ALLOCATABLE  :: value(:,:)

    !=================================================================
    ! CreateTrees begins here!
    !=================================================================

    ! ----------------------------
    ! Read left leaf indeces
    srcFile = TRIM(srcFileRoot)//TRIM(spc)//'/forest_'//TRIM(spc)//'_lefts.txt' 
    ALLOCATE(left(nnd,ntr))
    left = 0
    CALL ReadTxtArrI( am_I_Root, srcfile, nnd, ntr, left, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ! ----------------------------
    ! Read right leaf indeces
    srcFile = TRIM(srcFileRoot)//TRIM(spc)//'/forest_'//TRIM(spc)//'_rights.txt' 
    ALLOCATE(right(nnd,ntr))
    right = 0
    CALL ReadTxtArrI( am_I_Root, srcfile, nnd, ntr, right, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ! ----------------------------
    ! Read feature indeces 
    srcFile = TRIM(srcFileRoot)//TRIM(spc)//'/forest_'//TRIM(spc)//'_features.txt' 
    ALLOCATE(feature(nnd,ntr))
    feature = 0
    CALL ReadTxtArrI( am_I_Root, srcfile, nnd, ntr, feature, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ! ----------------------------
    ! Read thresholds 
    srcFile = TRIM(srcFileRoot)//TRIM(spc)//'/forest_'//TRIM(spc)//'_thresholds.txt' 
    ALLOCATE(threshold(nnd,ntr))
    threshold = 0.0_fp
    CALL ReadTxtArrF( am_I_Root, srcfile, nnd, ntr, threshold, RC=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ! ----------------------------
    ! Read values 
    srcFile = TRIM(srcFileRoot)//TRIM(spc)//'/forest_'//TRIM(spc)//'_values.txt' 
    ALLOCATE(value(nnd,ntr))
    value = 0.0_fp
    CALL ReadTxtArrF( am_I_Root, srcfile, nnd, ntr, value, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ! Allocate forest vector
    ALLOCATE( MyRF%SpcRF(SpcID)%forest(ntr) )
    ALLOCATE( newid(nnd),isnode(nnd)        )

    ! Do for all trees in forest
    !DO itree = 1, ntr
    IF ( max_trees > 0 ) THEN
       ntruse = min(max_trees, ntr )
    ELSE
       ntruse = ntr
    ENDIF
    DO itree = 1, ntruse

       ! Initialize values
       MyRF%SpcRF(SpcID)%forest(itree)%node_left      => NULL()
       MyRF%SpcRF(SpcID)%forest(itree)%node_right     => NULL()
       MyRF%SpcRF(SpcID)%forest(itree)%node_feature   => NULL()
       MyRF%SpcRF(SpcID)%forest(itree)%node_threshold => NULL()
       MyRF%SpcRF(SpcID)%forest(itree)%leaf_value     => NULL()
 
       ! Helper arrays. Put condensed vectors into those, then truncate
       ! when storing
       newid(:)  = -1
       isnode(:) = .false.
   
       ! Initialize indeces
       nid = 0
       lid = 0
  
       ! Walk through all values and update features. Remember that all python
       ! indeces are zero-indexed!
       DO I=1,nnd
          ! Is this a node?      
          IF ( left(I,itree) /= right(I,itree) ) THEN
             nid       = nid + 1
             newid(I)  = nid
             isnode(I) = .true.
   
          ! Is this a leaf?
          ELSEIF ( left(I,itree)==-1 .AND. right(I,itree)==-1 ) THEN
             lid       = lid + 1
             newid(I)  = lid
          ENDIF
       ENDDO
   
       ! Set number of nodes and allocate vectors 
       MyRF%SpcRF(SpcID)%forest(itree)%n_nodes = nid
       MyRF%SpcRF(SpcID)%forest(itree)%n_leafs = lid
       ALLOCATE( MyRF%SpcRF(SpcID)%forest(itree)%node_left(nid)      )
       ALLOCATE( MyRF%SpcRF(SpcID)%forest(itree)%node_right(nid)     )
       ALLOCATE( MyRF%SpcRF(SpcID)%forest(itree)%node_feature(nid)   )
       ALLOCATE( MyRF%SpcRF(SpcID)%forest(itree)%node_threshold(nid) )
       ALLOCATE( MyRF%SpcRF(SpcID)%forest(itree)%leaf_value(lid)     )

       ! Fill condensed vectors
       DO I=1,nnd
           IF ( newid(I) < 0 ) CYCLE
   
           ! Node
           IF ( isnode(I) ) THEN
   
              ! Get original indeces. Shift by one because python is zero-indexed.
              oleft  = left(I,itree) + 1
              oright = right(I,itree) + 1
   
              ! Get new indeces
              nleft  = newid(oleft)
              nright = newid(oright)
   
              ! Adjust index if it belongs to a leaf. Leaf indeces are adjusted by
              ! number of nodes
              IF ( .not. isnode(oleft ) ) nleft  = nleft  + nid
              IF ( .not. isnode(oright) ) nright = nright + nid
   
              ! Fill vectors 
              MyRF%SpcRF(SpcID)%forest(itree)%node_left(     newid(I)) = nleft
              MyRF%SpcRF(SpcID)%forest(itree)%node_right(    newid(I)) = nright
              MyRF%SpcRF(SpcID)%forest(itree)%node_feature(  newid(I)) = feature(  I,itree) + 1
              MyRF%SpcRF(SpcID)%forest(itree)%node_threshold(newid(I)) = threshold(I,itree)
   
           ! Leaf
           ELSE
              MyRF%SpcRF(SpcID)%forest(itree)%leaf_value(newid(I))     = value(I,itree)
   
           ENDIF
       ENDDO
   
       ! Some verbose for testing
       IF ( am_I_Root ) THEN
          write(*,*) 'Created new tree: ',SpcID, itree, ntruse, &
             MyRF%SpcRF(SpcID)%forest(itree)%n_nodes, MyRF%SpcRF(SpcID)%forest(itree)%n_leafs
       ENDIF
    ENDDO ! itree  
 
    ! Cleanup
    DEALLOCATE(left,right,feature,threshold,value)
    DEALLOCATE(newid,isnode)
   
    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE CreateTrees
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadTxtArr
!
! !DESCRIPTION: Helper routine to read txt array 
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE ReadTxtArrI( am_I_Root, ifile, nnd, ntr, oarr_i, RC )  
!
! !USES:
!
    USE inquireMod,       ONLY : findFreeLUN
!
! !INPUT PARAMETERS:
!    
    LOGICAL,          INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
    CHARACTER(LEN=*), INTENT(IN)  :: ifile
    INTEGER,          INTENT(IN)  :: nnd 
    INTEGER,          INTENT(IN)  :: ntr 
!
! !OUTPUT PARAMETERS:
!    
    INTEGER,  INTENT(OUT)           :: oarr_i(nnd,ntr)
    INTEGER,  INTENT(OUT)           :: RC          ! Success or failure?
!    
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                :: I, IOS, lun, tmp
    CHARACTER(LEN=8191)    :: line
    LOGICAL                :: exists 

    !=================================================================
    ! ReadTxtArr begins here!
    !=================================================================

    ! Assume failure
    RC = GIGC_FAILURE

    ! ----------------------------
    ! Read features 
    INQUIRE(file=trim(ifile),exist=exists)
    IF ( .NOT. exists ) THEN
       IF ( am_I_Root ) WRITE(*,*) 'File not found: ',trim(ifile) 
       RETURN
    ENDIF
    ! Initialize array to be filled 
    oarr_i = 0

    ! Open file
    IF ( am_I_Root ) write(*,*) 'Reading ',TRIM(ifile)
    lun = FindFreeLun()
    OPEN(unit=lun,file=TRIM(ifile))
    ! Read header
    READ(lun,'(a)',IOSTAT=IOS) line
    IF ( IOS < 0 ) THEN
       IF ( am_I_Root ) write(*,*) 'Error reading header (A): ',TRIM(ifile)
       CLOSE(lun)
       RETURN 
    ENDIF
    READ(line(9:len(line)),*) tmp
    IF ( tmp /= ntr ) THEN
     IF (am_I_Root) write(*,*) 'ntrees /= n_trees: ',TRIM(ifile),ntr,tmp
       CLOSE(lun)
       RETURN 
    ENDIF
    READ(lun,'(a)',IOSTAT=IOS) line
    IF ( IOS < 0 ) THEN
       IF ( am_I_Root ) write(*,*) 'Error reading header (B): ',TRIM(ifile)
       CLOSE(lun)
       RETURN 
    ENDIF
    READ(line(9:len(line)),*) tmp
    IF ( tmp /= nnd ) THEN
       IF (am_I_Root) write(*,*) 'nnodes /= n_nodes: ',TRIM(ifile),nnd,tmp
       CLOSE(lun)
       RETURN
    ENDIF
    ! Read line by line
    DO I =1,nnd 
       READ(lun,'(a)',IOSTAT=IOS) line
       IF ( IOS < 0 ) EXIT
       READ(line,*) oarr_i(I,:)
    ENDDO
    IF ( IOS < 0 ) THEN
       IF(am_I_Root) write(*,*) 'Encountered eof in ',TRIM(ifile)
       CLOSE(lun)
       RETURN
    ENDIF

    ! Close
    CLOSE(lun)

    ! Return 
    RC = GIGC_SUCCESS

  END SUBROUTINE ReadTxtArrI
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadTxtArr
!
! !DESCRIPTION: Helper routine to read txt array 
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE ReadTxtArrF( am_I_Root, ifile, nnd, ntr, oarr_f, RC )  
!
! !USES:
!
    USE inquireMod,       ONLY : findFreeLUN
!
! !INPUT PARAMETERS:
!    
    LOGICAL,          INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
    CHARACTER(LEN=*), INTENT(IN)  :: ifile
    INTEGER,          INTENT(IN)  :: nnd 
    INTEGER,          INTENT(IN)  :: ntr 
!
! !OUTPUT PARAMETERS:
!    
    REAL(fp), INTENT(OUT)           :: oarr_f(nnd,ntr)
    INTEGER,  INTENT(OUT)           :: RC          ! Success or failure?
!    
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                :: I, IOS, lun, tmp
    CHARACTER(LEN=8191)    :: line
    LOGICAL                :: exists 

    !=================================================================
    ! ReadTxtArr begins here!
    !=================================================================

    ! Assume failure
    RC = GIGC_FAILURE

    ! ----------------------------
    ! Read features 
    INQUIRE(file=trim(ifile),exist=exists)
    IF ( .NOT. exists ) THEN
       IF ( am_I_Root ) WRITE(*,*) 'File not found: ',trim(ifile) 
       RETURN
    ENDIF
    ! Initialize array to be filled 
    oarr_f = 0.0_fp

    ! Open file
    IF ( am_I_Root ) write(*,*) 'Reading ',TRIM(ifile)
    lun = FindFreeLun()
    OPEN(unit=lun,file=TRIM(ifile))
    ! Read header
    READ(lun,'(a)',IOSTAT=IOS) line
    IF ( IOS < 0 ) THEN
       IF ( am_I_Root ) write(*,*) 'Error reading header (A): ',TRIM(ifile)
       CLOSE(lun)
       RETURN 
    ENDIF
    READ(line(9:len(line)),*) tmp
    IF ( tmp /= ntr ) THEN
     IF (am_I_Root) write(*,*) 'ntrees /= n_trees: ',TRIM(ifile),ntr,tmp
       CLOSE(lun)
       RETURN 
    ENDIF
    READ(lun,'(a)',IOSTAT=IOS) line
    IF ( IOS < 0 ) THEN
       IF ( am_I_Root ) write(*,*) 'Error reading header (B): ',TRIM(ifile)
       CLOSE(lun)
       RETURN 
    ENDIF
    READ(line(9:len(line)),*) tmp
    IF ( tmp /= nnd ) THEN
       IF (am_I_Root) write(*,*) 'nnodes /= n_nodes: ',TRIM(ifile),nnd,tmp
       CLOSE(lun)
       RETURN
    ENDIF
    ! Read line by line
    DO I =1,nnd 
       READ(lun,'(a)',IOSTAT=IOS) line
       IF ( IOS < 0 ) EXIT
       READ(line,*) oarr_f(I,:)
    ENDDO
    IF ( IOS < 0 ) THEN
       IF(am_I_Root) write(*,*) 'Encountered eof in ',TRIM(ifile)
       CLOSE(lun)
       RETURN
    ENDIF

    ! Close
    CLOSE(lun)

    ! Return 
    RC = GIGC_SUCCESS

  END SUBROUTINE ReadTxtArrF 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_chml
!
! !DESCRIPTION: Subroutine Cleanup\_Chml deallocates module variables.
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE Cleanup_Chml( am_I_Root, RC )  
!
! !USES:
!
!
! !INPUT PARAMETERS:
!    
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!    
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!    
! !REVISION HISTORY:
!  18 Feb 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: I, N

    !=================================================================
    ! Cleanup_Chml begins here!
    !=================================================================

    IF ( .NOT. ASSOCIATED(MyRF) ) RETURN

    IF ( MyRF%NSPEC > 0 ) THEN 
       DO N=1,MyRF%NSPEC
          IF ( MyRF%SpcRF(N)%hasRF ) THEN
             DO I=1,MyRF%SpcRF(N)%n_trees
                IF ( ASSOCIATED(MyRF%SpcRF(N)%forest(I)%node_left      ) ) &
                   DEALLOCATE(MyRF%SpcRF(N)%forest(I)%node_left      )
                IF ( ASSOCIATED(MyRF%SpcRF(N)%forest(I)%node_right     ) ) &
                   DEALLOCATE(MyRF%SpcRF(N)%forest(I)%node_right     )
                IF ( ASSOCIATED(MyRF%SpcRF(N)%forest(I)%node_feature   ) ) &
                   DEALLOCATE(MyRF%SpcRF(N)%forest(I)%node_feature   )
                IF ( ASSOCIATED(MyRF%SpcRF(N)%forest(I)%node_threshold ) ) &
                   DEALLOCATE(MyRF%SpcRF(N)%forest(I)%node_threshold )
                IF ( ASSOCIATED(MyRF%SpcRF(N)%forest(I)%leaf_value     ) ) &
                   DEALLOCATE(MyRF%SpcRF(N)%forest(I)%leaf_value     )
             ENDDO
             DEALLOCATE(MyRF%SpcRF(N)%forest)
          ENDIF
       ENDDO
       IF ( ASSOCIATED( MyRF%JVALid        ) ) DEALLOCATE( MyRF%JVALid        )
       IF ( ASSOCIATED( MyRF%TRCid         ) ) DEALLOCATE( MyRF%TRCid         )
       IF ( ASSOCIATED( MyRF%feature_names ) ) DEALLOCATE( MyRF%feature_names )
    ENDIF
    DEALLOCATE(MyRF)

    ! INitialize
    RC = GIGC_SUCCESS

  END SUBROUTINE Cleanup_Chml
!EOC
END MODULE Chml_Mod
