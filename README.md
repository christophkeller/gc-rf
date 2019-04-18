# gc-rf
Random forest emulator for GEOS-Chem chemistry mechanism, as described in Keller and Evans (2019).

The workflow is as follows:
1. Run GEOS-Chem model with full chemistry, output all concentrations before and after chemistry into netCDF files (the modified chemistry module in Fortran/chemistry_mod.F can be used as guide for this). Some sample netCDF output is provided in data/nc4.
2. Generate pickled files from GEOS-Chem netCDF output using nc2pkl.py:
    python nc2pkl.py
3. Separate into training and validation data 
    python split_data.py
4. Train random forests for every species using train_rf.py, e.g.:
    python train_rf.py -n 41 -o 'tend'    
5. Run GEOS-Chem model with random forest. Module Fortran/chml_mod.F90 contains the Fortran code to read the random forest information created in step 4 and make invdividual predictions. An example of how this can be integrated in GEOS-Chem is provided in Fortran/chemistry_mod.F.

Notes:
This code is still experimental and simplifications to the workflow / embedding into GEOS-Chem are expected. Also note that this version is not faster than the numerical solution of GEOS-Chem. 

References:
Keller, C. A. and Evans, M. J.: Application of random forest regression to the calculation of gas-phase chemistry within the GEOS-Chem chemistry model v10, Geosci. Model Dev., 12, 1209-1225, https://doi.org/10.5194/gmd-12-1209-2019, 2019.

Contact:
christoph.a.keller@nasa.gov
