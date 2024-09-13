# specter_optimization

## Requirements

- [sd_foregrounds_optimize](https://github.com/asabyr/sd_foregrounds_optimize) (specter branch, but can be easily modified to work with the main)
- [bolocalc-space](https://github.com/csierra2/bolocalc-space)
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)
- [pandas](https://pandas.pydata.org/)

## Description 
This code was used to compute the optimized bands/detectors/calibration requirements for SPECTER and test different sky models: [/ini_files](/ini_files) contains an example configuration file that can be used to find optimized bands, SNRs for different detector combinations, calibration requirements and test different sky models; [/job_files](/job_files) includes example job scripts to run the calculations. Example sensitivity/band files are in [/files](/files). 

Before running any optimizations, it may be helpful to clean bolocalc-space directory with:

```python code/clean_bolocalc_dir.py {DIR_TO_BOLOCALC-SPACE}/Experiments/specter_v1```


1. To run band optimization -- 

      ```python run_band_optimization.py {FILE_NAME}.ini```

   Make sure to set the "lowest frequency edge (bolo)=0" if only using HEMTS and to set "hemt amps=False" in the .ini file if using only bolometers).
   The run will output files {FILE_PREFIX}_optimized_bands.txt and {FILE_PREFIX}_combined_bands.txt in [/files](/files).

2. To run detector optimization -- 

   I. First, pre-compute the instantenous detector noise and save a file with detector combinations by running:

      ```python setup_detector_optimization.py {FILE_NAME}.ini {RUN#}``` (i.e. detector optimization section number)

   This will create files: {FILE_PREFIX}_raw_sensitivity.txt and {FILE_PREFIX}_ndet_combinations_run{RUN#}.txt in [/files](/files).
   
   II. Now you can compute SNRs for each detector combination via

      ```python run_detector_optimization.py {FILE_NAME}.ini {RUN#}```

   The run will output: {FILE_PREFIX}_SNRs_run{RUN#}_CORE#.txt file(s) in [/files](/files). 

2. To run calibration/bias calculations --

   I. First, pre-compute all the combinations of calibration errors using

      ```python setup_calibration_optimization.py {FILE_NAME}.ini {RUN#}```

   This will create files: {FILE_PREFIX}_calib_combinations_run{RUN#}.txt in [/calib_files](/calib_files).

   II. Now you can compute bias based on calibration errors via 

      ```python run_calibration_optimization.py {FILE_NAME}.ini {RUN#}``` 

   The run will output: {FILE_PREFIX}Bs_run{RUN#}.txt in [/calib_files](/calib_files).

4. To compute SNRs for different sky models --

   I. Again, first precompute all possible sky variations via 

      ```python setup_skymodel_bias.py {FILE_NAME}.ini {RUN#}```
   
   This will create files: {FILE_PREFIX}_sky_combinations_run{RUN#}.txt in [/grid_files](/grid_files).

   II. Then run the calculations via

      ```python run_skymodel_bias.py {FILE_NAME}.ini {RUN#}``` 

   The run will output: {FILE_PREFIX}sky_SNR_bias_run{RUN#}.txt in [/grid_files](/grid_files). 

## Other notes:
For the calculations (2),(3),(4), it is best to split the calculations across cores/into some number of input combinations/output files N. This can be specified via the .ini file too. In this case, command line argument {RUN#} should be followed by another argument {n}, where n=0:N (see examples in [/job_files](/job_files)). Note that the project right now assumes a structure where `sd_foregrounds_optimize` and `bolocalc-space` are in a folder titled "software" and that you will need to modify the paths at the top of the scripts. For some examples on how to read the outputs, see [notebooks/examples.ipynb](notebooks/examples.ipynb).


Please feel free to contact [Alina Sabyr](as6131@columbia.edu) with any questions about the code or open a pull request.



