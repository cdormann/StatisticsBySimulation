Readme to the data and programming scripts for reproducing the results from the Nature publication titled:

"New land-use-change emissions indicate a declining CO2 airborne fraction".

Authors: Margreet J. E. van Marle*, Dave van Wees*, Richard A. Houghton, Robert D. Field, Jan Verbesselt, and Guido R. van der Werf
* These authors contributed equally.

DOI: https://doi.org/10.1038/s41586-021-04376-4



This dataset includes the following (All files are preceded by "Marle_et_al_Nature_AirborneFraction_"):

- "Datasheet.xlsx": Excel dataset containing all annual and monthly emissions and CO2 time series used for the analysis, and the resulting airborne fraction time series.

- "Script.py": 
BEFORE RUNNING THE SCRIPT: change the 'wdir' variable to the directory containing the provided script and files.
NOTE: This script requires the Python module: 'pymannkendall'
Python script used for reproducing the results and figures from the paper. The provided Datasheet.xlsx file and the .zip and .npz files are required for this program. In case all these files are found by the script, it should run within several seconds. Successful execution of the script will save Figures 1-4 from the main text and print the data from Table 1. In case script execution takes longer, please check if the .xlsx, .zip and .npz files are correctly present in the assigned 'wdir' directory. Otherwise the script will start recalculating these files, which might take a while (see notes below).

- "MC10000_MK_ts_TRENDabs.zip": .zip file containing all results from the Monte-Carlo simulation for trend estimation for Figure 3 (calculated using Python function 'calc_AF_MonteCarlo()'). This .zip file contains multiple .npz files for different emission scenarios and data treatments. This .zip file is managed by the Python script function 'calc_AF_MonteCarlo_filemanager()', there is no need to unzip the file manually. In case the .zip file is not found by the Python script (e.g. because the .zip file was unpacked manually and deleted), the program will start recalculating and save a new .zip file. This can take several minutes dependent on the computer used. Recalculated results could differ very slightly due to the random factor in the Monte-Carlo approach, even though the 10,000 iterations bring this variation to a minimum.

- "MC1000_MK_run50x50_TRENDabs.npz": .npz file containing the Monte-Carlo results used for producing Figure 4 (calculated using Python function 'calc_AF_MonteCarlo_ARR()'). In case the .npz file is not found by the Python script (e.g. because it was deleted or not downloaded), the program will start recalculating and save a new file. This can take around 30 hours(!) dependent on the computer used. Recalculated results could differ slightly due to the random factor in the Monte-Carlo approach.

- "tol_colors.py": Additional Python module used in script.py, required for producing the colors used in the Main text figures. Source: https://personal.sron.nl/~pault/

- Figure files: Figures 1-4 from the Main text saved as .pdf files. Figure 3 is saved as three independent panels. The Figures are also reproduced by script.py if executed successfully.
