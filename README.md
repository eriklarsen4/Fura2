# Fura2
Collection of the scripts used for calcium imaging on TMEM184B/IL31RA paper


Experiments were divided by imaging day (mouse). These experiments were given their own branches on Github in the same repository.
Each imaging day contains variable numbers of imaging runs. Each imaging run contains a script that processes the fluorescence data and its paired background fluorescence data (1 R script will process two CSVs and yield a variable objects to call in the master analysis script "Ca Img Response Quant Script").

All raw data (CSVs) and scripts not in the master branch should be stored all in one directory.

Each script can be called from the master script "Ca Img Response Quant Script" when re-defining the path to the directory where all the scripts and CSVs are stored.
All imaging scripts are pre-set to be run in their entireity in one call. However, each script should be run one at a time to acquire the necessary variables for analysis in the "Ca Img Response Quant Script". It is recommended to collapse functions and sections for organization within the analysis script. In brief, once an imaging run script has been called (e.g. 2019_10_30 A4 script), and all appropriate "Ca Img Response Quant Script" functions have been called into the environment, a subsequent call from line 124-144 will load the data to be analyzed. Running the "Analysis" function (line 381) will return putative responders of agonists for the given run. These were individually evaluated with the GGplot call from lines 391-413.

The ROIs (columns from the dataframe, FLUO) determined to be "responses" are included at the end of the "Ca Img Response Quant Script," sectioned by agonist and genotype, for corresponding magnitude acquisition. Magnitudes from each script were then copied and put into the "Ca Img Mag Export" script, from which these values are exported into a CSV for each agonist (two columns, one for each genotype). These CSVs are included in the main branch of the repository as "[agonist]_mags".

Statistical evaluation of the data is performed in the "Ca Img Response Stats Script", where weighted linear regression is performed on each agonist's percent response estimates and magnitudes.

The aggregate analysis global environment was saved as an R Global Environment file, "Ca Img Response Env." The aggregate statistics global environment was saved as an R Global Environment file, "Ca Img Response Stats Env." Both environments are included in the master branch of the repository.
