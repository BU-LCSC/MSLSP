# MultiSource Land Surface Phenology (MSLSP)
Algorithms to derive an operational multisource land surface phenology product from Harmonized Landsat and Sentinel-2 (HLS). Code builds upon algorithms written by Dr. Joshua Gray for the MODIS land surface phenology product. These scripts included added synthetic spectra compositing functionality derived from the spline fits.

Publications:
* Bolton, D. K., Gray, J. M., Melaas, E. K., Moon, M., Eklundh, L., & Friedl, M. A. (2020). Continental-scale land surface phenology from harmonized Landsat 8 and Sentinel-2 imagery. Remote Sensing of Environment, 240, 111685. https://doi.org/10.1016/j.rse.2020.111685
* Friedl, M. (2021). MuSLI Multi-Source Land Surface Phenology Yearly North America 30 m V011 [Data set]. NASA EOSDIS Land Processes Distributed Active Archive Center. Accessed 2024-03-06 from https://doi.org/10.5067/Community/MuSLI/MSLSP30NA.011


## Running MSLSP
### Step 1: Parameterization - [MSLSP_Parameters.json](MSLSP_Parameters.json)
Parameters for MSLSP runs can be adjusted in the parameters json file. Under "setup", we can set the imgStartYr, imgEndYr, phenStartYr, and phenEndYr to select the range over which to collect HLS data and run LSP. The remaining options should all be set as true for the first run. Options include whether or not to download imagery, process imagery, run phenology, use Fmask on S2 imagery, include Landsat, and include Sentinel-2 imagery.

Next, under "AWS" are file paths and job submission settings for submitting jobs via Amazon Web Services. Under "SCC" are file paths and job submission settings for submitting jobs via the BU Shared Compute Cluster (SCC). The SCC option is the submission method that should be used for smaller regions (<50 tiles) and for any job other than the final run. Make sure to change all file paths to match where your scripts are located and where you want to store the data. Job submission settings (i.e. numCores, numChunks, runS10) should not be changed.

Under "phenology_parameters", you can choose the vegetation_index to run with. The default is EVI2. Toggle doComposites to set whether or not the algorithm produces synthetic spectra composite images for each phenometric. All other settings within this section and the rest of the parameters file do not need to be changed. 

### Step 2: Tile list - [SCC/tileLists/](SCC/tileLists/)
Make a tileLists directory. Then, create a text file that stores the five-character HLS tile name for the tiles you want (##ABC), separating each tile name with a new line.

### Step 3: Setup of MSLSP_submitTiles_SCC.sh - [MSLSP_submitTiles_SCC.sh](SCC/MSLSP_submitTiles_SCC.sh)
Make sure the parameters variable is set to the file path where your MSLSP_Parameters.json file is located. The tileList variable is the path to the txt file created in Step 2.

### Step 4: Running MSLSP - [MSLSP_submitTiles_SCC.sh](SCC/MSLSP_submitTiles_SCC.sh)
At this point, you are ready to run MSLSP. Run MSLSP_submitTiles_SCC.sh to submit all download and LSP tasks. Keep track of tasks with "qstat -u *yourusername* " and delete tasks with "qdel *enter job-ID number*". View log files to track progress and any potential errors. The final products are output as netcdf files containing the product layers for each year.
   

## Additional Notes
Functions are contained in [MSLSP_Functions.r](MSLSP_Functions.r). The MSLSP_nonvegcomp branch includes the option to create composite imagery for non-vegetated areas without LSP. This is done by taking the average phenometrics for each image chunk and creating a +/-14 days range to gather images for mean compositing. Using this method, synthetic imagery will be gap-free where imagery is available. Vegetated areas will have synthetic imagery modeled from LSP while non-vegetated areas will have composite imagery from date rangse close to the surrounding area's LSP metrics. The compositing average method and date range can be changed within the DoNonvegComp function in MSLSP_Functions.r.
