This repo is for code for processing wind-wave tank data from Fabrice in Deleware.

# Scripts and what they do
## M-Files_FabMarcNovDec2014
These files are for processing the Langmuir Circulation experiments from Delaware. Many of them were written by Fabrice Veron and Marc Buckley. New files that I have added are listed below. For using files in this directory, you should add this directory and all subdirectories to the path in MATLAB.

**Main_LC_2023_Plotting.m:** Plots PIV output, allows you to flip through frames of surface images. Also draws surface following coordinate system on the PIV velocity fields. 

**CrapperOptimizedFindSurface:** More sensitive surface detection. Loops through as many filter sized as you choose.

**PIVLaserSheetPositionInIR_LON.m:** Script for determining the location of the PIV laser sheet in the transformed IR images for the longitudinal PIV cases.
