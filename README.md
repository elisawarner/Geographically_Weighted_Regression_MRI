# README

**Author:** Elisa Warner  
**Date:** 01/03/2023  
**Questions?** Please email elisawa@umich.edu  

# Description
This folder is part of a project to detect Pseudoprogression v True Progression

# Workflow
1. Download the raw data from Nick's Dropbox. Download the `\_rs\_norm` file. This should be the same as his matlab version but will include t2, t1post, t1pre. These files are additionally whitestriped already, so no whitestripe needed.   
2. Run `MoveFiles.ipynb` to copy over the psuedoprogression cases 5 times.  
3. Run `Edit_RUNFILE.ipynb` in order to create a proper `ws_success.txt` file. This file is created in the original code to document successful whitestripe. It will break in this dataset where there are too many 0s in the FLAIR file. The file `ws_success_ref.txt` was created based on an analysis of the successful files already, so that is used as a reference to create thh ws_success.txt file.  
4. Run `runSliceGWR.R` in order to create the slice file. It specially processes the copied folders and takes either up to 2 slices above or below the original best slice. Once this is saved you don't need to try to figure out the best slice anymore. So you don't need to run the Matlab code.  
5. Run `Create_NIFTI.R` to create (t1post - t1pre) using the `rs_norm` files.  
6. Run `EXTRACT_SLICES.R` to create the spdf files which are required for the GWR run.  
7. Run `RUN_RESD.R` to run GWR and extract the residuals from the slices. Note that this cannot be done on GreatLakes because of GWmodel package requirements.  
8. Run `RUN_LOOCV.R` to run the leave-one-out cross validation procedure for all patients.  

# Contents
1. `Create_NIFTI.R` [R]: This code is supposed to take the raw .nii files and create T1post-T1pre and stuff 
2. `Data` [folder] : (not in github code) contains the MRI data. Also contains the `ws_success.txt` document.  
3. `Data_Full42_2023`[folder] : contains the final residuals for this project.  
4. `Edit_RUNFILE.ipynb` [Jupyter NB]:  The purpose of this file is the edit `Data/ws_success.txt` to define which patients will be included and therefore from which patients the largest mask slice will be extracted.    
5. `EXTRACT_SLICES.R` [R] : This takes the saved slice numbers and then creates an spdf file of the slice.  
6. `density_features.R` [R]: a type of geomPCA with a train/test.  
7. `density_featuresOG.R` [R] : used to calculate the geometric PCs.   
8. `gwr.R` [R]: contains the code for GWR. If you have GWmodel, this shouldn't be used.    
9. `LargestCC.m` [Matlab]: contains the Matlab code to find the mask slice with the largest area.  
10. `RUN_RESD.R` [R] : contains only the code which runs GWR.  
11. `Run_LOOCV.R` [R] : The code for running patient-wise leave-one-out cross validation.  
12. `Run_LOOCV_ALL.R` [R] : The code for running patient-wise LOOCV but combining modality pairs together.  
13. `runSliceGWR.R` [R] : This code is designed to only extract the slice numbers for each file. It has been adapted to the -1/+1 extra folders that I added for the psuedoprogression cases.    
14. `slice_number.RData` [R]: This is a saved file for the proper slice number. An easy solution is just choose the slices 2 above and 2 below this number    
15. `whitestripe.R` [R]: contains only the code for conducting whitestripe according to Shariq's methods  
16. `ws_norm.R` [R]: Same as whitestripe.R code above but newer. I believe this was more for debugging/testing purposes because it has more commentary in it. It also does not normalize them all at once.  
17. `Data/ws_success.txt` [Txt] : This is an example text of how the ws_success.txt document should appear. This document is a list of files which have successfully passed whitestripe and therefore won't break later on in the code.  
18. `Naive Bayes.ipynb` [Jupyter NB]: NOTE THIS CODE IS AN EXAMPLE, DON'T TRY TO RUN OFF THE BAT. You may have a lot of errors because I ran a lot of different ideas here. This code should be used as a reference. Take useful elements and throw out what you don't like. Don't use it trying to run from start to finish. Or ask me what you want to do and I'll try to help.  

# MRI files
1. `<id>_<modality>.nii.gz` : this is the raw file.  
2. `<id>_<modality>_rs.nii.gz` : this is the resized file.  
3. `<id>_<modality>_rs_norm.nii.gz` : this is the normalized resized file.  
4. `<id>_<modality>_rs_ws.nii.gz` : this is the whitestripe mask.  

# How to structure your data
1. Create an MRI folder.  
2. Inside the folder should be folders which contain the ID of every patient.  
3. Inside those folders should be the MRI files listed above.  