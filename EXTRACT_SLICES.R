##################################################################
# Restricted usage under the Research End-User License Agreement #
# OTT 2021-042 by The University of Michigan.                    #
# License text provided in the README.pdf file                   #
##################################################################

################################################################################
# Author and maintainer: Elisa Warner < elisawa@umich.edu > 3/3/2023      #
# For any queries please contact the author or Arvind Rao <ukarvind@umich.edu> #
################################################################################

## SET UP
setwd("Z:/.../1Pseudoprogression/UPenn-UMich Code Sharing")
#setwd("/nfs/turbo/umms-ukarvind/elisawa/1Pseudoprogression/UPenn-UMich Code Sharing")
folder.name = getwd()

packages = c("WhiteStripe", "oro.nifti", "doParallel", "stringr", "R.utils", "matlabr",
             "R.matlab", "sp", "RNifti")
package.check = lapply(packages,
                       function(x){
                         if(!x %in% rownames(installed.packages())) {
                           install.packages(x, dependencies = TRUE)
                         }
                         library(x, character.only = TRUE)
                       })
rm(package.check,packages)

subfolders = list.files(file.path(folder.name, 'MRI_full'))

#######################

pairList = list(c("t1","flair"),
                c("t1post","flair"),
                c("t2","flair"),
                c("t1","t1post"),
                c("t1","t2"),
                c("t1","t1postpre"),
                c("t1post","t2"),
                c("t1post","t1postpre"),
                c("t1","t1post"),
                c("t1postpre","flair"),
                c("t1postpre","t2"),
                c("adc","flair"),
                c("adc","t1"),
                c("adc","t1post"),
                c("adc","t1postpre"),
                c("adc","t2"))

for (i in 1:length(pairList)){
  x_pred = pairList[[i]][1]
  y_pred = pairList[[i]][2]
  print(x_pred)
  print(y_pred)

  # cycle through patients...if the file exists then delete it???
  for(s in subfolders){
    fld.name = file.path(folder.name, 'MRI_full', s)
    sub.files = list.files(fld.name)
  
    if('flair.nii.gz' %in% sub.files) file.remove(paste0(fld.name,'/flair.nii.gz'))
    if('t2.nii.gz' %in% sub.files) file.remove(paste0(fld.name,'/t2.nii.gz'))
    if('mask.nii.gz' %in% sub.files) file.remove(paste0(fld.name,'/mask.nii.gz'))
    if('spdf.RData' %in% sub.files) file.remove(paste0(fld.name,'/spdf.RData'))
    if('temp.mat' %in% sub.files) file.remove(paste0(fld.name,'/temp.mat'))
    rm(fld.name, sub.files)
  }
  rm(s)

  ############################# EXTRACT SLICE
  # extracts the slice from slicenumber file
  load("slice_number.RData") # loads slice.number, sub.folders
  subfolders = sub.folders # get folders that ran through spdf
  rm(sub.folders)

  folder.name = getwd() # overwrite slice_number.RData

  # unique names
  unique = str_count(subfolders,"[-+]")
  names(unique) = subfolders

  for(s in subfolders){
    print(s)
  
    # s_base so I can load the mask
    if (unique[s] == 2){
      s_base = substr(s, start = 1, stop = nchar(s)-2)
    }
    else{
      s_base = s
    }
  
    save_fname = file.path(folder.name, 'MRI_full', s)
    mask = readNIfTI(paste0(save_fname, '/', 'label_rs.nii.gz'), reorient = FALSE)@.Data
  
    l.final = slice.number[s]
    mask = mask[,,l.final]
  
    if (sum(mask) == 0) {next}
  
    # isloate the t2 and flair of affected area
    y_pred_img = readNIfTI(paste0(save_fname, '/', y_pred, '_rs_norm.nii.gz'), reorient = FALSE)@.Data[,,l.final]*mask
    if (x_pred == "adc") {
      x_pred_img = readNIfTI(paste0(save_fname, '/', x_pred,'_rs.nii.gz'), reorient = FALSE)@.Data[,,l.final]*mask
    }
    else {
      x_pred_img = readNIfTI(paste0(save_fname, '/', x_pred,'_rs_norm.nii.gz'), reorient = FALSE)@.Data[,,l.final]*mask
    }

    # CHECK WHAT THESE DO
    print("saving...")
    sp_coords = coordinates(which(mask!=0, arr.ind = T))
    sp_data_t2 = x_pred_img[mask!=0]
    sp_data_flair = y_pred_img[mask!=0]
  
    sp_df = data.frame(x_pred_img = sp_data_t2, y_pred_img = sp_data_flair)
    spdf = SpatialPointsDataFrame(sp_coords, sp_df)
    dims = dim(mask)
    rm(sp_coords, sp_df, sp_data_flair, sp_data_t2, y_pred_img, l.final, mask, x_pred_img)
    save(spdf, dims, file = paste0(save_fname,'/spdf','_',x_pred,'_',y_pred,'.RData'))
    rm(spdf, dims, save_fname)
  }
  rm(s)
}
