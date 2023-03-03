################################################################################
## Author: Elisa Warner
## Adapted and improved upon from the work of: Shariq Mohammed
## Date: 05/11/2022
## GET SLICE NUMBERS
## Description: This code will find the slice of the MRI which corresponds to the
##  largest mask and save it all iin a file called slice_number.RData
################################################################################

# Set the working directory
setwd("Z:/.../1Pseudoprogression/UPenn-UMich Code Sharing")
folder.name = getwd()

#### LOAD PACKAGES
packages = c("WhiteStripe", "oro.nifti", "doParallel", "R.utils", "stringr", "matlabr",
             "R.matlab", "sp", "GWmodel", "RNifti")
package.check = lapply(packages,
                       function(x){
                         if(!x %in% rownames(installed.packages())) {
                           install.packages(x, dependencies = TRUE)
                         }
                         library(x, character.only = TRUE)
                       })

rm(package.check,packages)

############################# FIND BEST SLICES
# this opens up the mask and takes the largest blob

##### Location of successful files document
sub.folders = read.table('./Data/ws_success.txt',row.names = NULL, sep=',')$x
slice.number = rep(NA, length(sub.folders))
names(slice.number) = sub.folders

# unique names
unique = str_count(sub.folders,"[-+]")
names(unique) = sub.folders

# loop to find slice
for(s in sub.folders){
  print(s)
  
  if (unique[s] == 2){
    s_base = substr(s, start = 1, stop = nchar(s)-2)
  }
  else{
    s_base = s
  }
  
  # folder for MRI data
  save_fname = file.path(folder.name, 'MRI_full', s)

  mask = readNIfTI(paste0(save_fname, '/label_rs.nii.gz'), reorient = FALSE)@.Data
  writeMat(con = paste0(save_fname,'/temp.mat'), mask = mask, wd = save_fname)
  rm(mask)
  run_matlab_code(c(paste0('load(',"'", save_fname, '/temp.mat',"'",');'),
                    "nslice = size(mask,3);",
                    "largestblobsize = [];",
                    "for l = 1:nslice",
                    "BW = mask(:,:,l)~=0;",
                    "largestBW = bwareafilt(BW, 1);",
                    "largestblobsize(l) = sum(sum(largestBW));",
                    "end",
                    "l_final = find(largestblobsize == max(largestblobsize));",
                    "if size(l_final,2)>1",
                    "l_final = l_final(1);",
                    "end",
                    "BW = mask(:,:,l_final)~=0;",
                    "largestBW = bwareafilt(BW, 1);",
                    paste0('save(', "'", save_fname, '/temp.mat', "'",
                           ", 'largestBW','l_final');")),
                  endlines = F, verbose = F)
  mask.d = readMat(con = paste0(save_fname,'/temp.mat'))
  mask = mask.d$largestBW
  l.final = mask.d$l.final
  rm(mask.d)
  
  ################# 01/04/23 for other slices
  if (unique[s] == 2){
    minus_sign = (substr(s, start = nchar(s)-1, stop = nchar(s)-1) == "-")
    add_amount = strtoi(substr(s, start=nchar(s), stop=nchar(s)))
    
    l.final = l.final + (-2 * minus_sign + 1) * add_amount
  }
  
  # add to slice number list
  slice.number[s] = l.final
}
rm(s)


# save slice
############
save(slice.number, sub.folders, folder.name, file = paste0(folder.name, "/slice_number.RData"))
load(paste0(folder.name, "/slice_number.RData"))

