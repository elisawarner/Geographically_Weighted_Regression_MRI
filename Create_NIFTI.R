###########################################################################
##                      Create_NIFTI
##
##  Author: Elisa Warner
##  Date: 03/21/2022
##  Description: Will cycle through folders according to Shariq's preferred
##    folder structure and create a file called ******_t1pp.nii.gz, which is in
##    fact raw t1post - raw t1pre.
##  Requires: An MRI folder. And within MRI folder, folders labeled after the 
##    patID. Within the patID folder should be a file containing the text
##    t1post and t1pre, respectively.
########################################################################

library("oro.nifti")

minuend = 't1post_rs_norm'
subtrahend = 't1_rs_norm'

folderName = getwd()
subfolders = list.files(file.path(folderName, 'MRI_full'))

for (s in subfolders){
  
  subFiles= list.files(file.path(folderName, 'MRI_full', s))
  for (chStr in subFiles){
    if (grepl(minuend, chStr, fixed = TRUE)){
      x1 = readNIfTI(file.path(folderName, 'MRI_full', s, chStr), reorient = FALSE)
    }
    else if (grepl(subtrahend, chStr, fixed = TRUE)){
      x2 = readNIfTI(file.path(folderName, 'MRI_full', s, chStr), reorient = FALSE)
    }
  }
  
  x3 = x1 - x2
  outName = 't1postpre_rs_norm' # paste0(s, '_t1postpre')
  outPath = file.path(folderName, 'MRI_full', s, outName)
  writeNIfTI(x3, filename=outPath)
  rm(x1, x2, x3)
}
