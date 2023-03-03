###################################################
### RUN WHITESTRIPE
###################################################


folder.name = getwd()

x_pred = 't1postpre'
y_pred = 'flair'

packages = c("WhiteStripe", "oro.nifti", "doParallel", "R.utils", "matlabr",
             "R.matlab", "sp", "GWmodel", "RNifti")
package.check = lapply(packages,
                       function(x){
                         if(!x %in% rownames(installed.packages())) {
                           install.packages(x, dependencies = TRUE)
                         }
                         library(x, character.only = TRUE)
                       })
rm(package.check,packages)

sub.folders = list.files(file.path(folder.name, 'MRI'))

# cycle through patients...if the file exists then delete it???
for(s in sub.folders){
  fld.name = file.path(folder.name, 'MRI', s)
  sub.files = list.files(fld.name)
  
  if('flair.nii.gz' %in% sub.files) file.remove(paste0(fld.name,'/flair.nii.gz'))
  if('t2.nii.gz' %in% sub.files) file.remove(paste0(fld.name,'/t2.nii.gz'))
  if('mask.nii.gz' %in% sub.files) file.remove(paste0(fld.name,'/mask.nii.gz'))
  if('spdf.RData' %in% sub.files) file.remove(paste0(fld.name,'/spdf.RData'))
  if('temp.mat' %in% sub.files) file.remove(paste0(fld.name,'/temp.mat'))
  rm(fld.name, sub.files)
}
rm(s)

# preprocessing function
ws_norm = function(in_fname, tumor_mask_fname, modality, out_fname_int){
  img = readNIfTI(in_fname, reorient = FALSE)
  t_mask = readNIfTI(tumor_mask_fname, reorient = FALSE)
  ws = whitestripe(img = img, type = modality, stripped = TRUE)
  tmr.ind = which(t_mask!=0)
  ws_mask = setdiff(ws$whitestripe.ind, tmr.ind)
  ws$mask.img[tmr.ind] = 0
  ws.norm = whitestripe_norm(img = img, indices = ws_mask)
  writeNIfTI(ws.norm, filename=out_fname_int)
}

########################### PREPROCESSING
# it will only pay attention to these files here
mods = c('flair','t2','t1','t1post')

write.table(numeric(), file = './Data/ws_success.txt')
for(s in sub.folders){
  print(s)
  sub.files = list.files(file.path(folder.name, 'MRI', s)) # get files in each folder again
  
  # find the file called flair and use it to get root name
  s.name = which(sapply(1:length(sub.files),
                        function(sf){
                          x = sub.files[sf]
                          substr(x, start = nchar(x)-11, stop = nchar(x))
                        }) == 'flair.nii.gz') # get idx
  chr.str = strsplit(sub.files[s.name], split = 'flair.nii.gz')[[1]]
  
  file_init = file.path(folder.name, 'MRI', s, chr.str) # the root of the filenames
  save_fname = file.path(folder.name, 'MRI', s) # path of this folder
  
  # call preprocessing
  cl = makeCluster(2) # parallel
  registerDoParallel(cl, cores = 2)
  registerDoParallel(cl, cores = 2)
  tryCatch({
    withTimeout({
      Nsim = foreach(m = 1:length(mods),
                     .packages = c("WhiteStripe", "oro.nifti")) %dopar% {
                       ws_norm(paste0(file_init, mods[m], '.nii.gz'), 
                               paste0(file_init, 'mask.nii.gz'), 'T1',
                               file.path(save_fname, mods[m]))
                     }
      
      mask = readNIfTI(paste0(file_init, 'mask.nii.gz'), reorient = FALSE)
      writeNIfTI(mask, file.path(save_fname, 'mask'))
      write.table(s, "./Data/ws_success.txt", sep = ",",
                  col.names = !file.exists("./Data/ws_success.txt"), append = T)
      rm(Nsim, mask) # del fro mem
    }, timeout = 6000)
  }, error = function(e) cat('\n', s, ': WhiteStripe Failed.', '\n', sep = ''))
  stopCluster
  rm(cl, save_fname, file_init, chr.str, s.name, sub.files)(cl)
}

rm(s, mods, ws_norm, sub.folders)