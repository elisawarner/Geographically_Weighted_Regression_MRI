###################################################################################
### RUN Residual Densities (RUN_RESD.R)
### Author and Maintainer: Elisa Warner (elisawa@umich.edu)
### Last Update : 3/3/2023
###
### Description: This code takes in the Slice numbers obtained from runSliceGWR.R
###  as well as .spdf files from the previous program and runs Geographically Weighted
###  Regression (GWR) to obtain predictions from MRI modality X -> MRI modality Y.
###  In other words, we fit a model s.t. ideally (Y = pred(X)). Then, residuals are
###  calculated pred(X) - Y and flattened.
###################################################################################

setwd("Z:/.../1Pseudoprogression/UPenn-UMich Code Sharing")
folder.name = getwd()

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


  folder = file.path(folder.name,"Data_Full42_2023",paste0(x_pred, "_", y_pred))
  dir.create(folder)


  packages = c("WhiteStripe", "oro.nifti", "doParallel", "stringr", "R.utils", "matlabr",
             "R.matlab", "sp", "GWmodel", "RNifti")
  package.check = lapply(packages,
                       function(x){
                         if(!x %in% rownames(installed.packages())) {
                           install.packages(x, dependencies = TRUE)
                         }
                         library(x, character.only = TRUE)
                       })
  rm(package.check,packages)

  subfolders = list.files(file.path(folder.name, 'MRI_full'))

  ############################## check for spdf
  # delete the subfolders which have no valid spdf
  temp = c()
  for (i in 1:length(subfolders)){
    save_fname = file.path(folder.name, 'MRI_full', subfolders[i])
    filename = paste0(save_fname,'/spdf','_',x_pred,'_',y_pred,'.RData')
    if (file.exists(filename)){
      temp = c(temp, subfolders[i])
    }
  }
  subfolders = temp
  rm(filename, temp)

  ############################### GET RESD
  nCores = 6
  cl = makeCluster(nCores)
  registerDoParallel(cl, cores = nCores)

  # run GWR with parallel processing
  resd = foreach(s = subfolders, .packages = 'GWmodel') %dopar% {
    #for(s in subfolders) {
      save_fname = file.path(folder.name, 'MRI_full', s)
      f.path = file.path(folder.name, 'MRI_full', s)
      load(paste0(save_fname,'/spdf','_',x_pred,'_',y_pred,'.RData')) # open spatial data
  
      #file.remove(paste0(f.path,'/spdf.RData'))
      dist.gw = gw.dist(spdf@coords)
      gwr.bw = bw.gwr(x_pred_img ~ y_pred_img, data = spdf, dMat = dist.gw, kernel = "gaussian",
                  adaptive = TRUE, approach = 'AICc')
      gwr_res = gwr.basic(x_pred_img ~ y_pred_img, data = spdf, dMat = dist.gw,
                      kernel = "gaussian", adaptive = TRUE, bw = gwr.bw)
      r = matrix(nrow = dims[1], ncol = dims[2])
      for(k in 1:length(spdf)){
        r[spdf@coords[k,1],spdf@coords[k,2]] = gwr_res$SDF$residual[k]
      }
      rm(f.path, spdf, dims, dist.gw, gwr.bw, gwr_res, k)
      r
    }
  stopCluster(cl)
  names(resd) = subfolders
  rm(s, cl, nCores)

  ########## density
  den_npoints = 1000

  n = length(resd)
  den = matrix(nrow = n, ncol = den_npoints)
  mn = -10
  mx = 10
  for(i in 1:n){
    temp_vec = c(resd[[i]][!is.na(resd[[i]])])
    temp_vec[temp_vec<mn] = mn
    temp_vec[temp_vec>mx] = mx
    den[i,] = density(temp_vec, from = mn, to = mx,
                    n = den_npoints, bw='nrd')$y
    den[i,] = den[i,]/(mx - mn)
    den[i,] = den[i,]/sum(diff(seq(0,1,length = den_npoints))*
                        (den[i,-den_npoints]+ den[i,-1])/2)
  }
  rownames(den) = subfolders
  #rm(temp_vec, n, den_npoints, mn, mx, i, folder.name)
  save(resd, den, subfolders, file = paste0(folder, '/resd_', x_pred, '_', y_pred, '.RData'))
}
