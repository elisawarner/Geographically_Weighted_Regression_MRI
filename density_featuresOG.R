#############################################################################################
### DENSITY FEATURES CODE
### Description : function which inputs a training list (x.trn), a test list (x.tst),
### percent variance you want covered by the principal components (perc.var), and 
### the size of the density curves (1 x m)
### 
### Author: Shariq Mohammed, maintained by Elisa Warner (elisawa@umich.edu)
### Date: 3/3/2023
#############################################################################################

density_features = function(x.trn, # list of length n with x-values for n training cases
                            x.tst = NULL, # list of length n.new with x-values for n.new testing cases
                            perc.var = 99, # % variance explained by PCs; a number between 0 to 100
                            m = 512 # number of values to evaluate PDF on x-axis
                            ){
  perc.var = perc.var/100
  stopifnot(perc.var > 0)
  if(perc.var>1) perc.var = 1
  
  n = length(x.trn)
  
  # construct densities for training data
  min_coef = min(unlist(x.trn), na.rm=T)
  max_coef = max(unlist(x.trn), na.rm=T)
  den = matrix(nrow = n, ncol = m)
  for(i in 1:n){
    temp_vec = c(x.trn[[i]][!is.na(x.trn[[i]])])
    den[i,] = density(temp_vec, from = min_coef, to = max_coef, 
                      n = m, bw='nrd')$y
    den[i,] = den[i,]/(max_coef - min_coef)
    den[i,] = den[i,]/sum(diff(seq(0,1,length = m))*(den[i,-m]+ den[i,-1])/2)
  }
  
  # replace testing data extreme values with training data extremes
  x.tst_temp = lapply(1:length(x.tst),
                 function(l){
                   x = x.tst[[l]]
                   x[x>max_coef] = max_coef
                   x[x<min_coef] = min_coef
                   x
                 })
  x.tst = x.tst_temp
  rm(x.tst_temp)
  
  # square-root transform of the densities for training data
  sqrt_den = sqrt(den)
  # compute Karcher mean on the sqrt space for training data
  eps = .5
  nmv = 100
  kmean_sqrt = rep(1, m)
  
  iter = 1
  while(nmv[iter] > 1e-10){
    vbar = rep(0, m)
    v = matrix(nrow = n, ncol = m)
    for(i in 1:n){
      if(abs(sum(kmean_sqrt-(sqrt_den[i,])))<1e-10){
        v[i,] = rep(0,m)
      } else{
        th = acos(sum(diff(seq(0,1,length=m))*
                        ((kmean_sqrt*(sqrt_den[i,]))[-m]+
                           (kmean_sqrt*(sqrt_den[i,]))[-1])/2))
        v[i,] = (th/sin(th))*((sqrt_den[i,])-kmean_sqrt*cos(th))
      }
      
      vbar = vbar + v[i,]/n
    } 
    
    nv = sqrt(sum((((eps*vbar)*(eps*vbar))[-m]+
                     ((eps*vbar)*(eps*vbar))[-1])/2)/(m-1))
    kmean_sqrt = cos(nv)*kmean_sqrt + sin(nv)*(eps*vbar)/nv
    iter = iter+1
    nmv = append(nmv, sqrt(sum(((vbar*vbar)[-m]+(vbar*vbar)[-1])/2)/(m-1)))
  }
  
  # project sqrt transforms onto the tangent space of Karcher mean 
  ieden = matrix(0, nrow = n, ncol = m)
  for(i in 1:n){
    if(abs(sum(kmean_sqrt-(sqrt_den[i,])))<1e-10){
      ieden[i,] = rep(0,m)
    } else{
      th = acos(sum(diff(seq(0,1,length=m))*
                      ((kmean_sqrt*(sqrt_den[i,]))[-m]+ 
                         (kmean_sqrt*(sqrt_den[i,]))[-1])/2))
      ieden[i,] = (th/sin(th))*((sqrt_den[i,])-kmean_sqrt*cos(th))
    }
  }
  
  # compute the sample covariance matrix
  K = matrix(0, nrow = m, ncol = m)
  for(i in 1:n) K = K+tcrossprod(ieden[i,])
  K = K/(n-1)
  
  # Compute PC Scores for training data
  svdecomp = svd(K)
  p = sum(cumsum((svdecomp$d^2)/sum(svdecomp$d^2))<perc.var)+1
  U = svdecomp$u[,1:p]
  pScore = ieden%*%U
  
  # standardization for training data
  cl.mean = colMeans(pScore)
  cl.sd = apply(pScore, 2, sd)
  
  pScore = t(t(pScore)-cl.mean)
  pScore = t(t(pScore)/cl.sd)
  
  if(is.null(x.tst)){
    pScore.new = NULL
  }else{
    n.new = length(x.tst)
    
    # construct densities for testing data
    den_new = matrix(nrow = n.new, ncol = m)
    for(i in 1:n.new){
      temp_vec = c(x.tst[[i]][!is.na(x.tst[[i]])])
      den_new[i,] = density(temp_vec, from = min_coef, to = max_coef,
                            n = m, bw='nrd')$y
      den_new[i,] = den_new[i,]/(max_coef - min_coef)
      den_new[i,] = den_new[i,]/sum(diff(seq(0,1,length = m))*
                                      (den_new[i,-m]+ den_new[i,-1])/2)
    }
    
    # square-root transform of the densities for testing data
    sqrt_den_new = sqrt(den_new)
    if(length(den_new) == m){
      n.new = 1
      sqrt_den_new = as.matrix(sqrt_den_new, nrow = n.new)
    }
    
    # Compute PC Scores for testing data
    ieden_new = matrix(0, nrow = n.new, ncol = m)
    for(i in 1:n.new){
      if(abs(sum(kmean_sqrt-(sqrt_den_new[i,])))<1e-10){
        ieden_new[i,] = rep(0,m)
      } else{
        th = acos(sum(diff(seq(0,1,length=m))*
                        ((kmean_sqrt*(sqrt_den_new[i,]))[-m]+ 
                           (kmean_sqrt*(sqrt_den_new[i,]))[-1])/2))
        ieden_new[i,] = (th/sin(th))*((sqrt_den_new[i,])-kmean_sqrt*cos(th))
      }
    }
    pScore.new = ieden_new%*%U
    
    # standardization for testing data
    pScore.new = t(t(pScore.new)-cl.mean)
    pScore.new = t(t(pScore.new)/cl.sd)
  }
  
  list(pScore = pScore,
       pScore.new = pScore.new)
}