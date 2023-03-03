# Shariq Mohammed <shariqm@bu.edu>
# Last edited: Feb 7, 2022
# 
# Input
# sp_df: Should be a data frame with four columns with column names as
#        c("x","y",'im.resp','im.pred'). The first two columns are the 
#        x-y coordinates of the tumor pixels. The third and fourth columns 
#        are the corresponding values for the response and predictor images.
#        Number of rows should be same as number of pixels in the largest 
#        connected component of the tumor.
#
# Output
# sp_df: Will be the same data frame as the input but with two additional 
#        columns. That is, the fifth and sixth columns are named as 
#        "coef" and "resid", which correspond to the GWR model coefficients
#        and the residuals.
gwr <- function(sp_df){
  # Check for packages sp and GWmodel
  # install if they are not installed
  if(!"sp" %in% rownames(installed.packages())){
    install.packages("sp")
    library(sp)
  }
  if(!"GWmodel" %in% rownames(installed.packages())){
    install.packages("GWmodel")
    library(GWmodel)
  }
  spdf <- SpatialPointsDataFrame(sp_df[,c("x","y")],
                                 sp_df[,c("im.resp","im.pred")])
  # compute distance matrix
  dist.gw <- gw.dist(spdf@coords)
  # compute adaptive bandwidth
  gwr.bw <- bw.gwr(im.resp ~ im.pred, data = spdf, dMat = dist.gw,
                  kernel = "gaussian", adaptive = TRUE, approach = 'AICc')
  # fit GWR model
  gwr_res <- gwr.basic(im.resp ~ im.pred, data = spdf, dMat = dist.gw,
                      kernel = "gaussian", adaptive = TRUE, bw = gwr.bw)

  # extract GWR model coefficients and residuals and append to sp_df
  sp_df$coef <- gwr_res$SDF$im.pred
  sp_df$resid <- gwr_res$SDF$residual
  
  # return sp_df
  sp_df
}