# This script is for performing mantel regression (Distance matrix regression)
# 
# written by Suhwan Gim 
# 2020. Nov. 13
# References 
# 1. "Multiple regression on distance matrices: a multivariate spatial analysis tool (2008)"
#                            https://link.springer.com/article/10.1007/s11258-006-9126-3
# 2. "Package 'ecodist' https://cran.r-project.org/web/packages/ecodist/ecodist.pdf

# ==================================    
# Load packages
library(ecodist)                        # For distance matrix regression 
library(ggplot2)                        # For visualization 
library(ggcorrplot)
library(psych)  
library(R.matlab)                       # library to read matlab data formats into R
library(corrr) # network 
library(RColorBrewer)
library(doParallel) # for Parallel processing
library(foreach)    # for Parallel processing
library(itertools) 
library(doSNOW)
# ==================================


## Load trajectory
temp_ISRSA <- readMat("/Users/suhwan/Dropbox/Projects/SEMIC2/data/ISRSA_structures/trajectory_N58.mat");
traj_raw<-temp_ISRSA$traj

cor_trj<-lower(cor(traj_raw))
ggcorrplot(full(cor_trj),method="square",outline.col = "white",type = "lower", 
           ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"))

## Load whole-brain ISRSA 
cue_vox<-read.csv('/Users/suhwan/Dropbox/Projects/SEMIC2/data/ISRSA_structures/cue_voxel.csv',header = FALSE)
nc_vox<-read.csv('/Users/suhwan/Dropbox/Projects/SEMIC2/data/ISRSA_structures/nc_voxel.csv',header = FALSE)

## Load individual differnces measurement (Lower triangle)
IND_matrix<-read.csv('/Users/suhwan/Dropbox/Projects/SEMIC2/data/IND_matrix.csv',header = FALSE)
Col_matrix<-read.csv('/Users/suhwan/Dropbox/Projects/SEMIC2/data/COL_matrix.csv',header = FALSE)
IRI_matrix<-read.csv('/Users/suhwan/Dropbox/Projects/SEMIC2/data/IRI_matrix.csv',header = FALSE)
## load first_level beta trajectory
first_levelN59<-readMat('/Users/suhwan/Dropbox/Projects/SEMIC2/data/ISRSA_structures/firstlevel_beta_N59.mat')
cue_col <- lower(cor(first_levelN59$firstlevel.cue[,c(1:16,18:59)]))
stim_col <- lower(cor(first_levelN59$firstlevel.stim[,c(1:16,18:59)]))

## Do Parallel processing 
# 1.set function 
multiple_MRM<-function(voxel_i){
  #res<-MRM_one(cue_vox[,voxel_i] ~ IND_matrix$V1 + Col_matrix$V1 + IRI_matrix$V1 + cor_trj + cue_col +  stim_col, nperm=10000, mrank=TRUE)
  res<-MRM_one(cue_vox[,voxel_i] ~ IND_matrix$V1 + Col_matrix$V1 + cor_trj, nperm=10000, mrank=TRUE)
  
  
  new_res = data.frame(interr = res$coef[1,1], interpp1 = res$coef[1,2], interpp2= res$coef[1,3], interpp3 = res$coef[1,4],
                       indr=res$coef[2,1],indp1=res$coef[2,2],indp2=res$coef[2,3],indp3=res$coef[2,4],
                       colr=res$coef[3,1],colp1=res$coef[3,2],colp2=res$coef[3,3],colp3=res$coef[3,4],
                       trjr=res$coef[4,1],trjp1=res$coef[4,2],trjp2=res$coef[4,3],trjp3=res$coef[4,4],
                       idx = voxel_i,
                       r2=res$r.squared[1],r2_p = res$r.squared[2],Fscore = res$F.test[1])
  if (voxel_i %% 1000 == 0) {
    print(voxel_i)
  }
  return(new_res)
}

multiple_MRM2<-function(voxel_i){
  res<-MRM_one(nc_vox[,voxel_i] ~ IND_matrix$V1 + Col_matrix$V1 + IRI_matrix$V1 + cor_trj + cue_col +  stim_col, nperm=10000, mrank=TRUE)
  
  
  
  new_res = data.frame(interr = res$coef[1,1], interpp1 = res$coef[1,2], interpp2= res$coef[1,3], interpp3 = res$coef[1,4],
                       indr=res$coef[2,1],indp1=res$coef[2,2],indp2=res$coef[2,3],indp3=res$coef[2,4],
                       colr=res$coef[3,1],colp1=res$coef[3,2],colp2=res$coef[3,3],colp3=res$coef[3,4],
                       trjr=res$coef[4,1],trjp1=res$coef[4,2],trjp2=res$coef[4,3],trjp3=res$coef[4,4],
                       idx = voxel_i,
                       r2=res$r.squared[1],r2_p = res$r.squared[2],Fscore = res$F.test[1])
  if (voxel_i %% 1000 == 0) {
    print(voxel_i)
  }
  return(new_res)
}

set.seed(150421)

# 3. run Parallel processing 
numCores <- 30 #detectCores() -1
myCluster <- makeCluster(numCores,outfile="")
registerDoParallel(myCluster)
record<-numeric(0)
parallel::clusterExport(myCluster,varlist = c("cue_vox","nc_vox", "IND_matrix", "Col_matrix","record","IRI_matrix","cor_trj","cue_col","stim_col"))


system.time({
  total_result_cue<-foreach(n = 1:199209, .combine = rbind, .packages = "ecodist", .noexport = c('temp_ISRSA') ) %dopar% {
    #total_res<-multiple_MRM(n)}
    total_res <- multiple_MRM(n)}
  write.csv(total_result_cue,file = "/Users/suhwan/Dropbox/Projects/SEMIC2/data/ISRSA_results_resting/heat_cue_results_regrs_reduced.csv")
})

system.time({
  total_result_nc<-foreach(n = 1:199209, .combine = rbind, .packages = "ecodist", .noexport = c('temp_ISRSA') ) %dopar% {
    #total_res<-multiple_MRM(n)}
    total_res <- multiple_MRM2(n)}
  write.csv(total_result_nc,file = "/Users/suhwan/Dropbox/Projects/SEMIC2/data/ISRSA_results_resting/heat_nc_results_regrs_reduced.csv")
})


stopCluster(myCluster)



# Define function
#####
MRM_one <- function(formula = formula(data), data, nperm = 1000, method="linear", mrank = FALSE)
{
  # MRM: Multiple Regression on distance Matrices
  # Sarah Goslee 2008-07-18
  # tests R^2 and regression coefficients using a
  # permutation test
  
  
  # added method argument
  # linear is the default, standard option
  # logistic adds capability to do logistic regression
  # - this will be much slower
  # mrank is ignored if method is not "linear"
  
  # Stuff R needs to be able to use a formula
  m <- match.call(expand.dots = FALSE)
  m2 <- match(c("formula", "data"), names(m), nomatch=0)
  m <- m[c(1, m2)]
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  m <- as.matrix(m)
  
  
  # End of R stuff. m is now the data for the MRM test as
  # columns y, x, n1, n2, n3, ...
  # Determine the size of the matrices & do some error checking.
  n <- (1 + sqrt(1 + 8 * nrow(m)))/2
  if(abs(n - round(n)) > 0.0000001)
    stop("Matrix not square.\n")
  n <- round(n)
  if(ncol(m) < 2) stop("Not enough data. \n")
  
  
  if(method == "linear") {
    if(mrank) {
      m <- apply(m, 2, rank)
    }
    
    # convert matrices to column order to ensure compatibility with C
    for(thiscol in seq_len(ncol(m))) {
      tempmat <- full(m[,thiscol])
      m[,thiscol] <- tempmat[col(tempmat) > row(tempmat)]
    }
    
    # use matrix form to speed up calculations
    X <- m[ ,2:ncol(m), drop=FALSE]
    X <- cbind(rep(1, nrow(X)), X)
    Y <- m[ ,1, drop=FALSE]
    
    nd <- nrow(X)
    
    # only need to calculate (X'X)^-1 once
    XX <- crossprod(X)
    XX <- solve(XX)
    
    # will need to calculate Xy for each permutation
    XY <- crossprod(X, Y)
    YY <- crossprod(Y)
    
    # regression coefficients
    b <- XX %*% XY
    rownames(b) <- c("Int", colnames(X)[2:ncol(X)])
    
    bXY <- crossprod(b, XY)
    SSE <- YY - bXY
    
    SSTO <- YY - sum(Y)^2/nd
    SSR = SSTO - SSE
    
    # R2 = 1 - SSE/SSTO
    R2 <- 1 - SSE/SSTO
    R2 <- as.vector(R2)
    
    # F* = MSR / MSE
    # MSR = SSR / (p - 1) 
    # MSE = SSE / (n - p)
    p <- ncol(X) # number of parameters estimated
    F <- (SSR / (p - 1)) / (SSE / (nd - p))
    
    R2.pval <- NA
    b.pval <- rep(NA, ncol(X))
    F.pval <- NA
    
    if(nperm > 0) {
      
      R2.all <- numeric(nperm)
      
      # for regression coefficients, use pseudo-t of Legendre et al. 1994
      b.all <- numeric(nperm*p)
      
      # try out an overall F-test for lack of fit
      F.all <- numeric(nperm)
      
      cresults <- .C("mrmperm", as.double(as.vector(X)), as.double(as.vector(Y)), 
                     as.integer(p), as.integer(nd), as.integer(n), as.integer(nperm), 
                     R2.all = as.double(R2.all), b.all = as.double(b.all), 
                     F.all = as.double(F.all), as.double(numeric(n*n)), 
                     as.integer(numeric(n)), as.double(as.vector(XX)), as.double(numeric(p)),
                     as.double(0), as.double(numeric(p)), PACKAGE = "ecodist")
      
      R2.all <- cresults$R2.all
      R2.pval <- length(R2.all[R2.all >= R2.all[1]])/nperm
      
      F.all <- cresults$F.all
      F.pval <- length(F.all[F.all >= F.all[1]])/nperm
      
      # b.all contains pseudo-t of Legendre et al. 1994
      b.all <- matrix(cresults$b.all, nrow=nperm, ncol=p, byrow=TRUE)
      #b.pval <- apply(b.all, 2, function(x)length(x[abs(x) >= abs(x[1])])/nperm)
      b.pval1 <- apply(b.all, 2, function(x)length(x[x >= x[1]])/nperm)
      b.pval2 <- apply(b.all, 2, function(x)length(x[x <= x[1]])/nperm)
      b.pval3 <- apply(b.all, 2, function(x)length(x[abs(x) >= abs(x[1])])/nperm)
    }
    results <- list(coef=cbind(b, pval1=b.pval1, pval2=b.pval2, pval3=b.pval3), r.squared=c(R2=R2, pval = R2.pval),F.test=c(F=F, F.pval = F.pval))
  } else {
    if(method == "logistic") {
      # extract data from formula object
      X <- m[ ,2:ncol(m), drop=FALSE]
      Y <- m[ ,1, drop=FALSE]
      
      colnames(Y) <- "Y"
      newdata <- data.frame(Y=Y, X)
      fit1 <- glm(Y ~ ., data=newdata, family=binomial(link = "logit"))
      
      # want to save coefficients, deviance & df
      b <- coefficients(fit1)
      dev <- summary(fit1)$deviance
      dev.df <- summary(fit1)$df.residual
      
      b.pval <- NA
      dev.pval <- NA
      
      if(nperm > 0) {
        b.all <- matrix(NA, nrow=nperm, ncol=length(b))
        b.all[1,] <- b
        dev.all <- rep(NA, nperm)
        dev.all[1] <- dev
        
        for(i in 2:nperm) {
          newSample <- sample(n)
          newY <- full(Y)
          newY <- newY[newSample, newSample]
          newY <- lower(newY)
          newdata <- data.frame(Y=newY, X=X)
          newfit <- glm(Y ~ ., data=newdata, family=binomial(link = "logit"))
          b.all[i,] <- coefficients(newfit)
          dev.all[i] <- summary(newfit)$deviance
        }
        b.pval <- apply(b.all, 2, function(x) length(x[abs(x) >= abs(x[1])])/nperm)
        dev.pval <- length(dev.all[dev.all >= dev.all[1]])/nperm
      }
      results <- list(coef = cbind(b, pval = b.pval), dev = c(resid.dev = dev, resid.df = dev.df, dev.pval = dev.pval))
    } else {
      stop("method must be 'linear' or 'logistic'\n")
    }
  }
  results
}
