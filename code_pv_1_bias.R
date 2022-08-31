## Create a function called "DGM" to  generate standard data
DGM <- function(N,rho){
  p <- 10
  mu <- rep(0,p)
  cov_mat <- diag(1,p)
  cov_mat[2,1] <- rho
  cov_mat[1,2] <- rho
  x <- MASS::mvrnorm(n = N, mu = mu, Sigma = cov_mat)
  y0 <- rnorm(n = N, mean = (1*x[,1] + 2*x[,7] + 3*x[,10]))
  y <- rbinom(n = N, 1, prob = exp(y0)/(1+exp(y0)))
  mydata <- data.frame(cbind(x,y))
  colnames(mydata)[1:10] <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10")
  return(mydata)
}

## Create a function called "IMD" to introduce the missing data
ITM <- function(data,proportion){
  y <- data$y
  mydata_x <- data %>% dplyr::select(-y)
  freq_0 <- (1 - 2*proportion)/8
  myfreq <- c(proportion, proportion, freq_0, freq_0, freq_0, freq_0, freq_0, freq_0, freq_0, freq_0)
  mads_mydata <- ampute(mydata_x, mech = "MAR", prop = 0.7, freq = myfreq) 
  mydata_mis <- cbind(mads_mydata$amp,y)
  return(mydata_mis)
}

## Create a function called "DMI" to do the Multiple Imputation
DMI <- function(data,imps){
  imp_mydata_mis <- mice(data, m=imps, pri = FALSE)
  return(imp_mydata_mis)
}

#### Do the summarise of the results(What is bias of each variable) 
## Create a function to use "Rubin" rules to pull the result of coefficient (bias)
Rubin <- function(data){
  mydata_mi <- complete(data, action = "long", include = TRUE)
  imp <- as.mids(mydata_mi)
  analysis_models <- with(imp, glm(y ~ x1 + x2 + x3 + x4 + x5 
                                  + x6 + x7 + x8 + x9 + x10,family = "binomial"))
  a <- as.data.frame(summary(pool(analysis_models)))
  results <- a[2:11,] %>%
    select(c(term, estimate, p.value)) %>%
    filter(p.value <= 0.05) %>%
    select(c(term, estimate)) %>%
    tidyr::pivot_wider(names_from = "term",
                       values_from = "estimate")
  return(results)
}

## Create a function to do the whole process of simulation one time
simulation_singlerun_C <- function(N,C,P,M){
  mydata <- DGM(N = N,
                rho = C) ## Generate dataset
  mydata_mis <- ITM(data = mydata,
                    proportion = P) ## Introduce the missing data
  imp_mydata_mis <- DMI(data = mydata_mis,
                        imps = M)
  var_names <- names(mydata %>% select(starts_with("x")))
  results <- data.frame(matrix(NA, ncol = length(var_names), nrow = 1))
  names(results) <- var_names
  results <- results %>% 
    bind_rows(Rubin(imp_mydata_mis)) %>%
    slice(2) %>%
    replace(is.na(.), 0)
  return(results) 
}

## Do the whole process of simulation lots of times
library(MASS)
library(mice)
library(tidyverse)
library(foreach)
library(doParallel)
library(parallel)
## standard
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results0 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0,0,20)
  results
}
total_row <- dim(Results0)[1]
bias_x1 <- abs((sum(abs(Results0$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results0$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results0$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results0$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results0$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results0$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results0$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results0$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results0$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results0$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
                bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia0_x <- data.frame(t(data.frame(bia_x)))
## c=0 p=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results1 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0,0.1,20)
  results
}
bias_x1 <- abs((sum(abs(Results1$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results1$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results1$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results1$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results1$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results1$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results1$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results1$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results1$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results1$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia1_x <- data.frame(t(data.frame(bia_x)))
## c=0 p=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results2 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0,0.2,20)
  results
}
bias_x1 <- abs((sum(abs(Results2$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results2$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results2$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results2$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results2$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results2$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results2$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results2$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results2$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results2$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia2_x <- data.frame(t(data.frame(bia_x)))
## c=0 p=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results3 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0,0.3,20)
  results
}
bias_x1 <- abs((sum(abs(Results3$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results3$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results3$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results3$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results3$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results3$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results3$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results3$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results3$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results3$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia3_x <- data.frame(t(data.frame(bia_x)))
## c=0 p=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results4 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0,0.5,20)
  results
}
bias_x1 <- abs((sum(abs(Results4$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results4$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results4$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results4$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results4$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results4$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results4$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results4$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results4$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results4$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia4_x <- data.frame(t(data.frame(bia_x)))
###################################################################################
## c=0.2 p=0
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results5 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.2,0,20)
  results
}
bias_x1 <- abs((sum(abs(Results5$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results5$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results5$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results5$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results5$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results5$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results5$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results5$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results5$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results5$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia5_x <- data.frame(t(data.frame(bia_x)))
## c=0.2 p=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results6 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.2,0.1,20)
  results
}
bias_x1 <- abs((sum(abs(Results6$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results6$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results6$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results6$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results6$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results6$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results6$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results6$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results6$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results6$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia6_x <- data.frame(t(data.frame(bia_x)))
## c=0.2 p=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results7 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.2,0.2,20)
  results
}
bias_x1 <- abs((sum(abs(Results7$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results7$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results7$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results7$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results7$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results7$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results7$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results7$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results7$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results7$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia7_x <- data.frame(t(data.frame(bia_x)))
## c=0.2 p=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results8 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.2,0.3,20)
  results
}
bias_x1 <- abs((sum(abs(Results8$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results8$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results8$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results8$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results8$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results8$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results8$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results8$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results8$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results8$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia8_x <- data.frame(t(data.frame(bia_x)))
## c=0.2 p=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results9 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.2,0.5,20)
  results
}
bias_x1 <- abs((sum(abs(Results9$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results9$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results9$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results9$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results9$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results9$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results9$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results9$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results9$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results9$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia9_x <- data.frame(t(data.frame(bia_x)))
###############################################################################
## c=0.5 p=0
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results10 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.5,0,20)
  results
}
bias_x1 <- abs((sum(abs(Results10$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results10$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results10$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results10$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results10$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results10$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results10$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results10$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results10$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results10$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia10_x <- data.frame(t(data.frame(bia_x)))
## c=0.5 p=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results11 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.5,0.1,20)
  results
}
bias_x1 <- abs((sum(abs(Results11$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results11$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results11$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results11$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results11$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results11$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results11$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results11$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results11$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results11$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia11_x <- data.frame(t(data.frame(bia_x)))
## c=0.5 p=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results12 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.5,0.2,20)
  results
}
bias_x1 <- abs((sum(abs(Results12$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results12$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results12$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results12$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results12$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results12$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results12$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results12$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results12$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results12$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia12_x <- data.frame(t(data.frame(bia_x)))
## c=0.5 p=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results13 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.5,0.3,20)
  results
}
bias_x1 <- abs((sum(abs(Results13$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results13$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results13$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results13$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results13$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results13$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results13$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results13$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results13$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results13$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia13_x <- data.frame(t(data.frame(bia_x)))
## c=0.5 p=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results14 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.5,0.5,20)
  results
}
bias_x1 <- abs((sum(abs(Results14$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results14$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results14$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results14$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results14$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results14$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results14$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results14$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results14$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results14$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia14_x <- data.frame(t(data.frame(bia_x)))
###################################################################################################
## c=0.7 p=0
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results15 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.7,0,20)
  results
}
bias_x1 <- abs((sum(abs(Results15$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results15$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results15$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results15$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results15$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results15$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results15$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results15$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results15$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results15$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia15_x <- data.frame(t(data.frame(bia_x)))
## c=0.7 p=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results16 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.7,0.1,20)
  results
}
bias_x1 <- abs((sum(abs(Results16$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results16$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results16$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results16$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results16$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results16$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results16$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results16$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results16$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results16$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia16_x <- data.frame(t(data.frame(bia_x)))
## c=0.7 p=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results17 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.7,0.2,20)
  results
}
bias_x1 <- abs((sum(abs(Results17$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results17$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results17$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results17$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results17$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results17$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results17$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results17$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results17$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results17$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia17_x <- data.frame(t(data.frame(bia_x)))
## c=0.7 p=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results18 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.7,0.3,20)
  results
}
bias_x1 <- abs((sum(abs(Results18$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results18$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results18$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results18$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results18$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results18$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results18$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results18$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results18$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results18$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia18_x <- data.frame(t(data.frame(bia_x)))
## c=0.7 p=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results19 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.7,0.5,20)
  results
}
bias_x1 <- abs((sum(abs(Results19$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results19$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results19$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results19$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results19$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results19$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results19$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results19$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results19$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results19$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia19_x <- data.frame(t(data.frame(bia_x)))
##########################################################################
## c=0.95 p=0
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results20 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.95,0,20)
  results
}
bias_x1 <- abs((sum(abs(Results20$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results20$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results20$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results20$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results20$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results20$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results20$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results20$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results20$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results20$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia20_x <- data.frame(t(data.frame(bia_x)))
## c=0.95 p=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results21 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.95,0.1,20)
  results
}
bias_x1 <- abs((sum(abs(Results21$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results21$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results21$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results21$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results21$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results21$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results21$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results21$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results21$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results21$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia21_x <- data.frame(t(data.frame(bia_x)))
## c=0.95 p=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results22 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.95,0.2,20)
  results
}
bias_x1 <- abs((sum(abs(Results22$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results22$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results22$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results22$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results22$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results22$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results22$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results22$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results22$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results22$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia22_x <- data.frame(t(data.frame(bia_x)))
## c=0.95 p=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results23 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.95,0.3,20)
  results
}
bias_x1 <- abs((sum(abs(Results23$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results23$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results23$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results23$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results23$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results23$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results23$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results23$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results23$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results23$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia23_x <- data.frame(t(data.frame(bia_x)))
## c=0.95 p=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results24 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_C(10000,0.95,0.5,20)
  results
}
bias_x1 <- abs((sum(abs(Results24$x1))/total_row) - 1)
bias_x2 <- abs((sum(abs(Results24$x2))/total_row) - 0)
bias_x3 <- abs((sum(abs(Results24$x3))/total_row) - 0)
bias_x4 <- abs((sum(abs(Results24$x4))/total_row) - 0)
bias_x5 <- abs((sum(abs(Results24$x5))/total_row) - 0)
bias_x6 <- abs((sum(abs(Results24$x6))/total_row) - 0)
bias_x7 <- abs((sum(abs(Results24$x7))/total_row) - 2)
bias_x8 <- abs((sum(abs(Results24$x8))/total_row) - 0)
bias_x9 <- abs((sum(abs(Results24$x9))/total_row) - 0)
bias_x10 <- abs((sum(abs(Results24$x10))/total_row) - 3)
bia_x <- c(bias_x1,bias_x2,bias_x3,bias_x4,bias_x5,
           bias_x6,bias_x7,bias_x8,bias_x9,bias_x10)
bia24_x <- data.frame(t(data.frame(bia_x)))