## Create a function called "DGM" to  generate standard data
DGM <- function(N,rho){
  p <- 10
  mu <- rep(0,p)
  cov_mat <- diag(1,p)
  cov_mat[2,1] <- rho
  cov_mat[1,2] <- rho
  x <- MASS::mvrnorm(n = N, mu = mu, Sigma = cov_mat)
  y0 <- rnorm(n = N, mean = (0.1*x[,1] + 2*x[,7] + 3*x[,10]))
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
## create a function to do the "Majority Method" selection
DMM <- function(data){
  scope0 <- list(upper = ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, 
                 lower = ~ 1)
  expr <- expression(f1 <- glm(y ~ x1 + x2 + x3 + x4 
                               + x5 + x6 + x7 + x8 + x9 + x10,family = "binomial"),
                     f2 <- step(f1, scope = scope0, trace = FALSE))
  
  fit <- with(data, expr)
  ## Apply stepwise on each of the imputed dataset separately
  formulas <- lapply(fit$analyses,formula)
  ## Return the selection result
  terms <- lapply(formulas, terms)
  votes <- unlist(lapply(terms, labels))
  ## Check the terms on each models
  return((as.data.frame(table(votes))) %>%
           tidyr::pivot_wider(names_from = "votes",
                              values_from = "Freq"))
}
#### Do the summarise of the results(Which variable are selected) 
## Create a function to do the whole process of simulation one time
simulation_singlerun_W <- function(N,C,P,M){
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
    bind_rows(DMM(imp_mydata_mis)) %>%
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
  results <- simulation_singlerun_W(10000,0,0,20)
  results
}
Results0 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0 P=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results1 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0,0.1,20)
  results
}
Results1 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0 P=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results2 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0,0.2,20)
  results
}
Results2 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0 P=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results3 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0,0.3,20)
  results
}
Results3 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0 P=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results4 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0,0.5,20)
  results
}
Results4 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
######################################################################################################
## c=0.2 p=0
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results5 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.2,0,20)
  results
}
Results5 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.2 P=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results6 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.2,0.1,20)
  results
}
Results6 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.2 P=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results7 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.2,0.2,20)
  results
}
Results7 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.2 P=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results8 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.2,0.3,20)
  results
}
Results8 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.2 P=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results9 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.2,0.5,20)
  results
}
Results9 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
######################################################################################################
## c=0.5 p=0
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results10 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.5,0,20)
  results
}
Results10 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.5 P=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results11 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.5,0.1,20)
  results
}
Results11 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.5 P=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results12 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.5,0.2,20)
  results
}
Results12 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.5 P=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results13 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.5,0.3,20)
  results
}
Results13 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.5 P=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results14 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.5,0.5,20)
  results
}
Results14 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
######################################################################################################
## c=0.7 p=0
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results15 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.7,0,20)
  results
}
Results15 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.7 P=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results16 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.7,0.1,20)
  results
}
Results16 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.7 P=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results17 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.7,0.2,20)
  results
}
Results17 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.7 P=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results18 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.7,0.3,20)
  results
}
Results18 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.7 P=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results19 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.7,0.5,20)
  results
}
Results19 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
######################################################################################################
## c=0.95 p=0
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results20 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.95,0,20)
  results
}
Results20 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.95 P=0.1
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results21 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.95,0.1,20)
  results
}
Results21 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.95 P=0.2
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results22 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.95,0.2,20)
  results
}
Results22 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.95 P=0.3
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results23 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.95,0.3,20)
  results
}
Results23 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))
## C=0.95 P=0.5
a <- detectCores()
cl <- makeCluster(a[1]-1)
registerDoParallel(cl)
Results24 = foreach(i = 1:200, .combine=rbind, .packages=c("MASS", "mice","tidyverse")) %dopar% {
  results <- simulation_singlerun_W(10000,0.95,0.5,20)
  results
}
Results24 %>%
  summarise_all(.funs = list("min" = min,
                             "mean" = mean,
                             "max" = max))






















