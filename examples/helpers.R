
### Required R packages: #############################################

library("archetypes")
library("colorspace")
library("classInt")
library("xtable")
library("R.matlab")
library("TSP")



### Deviances: #######################################################

gaussian_deviance <- function(object, data) {
  t <- t(object$Z %*% object$H)
  sqrt(rowSums((t - data)^2))  
}


poission_deviance <- function(object, data) {
  t <- t(object$Z %*% object$H)
  rowSums(2 * (data * log((data +.Machine$double.eps)/t) - data + t))
}


bernoulli_deviance <- function(object, data) {
  t <- t(object$Z %*% object$H)
  t[t > 1] <- 1.0
  rowSums(2 * (data * log((data +.Machine$double.eps)/(t+.Machine$double.eps)) + 
                 (1-data) * log(((1-data) +.Machine$double.eps)/(1-t+.Machine$double.eps))))
}



### Simplex plot wrapper for PAA: ####################################

simplexplot_paa <- function(paa, ...) {
  a <- as.archetypes(t(paa$Z), k = ncol(paa$Z), alphas = t(paa$H), rss = NA)
  archetypes::simplexplot(a, ...)
}



### Loader for pre-computed PAA: #####################################

load_paa <- function(matfile, dat) {
  a <- readMat(matfile)

  r <- list()
  r$W <- a$matSamLat
  r$H <- a$matLatSam
  r$Z <- t(dat) %*% r$W
  r$obj <- a$obj
  
  r
}



### Greedy match: ####################################################

greedy_match <- function(aa, paa) {
  stopifnot(require("flexclust"))
  
  d <- distEuclidean(aa, paa)
  
  match <- list()
  for ( i in 1:ncol(d) ) {
    match[[i]] <- arrayInd(which.min(d), dim(d))
    d[match[[i]][1], ] <- Inf
    d[, match[[i]][2]] <- Inf
  }
  
  match <- do.call(rbind, match)
  colnames(match) <- c("AA", "PAA")
  match <- match[order(match[, "PAA"]), ]

  match
}



### Percentiles: #####################################################

perc <- function(aa, dat) {
  l <- lapply(1:ncol(dat),
              function(i) {
                aa[i, ] / max(dat[, i])})
  do.call(rbind, l)
}



