#!usr/bin/env Rscript
# Supporting function for classic_aa_test.m
#
# Copyright Sohan Seth and Manuel Eugster

set.seed(1)

library("archetypes")
library(R.matlab)
load("bestModel.RData")

predict.archetypes <- 
function (object, newdata, ...)
{
   stopifnot(object$family$which == "original")
   scale <- object$scaling
   x <- t(newdata)
   x <- x - scale$mean
   x <- x/scale$sd
   x <- object$family$dummyfn(x, ...)
   zs <- t(parameters(object))
   zs <- zs - scale$mean
   zs <- zs/scale$sd
   zs <- rbind(zs, 200)
   alphas <- matrix(NA, ncol = ncol(x), nrow = ncol(coef(object)))
   alphas <- object$family$alphasfn(alphas, zs, x)
   t(alphas)
}

inp <- readMat('classic_aa_test_Input.mat')

op <- predict.archetypes(a, t(inp$matFeatSam))
writeMat('classic_aa_test_Output.mat', matLatSam = t(op))