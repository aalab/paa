#!usr/bin/env Rscript
# Supporting function for classic_AA.m
#
# Copyright Sohan Seth and Manuel Eugster

# set.seed(1)

library("archetypes")
library(R.matlab)
inp <- readMat('classic_aa_Input.mat')

# matFeatSam available in the workspace (from Matlab)
# data("toy")
# matFeatSam <- t(toy)

as <- stepArchetypes(t(inp$matFeatSam), k = inp$nLat, nrep = 1)
a <- bestModel(as)

s <- list(matSamLat = a$betas,  # W
          matLatSam = a$alpha,  # H
          obj = 
            sapply(seq(length(a$history$states())-1), 
                   function(i) a$history$get(i)$archetypes$rss))

writeMat('classic_aa_Output.mat',matSamLat = t(s$matSamLat), matLatSam = t(s$matLatSam), obj = s$obj)
save(a, file = "bestModel.RData")