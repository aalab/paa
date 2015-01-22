
## ----echo=FALSE, message=FALSE, warnings=FALSE---------------------------
knitr::opts_chunk$set(echo = FALSE, messages = FALSE, 
                      fig.height = 6, fig.width = 6)


## ----results='hide', message=FALSE, warning=FALSE, error=FALSE-----------
library("SportsAnalytics")
source("helpers.R")


## ----results='asis'------------------------------------------------------
cat("* Load soccer data from package `SportsAnalytics`\n")

data("EURO4PlayerSkillsSep11")

dat <- subset(EURO4PlayerSkillsSep11,
              Position != "Goalkeeper",
              select = -c(Birthday, Positions))
dat <- subset(dat,  Attack != 6 & TopSpeed > 0)

mat <- as.matrix(subset(dat, select = -c(Team, Name, Number, Nationality,
                                         Age, InjuryTolerance, Foot, Side,
                                         Position, League, KeeperSkills,
                                         Height, Weight, ConditionFitness,
                                         WeakFootAccuracy, WeakFootFrequency)))
rownames(mat) <- NULL


cat("* Load precomputed PAA solution from file `Soccer_paa.mat`\n")
paa <- load_paa("data_Soccer/Soccer_paa.mat", mat)


## ----fig.height=4--------------------------------------------------------
op <- par(mar = c(4, 4, 1, 1))
plot(paa$obj, type = "l", xlab = "Iterations", ylab = "Cost")
par(op)


## ------------------------------------------------------------------------
a <- as.archetypes(t(paa$Z), k = ncol(paa$Z), alpha = paa$H, rss = NA)

## Change order to match Eugster(2012):
a$archetypes <- a$archetypes[c(3, 4, 1, 2), ]
a$alphas <- a$alphas[c(3, 4, 1, 2), ]

barplot(a, mat, percentiles = TRUE, below.compressed.srt = 90, 
        below.compressed.height = 1.1, border = "white")


