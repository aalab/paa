
## ----echo=FALSE, message=FALSE, warnings=FALSE---------------------------
opts_chunk$set(echo = FALSE, messages = FALSE, 
               fig.height = 6, fig.width = 6)


## ----results='hide', message=FALSE, warning=FALSE, error=FALSE-----------
source("helpers.R")


## ----results='asis', message=FALSE, warning=FALSE------------------------
cat("* Read tourist data from file `GSAW97_bin.csv`\n")
GSAW97bin <- read.csv2("data_GSAW97/GSAW97_bin.csv")

cat("* Load precomputed PAA solution from file `GSAW97_paa.mat`\n")
paa <- load_paa("data_GSAW97/GSAW97_paa.mat", GSAW97bin)

cat("* Load precomputed AA solution from file `GSAW97_aa.Rds`\n")
aa <- readRDS("data_GSAW97/GSAW97_aa.Rds")

cat("* Match AA solution with PAA solution\n")
match <- greedy_match(aa$archetypes, t(paa$Z))
match <- match[c(1, 3, 5, 4, 6, 2), ]


## ------------------------------------------------------------------------
thepaa <- paa$Z
colnames(thepaa) <- sprintf("PAA%s", 1:6)
round(thepaa[, match[, "PAA"]], 1)


## ------------------------------------------------------------------------
theaa <- t(aa$archetypes)
colnames(theaa) <- sprintf(" AA%s", 1:6)
round(theaa[, match[, "AA"]], 1)

