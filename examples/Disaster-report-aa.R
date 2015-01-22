
## ----echo=FALSE, message=FALSE, warnings=FALSE---------------------------
opts_chunk$set(echo = FALSE, messages = FALSE, 
               fig.height = 6, fig.width = 6)


## ----results='hide', message=FALSE, warning=FALSE, error=FALSE-----------
source("helpers.R")


## ----results='asis', message=FALSE, warning=FALSE------------------------
cat("* Read disaster data from file `emdata.csv`\n")
emdata <- read.table("data_Disaster/emdata.csv", stringsAsFactors = FALSE, sep = ";", header = TRUE)
emdata_mat <- as.matrix(emdata[, -1])

cat("* Load precomputed PAA solution from file `Disaster_paa.mat`\n")
paa <- load_paa("data_Disaster/Disaster_paa.mat", emdata_mat)

cat("* Load precomputed AA solution from file `Disaster_aa.Rds`\n")
aa <- readRDS("data_Disaster/Disaster_aa.Rds")

cat("* Match AA solution with PAA solution\n")
match <- greedy_match(aa$archetypes, t(paa$Z))


## ------------------------------------------------------------------------
thepaa <- paa$Z
thepaa <- perc(thepaa, emdata_mat)
colnames(thepaa) <- sprintf("PAA%s", 1:7)
thepaa <- thepaa[, match[, "PAA"]]
round(thepaa, 2)


## ------------------------------------------------------------------------
paa_perc <- paa
paa_perc$Z <- perc(paa$Z, emdata_mat)


## ------------------------------------------------------------------------
op <- par(mar = c(1.5, 1.5, 1.5, 1.5), xpd = TRUE)
set.seed(1239)
s <- simplexplot_paa(paa_perc, show_annotation = TRUE, show_circle = FALSE, 
                 show_points = TRUE, points_col = "black",
                 projection = tspsimplex_projection,
                 labels = colnames(thepaa))
par(op)


## ------------------------------------------------------------------------
theaa <- t(aa$archetypes)
theaa <- perc(theaa, emdata_mat)
colnames(theaa) <- sprintf(" AA%s", 1:7)
theaa <- theaa[, match[, "AA"]]
round(theaa, 2)


## ------------------------------------------------------------------------
aa$archetypes <- aa$archetypes[, match[, "AA"]]
aa$alphas <- aa$alphas[, match[, "AA"]]
aa$betas <- aa$betas[match[, "AA"], ]


## ------------------------------------------------------------------------
op <- par(mar = c(1.5, 1.5, 1.5, 1.5), xpd = TRUE)
set.seed(1239)
s <- simplexplot(aa, show_annotation = TRUE, show_circle = FALSE, 
                 show_points = TRUE, points_col = "black",
                 projection = tspsimplex_projection, 
                 labels = colnames(theaa))
par(op)

