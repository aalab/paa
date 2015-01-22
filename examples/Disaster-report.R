
## ----echo=FALSE, message=FALSE, warnings=FALSE---------------------------
opts_chunk$set(echo = FALSE, messages = FALSE, 
               fig.height = 6, fig.width = 6)


## ----results='hide', message=FALSE, warning=FALSE, error=FALSE-----------
library("ggplot2")
library("reshape")
library("rworldmap")
library("colorspace")
library("maptools")
library("plotrix")

data("wrld_simpl")

source("mapCountryData.R")
source("helpers.R")


## ----results='asis'------------------------------------------------------
cat("* Read tourist data from file `emdata.csv`\n")
emdata <- read.table("data_Disaster/emdata.csv", stringsAsFactors = FALSE, sep = ";", header = TRUE)
emdata_mat <- as.matrix(emdata[, -1])

cat("* Load precomputed PAA solution from file `Disaster_paa.mat`\n")
paa <- load_paa("data_Disaster/Disaster_paa.mat", emdata_mat)


## ----fig.height=4--------------------------------------------------------
op <- par(mar = c(4, 4, 1, 1))
plot(paa$obj, type = "l", xlab = "Iterations", ylab = "Cost")
par(op)


## ----results='asis'------------------------------------------------------
colnames(paa$Z) <- sprintf("A%s", 1:7)

print(xtable(paa$Z), type = "html")


## ------------------------------------------------------------------------
aa_perc <- t(perc(paa$Z, emdata_mat))

colnames(aa_perc) <- colnames(emdata_mat)
rownames(aa_perc) <- sprintf("A%s", seq(ncol(paa$Z)))
colnames(aa_perc) <- c("CD", "DR", "EQ", "EP", "ET", "FL",
                       "IA", "II", "MD", "MW", "MA", "ST",
                       "TA", "VO", "WF")

aa_perc[aa_perc > 1] <- 1.0
aa_perc <- aa_perc * 100

aa_perc <- melt(aa_perc)
aa_perc <- rename(aa_perc, c(X2 = "Disaster", X1 = "Archetype"))

ot <- theme_set(theme_bw(base_size = 9))

ggplot(aa_perc, aes(Disaster, value)) +
    geom_bar(stat = "identity") +
    facet_grid(Archetype ~ ., scales = "fixed") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_continuous(name = "Archetype profiles (%)")


## ------------------------------------------------------------------------
plot_map <- function(H, i, borderLwd = 0.01, borderCol = "white") {
  dat <- data.frame(H3 = H[i, ], Country = emdata$Country, stringsAsFactors = FALSE)
  map <- joinCountryData2Map(dat, joinCode = "NAME", nameJoinColumn = "Country")

  cols <- rev(heat_hcl(51, c = c(80, 30), l = c(30, 90), power = c(1/5, 2)))
  cols[1] <- "lightgray"

  op <- par(mai = c(0,0,0,0), xaxs = "i", yaxs = "i")
  mp <- mapCountryData(map, nameColumnToPlot = "H3",
                       catMethod = "pretty", numCat = 50,
                       colourPalette = cols, addLegend = FALSE,
                       mapTitle = "", borderCol = borderCol,
                       borderLwd = borderLwd,
                       missingCountryCol = "lightgray")
  text(-152, 40, sprintf("A%s", i))

  invisible(mp)
}


## ----message=FALSE, warning=FALSE, results='hide', fig.height=3, fig.width=6----
plot_map(paa$H, 1)
plot_map(paa$H, 2)
plot_map(paa$H, 3)
plot_map(paa$H, 4)
plot_map(paa$H, 5)
plot_map(paa$H, 6)
plot_map(paa$H, 7)


## ------------------------------------------------------------------------
paa_perc <- paa
paa_perc$Z <- perc(paa$Z, emdata_mat)


## ------------------------------------------------------------------------
dev <- poission_deviance(paa_perc, emdata_mat)
val <- (dev / max(dev))
val <- (1 - (exp(-val * 10))) / (1 - exp(-10))
ci <- classIntervals(val, style = "pretty", n = 10)
point_cols <- findColours(ci, sequential_hcl(11))


## ------------------------------------------------------------------------
op <- par(mar = c(1.5, 1.5, 1.5, 1.5), xpd = TRUE)
set.seed(1239)
simplexplot_paa(paa_perc, show_circle = FALSE, 
             show_points = TRUE, points_col = point_cols,
             projection = tspsimplex_projection, show_edges = TRUE, labels_cex = 1)

t <- names(attr(point_cols, "table"))

legend("bottomleft", fill = attr(point_cols, "palette"), 
       border = attr(point_cols, "palette"),
       legend = t, bg = "white",
       horiz = FALSE, cex = 0.8, box.lwd = NA)

par(op)


## ------------------------------------------------------------------------
## Variable of interest:
i <- 8

## Color coding:
point_cols <- ifelse(emdata_mat[, i] > 0, "black", "darkgray")
dir_cols <- ifelse(emdata_mat[, i] > 0, "black", NA)
iicntry <- emdata_mat[, i]/max(emdata_mat[, i])
point_cexs <- ifelse(iicntry == 0, 0.4, iicntry)


## Simplex plot:
op <- par(mar = c(1.5, 1.5, 1.5, 1.5), xpd = TRUE)
set.seed(1239)
s <- simplexplot_paa(paa_perc, show_circle = FALSE, points_col = point_cols, show_direction=TRUE,
                  show_points = TRUE, directions_col = dir_cols, points_cex = point_cexs * 2,
                  projection = tspsimplex_projection, show_edges = TRUE, labels_cex = 1)

iso2 <- match(emdata$Country, wrld_simpl$NAME)
iso2 <- wrld_simpl$ISO2[iso2]
 
w <- which(iicntry > 0.5)
w <- w[order(iicntry[w], decreasing = TRUE)]

pc <- function(i, xshift = 0, yshift = 0) {
  boxed.labels(s$proj_h[i, 1, drop = FALSE] + xshift, s$proj_h[i, 2, drop = FALSE] + yshift,
               labels = iso2[i], bg = "white", col = "black", border = FALSE, cex = 1)
}

pc(w[1], yshift = 0.8)
pc(w[2], yshift = 0.8)
pc(w[3], yshift = -0.8)
pc(w[4], yshift = -0.8)
pc(w[5], yshift = 0.8)
pc(w[6], xshift = -0.8)
pc(w[7], yshift = -0.8)
pc(w[8], yshift = -0.8)
pc(w[9], yshift = 0.8)

p <- seq(0, 1, length.out = 5)

legend("bottomleft", pch = 19, col = "black", 
       title = "Insect\ninfestation",
       pt.cex = p * 2,
       legend = paste(p*100, "%", sep = ""), bg = "white",
       horiz = FALSE, cex = 1, box.lwd = NA)

par(op)


