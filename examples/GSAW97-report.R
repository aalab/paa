
## ----echo=FALSE, message=FALSE, warnings=FALSE---------------------------
opts_chunk$set(echo = FALSE, messages = FALSE, 
               fig.height = 6, fig.width = 6)


## ----results='hide', message=FALSE, warning=FALSE, error=FALSE-----------
source("helpers.R")


## ----results='asis'------------------------------------------------------
cat("* Read tourist data from file `GSAW97_bin.csv`\n")
GSAW97bin <- read.csv2("data_GSAW97/GSAW97_bin.csv")

cat("* Load precomputed PAA solution from file `GSAW97_paa.mat`\n")
paa <- load_paa("data_GSAW97/GSAW97_paa.mat", GSAW97bin)


## ----fig.height=4--------------------------------------------------------
op <- par(mar = c(4, 4, 1, 1))
plot(paa$obj, type = "l", xlab = "Iterations", ylab = "Cost")
par(op)


## ----results='asis'------------------------------------------------------
colnames(paa$Z) <- c("Maximal", "Cultural", "Minimal", "Modern", "Traditional", "Wellness")

print(xtable(paa$Z[, c(1, 3, 5, 4, 6, 2)]), type = "html")


## ------------------------------------------------------------------------
## Color coding:
dev <- bernoulli_deviance(paa, GSAW97bin)

val <- (dev / max(dev))
val <- (1 - (exp(-val*2)))/(1-exp(-2))
ci <- classIntervals(val, style = "pretty", n = 10)
point_cols <- findColours(ci, sequential_hcl(11))


## Simplex plot:
op <- par(mar = c(1.5, 1.5, 1.5, 1.5), xpd = TRUE)

set.seed(1239)
s <- simplexplot_paa(paa, show_direction = TRUE, direction_length = 1,
  show_circle = FALSE, show_points = FALSE, directions_col = point_cols,
  projection = tspsimplex_projection, show_edges = FALSE, 
  labels = colnames(paa$Z))

legend("topleft",  fill = attr(point_cols, "palette"), 
       border = attr(point_cols, "palette"), bg = "white", horiz = FALSE,
       legend = names(attr(point_cols, "table")), cex = 1, box.lwd = NA)

par(op)


## ------------------------------------------------------------------------
var_cod <- c(Yes = "black", No = "lightgray")

## Color coding:
var <- "Shopping"
point_cols <- GSAW97bin[, grepl(tolower(var), colnames(GSAW97bin))]


## Simplex plot:
op <- par(mar = c(1.5, 1.5, 1.5, 1.5), xpd = TRUE)

set.seed(1239)
s <- simplexplot_paa(paa, show_direction = FALSE, 
  show_circle = FALSE, show_points = FALSE,
  projection = tspsimplex_projection, show_edges = FALSE, 
  labels = colnames(paa$Z))

points(s$proj_h[point_cols == 0, ], col = var_cod["No"], pch = 19, cex = 0.5)
points(s$proj_h[point_cols == 1, ], col = var_cod["Yes"], pch = 19, cex = 0.5)

legend("topleft", fill = var_cod, border = var_cod, title = var, 
       legend = names(var_cod), bg = "white", horiz = FALSE, cex = 1,
       box.lwd = NA)

par(op)


## ------------------------------------------------------------------------
## Color coding:
var <- "Snowboard"

point_cols <- GSAW97bin[, grepl(tolower(var), colnames(GSAW97bin))]


## Simplex plot:
op <- par(mar = c(1.5, 1.5, 1.5, 1.5), xpd = TRUE)

set.seed(1239)
s <- simplexplot_paa(paa, show_direction = FALSE, 
  show_circle = FALSE, show_points = FALSE,
  projection = tspsimplex_projection, show_edges = FALSE, 
  labels = colnames(paa$Z))

points(s$proj_h[point_cols == 0, ], col = var_cod["No"], pch = 19, cex = 0.5)
points(s$proj_h[point_cols == 1, ], col = var_cod["Yes"], pch = 19, cex = 0.5)

legend("topleft", fill = var_cod, border = var_cod, title = var, 
       legend = names(var_cod), bg = "white", horiz = FALSE, cex = 1,
       box.lwd = NA)

par(op)


## ------------------------------------------------------------------------
## Color coding:
var <- "Museum"
point_cols <- GSAW97bin[, grepl(tolower(var), colnames(GSAW97bin))]


## Simplex plot:
op <- par(mar = c(1.5, 1.5, 1.5, 1.5), xpd = TRUE)

set.seed(1239)
s <- simplexplot_paa(paa, show_direction = FALSE, 
  show_circle = FALSE, show_points = FALSE,
  projection = tspsimplex_projection, show_edges = FALSE, 
  labels = colnames(paa$Z))

points(s$proj_h[point_cols == 0, ], col = var_cod["No"], pch = 19, cex = 0.5)
points(s$proj_h[point_cols == 1, ], col = var_cod["Yes"], pch = 19, cex = 0.5)

legend("topleft", fill = var_cod, border = var_cod, title = var, 
       legend = names(var_cod), bg = "white", horiz = FALSE, cex = 1,
       box.lwd = NA)

par(op)


