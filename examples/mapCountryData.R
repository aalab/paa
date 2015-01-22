

mapCountryData <-
function (mapToPlot = "", nameColumnToPlot = "", numCats = 7,
    xlim = NA, ylim = NA, mapRegion = "world", catMethod = "quantiles",
    colourPalette = "heat", addLegend = TRUE, borderCol = "grey",
    mapTitle = "columnName", oceanCol = NA, aspect = 1, missingCountryCol = NA,
    add = FALSE, nameColumnToHatch = "", borderLwd = 1)
{

    cat("\nHuom! You are using a slightly modified version\nof the mapCountryData function from the rworldmap package.\n\n")

    functionName <- as.character(sys.call()[[1]])
    require(sp)
    if (class(mapToPlot) == "SpatialPolygonsDataFrame") {
        if (length(mapToPlot@data[, 1]) < 1) {
            stop("seems to be no data in your chosen file or dataframe in ",
                functionName)
            return(FALSE)
        }
    }
    else if (mapToPlot == "") {
        message(paste("using example data because no file specified in",
            functionName))
        mapToPlot <- getMap(resolution = "coarse")
        if (nameColumnToPlot == "")
            nameColumnToPlot <- "POP_EST"
    }
    else {
        stop(functionName, " requires a SpatialPolygonsDataFrame object created by the joinCountryData2Map() function \n")
        return(FALSE)
    }
    if (is.na(match(nameColumnToPlot, names(mapToPlot@data)))) {
        stop("your chosen nameColumnToPlot :'", nameColumnToPlot,
            "' seems not to exist in your data, columns = ",
            paste(names(mapToPlot@data), ""))
        return(FALSE)
    }
    dataCategorised <- mapToPlot@data[[nameColumnToPlot]]
    if (!is.numeric(dataCategorised) && catMethod != "categorical") {
        catMethod = "categorical"
        message(paste("using catMethod='categorical' for non numeric data in",
            functionName))
    }
    if (length(catMethod) == 1 && catMethod == "categorical") {
        dataCategorised <- as.factor(dataCategorised)
        cutVector <- levels(dataCategorised)
    }
    else {
        if (is.character(catMethod) == TRUE) {
            cutVector <- rworldmap:::rwmGetClassBreaks(dataCategorised, catMethod = catMethod,
                numCats = numCats, verbose = TRUE)
        }
        else if (is.numeric(catMethod) == TRUE) {
            cutVector <- catMethod
        }
        dataCategorised <- cut(dataCategorised, cutVector, include.lowest = TRUE)
    }
    colNameRaw <- nameColumnToPlot
    colNameCat <- paste(colNameRaw, "categorised", sep = "")
    mapToPlot@data[[colNameCat]] <- dataCategorised
    numColours <- length(levels(dataCategorised))
    colourVector <- rworldmap:::rwmGetColours(colourPalette, numColours)
    dataCatNums <- as.numeric(dataCategorised)
    if (!is.na(missingCountryCol)) {
        colourVector <- c(colourVector, missingCountryCol)
        dataCatNums[is.na(dataCatNums)] <- length(colourVector)
    }
    hatchVar = NULL
    if (nameColumnToHatch == "") {
        if (!add)
          rworldmap:::rwmNewMapPlot(mapToPlot, mapRegion = mapRegion, xlim = xlim,
                ylim = ylim, oceanCol = oceanCol, aspect = aspect)
        plot(mapToPlot, col = colourVector[dataCatNums], border = borderCol,
            add = TRUE, usePolypath = FALSE, lwd = borderLwd)
    }
    else {
        hatchVar = mapToPlot@data[[nameColumnToHatch]]
        hatchVar = (hatchVar - min(hatchVar, na.rm = TRUE))/max(hatchVar,
            na.rm = TRUE)
        hatchVar = 1 - hatchVar
        hatchVar = (hatchVar * 50) + 30
        hatchVar[hatchVar > 79] = -1
        if (!add)
          rworldmap:::rwmNewMapPlot(mapToPlot, mapRegion = mapRegion, xlim = xlim,
                ylim = ylim, oceanCol = oceanCol, aspect = aspect)
        plot(mapToPlot, col = colourVector[dataCatNums], border = borderCol,
            density = hatchVar, angle = 135, lty = 1, add = TRUE,
            usePolypath = FALSE, lwd = 10)
        plot(mapToPlot, col = colourVector[dataCatNums], border = borderCol,
            density = hatchVar, angle = 45, lty = 1, add = TRUE,
            usePolypath = FALSE, lwd = 10)
    }
    if (addLegend) {
        if ((length(catMethod) == 1 && catMethod == "categorical") ||
            !require("spam") || !require("fields")) {
            addMapLegendBoxes(colourVector = colourVector, cutVector = cutVector,
                catMethod = catMethod)
        }
        else {
            addMapLegend(cutVector = cutVector, colourVector = colourVector,
                catMethod = catMethod)
        }
    }
    if (mapTitle == "columnName") {
        title(nameColumnToPlot)
    }
    else {
        title(mapTitle)
    }
    invisible(list(colourVector = colourVector, cutVector = cutVector,
        plottedData = mapToPlot[[nameColumnToPlot]], catMethod = catMethod,
        colourPalette = colourPalette))
}




image.plot <-
function (..., add = FALSE, nlevel = 64, horizontal = FALSE,
    legend.shrink = 0.9, legend.width = 1.2, legend.mar = ifelse(horizontal,
        3.1, 5.1), legend.lab = NULL, graphics.reset = FALSE,
    bigplot = NULL, smallplot = NULL, legend.only = FALSE, col = tim.colors(nlevel),
    lab.breaks = NULL, axis.args = NULL, legend.args = NULL,
    midpoint = FALSE, border = NA, lwd = 1)
{

    print("hello")
    old.par <- par(no.readonly = TRUE)
    info <- image.plot.info(...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink,
        legend.width = legend.width, legend.mar = legend.mar,
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            image(..., add = add, col = col)
        }
        else {
            poly.image(..., add = add, col = col, midpoint = midpoint,
                border = border, lwd.poly = lwd)
        }
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!is.null(breaks) & !is.null(lab.breaks)) {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4),
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2),
            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4),
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
            axis.args)
    }
    if (!horizontal) {
        if (is.null(breaks)) {
            print("1")
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col)
        }
        else {
            print("2")
            image(ix, iy, iz, xaxt = "n", yaxt = "n", bty = "n", xlab = "",
                ylab = "", col = col, breaks = breaks)
        }
    }
    else {
        if (is.null(breaks)) {
            print("3")
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col)
        }
        else {
            print("4")
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col, breaks = breaks)
        }
    }
    do.call("axis", c(axis.args, lty = 0, hadj = 1.5))
    ##box()
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal,
            1, 4), line = legend.mar - 2)
    }
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}
