format.series <- function (scoreDataFrame, gapSize = Inf,
    bridge = TRUE, check.sort = TRUE) {

# Author: Guillaume Filion.
# Date: September 14, 2009.
# scoreDataFrame: [GATCfragment, chromosome, start, end, xhat]
# -----------------------: [probeID, chromosome, start, end, y, x1, x2, x3, ...]

    D <- ncol(scoreDataFrame)-5

    scoreDataFrame[scoreDataFrame[,5] == 0,6:ncol(scoreDataFrame)] <- NA

    n <- nrow(scoreDataFrame)
    nonvirtuals <- rep(TRUE, n)

    # Compute the mean distance between consecutive start sites.
    distances <- diff(scoreDataFrame[-1,3])
    meanDistance <- mean(distances[distances > 0], na.rm = TRUE)

    # Segments are delimited by breaks: chromosome ends and gaps larger than gapSize.
    gaps <- scoreDataFrame[-1,3] - scoreDataFrame[-n,4]
    breaks <- gaps > gapSize | scoreDataFrame[-1,2] != scoreDataFrame[-n,2]

    # Virtual bridging: gaps are filled with NAs.
    if (bridge) {
        
        bridges <- which(gaps > meanDistance)
        bridges <- bridges[!bridges %in% breaks]
        virtuals <- as.integer(ceiling(gaps[bridges] / meanDistance))

        # Make a selection vector by 'zipping'.
        nonvirtuals <- rep(rep(c(TRUE, FALSE), length(bridges)+1),
            times = t(cbind(diff(c(0, bridges, n)), c(virtuals, 0))))
        n <- length(nonvirtuals)

        # Build the virtually bridged table.
        bridgedDataFrame <- data.frame(matrix(NA_character_, ncol = 2, nrow = n),
            matrix(NA_integer_, ncol = 3, nrow = n), matrix(NA_real_, ncol = D, nrow = n))
        colnames(bridgedDataFrame) <- colnames(scoreDataFrame)
        bridgedDataFrame[nonvirtuals,] <- scoreDataFrame
        scoreDataFrame <- bridgedDataFrame

        # Update the break positions.
        updatedBreaks <- rep(FALSE, n)
        updatedBreaks[nonvirtuals] <- c(breaks, FALSE)
        breaks <- updatedBreaks

    }

    breaks <- which(breaks)

    # Format the data in a segment-indexed list.
    segmentIndicator <- rep(1:(length(breaks)+1), times = diff(c(0, breaks, n)))
    
    singletons <- ! (segmentIndicator %in%
        unique(segmentIndicator[duplicated(segmentIndicator)]))
    segmentIndicator <- segmentIndicator[!singletons]
    scoreDataFrame <- scoreDataFrame[!singletons,]
    
    x <- tapply(X = scoreDataFrame[,6], INDEX = segmentIndicator,
        FUN = function(x) x)
    y <- tapply(X = scoreDataFrame[,5], INDEX = segmentIndicator,
        FUN = function(x) x)
    # Remove attributes of x to turn it into a proper list.
    attributes(x) <- NULL
    attributes(y) <- NULL

    segmentNames <- unique(segmentIndicator)

    if (check.sort) {
        checksort <- list()
        for (i in 1:length(x)) {
            checksort[[i]] <- scoreDataFrame[segmentIndicator == segmentNames[i],3]
            x[[i]] <- as.matrix(scoreDataFrame[segmentIndicator == segmentNames[i],6:(D+5)])
        }
        for (i in 1:length(checksort)) {
            if (is.unsorted(checksort[[i]], na.rm = TRUE)) {
                if (is.unsorted(-checksort[[i]], na.rm = TRUE)) {
                    stop ("data not sorted")
                }
                else {
                    x[[i]] <- rev(x[[i]])
                }
            }
        }
    }

    return(list(x = x, y = unlist(y), nonvirtuals = nonvirtuals))

}

