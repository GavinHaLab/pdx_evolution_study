# Utils for create_matrix_updated.R (PDX CN analysis step 1)
# Written by Dr. Gavin Ha, updated by Anna Hoge
# Ha Lab
# Fred Hutchinson Cancer Research Center
# February 27th, 2020


##############################################
## FUNCTION TO GENERATE MATRIX BY REGIONS ####
##############################################
## regionAnnot must have at least 3 columns: chr, start, stop
createProbeMatrix <- function(x, regionAnnot, by="Coordinates"){
	require(IRanges)
	## if x is a character filename, then load the file
	if (class(x)=="character"){
		x <- read.delim(x, header=TRUE, stringsAsFactors=FALSE, sep="\t")
	}
	## convert x to a RangedData object ##
	x <- cbind(x,x[,4])
	colnames(x)[c(3,4,ncol(x))] <- c("space","start","end")
	x <- as(x, "RangedData")

	if (class(regionAnnot)=="character"){
		regionAnnot <- read.delim(regionAnnot, header=TRUE, stringsAsFactors=FALSE, sep="\t")
	}
	if (class(regionAnnot)=="data.frame"){
		colnames(regionAnnot) <- c("Ensembl.Gene.ID", "Associated.Gene.Name",
								"space","start","end", "Band","Status..gene.")
	}
	## convert regionAnnot to RangedData object ##
	if (class(regionAnnot) != "RangedData"){
		regionAnnot <- as(regionAnnot, "RangedData")
	}

	## create the list of column headers as the gene_chrom_position
	rowNames <- regionAnnot$Associated.Gene.Name
	colNames <- colnames(x)[-c(1:2)]
	## instantiate the matrix
	mat <- matrix(NA, nrow=length(rowNames), ncol=length(colNames),
				dimnames=list(rowNames,colNames))

	if (by=="Coordinates"){
		## find overlap
		hits <- findOverlaps(query=x, subject=regionAnnot)
		indAnnot <- subjectHits(hits)
		indX <- queryHits(hits)
		## assign each sample (column) to the matrix
		for (i in 1:length(colNames)){
			mat[indAnnot, colNames[i]] <- as.data.frame(x[indX,])[, colNames[i]]
		}
	}else if (by=="Genes"){
		indGene <- !duplicated(x$name) & x$name %in% rowNames & x$name != " "
		for (i in 1:length(colNames)){
			mat[x$name[indGene], colNames[i]] <- as.data.frame(values(x))[indGene, colNames[i]]
		}
	}

	return(mat)

}

createRegionMatrix <- function(segs, regionAnnot, annot.type="Band", seg.colName="Segment_Mean", fun = c("largest", "severe", "weighted_average")){
	require(IRanges)
	require(plyr)
	require(foreach)
	require(data.table)
	fun <- match.arg(fun)
	samples <- unique(segs$Sample)
	#mat <- matrix(NA, nrow=length(regionAnnot), ncol=length(samples),
	#							dimnames=list(values(regionAnnot)[,annot.type], samples))
	mat <- foreach(i = 1:length(samples), .combine = cbind)  %dopar% {
	#for (i in 1:length(samples)){
		id <- samples[i]
		message("Analyzing ", id)
		segSample <- segs[segs$Sample == id, ]
      ###added check.names = FALSE --Anna 2/10/2020
	  regionStates <- data.frame(matrix(NA, nrow=length(regionAnnot), ncol=1,
	  									dimnames=list(values(regionAnnot)[,annot.type], samples[i])), check.names = FALSE)
  	hits <- findOverlaps(subject = segSample, query = regionAnnot)
  	segInd <- subjectHits(hits)
		annotInd <- queryHits(hits)
	  annotNameOverlap <- as.data.frame(values(regionAnnot[annotInd, annot.type]))[, 1]
		# assign all CN calls to cnMat
		regionStates[annotNameOverlap, id] <- values(segSample[segInd, seg.colName])[,1]

	  ## resolve genes overlapping multiple segments ##
		dupInd <- annotInd[which(duplicated(annotInd))]
		if (length(dupInd) > 0){
			dupAnnotNames <- as.data.frame(values(regionAnnot[dupInd, annot.type]))[, 1]
			if (fun == "severe"){
                dupCN <- adply(dupInd, .margins=1, .fun=function(x){
                    segDupInd <-segInd[annotInd == x]
                    severityInd <- which.max(abs(segSample[segDupInd]$Segment_Mean))
                    dupResults <- segSample[segDupInd[severityInd], "Segment_Mean"]
                    as.data.frame(values(dupResults))
                })

			###added by Anna 2/10/2020
            }else if (fun == "weighted_average"){
                    dupCN <- adply(dupInd, .margins=1, .fun=function(x){
                        segDupInd <-segInd[annotInd == x]

                        numerator <- 0
                        for(i in range(1:length(segDupInd))){
                            segSampleIndex <- segDupInd[i]
                            segSampleScore <- segSample[segSampleIndex]$Segment_Mean
                            segSampleLength <- sum(width(reduce(segSample[segSampleIndex])))
                            numerator <- numerator + segSampleScore * segSampleLength
                        }
                        denominator <- sum(width(reduce(segSample[segDupInd])))
                        dupResults <- numerator/denominator
                        data.frame(data.frame("Segment_Mean" = dupResults))
                    })

            ###added 2/25/2020
			}else if (fun == "largest"){
                dupCN <- adply(dupInd, .margins=1, .fun=function(x){
                    segDupInd <- segInd[annotInd == x]
                    #dupHit <- findOverlaps(query=segSample[segDupInd,], subject=regionAnnot[x,])
                    lengths <- width(reduce(segSample[segDupInd, ]))
                    dupResults <- values(segSample[segDupInd[which.max(lengths)], seg.colName])[[1]]
                    data.frame("Segment_Mean" = dupResults)
                })
            }

			regionStates[dupAnnotNames, id] <- dupCN$Segment_Mean
		}
		regionStates
	}
	return(mat)
}

















#######other functions provided in original script Gavin wrote:


aggregateSegs <- function(x, mat, fun, regionAnnot){
	if (length(x) > 0){
		if (fun == "largest"){
			xMat <- as.data.frame(mat[x, ])
			state <- xMat$Segment_Mean
			lengths <- getOverlapLength(xMat, start(ranges(regionAnnot[i, ])), end(ranges(regionAnnot[i, ])))
			stateLength <- aggregate(lengths, by=list(state), FUN=sum)
			state <- stateLength[which.max(stateLength[,2]), 1]
			return(state)
		}
	}else{
		return(0)
	}
}


aggregateSegs <- function(i, hits, mat, fun, regionAnnot){
	x <- hits[[i]]
	if (length(x) > 0){
		if (fun == "largest"){
			xMat <- as.data.frame(mat[x, ])
			state <- xMat$Segment_Mean
			lengths <- getOverlapLength(xMat, start(ranges(regionAnnot[i, ])), end(ranges(regionAnnot[i, ])))
			stateLength <- aggregate(lengths, by=list(state), FUN=sum)
			state <- stateLength[which.max(stateLength[,2]), 1]
			return(state)
		}
	}else{
		return(0)
	}
}

## input: lohHits = loh rows that overlap region of interest
# start = start coordinate of region of interest
# end = end coordinate of region of interest
getOverlapLength <- function(lohHits, start, end){
	coords <- cbind(lohHits[, c("start","end")], as.numeric(start), as.numeric(end))
	coordsSort <- t(apply(coords, 1, sort))
	dist <- coordsSort[, 3] - coordsSort[, 2] + 1
	return(dist)
}


####################################################
###### FUNCTION TO COMPUTE FREQUENCY ###############
####################################################
computeFrequency <- function(mat, neutral=2, countNA=TRUE, ind = TRUE, numSamples = NULL){
	if (is.null(numSamples)){
		if (countNA){
			numSamples <- ncol(mat)
		}else{
			numSamples <- rowSums(!is.na(mat))
		}
	}
	indGain <- mat > neutral & ind
	gainFreq <- rowSums(indGain, na.rm=TRUE) / numSamples
	indLoss <- mat < neutral & ind
	lossFreq <- rowSums(indLoss, na.rm=TRUE) / numSamples
	return(list(gainFreq=gainFreq, lossFreq=lossFreq))
}



##############################################
##### FUNCTION TO PLOT GENE FREQUENCY ########
##############################################
## annot must have at least 3 columns: chr, start, stop
## annot must correspond to same rows as lossFreq and gainFreq
plotFrequency <- function(lossFreq, gainFreq, annot, type = "lines", chrs = c(1:22,"X"), xlim=NULL, yAxis=seq(-1,1,0.1), xAxis="chrs", build="hg19", ...){
	require(GenomeInfoDb)
	require(GenomicRanges)
	if (nrow(annot) != length(lossFreq) || nrow(annot) != length(gainFreq)){
		stop("annot does not have the same number of rows as lossFreq and gainFreq.")
	}

	dots <- list(...)
	#colnames(annot)[1] <- c("seqnames")
	#annot <- as(annot, "GRanges")
	indChrs <- annot$chrs %in% chrs
	annot <- annot[indChrs, ]

	# get chr info from USCS using GenomeInfoDB #
  chrInfo <- fetchExtendedChromInfoFromUCSC(build, quiet = TRUE)
  chrLens <- chrInfo$UCSC_seqlength[chrInfo$NCBI_seqlevel %in% chrs]

	## get plotting coords (X-axis) ##
	#midpt <- ((annot$stop[indChrs] - annot$start[indChrs] + 1) / 2) + annot$start[indChrs]
	coord <- getGenomeWidePositions(annot, chrLens)
	if (xAxis == "arms"){
		colnames(coord$regions)[c(1,4)] <- c("chrName", "chrs")
		coord$chrBkpt <- c(0, coord$regions$stop)
	}
	#coordStart <- getGenomeWidePositions(annot$chr[indChrs], annot$start[indChrs])
	#coordStop <- getGenomeWidePositions(annot$chr[indChrs], annot$stop[indChrs])

	# setup plot #
	if (is.null(xlim)){
		xlim <- as.numeric(c(1, tail(coord$regions$stop, 1)))
	}
	if (is.null(dots$cex.axis)){
		cex.axis <- 1
	}else{
		cex.axis <- dots$cex.axis
	}

	plot(0, type="n", bty="n", yaxt="n", xaxt="n", xlim=xlim, ...)
	axis(2, at=yAxis, label=yAxis, las=2, cex=cex.axis)
	if (type == "lines"){
		# plot loss frequency #
		lines(coord$posn, -1*lossFreq[indChrs], type = "h", col="blue")
		# plot gain frequency #
		lines(coord$posn, gainFreq[indChrs], type = "h", col="red")
	}else if (type == "bars"){
		for (j in 1:nrow(coord$regions)){
			# plot loss frequency
			rect(xleft=coord$regions$start[j], ybottom=-1*lossFreq[indChrs][j],
			     xright=coord$regions$stop[j], ytop=0, col="blue")
			# plot gain frequency
			rect(xleft=coord$regions$start[j], ybottom=0,
			     xright=coord$regions$stop[j], ytop=gainFreq[indChrs][j], col="red")
		}
	}
	plotChrLines(unique(coord$regions$chrs), coord$chrBkpt, c(-1, 1))

}

##############################################
##### FUNCTION TO PLOT pvalue ########
##############################################
## annot must have at least 3 columns: chr, start, stop
## annot must correspond to same rows as lossFreq and gainFreq
plotGenomeWidePval <- function(loss, gain, annot, signif=0.05, chrs = c(1:22,"X"), plot.new=TRUE, plot.type="points", ...){

	dots <- list(...)

	indChrs <- annot$chr %in% chrs
	## get plotting coords (X-axis) ##
	coord <- getGenomeWidePositions(annot)

	# setup plot #
	if (is.null(dots$xlim)){
		dots$xlim <- as.numeric(c(1, tail(coord$posns, 1)))
	}
	if (is.null(dots$ylim)){
		dots$ylim <- c(0, max(gain,loss, 1, na.rm=T))
	}
	if (is.null(dots$cex.axis)){
		dots$cex.axis <- 1
	}


	if (plot.new){
		do.call(plot, c(0, type="n", bty="n", yaxt="n", xaxt="n", dots))
	}
	axis(2, las=2, cex.axis=dots$cex.axis)

	sigLossInd <- loss[indChrs] < -log(signif)
	sigGainInd <- gain[indChrs] < -log(signif)
	if (plot.type == "points"){
		# plot data #
		points(coord$posn, loss[indChrs], pch=19, col=rgb(0,0,1,0.75))
		points(coord$posn, gain[indChrs], pch=19, col=rgb(1,0,0,0.75))

		# plot grey for non-significant data #
		points(coord$posn[sigLossInd], loss[indChrs][sigLossInd], pch=19, col="grey")
		points(coord$posn[sigGainInd], gain[indChrs][sigGainInd], pch=19, col="grey")
	}else if (plot.type == "lines"){
		# plot data #
		x.coords <- cbind(c(1, coord$posns[-length(coord$posns)]), coord$posns)
		tmp <- apply(cbind(x.coords, loss[indChrs]), 1, function(x){
			lines(x[c(1,2)], rep(x[3], 2), lwd=5, col=rgb(0,0,1,0.75))
		})
		tmp <- apply(cbind(x.coords, gain[indChrs]), 1, function(x){
			lines(x[c(1,2)], rep(x[3], 2), lwd=5, col=rgb(1,0,0,0.75))
		})

		tmp <- apply(cbind(x.coords[sigLossInd, ], loss[indChrs][sigLossInd]), 1, function(x){
			lines(x[c(1,2)], rep(x[3], 2), lwd=6, col="black")
		})
		tmp <- apply(cbind(x.coords[sigGainInd, ], gain[indChrs][sigGainInd]), 1, function(x){
			lines(x[c(1,2)], rep(x[3], 2), lwd=6, col="black")
		})
	}


    if (plot.new){
    	plotChrLines(unique(annot$chr[indChrs]), coord$chrBkpt, yrange=dots$ylim, cex = dots$cex.axis)
    }
}


##############################################
############# HELPER FUNCTION ################
##############################################
plotChrLines <- function(chrs, chrBkpt, yrange, cex = 1) {
    # plot vertical chromosome lines
    for (j in 1:length(chrBkpt)) {
        lines(rep(chrBkpt[j], 2), yrange, type = "l",
            lty = 2, col = "black", lwd = 0.75)
    }
    numLines <- length(chrBkpt)
    mid <- (chrBkpt[1:(numLines - 1)] + chrBkpt[2:numLines])/2
    #chrs[chrs == "X"] <- 23
    #chrs[chrs == "Y"] <- 24
    #chrsToShow <- sort(unique(as.numeric(chrs)))
    #chrsToShow[chrsToShow == 23] <- "X"
    #chrsToShow[chrsToShow == 24] <- "Y"
    axis(side = 1, at = mid, line = -1, labels = chrs, cex.axis = cex, tick = FALSE)
    #text(y = yrange[1] * 1.1, x = mid, labels = chrs, cex = cex, xpd = TRUE)
}

##############################################
############# HELPER FUNCTION ################
##############################################
# regions is data frame with "chrs", "start", "stop"
getGenomeWidePositions <- function(regions, chrLens = NULL) {
    # create genome coordinate scaffold
    #positions <- as.numeric(posns)
    if (class(regions) == "GRanges"){
    	regions <- as.data.frame(regions)
    	colnames(regions)[1:3] <- c("chrs", "start", "stop")
    }
    positions <- regions
    chrsNum <- unique(regions$chrs)
    chrBkpt <- rep(0, length(chrsNum) + 1)
    if (!is.null(chrLens)){
			prevChrEnd <- as.numeric(chrLens[1])
		}else{
			prevChrEnd <- as.numeric(max(regions[regions$chrs == 1, c("start","stop")]))
		}
    for (i in 2:length(chrsNum)) {
        chrInd <- which(regions$chrs == chrsNum[i])
        if (!is.null(chrLens)){
        	chrMaxEnd <- as.numeric(chrLens[i])
        }else{
        	chrMaxEnd <- as.numeric(max(regions[regions$chrs == chrsNum[i], c("start","stop")]))
        }
        chrBkpt[i] = prevChrEnd
        positions$start[chrInd] = positions$start[chrInd] + prevChrEnd
        positions$stop[chrInd] = positions$stop[chrInd] + prevChrEnd
      	prevChrEnd <- prevChrEnd + chrMaxEnd
    }
    chrBkpt[i + 1] <- prevChrEnd
    posns <- (positions$stop - positions$start) / 2 + positions$start
    return(list(posns = posns, chrBkpt = chrBkpt, regions = positions))
}
