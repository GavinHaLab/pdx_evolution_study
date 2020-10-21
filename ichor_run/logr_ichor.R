###imports
library(HMMcopy)
library(GenomicRanges)
library(GenomeInfoDb)
library(tidyverse)
library(foreach)
library(doMC)
library(optparse)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')
options(scipen=0, stringsAsFactors=F)

option_list <- list(
  make_option(c("-i", "--inputMatrixFile"), type = "character", help = "Path to file containing rds file of the cohort matrix. Required."),
  make_option(c("-s", "--sampleFile"), type = "character", default = NULL, help = "Path to file containing text file with list of samples to run. Optional."),
  make_option(c("--ploidy"), type="character", default="c(2,3)", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
  make_option(c("--cores"), type="numeric", default = 1, help = "Number of cores to use for EM. Default: [%default]"),
  make_option(c("-o", "--outDir"), type="character", default = "./", help="Output directory")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)

###load in ichor functions
source("utils.R")
source("/fh/fast/ha_g/user/ahoge/software/R/ichorCNA_ghlab/R/utils.R")
#source("/fh/fast/ha_g/user/ahoge/software/R/ichorCNA_ghlab/R/segmentation.R")
source("/fh/fast/ha_g/user/gha/software/git/GavinHaLab/ichorCNA/R/segmentation.R")
source("/fh/fast/ha_g/user/gha/software/git/GavinHaLab/ichorCNA/R/EM.R")
source("/fh/fast/ha_g/user/gha/software/git/GavinHaLab/ichorCNA/R/output.R")
source("plotting.R")

inputMatrixFile <- opt$inputMatrixFile
sampleFile <- opt$sampleFile
ploidy <- eval(parse(text = opt$ploidy))
outDirPrefix <- opt$outDir

registerDoMC()
options(cores = opt$cores)
message("ichorCNA: Using ", getDoParWorkers(), " core(s).")

###get logR data from each sample into one big matrix of granges called logR.data
#in /fh/fast/ha_g/projects/PDX/ichor
args = commandArgs(trailingOnly = TRUE)
print(args)
input_data <- readRDS(inputMatrixFile)
#input_data <- readRDS("/fh/fast/ha_g/projects/PDX/matrices_full_set_1mb_severe/PDX_armMatrix_bcm_snp_breast.rds")


chrs_with_data <- as.character(c(1:22))
cohort_name <- names(input_data)
print(cohort_name)
#dir.create(cohort_name)
data <- input_data[[cohort_name]]
samples <- names(data)

if (!is.null(sampleFile)){
  sampleList <- read.delim(sampleFile, header = TRUE, row.names = 1)
  samples <- samples[samples %in% rownames(sampleList)]

}


###configure ichor settings
gender <- "female" #this shouldn't matter because we're only evaluating autosomes
chrTrain <- as.character(c(1:22))
estimateNormal <- TRUE
estimatePloidy <- TRUE
estimateScPrevalence <- TRUE
maxCN <- 8
includeHOMD <- FALSE
scStates <- c(1, 3)
normal2IgnoreSC <- 1
#ploidy <- input_ploidy #c(2, 3) ###
txnE <- 0.99
subclone.penalty <- 1
txnStrength <- 100
likModel <- "t"
minSegmentBins <- 50
altFracThreshold <- 0.05
plotFileType <- "pdf"
S <- 1
maxiter <- 50
#mainName <- vector('list', S) #rep(NA, nrow(normal.restarts) * nrow(ploidy.restarts))
genomeStyle <- "NCBI"
genomeBuild <- "hg38"
seqinfo <- getSeqInfo(genomeBuild, genomeStyle, chrs_with_data)
plotYLim <- c(-2, 4)
coverage <- NULL
s <- 1
maxFracCNASubclone <- 0.75
maxFracGenomeSubclone <- 0.75
numSamples <- 1


num_samples_in_cohort <- length(samples)
tmp <- foreach(k = 1:num_samples_in_cohort) %dopar% {
#for (k in 1:num_samples_in_cohort){
    id <- samples[[k]]
    message("ichorCNA: Running analysis for ", id)
    logR.data <- GRangesList()

    sample_name <- samples[[k]]
    outDir <- paste0(outDirPrefix, "/", cohort_name, "/", sample_name)
    dir.create(outDir, recursive = TRUE)
    dir.create(paste0(outDir, "/", id, "/"), recursive = TRUE)
    
    sample_data <- data[sample_name]
    patientID <- sample_name
    message("Loading normalized log ratios: ", sample_name)

    #if the sample is a PT, use different normals than otherwise
    if(grepl("PT", sample_name, fixed = TRUE)){
        normal <- c(0.25, 0.5, 0.75)
        #normal <- c(0.5)
        sp_hyperA <- 200; sp_hyperB <- 200
        n_hyperA <- 1; n_hyperB <- 1
    }else{
        normal <- c(0.05)
        sp_hyperA <- 200; sp_hyperB <- 200
        n_hyperA <- 1; n_hyperB <- 1000
    }
    
    # get ploidy
    if (is.null(ploidy)){
      # ensure it is never smaller than 2 
      ploidy <- max(round(sampleList[sample_name, 1]), 2)
    }
    
    sample_data_df <- sample_data %>%
      rownames_to_column() %>%
      separate(rowname, "[[:punct:]]", into = c("chr", "start", "end"), remove = TRUE) %>%
      add_column(valid = TRUE) %>%
      filter(chr %in% chrs_with_data)
    
    names(sample_data_df)[4] <- "copy"
    sample_data_granges <- makeGRangesFromDataFrame(sample_data_df, keep.extra.columns = TRUE)
    
    logR.in <- sample_data_granges
    logR.data[[sample_name]] <- logR.in
        
    valid <- logR.data[[1]]$valid & !is.na(logR.data[[1]]$copy)
    chrInd <- as.character(seqnames(logR.data[[1]])) %in% chrTrain


    ###storing results for different runs
    ###
    counter <- 1
    ptmTotalSolutions <- proc.time() # start total timer
    results <- list()
    numCombinations <- (length(normal) * length(ploidy)) ^ S
    loglik <- as.data.frame(matrix(NA, nrow = numCombinations, ncol = 7, 
                     dimnames = list(c(), c("n_0", "phi_0", "n_est", "phi_est", 
                     												"Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
    fracGenomeSub <- as.data.frame(matrix(NA, nrow = numCombinations, ncol = S))
    fracAltSub <- as.data.frame(matrix(NA, nrow = numCombinations, ncol = S))
    normal.restarts <- expand.grid(rep(list(normal), S))
    compNames <- rep(NA, nrow(loglik))
    mainName <- vector('list', S) #rep(NA, nrow(normal.restarts) * nrow(ploidy.restarts))

    for (i in 1:length(ploidy)){
      p <- rep(ploidy[i], S)
      for (j in 1:nrow(normal.restarts)){
        n <- as.numeric(normal.restarts[j, ])
        ## skip restarts where normal=0.95 and ploidy not diploid (2)
        if (sum(n == 0.95 & p != 2) > 0) {
          next
        }
        message("Running EM for normal=", paste(n, collapse=","), ", ploidy=", paste(p, collapse=","))

        ###run ichor
        logR <- as.data.frame(lapply(logR.data, function(x) { x$copy })) # NEED TO EXCLUDE CHR X #
        param <- getDefaultParameters(logR[valid & chrInd, , drop=F], n_0 = n, maxCN = maxCN, includeHOMD = includeHOMD, 
                    ct.sc=scStates, normal2IgnoreSC = normal2IgnoreSC, ploidy_0 = floor(p), 
                    e=txnE, e.subclone = subclone.penalty, e.sameState = 50, strength=txnStrength, likModel = likModel)
        K <- length(param$ct)
        logR.lambda <- (1 / ( apply(logR, 2, sd, na.rm = TRUE) / K ) ^ 2) / 100
        param$lambda <- matrix(logR.lambda, K)
        param$lambda[param$ct %in% c(0)] <- logR.lambda / 5
        param$lambda[param$ct %in% c(1)] <- logR.lambda * 1
        param$lambda[param$ct %in% c(2)] <- logR.lambda * 3
        param$lambda[param$ct %in% c(3)] <- logR.lambda * 1.25
        param$lambda[param$ct >= 4] <- logR.lambda * 1
        param$lambda[param$ct == max(param$ct)] <- logR.lambda * 1
        param$lambda[param$ct.sc.status] <- logR.lambda * 1
        param$betaLambda <- (param$alphaLambda / param$lambda) / 1
        param$alphaSp <- sp_hyperA; param$betaSp <- sp_hyperB
        param$alphaN <- n_hyperA; param$betaN <- n_hyperB

        logR.var <- ( apply(logR, 2, sd, na.rm = TRUE) / K ) ^ 2 * 50
        param$var <- matrix(logR.var, ncol = S, nrow = K, byrow = TRUE)
        param$var[param$ct %in% c(2)] <- logR.var 
        param$var[param$ct %in% c(1,3)] <- logR.var 
        param$var[param$ct >= 4] <- logR.var
        param$var[param$ct == max(param$ct)] <- logR.var 
        param$var[param$ct.sc.status] <- logR.var * 2
        param$alphaVar <- param$betaVar / param$var

        hmmResults.cor <- HMMsegment(logR.data, valid, dataType = "copy", 
                                     param = param, chrTrain = chrTrain, maxiter = maxiter,
                                     estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, 
                                     estimatePrecision = TRUE, estimateVar = TRUE, 
                                     estimateSubclone = estimateScPrevalence, estimateTransition = TRUE,
                                     estimateInitDist = TRUE, verbose = TRUE)

        iter <- hmmResults.cor$results$iter
        id <- names(hmmResults.cor$cna)[s]

        correctedResults <- correctIntegerCN(cn = hmmResults.cor$cna[[s]],
           segs = hmmResults.cor$results$segs[[s]],
           purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
           cellPrev = 1 - hmmResults.cor$results$sp[s, iter],
           maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = 0.03,
           gender = gender, chrs = chrs_with_data, correctHOMD = includeHOMD)

        hmmResults.cor$results$segs[[s]] <- correctedResults$segs
        hmmResults.cor$cna[[s]] <- correctedResults$cn

        segsS <- hmmResults.cor$results$segs[[s]]
        segsS <- segsS[segsS$chr %in% chrTrain, ]
        segAltInd <- which(segsS$event != "NEUT")
        maxBinLength = -Inf
        if (length(segAltInd) > 0){
            maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
        maxSegRD <- GRanges(seqnames=segsS$chr[segAltInd[maxInd]], 
                  ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
        hits <- findOverlaps(query=maxSegRD, subject=logR.data[[s]][valid, ])
        maxBinLength <- length(subjectHits(hits))
        }
        cnaS <- hmmResults.cor$cna[[s]]
        altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
        altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
        if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
            hmmResults.cor$results$n[s, iter] <- 1.0
        }

        mainName[[s]] <- c(mainName[[s]], paste0("Solution: ", counter, ", Sample: ", id, ", n: ", n[s], ", p: ", p[s], 
          ", log likelihood: ", signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4)))
        ## plot individual samples if only single-sample analysis
        if (numSamples == 1){
          #dir.create(paste0(outDir, "/", id))
          #outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_", "n", n[s], "-p", p[s])      
          outPlotFile <- paste0(outDir, "/", id, "_genomeWide_", "n", n[s], "-p", p[s])      
          plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType=plotFileType, 
                logR.column = "logR", call.column = "event",
                       plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, seqinfo=seqinfo, main=mainName[[s]][counter])
        }
        iter <- hmmResults.cor$results$iter
        results[[counter]] <- hmmResults.cor
        loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
        subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
        fracGenomeSub[counter, ] <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
        fracAltSub[counter, ] <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
        fracAltSubVect <- lapply(fracAltSub[counter, ], function(x){if (is.na(x)){0}else{x}})
        loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub[counter, ], digits=2), collapse=",")
        loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSubVect), digits=2), collapse=",")
        loglik[counter, "n_0"] <- paste0(n, collapse = ",")
        loglik[counter, "phi_0"] <- paste0(p, collapse = ",")
        loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
        loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")
        
        counter <- counter + 1
      }
    }

    ## remove solutions witn a likelihood of NA (i.e. no solution)
    loglik <- loglik[!is.na(loglik$loglik), ]
    ## get total time for all solutions ##
    elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
    message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

    ### SAVE R IMAGE ###
    outImage <- paste0(outDir, "/", id, ".RData")
    save.image(outImage)
    #save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))

    ### SORT SOLUTIONS BY DECREASING LIKELIHOOD ###
    if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
        if (numSamples > 1){
            fracInd <- which(rowSums(fracAltSub <= maxFracCNASubclone) == S &
                            rowSums(fracGenomeSub <= maxFracGenomeSubclone) == S)
        }else{
            fracInd <- which(fracAltSub <= maxFracCNASubclone &
                            fracGenomeSub <= maxFracGenomeSubclone)
        }
    	if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
    		ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
    	}else{ # otherwise just take largest likelihood
    		ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
    	}
    }else{#sort by likelihood only
      ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
    }

    ## print combined PDF of all solutions
    #new loop by order of solutions (ind)
    for (s in 1:numSamples){
      id <- names(results[[1]]$cna)[s]
      #outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_all_sols")
      outPlotFile <- paste0(outDir, "/", id, "_genomeWide_all_sols")
      for (i in 1:length(ind)) {
        hmmResults.cor <- results[[ind[i]]]
        turnDevOff <- FALSE
        turnDevOn <- FALSE
        if (i == 1){
        	turnDevOn <- TRUE
        }
        if (i == length(ind)){
        	turnDevOff <- TRUE
        }
        plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
                           logR.column = "logR", call.column = "event",
                           plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
                           seqinfo = seqinfo,
                           turnDevOn = turnDevOn, turnDevOff = turnDevOff, main=mainName[[s]][ind[i]])
      }
    }


    hmmResults.cor <- results[[ind[1]]]
    hmmResults.cor$results$loglik <- as.data.frame(loglik)
    hmmResults.cor$results$gender <- gender
    hmmResults.cor$results$chrYCov <- NULL
    hmmResults.cor$results$chrXMedian <- NULL
    hmmResults.cor$results$coverage <- coverage

    outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
                          results = hmmResults.cor$results, patientID = patientID, outDir=outDir)
    outFile <- paste0(outDir, "/", patientID, ".params.txt")
    outputParametersToFile(hmmResults.cor, file = outFile)

    ## plot solutions for all samples 
    counts <- list(list(counts=logR.data[[1]]))
    counts[[1]]$counts$reads <- counts[[1]]$counts$reads.tumor
    counts[[1]]$counts$ideal <- counts[[1]]$counts$ideal.tumor
    plotSolutions(hmmResults.cor, logR.data, chrs_with_data, outDir, counts, numSamples=numSamples,
               logR.column = "logR", call.column = "Corrected_Call", likModel = likModel,
               plotFileType=plotFileType, plotYLim=plotYLim, seqinfo = seqinfo,
               estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)

}


