# Creating matrices from seg files (PDX CN analysis step 1)
# Run using create_matrix_commands.sh
# Written by Dr. Gavin Ha, updated by Anna Hoge
# Ha Lab
# Fred Hutchinson Cancer Research Center
# February 27th, 2020

#settings used for study:
#segFileList = segment_files_list_full_set.txt
#normalList = PDX_normal_ids.txt
#regionAnnotFile = hg38_1mb_bins.txt
#annotType = chrPosn
#duplicateFunction = severe

library(foreach)
library(plyr)
library(doMC)
library(GenomicRanges)
library(stringr)
registerDoMC()
options(cores = 4)
args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE)

scriptDir <- args[1]
segFileList <- args[2]
normalList <- args[3]
regionAnnotFile <- args[4]
outRoot <- args[5]
annotType <- args[6]
duplicateFunction <- args[7]
imageFile <- args[8]


source(paste0(scriptDir, "/utils_updated.R"))
chrs <- c(1:22,"X")
## load normal list ##
normals <- read.delim(normalList, header=F, as.is=T)[,1]

## load the mouse annotations ##
regionAnnot <- read.delim(regionAnnotFile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
colnames(regionAnnot)[1:3] <- c("Chromosome","Start","End")
regionAnnot <- regionAnnot[regionAnnot$Chromosome %in% chrs, ]
## order the chromosomes using factors ##
regionAnnot$Chromosome <- factor(regionAnnot$Chromosome, levels=chrs)
regionAnnot <- regionAnnot[order(regionAnnot$Chromosome), ]
#regionAnnot$Band <- paste0(regionAnnot$Chromosome, regionAnnot$Band)
#regionAnnot <- with(regionAnnot, GRanges(seqnames = Chromosome, ranges = IRanges(start = Start, end = End), Band = Band))
regionAnnot <- as(regionAnnot, "GRanges")
save.image(imageFile)

segFiles <- read.delim(segFileList, header=F, stringsAsFactors=F, sep="\t")[,1]
mat <- list()
for (i in 1:length(segFiles)){
    file_name <- basename(segFiles[i])
	id <- substr(file_name, 1, nchar(file_name) - 4)
	segs <- read.delim(segFiles[i], header = TRUE, stringsAsFactors=F, sep="\t", skipNul = TRUE)
	segs <- segs[, c(2:4, 1, 5:ncol(segs))]
	segs[,1] <- factor(segs[,1], levels=c(1:22,"X","Y"))
	segs <- with(segs, GRanges(seqnames = chromosome, ranges = IRanges(start = start, end = end), Sample = ID, Num_Probes = N, Segment_Mean = log2.ratio))

    ###adjusted annot.type and fun for 1mb bin analysis
	mat.tmp <- createRegionMatrix(segs, regionAnnot, annot.type=annotType, fun = duplicateFunction)

	matOut <- list(); matOut[[id]] <- mat.tmp
	outFile <- paste0(outRoot, "_", id, ".rds")
	saveRDS(matOut, file=outFile)
	#write.table(mat, file = outFile, row.names=F, col.names=T, quote=F, sep="\t")
	mat[[id]] <- mat.tmp
}

outFile <- paste0(outRoot, "_allTissue.rds")
saveRDS(mat, file=outFile)
