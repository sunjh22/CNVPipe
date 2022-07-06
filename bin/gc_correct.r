#! /usr/bin/env Rscript
# This Rscript is written for performing cbs segmentation and MergeLevel on preliminary CNV profile
# Usage: Rscript 02-cbs.plot.r SRR.varbin.txt SRR SRR.mergeLevel.txt

suppressMessages(library(aCGH))
suppressMessages(library(DNAcopy))

args=commandArgs(T)

# gc normalization
lowess.gc <- function(jtkx, jtky) {
    jtklow <- lowess(jtkx, log(jtky), f=0.05)
	jtkz <- approx(jtklow$x, jtklow$y, jtkx)
	return(exp(log(jtky) - jtkz$y))
}

cbs.segment01 <- function(bad.bins, varbin.gc, varbin.data, sample.name, outputfile, 
                          alpha, nperm, undo.prune, min.width) {
	gc <- read.table(varbin.gc, header=T)
	bad <- read.table(bad.bins, header=F)

	chrom.numeric <- sub('chr', '', gc$bin.chrom)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisRatio <- read.table(varbin.data, header=F) 
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
	a <- thisRatio$bincount + 1
	thisRatio$ratio <- a / mean(a)
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)

	a <- quantile(gc$bin.length, 0.985)
	thisRatioNobad <- thisRatio[which(bad[, 1] == 0), ]

	# segment and smooth bin, here lowratio was logged by 2
	set.seed(25) 
	CNA.object <- CNA(log(thisRatioNobad$lowratio, base=2), thisRatioNobad$chrom, thisRatioNobad$chrompos, 
	                  data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="prune", 
	                                       undo.prune = undo.prune, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatioNobad), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}

	thisRatioNobad$seg.mean.LOWESS <- m[, 1]

	merge.ratio <- mergeLevels(thisRatioNobad$lowratio, thisRatioNobad$seg.mean.LOWESS)
	thisRatioNobad$mergeRatio <- merge.ratio$vecMerged

	if(FALSE){
		# plot the segmented bins
		chr <- thisRatioNobad$chrom
		chr.shift <- c(chr[-1], chr[length(chr)])
		vlines <- c(1, thisRatioNobad$abspos[which(chr != chr.shift) + 1], thisRatioNobad$abspos[nrow(thisRatioNobad)])
		hlines <- c(0.5, 1.0, 1.5, 2.0)
		chr.text <- c(1:22, "X", "Y")
		chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
		vlines.shift <- c(vlines[-1], 4*10^9)
		chr.at <- vlines + (vlines.shift - vlines) / 2
		chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]
		x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
		x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
		y.at <- c(0.005, 0.020, 0.100, 0.500, 2.000)
		y.labels <- c("0.005", "0.020", "0.100", "0.500", "2.000")


		postscript(paste(sample.name, ".wg.nobad.ps", sep=""), height=400, width=600)
		par(mar=c(5.1,4.1,4.1,4.1))
		plot(x=thisRatioNobad$abspos, y=thisRatioNobad$lowratio, log="y", main=sample.name, xaxt="n", xlab="Genome Position Gb", yaxt="n", ylab="Ratio", col="#CCCCCC", cex=0.5)
		axis(1, at=x.at, labels=x.labels)
		axis(2, at=y.at, labels=y.labels)
		lines(x=thisRatioNobad$abspos, y=thisRatioNobad$lowratio, col="#CCCCCC", cex=0.5)
		points(x=thisRatioNobad$abspos, y=thisRatioNobad$mergeRatio, col="#0000AA", cex=0.5)
		lines(x=thisRatioNobad$abspos, y=thisRatioNobad$mergeRatio, col="#0000AA", cex=0.5)
		abline(h=hlines)
		abline(v=vlines)
		mtext(chr.text, at = chr.at)
		dev.off()
	}
	
	# write.table(thisShort, sep="\t", file=paste(sample.name, ".varbinShort.txt", sep=""), quote=F, row.names=F)
	# write.table(thisRatioNobad, sep="\t", file=paste(sample.name, ".varbin.cbsSegment.mergeLevel.txt", sep=""), quote=F, row.names=F)
	write.table(thisRatioNobad, file=outputfile, sep="\t", quote=F, row.names=F)
}

cbs.segment01(bad.bins="/home/sunjh/data3/project/scDNA-CNV/refs/varbin-boundary/hg38_bwt2_12k_k76.badbins.txt", 
varbin.gc="/home/sunjh/data3/project/scDNA-CNV/refs/varbin-boundary/hg38_bwt2_12k_k76.gc.content.txt", 
varbin.data=args[1], sample.name=args[2], outputfile=args[3], alpha=0.0001, nperm=1000, undo.prune=0.05, min.width=5)
