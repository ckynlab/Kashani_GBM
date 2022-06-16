######################################################################################################
# plot_frequency function plots the frequency of copy number gains and losses (type="GL") and/or
# amplifications and homozygous deletions (type="AD"), between any number of groups (>=1). 
# Written to compare gene-level copy number alterations, but could be used as long as the input data
# is a data.frame of copy number states.
# If the number of groups is 2, optionally do a multi-Fisher's exact test to compare
# if there are significant differences between the two groups
#
# facets_tab : a data.frame of gene-level copy number states, see details below.
# sample_names: a list of sample names, where each element of the list consists of a character vector
# 	of sample names within each group. See details below.
# plot_file_prefix: prefix for output files
# width : width of plots (default: 7)
# height : height of plots (default: 2.5 per panel to be plotted)
# chr : vector of chromosomes to include (default: c(1:22,"X"))
# doFishers : should it do the Fisher's exact test between 2 groups? 
#	This would fail if the number of groups is not two. (default: F)
# type : one or both of "GL" and "AD", as a vector. 
# remove_genes_with_incomplete_data : should it remove genes with missing data? (default: T)
# anno_col : which columns or column names are annotations data and not sample data? 
#	These columns are included in the text report. (default: 1:5)
#
#
# Details
# facets_tab : a data frame of gene-level copy number data, encoded as 
#	amplification = 2, gain = 1, neutral = 0, loss = -1, deletion = -2
# 	Columns should consist of a minimum of chromosomes (named "chrom") and sample names
#	If sample_names is null, then sample_names will take the names from column 6 and treat all samples as a single group.
# 	Rows of facets_tab should be genes, sorted according to genomic positions
#	This function was originally written to plot output from FACETS, hence the first 5 columns would usually be
#	chrom, start, end, band, gene or some variant. Only chrom is used. The remaining annotation columns will be ignored.
# 	Alternatively, specify which columns are annotation columns with anno_col.
#	However, ignoring them is not recommended as the output table/s include these annotation columns.
# 
# sample_names : a list of sample names, with each element of the list containing a character vector of samples names
# 	The sample names should match the column names of facets_tab. If there are names in sample_names not found
# 	in the header of facets_tab, the function will fall over.
# 	The names of the list sample_names will be written as the titles of the plots.
# 	If sample_names is null, then sample_names will take the names from column 6 and treat all samples as a single group.
# 	If sample_names is a vector, all samples would be treated as a single group.
#######################################################################################################

#######################
plot_heatmapgistic <- function(facets_tab, plot_file, sample_names=NULL, col=c("blue","lightblue", NA, "darksalmon", "red"), zlim=c(-2,2)) {
  mm <- facets_tab
  if (is.null(sample_names)) { sample_names <- list(colnames(mm)[-c(1:5)]) }
  chrsep <- cumsum(rle(mm$chrom)$lengths)
  chrmid <- c(0,chrsep[-length(chrsep)]) + (rle(mm$chrom)$lengths/2)
  
  par(mfrow=c(length(sample_names),1), mar=c(8,.7*(max(sapply(sample_names,nchar))),1,2))
  lapply(sample_names, function(x, mm) {
    mm2 <- mm[,rev(x)];
    for (i in 1:ncol(mm2)) { mm2[,i] <- as.numeric(mm2[,i]) }
    image(as.matrix(mm2), col=col, xaxt='n', yaxt='n', zlim=zlim)
    box()
    for (i in (chrsep*2)-1) { abline(v=i/((max(chrsep)-1)*2), col="grey") }
    for (i in seq(-1, max(((2*(ncol(mm2)-1))+1),1), 2)) { abline(h=i/(2*(ncol(mm2)-1)), col="white", lwd=2)}
    
    axis(1,at=(chrmid/(max(chrsep)-1))[seq(1,length(chrmid),by=2)], 
         label=rle(mm$chrom)$values[seq(1,length(chrmid),by=2)], cex.axis=0.8, tick=F, line=-0.8, cex.axis=0.6)
    axis(1,at=(chrmid/(max(chrsep)-1))[seq(2,length(chrmid),by=2)], 
         label=rle(mm$chrom)$values[seq(2,length(chrmid),by=2)], cex.axis=0.8, tick=F, line=0, cex.axis=0.6)

    axis(2,at=seq(0,1,1/max((ncol(mm2)-1),1)), label=colnames(mm2), las=2, cex.axis=0.9, tick=F)
  }, mm)
  legend("bottom", inset=c(0,-0.15), legend=c("Homozygous deletion", "Loss", "Gain", "Amplification"),
         fill=col[c(1,2,4,5)], xpd=T, ncol=2, bty='n')
}
#######################

#######################
plot_frequency <- function(facets_tab, sample_names=NULL, plot_file_prefix=NULL, 
	width=12, height=3.2*(as.numeric(doFishers)+length(sample_names)), chr=c(1:22,"X"), ylim=100,
	doFishers=F, type=c("GL", "AD"), remove_genes_with_incomplete_data=T, anno_col=1:5) {

	mm <- facets_tab

	if (!"chrom" %in% colnames(mm)) { stop ("Need chrom column") }
	mm <- mm[which(mm$chrom %in% chr), , drop=F]
	if (is.null(anno_col)){ 
		stop("anno_col should be a vector of column index. If nothing else is available, at least chrom should be there")}

	if (is.null(sample_names) & !is.null(anno_col)) { sample_names <- list(SAMPLES=colnames(mm)[-anno_col]) }
	if (is.character(sample_names)) { sample_names <- list(SAMPLES=sample_names) }

	if (remove_genes_with_incomplete_data) {
		idx_to_rm <- which(unlist(apply(mm[,unlist(sample_names),drop=F], 1, function(x) { all(!is.na(x))})))
		mm <- mm[idx_to_rm,,drop=F]
	}
	if (is.null(plot_file_prefix)) { plot_file_prefix <- toString(names(sample_names)) }
	if (doFishers & length(sample_names)!=2) { stop("Can only do Fisher's with 2 groups") }

	double_barplot <- function(top, bottom, ylim, chromosomes, chr_levels, main, ylab) {
			par(mar=c(0,5,2,1))
			barplot(top, names.arg=NA, ylim=c(0, ylim), 
				xlim=c(0,length(top)), border=F, col="darkorchid4", space=0, las=2)

			lapply(cumsum(table(factor(chromosomes, levels= chr_levels))), 
				function(z) { abline(v=z, col="grey", lty=2)})
				text(cumsum(table(factor(chromosomes, levels= chr_levels)))-
				table(factor(chromosomes, levels= chr_levels))/2, ylim*0.9, chr_levels)
			mtext(main, side=3, line=0.5)
			box()

			par(mar=c(2,5,0,1))
			barplot(bottom, names.arg=NA, ylim=c(ylim,0), 
				xlim=c(0,length(bottom)), border=F, col="orange", space=0, las=2)
			lapply(cumsum(table(factor(chromosomes, levels= chr_levels))), 
				function(z) { abline(v=z, col="grey", lty=2)})
			box()
			mtext(ylab, side=2, line=3, at=c(0,0))
	}

	for (type0 in type) {

		if (type0 =="GL") { 
			high=c(1,2); low=c(-1,-2); 
			file <- paste(plot_file_prefix, "_GL", sep="") 
		} else if (type0 =="AD") { 
			high=2; low=-2; 
			file <- paste(plot_file_prefix, "_AD", sep="") 
		} else { stop (paste("Don't recognise ", type0, sep=""))}

		pdf(paste(file, ".pdf", sep=""), width=width, height=height)

		if (doFishers) { mfrow=6 } else { mfrow=2*length(sample_names) }
		par(mfrow=c(mfrow,1), mar=c(1,5,1,1))

		if (ylim=="auto") { 

			m <- unlist(lapply(names(sample_names), function(x, mm) {
				top=100*unlist(apply(mm[,sample_names[[x]]], 1, 
					function(y) { length(which(y %in% high))}))/length(sample_names[[x]])
				bottom=100*unlist(apply(mm[,sample_names[[x]]], 1, 
					function(y) { length(which(y %in% low))}))/length(sample_names[[x]])
				c(top, bottom)
			}, mm))

			ylim0=round(max(m), digits= -1)
			cat("ylim0:",ylim0,"\n")
			if (ylim0<max(m)-5) { ylim0=ylim0+10 }
			cat("ylim0:",ylim0,"\n")
		} else { ylim0 = ylim }
	
		lapply(names(sample_names), function(x, mm) {
				top=100*unlist(apply(mm[,sample_names[[x]]], 1, 
					function(y) { length(which(y %in% high))}))/length(sample_names[[x]])
				bottom=100*unlist(apply(mm[,sample_names[[x]]], 1, 
					function(y) { length(which(y %in% low))}))/length(sample_names[[x]])
				double_barplot(top, bottom, ylim0, mm$chrom, chr, 
					paste(x, " (n=", length(sample_names[[x]]), ")", sep=""), "Frequency (%)")

		}, mm)

		if (doFishers) {
			tab <- do.call("cbind", lapply(names(sample_names), function(x, mm) {
					y = unlist(apply(mm[,sample_names[[x]]], 1, function(y) { 
						length(which(y %in% high))}
					))
				cbind(y, length(sample_names[[x]])-y)
			}, mm))
			p <- unlist(apply(tab, 1, function(x) { fisher.test(matrix(x,nrow=2))$p }));
			p.adj <- p.adjust(p, "BH"); 


			p.adj[which(p.adj>0.05)] <- NA
			p.adj.log <- -log10(p.adj)
			p.adj.log[which(p.adj.log>10)] <- 10
			top= p.adj.log

			yy <- cbind(tab,p,p.adj, p.adj.log); 
			rownames(yy) <- mm$hgnc; 
			if (type0=="GL") { high="GAIN" } else if (type0=="AD") { high="AMP"}
			colnames(yy) <- c(paste("NumCases", names(sample_names)[[1]], high, sep="_"),
				paste("NumCases", names(sample_names)[[1]], "Not", high, sep="_"),
				paste("NumCases", names(sample_names)[[2]], high, sep="_"),
				paste("NumCases", names(sample_names)[[2]], "Not", high, sep="_"),
				paste("p", high, sep="_"), paste("p.adj", high, sep="_"),paste("p.adj.log", high, sep="_"))

			tab <- do.call("cbind", lapply(names(sample_names), function(x, mm) {
				y = unlist(apply(mm[,sample_names[[x]]], 1, function(y) { 
					length(which(y %in% low))}))
				cbind(y, length(sample_names[[x]])-y)
			}, mm))
			p <- unlist(apply(tab, 1, function(x) { fisher.test(matrix(x,nrow=2))$p }))
			p.adj <- p.adjust(p, "BH")
			p.adj[which(p.adj>0.05)] <- NA
			p.adj.log <- -log10(p.adj)
			p.adj.log[which(p.adj.log>10)] <- 10
			bottom <- p.adj.log

			zz <- cbind(tab,p,p.adj, p.adj.log); 
			rownames(zz) <- mm$hgnc; 
			if (type0=="GL") { low="LOSS" } else if (type0=="AD") { low="HOMDEL"}
			colnames(zz) <- c(paste("NumCases", names(sample_names)[[1]], low, sep="_"),
				paste("NumCases", names(sample_names)[[1]], "Not", low, sep="_"),
				paste("NumCases", names(sample_names)[[2]], low, sep="_"),
				paste("NumCases", names(sample_names)[[2]], "Not", low, sep="_"),
				paste("p", low, sep="_"), paste("p.adj", low, sep="_"),paste("p.adj.log", low, sep="_"))

			res <- cbind(mm[,anno_col,drop=F],yy, zz)
			write.table(res, file=paste(file, ".txt", sep=""), 
				sep="\t", col.names=NA, quote=F, na="")

			double_barplot(top, bottom, 10, mm$chrom, chr, "Fisher's exact test", "-log10(p-value)")
		}
		dev.off()
	}
}
#######################
