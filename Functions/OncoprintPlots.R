library(ComplexHeatmap)
library(RColorBrewer)

oncoprint_alter_fun <- function(x, y, w, h, v) {
    n = sum(v)
    h = h*0.9
    if(n) {grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h,
            gp = gpar(fill = colors_oncoprint[names(which(v))], col = NA), just = "top")
} else { grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
    gp = gpar(fill = "grey90", col = NA))}}

effects_oncoprint <- list(
	Hotspot=c("Hotspot", "HOTSPOT"),
	Nonsense=c("STOP_GAINED", "Nonsense_Mutation", "stop_gained&splice_region_variant", "stop_gained", "nonsense",
		"stop_gained&disruptive_inframe_deletion", "stop_gained&disruptive_inframe_insertion"),
	Frameshift=c("FRAME_SHIFT", "Frame_Shift_Del", "Frame_Shift_Ins", "frameshift_variant", "frameshift_variant&stop_gained", 		"frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant",
		"frameshift_variant&splice_region_variant", "frame-shift", "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant",
		"frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant", "frameshift_variant&splice_donor_variant&intron_variant" ),
	Homozygous_deletion=c("Homozygous Deletion", "HOMDEL", "homozygous_deletion", "homozygous deletion", "Homozygous deletion"),
	Amplification=c("Amplification", "AMP", "amplification"),
	gain=c("gain","Gain"),
	loss=c("loss","Loss"),
	Promoter=c("upstream_gene_variant"),
	Missense=c("NON_SYNONYMOUS_CODING", "Missense_Mutation", "missense_variant", "missense_variant&splice_region_variant",
		"missense", "SPLICE_SITE_REGION|NON_SYNONYMOUS_CODING"),
	Inframe=c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_DELETION", "CODON_INSERTION", "In_Frame_Ins", "In_Frame_Del",
		"disruptive_inframe_deletion", "disruptive_inframe_insertion", "inframe_deletion", "inframe_insertion", "inframe_deletion&splice_region_variant",
		"splice_acceptor_variant&inframe_deletion&splice_region_variant&intron_variant", "disruptive_inframe_insertion&splice_region_variant",
		"splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant",
		"splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant", 
		"splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
		"disruptive_inframe_deletion&splice_region_variant", "SPLICE_SITE_REGION|CODON_CHANGE_PLUS_CODON_DELETION",
		"SPLICE_SITE_REGION|SPLICE_SITE_DONOR|CODON_CHANGE_PLUS_CODON_DELETION"),
	Splice_site=c("SPLICE_SITE_DONOR", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_REGION", "Splice_Site",
		"splice_donor_variant&intron_variant", "splice_acceptor_variant&intron_variant", "splicing", "splice_acceptor_variant&splice_donor_variant&intron_variant",
		"splice_donor_variant&splice_region_variant&intron_variant", "splice_acceptor_variant", "splice_donor_variant",
		"splice_acceptor_variant&splice_region_variant&intron_variant"),
	Other=c("STOP_LOST", "START_LOST", "START_GAINED", "UTR_5_PRIME", "start_lost", "Translation_Start_Site", "Nonstop_Mutation","stop_lost", "start_lost&splice_region_variant",
		"frameshift_variant&stop_lost", "stop_lost&splice_region_variant", "non_coding_exon_variant"),
	Fusion=c("fusion"),
	Exon_skipping=c("exon skipping"),
	rescued_hotspot=c("rescued hotspot"))

make_oncoprint_input <- function( muts=NULL, effects=effects_oncoprint, mutation_genes=NULL, sample_names=NULL,
    sample_name_col="TUMOR_SAMPLE", gene_name_col="GENE", eff_name_col ="EFFECT") {

        library(reshape2)
        if(any(!sample_name_col %in% colnames(muts), !gene_name_col %in% colnames(muts), !eff_name_col %in% colnames(muts))) {
            stop("At least one of the required columns not found in data\n")
        }

        effects <- melt(effects)
        mutation_genes <- unlist(mutation_genes) # included to avoid confusion with make_mutation_heatmap, which takes a list
        sample_names <- unlist(sample_names) # included to avoid confusion with make_mutation_heatmap, which takes a list
        
        y=do.call("cbind", lapply(sample_names, function(s) {
            smallmaf <- muts[which(muts[,sample_name_col]==s),,drop=F]
            x=tapply(smallmaf[,eff_name_col], smallmaf[,gene_name_col], function(eff) {
                toString(subset(effects, value %in% eff)$L1)
            })
            x[match(mutation_genes, names(x))]
        }))
        rownames(y) <- mutation_genes
        colnames(y) <- sample_names
        y
        
    }

# muttable needs 3 columns sample name, gene name, mutation effect. All other columns are ignored.
# make_oncoprint_input() looks for these 3 columns with names specified in the *_col arguments.
# the effects column can have any of the values specified in 'effects_oncoprint'. 
# The effects will be translated to the 10 simple consolidated terms


# the 3 arguments here are all required
###opi <- make_oncoprint_input(muttable, mutation_genes=unique(muttable$GENE), sample_names=unique(muttable$TUMOR_SAMPLE))

# Check out the ComplexHeatmap online documentation for endless customisation possibilities
###oncoPrint(opi, alter_fun = oncoprint_alter_fun, col = colors_oncoprint)
