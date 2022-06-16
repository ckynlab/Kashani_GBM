##############################################################################################################################################################################
##############################################################################################################################################################################
# Merge and filter mutation files; prepare mutation overview
# TMB paired boxplots -- Figure 1A
# Initial TMB vs Age at excision - Figure S3C
#####################################

#############################
# Merge VCFs
rm(list=ls())

allvcfs <- list.files("/cluster/variant_annotation/")[grep("ensemble-annotated.recall.som_ad_ft.pon.hotspot.pass.sufam.eff.gene_ann.nsfp.parsed.tab",list.files("/cluster/variant_annotation/"))]
allvcfs <- as.character(paste0("/cluster/variant_annotation/",allvcfs))
allvcfs <- setNames(as.data.frame(allvcfs),c("VCF"))


vcfs <- data.frame()
v<-1

for(v in 1:nrow(allvcfs)){
  f=try(read.table(as.character(allvcfs[v,"VCF"])))
  
  if(class(f)=='try-error') {
    print(paste0(allvcfs[v,"VCF"]," has no mutations"))
    next
  }else{
    
    thisvcf <- read.table(as.character(allvcfs[v,"VCF"]), header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")
    
    this.sample <- unlist(strsplit(basename(as.character(allvcfs[v,"VCF"])),split="-"))[1]
    
    
    colnames(thisvcf)[grep("^GEN.",colnames(thisvcf))] <- gsub("GEN.","",colnames(thisvcf)[grep("^GEN.",colnames(thisvcf))])
    colnames(thisvcf)[grep("_N\\.",colnames(thisvcf))] <- paste0("GL.",unlist(lapply(strsplit(colnames(thisvcf)[grep("_N\\.",colnames(thisvcf))],split="\\."), "[", 2)))
    colnames(thisvcf)[grep("_T\\.",colnames(thisvcf))] <- paste0("T.",unlist(lapply(strsplit(colnames(thisvcf)[grep("_T\\.",colnames(thisvcf))],split="\\."), "[", 2)))
    
    thisvcf <- cbind(thisvcf, setNames(as.data.frame(rep(this.sample, times=nrow(thisvcf))),c("sample")))
    
    vcfs <- rbind(vcfs,thisvcf)
  }
}
saveRDS(vcfs,"/VCFs/mergedVCF_all.RDS")
write.table(vcfs, "/VCFs/mergedVCF_all.txt", sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
#############################

#############################
# Filter VCFs
rm(list = ls())

allmuts <- read.table("/VCFs/mergedVCF_all.txt", header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")

library(tidyr)
allmuts <- allmuts %>% separate(T.AD, c("T.REFcount","T.ALTcount"), "[;]", extra = "merge")
allmuts <- allmuts %>% separate(GL.AD, c("GL.REFcount","GL.ALTcount"), "[;]", extra = "merge")
allmuts <- allmuts[,c("sample","CHROM","POS","ID","REF","ALT","FILTER","T.REFcount","T.ALTcount","T.AF","T.DP",colnames(allmuts)[grep("ANN",colnames(allmuts))],"CancerGeneSets","CALLERS", colnames(allmuts)[grep("HOTSPOT",colnames(allmuts))],"dbNSFP_rs_dbSNP151", "GL.REFcount","GL.ALTcount","GL.AF","GL.DP")]
write.table(allmuts, "/VCFs/mergedVCF.txt", sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
saveRDS(allmuts,"/VCFs/mergedVCF.RDS")

# remove the interrogation_Absent variants (interrogated but 0 alt reads)
detected <- allmuts[-which(allmuts$FILTER == "interrogation_Absent"),]
write.table(detected, "/VCFs/mergedVCF_withoutInterrogationAbsent.txt", sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
saveRDS(detected,"/VCFs/mergedVCF_withoutInterrogationAbsent.RDS")

# remove the interrogation_Weak variants (interrogated but weak evidence)
highconfidence <- allmuts[-which(allmuts$FILTER %in% c("interrogation_Absent","interrogation_Weak")),]
write.table(highconfidence, "/VCFs/mergedVCF_withoutInterrogationAbsentandWeak.txt", sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
saveRDS(highconfidence,"/VCFs/mergedVCF_withoutInterrogationAbsentandWeak.RDS")

excludedvariants <- highconfidence[which(highconfidence$ANN.IMPACT %in% c("LOW","MODIFIER","NA")),]
highconfidenceExonic <- highconfidence[-which(highconfidence$ANN.IMPACT %in% c("LOW","MODIFIER","NA")),]
table(highconfidenceExonic$ANN.EFFECT)

# change factor of samples
samples <- c("PreGBM002T", "PostGBM002T", "PreGBM005T", "PostGBM005T", "PreGBM006T", "PostGBM006T", "PreGBM009T", "PostGBM009T", 
             "PreGBM011T", "PostGBM011T", "PreGBM014T", "PostGBM014T", "PreGBM017T", "PostGBM017T", "PreGBM018T", "PostGBM018T", 
             "PreGBM026T", "PostGBM026T", "PreGBM027T", "PostGBM027T", "PreGBM036T", "PostGBM036T", "PreGBM040T", "PostGBM040T", 
             "PreGBM045T", "PostGBM045T", "PreGBM046T", "PostGBM046T", "PreGBM052T", "PostGBM052T", "PreGBM057T", "PostGBM057T", 
             "PreGBM059T", "PostGBM059T", "PreGBM064T", "PostGBM064T", "PreGBM066T", "PostGBM066T", "PreGBM073T", "PostGBM073T", 
             "PreGBM075T", "PostGBM075T", "PreGBM076T", "PostGBM076T", "PreGBM077T", "PostGBM077T", "PreGBM082T", "PostGBM082T", 
             "PreGBM083T", "PostGBM083T", "PreGBM084T", "PostGBM084T", "PreGBM086T", "PostGBM086T", "PreGBM088T", "PostGBM088T")
highconfidenceExonic$sample <- factor(highconfidenceExonic$sample, levels = samples)

saveRDS(highconfidenceExonic,"/VCFs/mergedVCF_withoutInterrogationAbsentandWeak_onlyexonic.RDS")
write.table(highconfidenceExonic, "/VCFs/mergedVCF_withoutInterrogationAbsentandWeak_onlyexonic.txt", sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
#############################

#############################
# nucelotide matrix exonic & intronic
rm(list = ls())

detected <- readRDS("/VCFs/mergedVCF_withoutInterrogationAbsentandWeak.RDS")
sample.mut.ref <- as.data.frame(detected)

library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)

sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,
                                sample.id = "sample",
                                chr = "CHROM",
                                pos = "POS",
                                ref = "REF",
                                alt = "ALT",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(sigs.input,"/VCFs/sigs.input_withoutInterrogationAbsentandWeak.RDS")
#############################

#############################
# nucelotide matrix exonic
rm(list = ls())

highconfidenceExonic <- readRDS("/VCFs/mergedVCF_withoutInterrogationAbsentandWeak_onlyexonic.RDS")
sample.mut.ref <- as.data.frame(highconfidenceExonic)

library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)

sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,
                                sample.id = "sample",
                                chr = "CHROM",
                                pos = "POS",
                                ref = "REF",
                                alt = "ALT",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(sigs.input,"/VCFs/sigs.input_withoutInterrogationAbsentandWeak_onlyexonic.RDS")
#############################

#############################
# summary on mutation types
rm(list = ls())

capturesize <- 46.346140 # Capture size = EXONIC base pairs covered with the capture. 46'346'140 bp = size of sureselect v7; 

allmuts <- readRDS("/VCFs/mergedVCF.RDS")
detected <- readRDS("/VCFs/mergedVCF_withoutInterrogationAbsent.RDS")
highconfidence <- readRDS("/VCFs/mergedVCF_withoutInterrogationAbsentandWeak.RDS")
highconfidenceExonic <- readRDS("/VCFs/mergedVCF_withoutInterrogationAbsentandWeak_onlyexonic.RDS")

sigs.input.all <- readRDS("/VCFs/sigs.input_withoutInterrogationAbsentandWeak.RDS")
sigs.input.exonic <- readRDS("/VCFs/sigs.input_withoutInterrogationAbsentandWeak_onlyexonic.RDS")

indels <- c("conservative_inframe_deletion", "conservative_inframe_insertion", "disruptive_inframe_deletion", "disruptive_inframe_deletion&splice_region_variant", 
            "disruptive_inframe_insertion", "frameshift_variant", "frameshift_variant&splice_region_variant", "frameshift_variant&start_lost", "frameshift_variant&stop_gained")
nssnvs <- c("missense_variant", "missense_variant&splice_region_variant", "splice_region_variant", "splice_region_variant&stop_retained_variant", 
            "splice_region_variant&synonymous_variant", "start_lost", "start_lost&splice_region_variant", "stop_gained", "stop_gained&splice_region_variant", 
            "stop_lost", "stop_lost&conservative_inframe_deletion")
synsnvs <- c("synonymous_variant","stop_retained_variant")

samples <- c("PreGBM002T", "PostGBM002T", "PreGBM005T", "PostGBM005T", "PreGBM006T", "PostGBM006T", "PreGBM009T", "PostGBM009T", 
             "PreGBM011T", "PostGBM011T", "PreGBM014T", "PostGBM014T", "PreGBM017T", "PostGBM017T", "PreGBM018T", "PostGBM018T", 
             "PreGBM026T", "PostGBM026T", "PreGBM027T", "PostGBM027T", "PreGBM036T", "PostGBM036T", "PreGBM040T", "PostGBM040T", 
             "PreGBM045T", "PostGBM045T", "PreGBM046T", "PostGBM046T", "PreGBM052T", "PostGBM052T", "PreGBM057T", "PostGBM057T", 
             "PreGBM059T", "PostGBM059T", "PreGBM064T", "PostGBM064T", "PreGBM066T", "PostGBM066T", "PreGBM073T", "PostGBM073T", 
             "PreGBM075T", "PostGBM075T", "PreGBM076T", "PostGBM076T", "PreGBM077T", "PostGBM077T", "PreGBM082T", "PostGBM082T", 
             "PreGBM083T", "PostGBM083T", "PreGBM084T", "PostGBM084T", "PreGBM086T", "PostGBM086T", "PreGBM088T", "PostGBM088T")

overview <- setNames(data.frame(c(samples),
                                matrix(nrow = length(samples),ncol=22)),
                     c("sample","all","detected","highconfidence", "highconfidence-NS-Exonic",
                       "nsSNVs","INDELs","TMB","TMBperMb",
                       "C>A-exonic","C>G-exonic","C>T-exonic","T>A-exonic","T>C-exonic","T>G-exonic","C>Tperc-exonic",
                       "C>A-all","C>G-all","C>T-all","T>A-all","T>C-all","T>G-all"))

r<-1
for(r in 1:nrow(overview)){
  this.sample <- overview[r,"sample"]
  
  overview[r,"all"] <- nrow(allmuts[which(allmuts$sample == this.sample),])
  overview[r,"detected"] <- nrow(detected[which(detected$sample == this.sample),])
  overview[r,"highconfidence"] <- nrow(highconfidence[which(highconfidence$sample == this.sample),])
  overview[r,"highconfidence-NS-Exonic"] <- nrow(highconfidenceExonic[which(highconfidenceExonic$sample == this.sample),])
  
  # make a mutation summary with the high confidence exonic variants and use them for TMB calculation
  vcf <- highconfidenceExonic[which(highconfidenceExonic$sample == this.sample),]
  
  overview[r,"nsSNVs"] <- nrow(vcf[which(vcf$ANN.EFFECT %in% nssnvs),])
  overview[r,"INDELs"] <- nrow(vcf[which(vcf$ANN.EFFECT %in% indels),])
  overview[r,"TMB"] <- sum(overview[r,"nsSNVs"] ,overview[r,"INDELs"] )
  overview[r,"TMBperMb"] <- (sum(overview[r,"nsSNVs"] ,overview[r,"INDELs"])/capturesize)
  overview[r,"C>A-exonic"] <- sum(sigs.input.exonic[this.sample,grep("\\[C>A\\]",colnames(sigs.input.exonic))])
  overview[r,"C>G-exonic"] <- sum(sigs.input.exonic[this.sample,grep("\\[C>G\\]",colnames(sigs.input.exonic))])
  overview[r,"C>T-exonic"] <- sum(sigs.input.exonic[this.sample,grep("\\[C>T\\]",colnames(sigs.input.exonic))])
  overview[r,"T>A-exonic"] <- sum(sigs.input.exonic[this.sample,grep("\\[T>A\\]",colnames(sigs.input.exonic))])
  overview[r,"T>C-exonic"] <- sum(sigs.input.exonic[this.sample,grep("\\[T>C\\]",colnames(sigs.input.exonic))])
  overview[r,"T>G-exonic"] <- sum(sigs.input.exonic[this.sample,grep("\\[T>G\\]",colnames(sigs.input.exonic))])
  overview[r,"C>Tperc-exonic"] <- sum(sigs.input.exonic[this.sample,grep("\\[C>T\\]",colnames(sigs.input.exonic))])/sum(overview[r,"synSNVs"],overview[r,"nsSNVs"])
  overview[r,"C>A-all"] <- sum(sigs.input.all[this.sample,grep("\\[C>A\\]",colnames(sigs.input.all))])
  overview[r,"C>G-all"] <- sum(sigs.input.all[this.sample,grep("\\[C>G\\]",colnames(sigs.input.all))])
  overview[r,"C>T-all"] <- sum(sigs.input.all[this.sample,grep("\\[C>T\\]",colnames(sigs.input.all))])
  overview[r,"T>A-all"] <- sum(sigs.input.all[this.sample,grep("\\[T>A\\]",colnames(sigs.input.all))])
  overview[r,"T>C-all"] <- sum(sigs.input.all[this.sample,grep("\\[T>C\\]",colnames(sigs.input.all))])
  overview[r,"T>G-all"] <- sum(sigs.input.all[this.sample,grep("\\[T>G\\]",colnames(sigs.input.all))])
  
  fullsamplename <- gsub("PreGBM","Pre_GBM",this.sample)
  fullsamplename <- gsub("PostGBM","Post_GBM",fullsamplename)
  fullsamplename <- gsub("T","_T",fullsamplename,ignore.case = FALSE)
}
# P046-Pre is a recurrence sample and P046-Post got excluded
overview <- overview[-which(overview$sample == "PostGBM046T"),]
overview[which(overview$sample == "PreGBM046T"),"sample"] <- "PostGBM046T"

write.table(overview, "/VCFs/Mutation-type-Summary.txt", sep = "\t", quote=FALSE, row.names=TRUE,col.names = TRUE)
#############################

#############################
# TMB paired boxplots - Figure 1A
rm(list = ls())

data <- read.table(file = "/VCFs/Mutation-type-Summary.txt", check.names = FALSE, header=TRUE,sep = "\t",na.string=".", stringsAsFactors=F)

data$patient <- gsub("Pre|Post|T","",data$sample)
data$timepoint <- unlist(lapply(strsplit(data$sample,split="GBM"), "[", 1))
data$timepoint <- gsub("Pre","Initial",data$timepoint)
data$timepoint <- gsub("Post","Recurrent",data$timepoint)

data$log10TMBperMb <- log10(data$TMBperMb)

# remove GBM046 because it has no pre-treatment sample
data <- dplyr::filter(data,patient!="GBM046")

### Packages
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(pheatmap)
library(gridExtra)

tmbComp <-ggpaired(data, x = "timepoint", y = "log10TMBperMb", id="patient", alpha=0.1,
                   color = "timepoint", line.color = "gray", line.size = 0.5,
                   palette = c("deepskyblue3","deepskyblue4"))+
  scale_y_continuous(limits=c(-0.2,2.2),breaks=seq(-0, 2, by = 0.5))+
  theme_classic()+
  xlab("")+
  theme(plot.margin = unit(c(1,1,1,1), "lines"))+ 
  ylab("log10(Mutations/Mb)")+
  ggtitle("") +
  theme(legend.position="none")+
  theme(plot.margin = unit(c(1,1,1,1), "lines"))+ 
  theme(title=element_text(size=11),
        axis.text.x=element_text(size=11), axis.title.x=element_text(size=13), 
        axis.text.y =element_text(size=11), axis.title.y=element_text(size=13))+
  stat_compare_means(paired = TRUE,method = "wilcox.test",label.y=1.9,size=4, label.x=1.1)
plot(tmbComp)
#############################

#############################
# Intial TMB vs Age at excision - Figure S3C
rm(list = ls())

pre <- data[which(data$timepoint == "Initial"),]
summary(pre$TMBperMb)

post <- data[which(data$timepoint == "Recurrent"),]
summary(post$TMBperMb)

library(openxlsx)
clinical <- read.xlsx("/Tables/SupplTable1noBlockInfo.xlsx")
initials <- data[which(data$timepoint == "Initial"),]
initials$Patient.ID <- gsub("Pre","",initials$sample)
initials$Patient.ID <- gsub("T","",initials$Patient.ID)

clinical <- dplyr::left_join(initials,clinical)

a <- ggplot(data=clinical, aes(x = as.numeric(clinical$Age.at.Initial.Excision), y = as.numeric(clinical$TMBperMb)),na.rm = FALSE) + 
  geom_point(data=clinical, aes(x = as.numeric(clinical$Age.at.Initial.Excision), y = as.numeric(clinical$TMBperMb)), alpha=0.6,size = 2.5, color="navy") +
  geom_smooth(method ="lm", color="black", size=0.5) +
  stat_cor(label.x = 45,label.y =20, size=5, method="spearman") + ## "pearson" == default method
  xlab("Age at excision") +
  ylab("TMB [mut/Mb]") +
  theme_bw()+
  theme(aspect.ratio = 1)+
  theme(legend.title=element_blank(),strip.text.x = element_text(size = 12),axis.text.x=element_text(size=12), axis.text.y =element_text(size=12), 
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))
print(a)
##############################################################################################################################################################################
##############################################################################################################################################################################


##############################################################################################################################################################################
##############################################################################################################################################################################
# Signature barplots -- Figure 1B
# Paired boxplots significant signatures, Absolute Number of Mutations -- Figure 1C
# Paired boxplots remaining signatures, Absolute Number of Mutations -- Figure S3A
# Paired boxplots all signatures, Percent Contribution -- plot to confirm that significance holds if % signature contribution instead of absolute mutation counts are used
#####################################
rm(list=ls())

##########################
# Run DeconstructSigs
out.dir <- "/DeconstructSigs/"
sigs.input <- readRDS("/VCFs/sigs.input_withoutInterrogationAbsentandWeak.RDS")

# make sample names shorter
rownames(sigs.input)[grep("Pre",rownames(sigs.input))] <- gsub("T$","-Pre",rownames(sigs.input)[grep("Pre",rownames(sigs.input))])
rownames(sigs.input)[grep("Pre",rownames(sigs.input))] <- gsub("^PreGBM","P",rownames(sigs.input)[grep("Pre",rownames(sigs.input))])
rownames(sigs.input)[grep("Post",rownames(sigs.input))] <- gsub("T$","-Post",rownames(sigs.input)[grep("Post",rownames(sigs.input))])
rownames(sigs.input)[grep("Post",rownames(sigs.input))] <- gsub("^PostGBM","P",rownames(sigs.input)[grep("Post",rownames(sigs.input))])

samples <- rownames(sigs.input)

# relevel samples and change order of rows
levelorder <- c("P002-Pre", "P002-Post", "P005-Pre", "P005-Post", "P006-Pre", "P006-Post", "P009-Pre", "P009-Post", "P011-Pre", "P011-Post", 
                "P014-Pre", "P014-Post", "P017-Pre", "P017-Post", "P018-Pre", "P018-Post", "P026-Pre", "P026-Post", "P027-Pre", "P027-Post", 
                "P036-Pre", "P036-Post", "P040-Pre", "P040-Post", "P045-Pre", "P045-Post", "P046-Pre", "P046-Post", "P052-Pre", "P052-Post", 
                "P057-Pre", "P057-Post", "P059-Pre", "P059-Post", "P064-Pre", "P064-Post", "P066-Pre", "P066-Post", "P073-Pre", "P073-Post", 
                "P075-Pre", "P075-Post", "P076-Pre", "P076-Post", "P077-Pre", "P077-Post", "P082-Pre", "P082-Post", "P083-Pre", "P083-Post", 
                "P084-Pre", "P084-Post", "P086-Pre", "P086-Post", "P088-Pre", "P088-Post")
sigs.input <- sigs.input[levelorder,]
levelorder <- factor(levelorder, levels = levelorder)
factor(rownames(sigs.input), levels = factor(levelorder))

# Run deconstructSigs
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
source("/Scripts/Signatures_AddAetiologyToRownames_DefineColorsForEachSignature.R")

# make GBM specific reference set
load("/DeconstructSigs/signatures.exome.cosmic.v3.may2019.rda")
signatures.exome.cosmic.v3.may2019.GBMselection <- signatures.exome.cosmic.v3.may2019[c(1,3,5,11,15,19,20,35,45),]
save(signatures.exome.cosmic.v3.may2019.GBMselection, file = "/DeconstructSigs/signatures.exome.cosmic.v3.may2019.GBMselection.rda")

for (s in samples){
  # Determine the signatures contributing to an already normalized sample
  sig = whichSignatures( tumor.ref = sigs.input,
                         signatures.ref = signatures.exome.cosmic.v3.may2019.GBMselection,
                         sample.id = s,
                         signature.cutoff = 0.05,
                         contexts.needed = TRUE)
  saveRDS(sig, file=file.path(out.dir, paste(s,"_sig.RDS",sep="")))
  
  # Plot output
  pdf(file=file.path(out.dir, paste(s,"_signatures.pdf",sep="")), width=12, height=6)
  plotSignatures(sig, sub = s)
  makePie(sig)
  dev.off()
}
### Summary
fns <- list.files(out.dir, pattern="_sig.RDS$")
output.fn <- file.path(out.dir,"Signatures_summary.pdf")

sigs <- c()
for (f in fns){
  f.path <- file.path(out.dir, f)
  f.rds <- readRDS(f.path)
  sigs <- rbind(sigs, f.rds$weights)
}

# change order of rows
sigs <- sigs[rownames(sigs.input),]
saveRDS(sigs,"/DeconstructSigs/ALL_DeconstructSigs.RDS")
#############################

#############################
#  Make consensus signatures
rm(list=ls())

# DeconstructSigs output
sigs <- readRDS("/DeconstructSigs/ALL_DeconstructSigs.RDS")
# P046-Pre is a recurrence sample and P046-Post got excluded
post <- which(rownames(sigs) == "P046-Post")
sigs <- sigs[-post,]
pre <- which(rownames(sigs) == "P046-Pre")
rownames(sigs)[pre] <- "P046-Post"
saveRDS(sigs,"/ConsensusSignatures/DeconstructSigsSignatures.RDS")

# Mutational patterns output. This output was generated with the SOCIBP pipeline using the same reference matrix as above for DeconstructSigs
mutpat <- read.table("/MutationalPatterns/MutationalPatterns-GBMsigs.signatures.exome.cosmic.v3.may2019.RelativeContribution.txt", header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")
# harmonize deconstruct sigs and mutational patterns outputs
rownames(mutpat) <- mutpat[,"SampleName"]
mutpat <- mutpat[,-which(colnames(mutpat)=="SampleName")]
# make sample names shorter
rownames(mutpat)[grep("Pre",rownames(mutpat))] <- gsub("_T$","-Pre",rownames(mutpat)[grep("Pre",rownames(mutpat))])
rownames(mutpat)[grep("Pre",rownames(mutpat))] <- gsub("^Pre_GBM","P",rownames(mutpat)[grep("Pre",rownames(mutpat))])
rownames(mutpat)[grep("Post",rownames(mutpat))] <- gsub("_T$","-Post",rownames(mutpat)[grep("Post",rownames(mutpat))])
rownames(mutpat)[grep("Post",rownames(mutpat))] <- gsub("^Post_GBM","P",rownames(mutpat)[grep("Post",rownames(mutpat))])
# P046-Pre is a recurrence sample and P046-Post got excluded
post <- which(rownames(mutpat) == "P046-Post")
mutpat <- mutpat[-post,]
pre <- which(rownames(mutpat) == "P046-Pre")
rownames(mutpat)[pre] <- "P046-Post"
# reorder
mutpat <- mutpat[rownames(sigs),]
# set any value <5% to 0; already true for deconstruct sigs output
for(r in 1:nrow(mutpat)){
  for (c in 1:ncol(mutpat)) {
    if(mutpat[r,c]<0.05){
      mutpat[r,c] <- 0
    }
  }
}
# rename the columns to add aetiology of the signatures
source("/Scripts/Signatures_AddAetiologyToRownames_DefineColorsForEachSignature.R")
sigs <- renameSigs(sigs)
mutpat <- renameSigs(mutpat)
saveRDS(mutpat,"/ConsensusSignatures/MutationalPatternsSignatures.RDS")

### Consensus
# df with mean between deconstruct sigs and mutational patterns signatures
meandf <- matrix(nrow = nrow(sigs),ncol = ncol(sigs))
rownames(meandf) <- rownames(sigs)
colnames(meandf) <- colnames(sigs)
for(r in 1:nrow(meandf)){
  for (c in 1:ncol(meandf)) {
    meandf[r,c] <- mean(c(sigs[r,c],mutpat[r,c]))
  }
}
saveRDS(meandf,"/ConsensusSignatures/ConsensusSignatures.RDS")
write.table(meandf, "/ConsensusSignatures/ConsensusSignatures.txt", sep = "\t", quote=FALSE, row.names=TRUE,col.names = TRUE)
#############################


#############################
#  make plots
rm(list=ls())
dev.off()

source("/Scripts/Signatures_AddAetiologyToRownames_DefineColorsForEachSignature.R")
library(ggplot2)
library(gridExtra)
library(ggpubr)
require(reshape2)

meandf <- readRDS("/ConsensusSignatures/ConsensusSignatures.RDS")

# make signature names shorter
colnames(meandf)[which(colnames(meandf)=="SBS1-Deamination of 5-methylcytosine")] <- "SBS1 Deamination of 5MeC"
colnames(meandf)[which(colnames(meandf)=="SBS3-Defective HR repair (BRCA1/2 mut)")] <- "SBS3 Defective HR repair"
colnames(meandf)[which(colnames(meandf)=="SBS5-Unknown - clock-like")] <- "SBS5 Unknown (clock-like)"
colnames(meandf)[which(colnames(meandf)=="SBS8-Unknown")] <- "SBS8 Unknown"
colnames(meandf)[which(colnames(meandf)=="SBS11-Temozolomide treatment")] <- "SBS11 Temozolomide treatment"
colnames(meandf)[which(colnames(meandf)=="SBS15-Defective MMR")] <- "SBS15 Defective MMR"
colnames(meandf)[which(colnames(meandf)=="SBS16-Unknown")] <- "SBS16 Unknown"
colnames(meandf)[which(colnames(meandf)=="SBS30-Defective BER (inactivating mut. in NTHL1)")] <- "SBS30 Defective BER"
colnames(meandf)[which(colnames(meandf)=="SBS40-Unknown")] <- "SBS40 Unknown"

# use sigs.input from deconstruct sigs to get the mutation counts for each sample
sigs.input <- readRDS("/VCFs/sigs.input_withoutInterrogationAbsentandWeak.RDS")

# make sample names shorter
rownames(sigs.input)[grep("Pre",rownames(sigs.input))] <- gsub("T$","-Pre",rownames(sigs.input)[grep("Pre",rownames(sigs.input))])
rownames(sigs.input)[grep("Pre",rownames(sigs.input))] <- gsub("^PreGBM","P",rownames(sigs.input)[grep("Pre",rownames(sigs.input))])
rownames(sigs.input)[grep("Post",rownames(sigs.input))] <- gsub("T$","-Post",rownames(sigs.input)[grep("Post",rownames(sigs.input))])
rownames(sigs.input)[grep("Post",rownames(sigs.input))] <- gsub("^PostGBM","P",rownames(sigs.input)[grep("Post",rownames(sigs.input))])
# P046-Pre is a recurrence sample and P046-Post got excluded
post <- which(rownames(sigs.input) == "P046-Post")
sigs.input <- sigs.input[-post,]
pre <- which(rownames(sigs.input) == "P046-Pre")
rownames(sigs.input)[pre] <- "P046-Post"
# rename to "Initial" and "Recurrent"
rownames(sigs.input)[grep("Pre",rownames(sigs.input))] <- gsub("Pre","Initial",rownames(sigs.input)[grep("Pre",rownames(sigs.input))])
rownames(sigs.input)[grep("Post",rownames(sigs.input))] <- gsub("Post","Recurrent",rownames(sigs.input)[grep("Post",rownames(sigs.input))])
rownames(meandf)[grep("Pre",rownames(meandf))] <- gsub("Pre","Initial",rownames(meandf)[grep("Pre",rownames(meandf))])
rownames(meandf)[grep("Post",rownames(meandf))] <- gsub("Post","Recurrent",rownames(meandf)[grep("Post",rownames(meandf))])

rwsms <- rowSums(sigs.input)
rwsms <- setNames(cbind(names(rwsms),as.data.frame(rwsms)), c("SampleName","N mutations"))
rwsms <- rwsms[rownames(meandf),]

meandf <- as.data.frame(meandf)

# remove dashes in the sample names for the plotting step
rwsms$SampleName <- gsub("-"," ",rwsms$SampleName)
rownames(meandf) <- gsub("-"," ",rownames(meandf))

######
# make Signature barplots - Figure 1B (functions defined in Signatures_AddAetiologyToRownames_DefineColorsForEachSignature.R)
# define the plotting-colors for the signatures
color <- defineColorsGBMsubset(meandf)
# percent instead of fraction
percent <- meandf*100
pmean <- makeSigBarplot(percent, thislabel = "", mutationscounts=rwsms)
pmean

saveRDS(meandf,"/ConsensusSignatures/ConsensusSignatures-renamed.RDS")
saveRDS(rwsms,"/ConsensusSignatures/MutationCountForSignatures-renamed.RDS")
#############################


#############################
# paired boxplot signatures
rm(list=ls())
dev.off()

source("/Scripts/Signatures_AddAetiologyToRownames_DefineColorsForEachSignature.R")
library(ggplot2)
library(gridExtra)
library(ggpubr)
require(reshape2)

meandf <- readRDS("/ConsensusSignatures/ConsensusSignatures-renamed.RDS")
rwsms <- readRDS("/ConsensusSignatures/MutationCountForSignatures-renamed.RDS")

meandf$SampleName <- rownames(meandf)
combined <- dplyr::left_join(meandf,rwsms)
rownames(combined) <- combined$SampleName

absolute <- combined
absolute[,-which(colnames(absolute) %in% c("SampleName","N mutations"))] <- round(absolute[,-which(colnames(absolute) %in% c("SampleName","N mutations"))]*absolute$`N mutations`)

prepareforplotting <- function(df){
  # prepare for plotting
  df <- df[-which(colnames(df) == "N mutations")]
  df <- melt(df, id.vars = c("SampleName"))
  df$patient <- unlist(lapply(strsplit(as.character(df$SampleName),split=" "), "[", 1))
  # remove P046 for this analysis because it's not a matched case
  df <- df[-which(df$patient == "P046"),]
  df$timepoint <- unlist(lapply(strsplit(as.character(df$SampleName),split=" "), "[", 2))
  df$variable <- as.character(df$variable)
  return(df)
}

fractions <- prepareforplotting(combined)
saveRDS(fractions,"/ConsensusSignatures/Signatures-Data-ForPlotting-Fractions.RDS")

absolute <- prepareforplotting(absolute)
saveRDS(absolute,"/ConsensusSignatures/Signatures-Data-ForPlotting-Absolute.RDS")

# choose if you want to plot the absolute count of mutations attributed to each signature or the signature fractions 
absoluteCounts <- TRUE

if(absoluteCounts == TRUE){
  combined <- absolute
}else{
  combined <- fractions
  combined$value <- combined$value*100
}

### run this code block manually and comment/uncomment the relevant bits
# To plot selected signatures use 'meandf[,selectedSignatures]' and 'combined <- dplyr::filter(combined, variable %in% selectedSignatures)';
# to plot all other signatures use 'meandf[,othersignatures]' and 'combined <- dplyr::filter(combined, variable %in% othersignatures)' and change the y-axis of the plot
allsignatures <- colnames(meandf)[-which(colnames(meandf)=="SampleName")]
# subset to selected signatures if required
selectedSignatures <- c("SBS1 Deamination of 5MeC","SBS11 Temozolomide treatment")
meandf <- meandf[,selectedSignatures]
combined <- dplyr::filter(combined, variable %in% selectedSignatures)
# subset to 'other' signatures if required
# othersignatures <- allsignatures[-which(allsignatures %in% selectedSignatures)]
# meandf <- meandf[,othersignatures]
# combined <- dplyr::filter(combined, variable %in% othersignatures)
###

# colors for plotting
color <- defineColorsGBMsubset(meandf)

# Compute statistical test
stat.test <- compare_means(
  value ~ timepoint, data = combined, 
  group.by = "variable", method = "wilcox.test", paired = TRUE
)

if(absoluteCounts==TRUE){
  ### transform the data -- here use log10(x+ 1) to avoid zeros to become -INF
  # Note: the stats are still based on the actual values, only the values for plotting are transformed
  combined$value <- combined$value + 1
  combined$value <- log10(combined$value)
}

# Plot without stats
if(absoluteCounts==TRUE){
  Comp <-ggpaired(combined, x = "timepoint", y = "value", id="patient",
                  line.color = "gray", line.size = 0.2,
                  color = "variable",palette = color)+
    scale_y_continuous(limits = c(0,round(max(combined$value), digits = 0)))+#, breaks=seq(0,round(max(combined$value)),by=1)))+
    theme_classic()+
    xlab("")+
    ylab("log10(N mutations+1)")+
    ggtitle("") +
    theme(legend.position="none")+
    theme(plot.margin = unit(c(1,1,1,1), "lines"))+ 
    facet_wrap(~ variable, ncol=4, scales = "free")+
    theme(title=element_text(size=10),strip.text = element_text(size=9),
          axis.text.x=element_text(size=10), axis.title.x=element_text(size=10), 
          axis.text.y =element_text(size=10), axis.title.y=element_text(size=10))
  plot(Comp)
}else{
  Comp <-ggpaired(combined, x = "timepoint", y = "value", id="patient",
                  line.color = "gray", line.size = 0.2,
                  color = "variable",palette = color)+
    scale_y_continuous(limits = c(0,100), breaks=seq(0,100,by=20))+
    theme_classic()+
    xlab("")+
    ylab("% signature contribution")+
    ggtitle("") +
    theme(legend.position="none")+
    theme(plot.margin = unit(c(1,1,1,1), "lines"))+ 
    facet_wrap(~ variable, ncol=5, scales = "free")+
    theme(title=element_text(size=10),strip.text = element_text(size=8.3),
          axis.text.x=element_text(size=10), axis.title.x=element_text(size=10), 
          axis.text.y =element_text(size=10), axis.title.y=element_text(size=10))
  plot(Comp)
}

# Add stats
stat.test$p <- round(stat.test$p,digits = 5)
Comp + stat_pvalue_manual(stat.test, label = "Wilcoxon, p = {p}", remove.bracket = TRUE, 
                          y.position =3.7, hjust=0.8,size=3.5)
##############################################################################################################################################################################
##############################################################################################################################################################################


##############################################################################################################################################################################
##############################################################################################################################################################################
# Mutation heatmaps
# MMR mutations -- Figure S3B
# Recurrence-only variants -- Figure S6
# GBM driver -- Mutations section of Figure 1D
#####################################

#########################
# MMR mutations -- Figure S3B
rm(list=ls())
dev.off()

library(ComplexHeatmap)
library(RColorBrewer)
source("/Scripts/oncoprint_examples.R")

# Prepare mutation data
muts <- readRDS("/VCFs/mergedVCF_withoutInterrogationAbsentandWeak_onlyexonic.RDS")

# remove PostGBM046T and rename PreGBM046T to PostGBM046T
pre <- which(muts$sample == "PreGBM046T")
post <- which(muts$sample == "PostGBM046T")
muts[pre,"sample"] <- "PostGBM046T"
muts <- muts[-post,]
samples <- c("PreGBM002T", "PostGBM002T", "PreGBM005T", "PostGBM005T", "PreGBM006T", "PostGBM006T", "PreGBM009T", "PostGBM009T", 
             "PreGBM011T", "PostGBM011T", "PreGBM014T", "PostGBM014T", "PreGBM017T", "PostGBM017T", "PreGBM018T", "PostGBM018T", 
             "PreGBM026T", "PostGBM026T", "PreGBM027T", "PostGBM027T", "PreGBM036T", "PostGBM036T", "PreGBM040T", "PostGBM040T", 
             "PreGBM045T", "PostGBM045T", "PostGBM046T", "PreGBM052T", "PostGBM052T", "PreGBM057T", "PostGBM057T", 
             "PreGBM059T", "PostGBM059T", "PreGBM064T", "PostGBM064T", "PreGBM066T", "PostGBM066T", "PreGBM073T", "PostGBM073T", 
             "PreGBM075T", "PostGBM075T", "PreGBM076T", "PostGBM076T", "PreGBM077T", "PostGBM077T", "PreGBM082T", "PostGBM082T", 
             "PreGBM083T", "PostGBM083T", "PreGBM084T", "PostGBM084T", "PreGBM086T", "PostGBM086T", "PreGBM088T", "PostGBM088T")#"PreGBM064T", 
samples <- factor(samples, levels = samples)
# change the effect to 'hotspot' for all hotspot variants
muts[which(muts$HOTSPOT == "true"),"ANN.EFFECT"] <- "Hotspot"
muttable <- setNames(muts[,c("sample", "ANN.GENE", "ANN.EFFECT")],c("TUMOR_SAMPLE", "GENE", "EFFECT"))

# Select genes to plot
mmrgenes <- c("EXO1", "LIG1", "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PCNA", "PMS2", "POLD1", "POLD2", "POLD3", "POLD4", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5", "RPA1", "RPA2", "RPA3", "RPA4", "SSBP1")
muttable <- muttable[which(muttable$GENE %in% mmrgenes),]
# prepare matrix for heatmap

opi <- make_oncoprint_input(muttable, mutation_genes=unique(muttable$GENE), sample_names=unique(muttable$TUMOR_SAMPLE))
# Add missing samples (samples with no mutation in any of the selected genes)
missingsamples <- as.character(samples[-which(samples %in% colnames(opi))])
missingsamples <- setNames(data.frame(matrix(nrow = nrow(opi),ncol=length(missingsamples))),c(missingsamples))
opi <- cbind(opi, missingsamples)
# reorder
opi <- opi[,levels(samples)]
opi <- as.matrix(opi)
# add column split by patient
columnsplit <- gsub("PreGBM0|PostGBM0|T","",colnames(opi))
colnames(opi) <- unlist(lapply(strsplit(colnames(opi),split="GBM"), "[", 1))
colnames(opi)[which(colnames(opi)=="Pre")] <- "Initial"
colnames(opi)[which(colnames(opi)=="Post")] <- "Recurrent"

# plot
anno_width = unit(3, "cm")
fontsizeforall <- 9
colors_oncoprint <- setNames(c("#88CCEE","#CC6677","#fcc5c0","forestgreen","#332288","#fee391","#9e9ac8"),
                             c("Frameshift","Hotspot","Hotspot, Missense","Inframe","Missense","Nonsense","Nonsense, Missense"))
colors_oncoprint <- colors_oncoprint[which(names(colors_oncoprint) %in% names(table(opi)))]
mutations_legend = Legend(labels = names(colors_oncoprint), title = "", legend_gp = gpar(fill = colors_oncoprint))
mutations_plot <- Heatmap(opi, col = colors_oncoprint, rect_gp = gpar(col = "grey"),
                          show_heatmap_legend = FALSE,
                          na_col = "white",name=" ",column_names_gp = gpar(fontsize = fontsizeforall),
                          row_names_gp = gpar(fontsize = fontsizeforall),column_split = columnsplit, column_title_gp = gpar(fontsize=fontsizeforall))

draw(mutations_plot, padding = unit(c(2, 2, 20, 2), "mm"), heatmap_legend_list = mutations_legend, merge_legend = FALSE, heatmap_legend_side = "right",annotation_legend_side = "bottom")


# TMB track
mutationcounts <- read.table(file = "/VCFs/Mutation-type-Summary.txt", check.names = FALSE, header=TRUE,sep = "\t",na.string=".", stringsAsFactors=F)
mutationcounts$sample
mutationcounts$patient <- gsub("Pre|Post|T","",mutationcounts$sample)
mutationcounts$timepoint <- unlist(lapply(strsplit(mutationcounts$sample,split="GBM"), "[", 1))
mutationcounts$timepoint <- gsub("Pre","Initial",mutationcounts$timepoint)
mutationcounts$timepoint <- gsub("Post","Recurrent",mutationcounts$timepoint)
mutationcounts$sample <- paste(mutationcounts$patient,mutationcounts$timepoint,sep=" ")
mutationcounts$TMBcategory <- "NA"
mutationcounts$TMBcategory[mutationcounts$TMBperMb < 6] <- "<6 mutations/Mb"
mutationcounts$TMBcategory[mutationcounts$TMBperMb >=6 & mutationcounts$TMBperMb <= 20] <- "6-20 mutations/Mb"
mutationcounts$TMBcategory[mutationcounts$TMBperMb >= 21 & mutationcounts$TMBperMb <= 50] <- "21-50 mutations/Mb"
mutationcounts$TMBcategory[mutationcounts$TMBperMb > 50] <- ">50 mutations/Mb"
# signature track
meandf <- readRDS("/ConsensusSignatures/ConsensusSignatures-renamed.RDS")
meandf$sample <- gsub("^P","GBM",rownames(meandf))

library(openxlsx)
relapsetime <- read.xlsx("/Tables/SupplTable1noBlockInfo.xlsx")
mutationcounts <- dplyr::left_join(mutationcounts,relapsetime,by=c("patient"="Patient.ID"))
clinical <- dplyr::left_join(mutationcounts,meandf)
rownames(clinical) <- clinical$sample
write.table(clinical, "/Tables/SupplTable1-additionalInfo.txt", sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
clinical <- clinical[,c("Days.to.relapse","SBS15 Defective MMR","SBS11 Temozolomide treatment","TMBcategory")]
clinical$`SBS11 Temozolomide treatment` <- as.numeric(as.character(clinical$`SBS11 Temozolomide treatment`))*100
clinical$`SBS15 Defective MMR` <- as.numeric(as.character(clinical$`SBS15 Defective MMR`))*100
clinicalannotation <- setNames(clinical,c("Days until relapse","% MMR signature","% Temozolomide signature","TMB category"))
clinical$TMBcategory <- factor(clinical$TMBcategory,levels = c(c("<6 mutations/Mb","6-20 mutations/Mb","21-50 mutations/Mb",">50 mutations/Mb")))

relapse_col <- circlize::colorRamp2(c(0, 2333), c("white", "#228B22"))
TMB_col <- c("<6 mutations/Mb"='#feedde',"6-20 mutations/Mb"='#fdbe85',"21-50 mutations/Mb"='#fd8d3c',">50 mutations/Mb"='#d94701')
temozolomide_col <- circlize::colorRamp2(c(0,100), c("white", "firebrick"))
MMR_col <- circlize::colorRamp2(c(0,100), c("white", "#35978f"))


fontsizeforall <- 9
clin_anno = HeatmapAnnotation(df = clinicalannotation,
                              col = list('TMB category' = TMB_col,
                                         '% MMR signature' = MMR_col,
                                         '% Temozolomide signature'=temozolomide_col,
                                         'Days until relapse' = relapse_col
                              ),
                              annotation_name_gp = gpar(fontsize = fontsizeforall),simple_anno_size= unit(4,"mm"),
                              show_legend = TRUE, gap = unit(c(1,1,1,1,1), "mm"))

bottom_anno = HeatmapAnnotation(text = anno_text(colnames(opi), gp = gpar(fontsize = fontsizeforall)))

ht_list = Heatmap(opi, col = colors_oncoprint,rect_gp = gpar(col = "grey"),
                            na_col = "white",name="Mutation type",show_column_names = FALSE,
                            row_names_gp = gpar(fontsize = fontsizeforall),column_split = columnsplit, column_title_gp = gpar(fontsize=fontsizeforall),
                            top_annotation = c(clin_anno),bottom_annotation = bottom_anno)


draw(ht_list, padding = unit(c(2, 2, 2, 2), "mm"),  merge_legend = TRUE, heatmap_legend_side = "right",annotation_legend_side = "top")
save.image("/Plots/HeatmapS3B.RData")
#################################

#################################
# Recurrence-only variants -- Figure S6
rm(list=ls())
dev.off()
 
library(ComplexHeatmap)
library(RColorBrewer)
source("/Scripts/oncoprint_examples.R")

muts <- readRDS("/VCFs/mergedVCF_withoutInterrogationAbsentandWeak_onlyexonic.RDS")

samples <- c("PreGBM002T", "PostGBM002T", "PreGBM005T", "PostGBM005T", "PreGBM006T", "PostGBM006T", "PreGBM009T", "PostGBM009T", 
             "PreGBM011T", "PostGBM011T", "PreGBM014T", "PostGBM014T", "PreGBM017T", "PostGBM017T", "PreGBM018T", "PostGBM018T", 
             "PreGBM026T", "PostGBM026T", "PreGBM027T", "PostGBM027T", "PreGBM036T", "PostGBM036T", "PreGBM040T", "PostGBM040T", 
             "PreGBM045T", "PostGBM045T", "PreGBM052T", "PostGBM052T", "PreGBM057T", "PostGBM057T", 
             "PreGBM059T", "PostGBM059T", "PreGBM064T", "PostGBM064T", "PreGBM066T", "PostGBM066T", "PreGBM073T", "PostGBM073T", 
             "PreGBM075T", "PostGBM075T", "PreGBM076T", "PostGBM076T", "PreGBM077T", "PostGBM077T", "PreGBM082T", "PostGBM082T", 
             "PreGBM083T", "PostGBM083T", "PreGBM084T", "PostGBM084T", "PreGBM086T", "PostGBM086T", "PreGBM088T", "PostGBM088T")
samples <- factor(samples, levels = samples)

# remove GBM046, it has no pre-treatment sample (only 2 post-treatment samples)
pre <- which(muts$sample == "PreGBM046T")
post <- which(muts$sample == "PostGBM046T")
muts <- muts[-c(pre,post),]

### change the Effect to hotspot for all hotspot variants
muts[which(muts$HOTSPOT == "true"),"ANN.EFFECT"] <- "Hotspot"
muts$matchID <- paste0(muts$CHROM,":",muts$POS,":",muts$REF,":",muts$ALT,"|",muts$ANN.GENE,"|",muts$ANN.EFFECT)

patients <- unique(gsub("Pre|Post|T","",muts$sample))
allpostonly <- setNames(data.frame(matrix(nrow = 0,ncol=4)),c("Patient","mutation","GENE","EFFECT"))

for(this.patient in patients){
  preMutations <- muts[which(muts$sample == paste0("Pre",this.patient,"T")),"matchID"]
  postMutations <- muts[which(muts$sample == paste0("Post",this.patient,"T")),"matchID"]
  overlap <- intersect(preMutations,postMutations)
  postonly <- postMutations[-which(postMutations %in% preMutations)]
  this.postonly <- setNames(data.frame(rep(this.patient, times=length(postonly)),c(postonly)),c("Patient","mutation"))
  this.postonly$GENE <- unlist(lapply(strsplit(as.character(this.postonly$mutation),split="\\|"), "[", 2))
  this.postonly$EFFECT <- unlist(lapply(strsplit(as.character(this.postonly$mutation),split="\\|"), "[", 3))
  allpostonly <- rbind(allpostonly,this.postonly)
}


muttable <- setNames(allpostonly[,c("Patient", "GENE", "EFFECT")],c("TUMOR_SAMPLE", "GENE", "EFFECT"))
# too many genes to plot --> reduce to genes which are post-only in at least 5 patients
allPatGeneCombos <- setNames(as.data.frame(unique(paste0(allpostonly$Patient,"-",allpostonly$GENE))),"code")
allPatGeneCombos$GENE <- unlist(lapply(strsplit(as.character(allPatGeneCombos$code),split="-"), "[", 2))
genestokeep <- as.character(names(table(allPatGeneCombos$GENE)[-which(table(allPatGeneCombos$GENE)<5)]))
muttable <- muttable[which(muttable$GENE %in% genestokeep),]

opi <- make_oncoprint_input(muttable, mutation_genes=unique(muttable$GENE), sample_names=unique(muttable$TUMOR_SAMPLE))

missingpatients <- as.character(patients[-which(patients %in% colnames(opi))])
missingpatients <- setNames(data.frame(matrix(nrow = nrow(opi),ncol=length(missingpatients))),missingpatients)
opi <- cbind(opi, missingpatients)

colnames(opi) <- gsub("GBM0","",colnames(opi))
opi <- as.matrix(opi)

anno_width = unit(3, "cm")
fontsizeforall <- 9
colors_oncoprint <- setNames(c("#88CCEE","#CC6677","#fcc5c0","forestgreen","#332288","#fee391","#9e9ac8"),
                             c("Frameshift","Hotspot","Hotspot, Missense","Inframe","Missense","Nonsense","Nonsense, Missense"))
colors_oncoprint <- colors_oncoprint[which(names(colors_oncoprint) %in% names(table(opi)))]
oncoPrint(opi, alter_fun = oncoprint_alter_fun, col = colors_oncoprint,show_column_names=TRUE)
#################################

#################################
# GBM driver -- Mutations section of Figure 1D
rm(list=ls())
dev.off()

library(ComplexHeatmap)
library(RColorBrewer)
source("/Scripts/oncoprint_examples.R")

# Prepare mutation data
muts <- readRDS("/VCFs/mergedVCF_withoutInterrogationAbsentandWeak_onlyexonic.RDS")

## remove PostGBM046T and rename PreGBM046T to PostGBM046T (both samples are post-treatment and we want to remove one of them)
pre <- which(muts$sample == "PreGBM046T")
post <- which(muts$sample == "PostGBM046T")
muts[pre,"sample"] <- "PostGBM046T"
muts <- muts[-post,]

samples <- c("PreGBM002T", "PostGBM002T", "PreGBM005T", "PostGBM005T", "PreGBM006T", "PostGBM006T", "PreGBM009T", "PostGBM009T", 
             "PreGBM011T", "PostGBM011T", "PreGBM014T", "PostGBM014T", "PreGBM017T", "PostGBM017T", "PreGBM018T", "PostGBM018T", 
             "PreGBM026T", "PostGBM026T", "PreGBM027T", "PostGBM027T", "PreGBM036T", "PostGBM036T", "PreGBM040T", "PostGBM040T", 
             "PreGBM045T", "PostGBM045T", "PostGBM046T", "PreGBM052T", "PostGBM052T", "PreGBM057T", "PostGBM057T", 
             "PreGBM059T", "PostGBM059T", "PreGBM064T", "PostGBM064T", "PreGBM066T", "PostGBM066T", "PreGBM073T", "PostGBM073T", 
             "PreGBM075T", "PostGBM075T", "PreGBM076T", "PostGBM076T", "PreGBM077T", "PostGBM077T", "PreGBM082T", "PostGBM082T", 
             "PreGBM083T", "PostGBM083T", "PreGBM084T", "PostGBM084T", "PreGBM086T", "PostGBM086T", "PreGBM088T", "PostGBM088T")#"PreGBM064T", 
samples <- factor(samples, levels = samples)

### change the Effect to hotspot for all hotspot variants
muts[which(muts$HOTSPOT == "true"),"ANN.EFFECT"] <- "Hotspot"

muttable <- setNames(muts[,c("sample", "ANN.GENE", "ANN.EFFECT")],c("TUMOR_SAMPLE", "GENE", "EFFECT"))

# Select genes to plotś
evgenes <- c("IDH1","EGFR","NF1","RB1","PTEN","PIK3CA","PIK3R1","TP53","ATRX","MTOR")
muttable <- muttable[which(muttable$GENE %in% evgenes),]

# prepare matrix for heatmap
opi <- make_oncoprint_input(muttable, mutation_genes=unique(muttable$GENE), sample_names=unique(muttable$TUMOR_SAMPLE))

# ADD missing samples (samples with no mutation in any of the selected genes)
missingsamples <- as.character(samples[-which(samples %in% colnames(opi))])
missingsamples <- setNames(data.frame(matrix(nrow = nrow(opi),ncol=length(missingsamples))),c(missingsamples))
opi <- cbind(opi, missingsamples)

# reorder
opi <- opi[,levels(samples)]
opi <- as.matrix(opi)
saveRDS(opi,"/Plots/MutationsForHeatmap-opiFormat.RDS")

# add column split by patient
columnsplit <- gsub("PreGBM0|PostGBM0|T","",colnames(opi))
colnames(opi) <- unlist(lapply(strsplit(colnames(opi),split="GBM"), "[", 1))

anno_width = unit(3, "cm")
fontsizeforall <- 9
colors_oncoprint <- setNames(c("#88CCEE","#CC6677","#fcc5c0","forestgreen","#332288","#fee391","#9e9ac8"),
                             c("Frameshift","Hotspot","Hotspot, Missense","Inframe","Missense","Nonsense","Nonsense, Missense"))

opi <- opi[c("IDH1","EGFR","NF1","RB1","PTEN","PIK3CA","PIK3R1","TP53","ATRX","MTOR"),]
colnames(opi)[which(colnames(opi)=="Pre")] <- "Initial"
colnames(opi)[which(colnames(opi)=="Post")] <- "Recurrent"

mutations_legend = Legend(labels = names(colors_oncoprint), title = "", legend_gp = gpar(fill = colors_oncoprint))

mutations_plot <- Heatmap(opi, col = colors_oncoprint, rect_gp = gpar(col = "grey"),
                          show_heatmap_legend = FALSE,
                          na_col = "white",name=" ",column_names_gp = gpar(fontsize = fontsizeforall),
                          row_names_gp = gpar(fontsize = fontsizeforall),column_split = columnsplit, column_title_gp = gpar(fontsize=fontsizeforall))

draw(mutations_plot, padding = unit(c(2, 2, 20, 2), "mm"), heatmap_legend_list = mutations_legend, merge_legend = FALSE, heatmap_legend_side = "right",annotation_legend_side = "bottom")

# save it to plot it together with the CN
saveRDS(mutations_legend,"/Plots/mutations_legend.RDS")
saveRDS(mutations_plot,"/Plots/mutations_plot.RDS")
##############################################################################################################################################################################
##############################################################################################################################################################################


##############################################################################################################################################################################
##############################################################################################################################################################################
# Copy number
# Gene-level amplifications, arm-level events section -- Figure 1D
# GISTIC analysis of tumors with late recurrence -- Figure S5A&S5B
#####################################

############
# preprocess data
rm(list = ls())

geneCN <- "/CNV/all_thresholded.by_genes.renamed.txt"

geneCN_tab <- read.table(geneCN, sep="\t", header=T, stringsAsFactors=F, check.names=F)

# P046-Pre is a recurrence sample and P046-Post got excluded
pre <- which(colnames(geneCN_tab) == "Pre_GBM046_T")
colnames(geneCN_tab)[pre] <- "Post_GBM046_T"


left <- geneCN_tab[,1:3]
left$chrom <- unlist(lapply(strsplit(geneCN_tab$Cytoband,split="p|q"), "[", 1))
left$arm <- unlist(lapply(strsplit(geneCN_tab$Cytoband,split="(?<=[/p|q])", perl = TRUE), "[", 1))
geneCN_tab <- geneCN_tab[,sort(colnames(geneCN_tab))]
pre <- geneCN_tab[,grep("Pre",colnames(geneCN_tab))]
post <- geneCN_tab[,grep("Post",colnames(geneCN_tab))]
geneCN_tab <- cbind(left,pre,post)

geneCN_tab <- geneCN_tab[geneCN_tab$chrom != "Y",]

write.table(geneCN_tab, "/CNV/all.geneCN.GL_ASCNA.txt", sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
geneCN_tab <- read.table("/CNV/all.geneCN.GL_ASCNA.txt", header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")
############

############
# print full CN profile and frequency plots
rm(list=ls())

source("/Scripts/process_facets.R")
geneCN <- "/CNV/all.geneCN.GL_ASCNA.txt"
outFile <- "/CNV/all.geneCN.GL_ASCNA.pdf"

geneCN_tab <- read.table(geneCN, header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")

colnames(geneCN_tab)[grep("Pre",colnames(geneCN_tab))] <- gsub("Pre_","",paste0(gsub("_T","",colnames(geneCN_tab)[grep("Pre",colnames(geneCN_tab))]),"-Initial"))
colnames(geneCN_tab)[grep("Post",colnames(geneCN_tab))] <- gsub("Post_","",paste0(gsub("_T","",colnames(geneCN_tab)[grep("Post",colnames(geneCN_tab))]),"-Recurrent"))

left <- geneCN_tab[,1:5]
right <- geneCN_tab[,6:ncol(geneCN_tab)]
right <- right[,order(colnames(right))]

geneCN_tab <- cbind(left,right)

plot_heatmapgistic(geneCN_tab, outFile)# plot with 14x10 inch

pre <- geneCN_tab[,c(1:5, grep("Initial",colnames(geneCN_tab)))]
post <- geneCN_tab[,c(1:5, grep("Recurrent",colnames(geneCN_tab)))]

plot_frequency(pre,plot_file_prefix="Pre")
plot_frequency(post,plot_file_prefix="Post")
plot_frequency(geneCN_tab,plot_file_prefix="Both")
##########

##########
# Armlevel calls
rm(list = ls())
library(ComplexHeatmap)
library(RColorBrewer)

armdata <- read.table("/CNV/all.geneCN.GL_ASCNA.txt", header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")

armcalls <- setNames(as.data.frame(matrix(nrow=length(unique(armdata$arm)),ncol = length(grep("Pre|Post",colnames(armdata))))),
                     colnames(armdata)[grep("Pre|Post",colnames(armdata))])
rownames(armcalls) <- unique(armdata$arm)

for(thisarm in rownames(armcalls)){
  for(thissample in colnames(armcalls)){
    thisarmdata <- armdata[which(armdata$arm==thisarm),thissample]
    if((length(which(thisarmdata>0))/length(thisarmdata))>0.6){armcalls[thisarm,thissample] <- "Gain"}
    if((length(which(thisarmdata<0))/length(thisarmdata))>0.6){armcalls[thisarm,thissample] <- "Loss"}
}}
arm <- as.data.frame(armcalls)

colnames(arm) <- gsub("_T","",colnames(arm))

# reorder and add dummy columns for the samples without CN data
samples <- c("Pre_GBM002", "Post_GBM002", "Pre_GBM005", "Post_GBM005", "Pre_GBM006", "Post_GBM006", "Pre_GBM009", "Post_GBM009", 
             "Pre_GBM011", "Post_GBM011", "Pre_GBM014", "Post_GBM014", "Pre_GBM017", "Post_GBM017", "Pre_GBM018", "Post_GBM018", 
             "Pre_GBM026", "Post_GBM026", "Pre_GBM027", "Post_GBM027", "Pre_GBM036", "Post_GBM036", "Pre_GBM040", "Post_GBM040", 
             "Pre_GBM045", "Post_GBM045", "Post_GBM046", "Pre_GBM052", "Post_GBM052", "Pre_GBM057", "Post_GBM057", 
             "Pre_GBM059", "Post_GBM059", "Pre_GBM064", "Post_GBM064", "Pre_GBM066", "Post_GBM066", "Pre_GBM073", "Post_GBM073", 
             "Pre_GBM075", "Post_GBM075", "Pre_GBM076", "Post_GBM076", "Pre_GBM077", "Post_GBM077", "Pre_GBM082", "Post_GBM082", 
             "Pre_GBM083", "Post_GBM083", "Pre_GBM084", "Post_GBM084", "Pre_GBM086", "Post_GBM086", "Pre_GBM088", "Post_GBM088")#"Pre_GBM046", 

missing <- samples[-which(samples %in% colnames(arm))]
dummycolumns <- setNames(data.frame(matrix(nrow = nrow(arm), ncol = length(missing),rep("Not available",each=nrow(arm)))),missing)
arm <- cbind(arm,dummycolumns)
samples <- factor(samples, levels = samples)
arm <- arm[,levels(samples)]
arm <- as.data.frame(arm)

saveRDS(arm,"/Plots/GisticAllSamples-GainLoss-AllArms-opiFormat.RDS")

# add column split by patient
columnsplit <- gsub("Pre_GBM0|Post_GBM0|T","",colnames(arm))
colnames(arm) <- unlist(lapply(strsplit(colnames(arm),split="_GBM"), "[", 1))
rownames(arm) <- unlist(lapply(strsplit(rownames(arm),split="--"), "[", 1))

colnames(arm)[which(colnames(arm)=="Pre")] <- "Initial"
colnames(arm)[which(colnames(arm)=="Post")] <- "Recurrent"

arm <- arm[c("1p","19q","7p","7q","10p","10q"),]
arm <- as.matrix(arm)

anno_width = unit(3, "cm")
fontsizeforall <- 9
colors_arm <- c("Gain"="darksalmon","Loss"="lightblue","Not available"="#D7D7D7")
arm_legend = Legend(labels = names(colors_arm), title = "", legend_gp = gpar(fill = colors_arm))
arm_plot <- Heatmap(arm, col = colors_arm,rect_gp = gpar(col = "grey"),
                    show_heatmap_legend = FALSE,column_split = columnsplit,
                    na_col = "white",name=" ",column_names_gp = gpar(fontsize = fontsizeforall),#column_names_gp = gpar(fontsize = fontsizeforall, fill=evolutioncols),
                    row_names_gp = gpar(fontsize = fontsizeforall),column_title_gp = gpar(fontsize=fontsizeforall))
arm_plot


# Gene-level calls: GisticOutput to oncoprint input
ad <- read.table("/CNV/all.geneCN.GL_ASCNA.txt", header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")

library(tidyr)
ad <- ad %>% tidyr::gather(TUMOR_SAMPLE, EFFECT, colnames(ad)[6:ncol(ad)])

effects_oncoprint <- list(
  Amplification=c("2"),
  'Homozygous deletion'=c("-2"),
  'NA'=c("0","1","-1"))

make_oncoprint_input <- function( muts=NULL, effects=effects_oncoprint, mutation_genes=NULL, sample_names=NULL,
                                  sample_name_col="TUMOR_SAMPLE", gene_name_col="hgnc", eff_name_col ="EFFECT") {
  
  library(reshape2)
  if(any(!sample_name_col %in% colnames(muts), !gene_name_col %in% colnames(muts), !eff_name_col %in% colnames(muts))) {
    stop("At least one of the required columns not found in data\n")
  }
  
  effects <- melt(effects)
  mutation_genes <- unlist(mutation_genes)
  sample_names <- unlist(sample_names)
  
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

ad$band <- unlist(lapply(strsplit(ad$Cytoband,split="p|q"), "[", 2))

ad$hgnc <- paste0(ad$Gene.Symbol,"--",ad$Cytoband)
geneCN <- make_oncoprint_input(ad, mutation_genes=unique(ad$hgnc), sample_names=unique(ad$TUMOR_SAMPLE))
geneCN <- apply(geneCN, MARGIN=c(1,2), function(x) gsub("NA",NA,x))

geneCN <- as.data.frame(geneCN)
colnames(geneCN) <- gsub("_T$","",colnames(geneCN))

geneCN <- apply(geneCN, MARGIN=c(1,2), function(x) gsub("^$|NA",NA,x))

# reorder and add dummy columns for the samples without CN data
samples <- c("Pre_GBM002", "Post_GBM002", "Pre_GBM005", "Post_GBM005", "Pre_GBM006", "Post_GBM006", "Pre_GBM009", "Post_GBM009", 
             "Pre_GBM011", "Post_GBM011", "Pre_GBM014", "Post_GBM014", "Pre_GBM017", "Post_GBM017", "Pre_GBM018", "Post_GBM018", 
             "Pre_GBM026", "Post_GBM026", "Pre_GBM027", "Post_GBM027", "Pre_GBM036", "Post_GBM036", "Pre_GBM040", "Post_GBM040", 
             "Pre_GBM045", "Post_GBM045", "Post_GBM046", "Pre_GBM052", "Post_GBM052", "Pre_GBM057", "Post_GBM057", 
             "Pre_GBM059", "Post_GBM059", "Pre_GBM064", "Post_GBM064", "Pre_GBM066", "Post_GBM066", "Pre_GBM073", "Post_GBM073", 
             "Pre_GBM075", "Post_GBM075", "Pre_GBM076", "Post_GBM076", "Pre_GBM077", "Post_GBM077", "Pre_GBM082", "Post_GBM082", 
             "Pre_GBM083", "Post_GBM083", "Pre_GBM084", "Post_GBM084", "Pre_GBM086", "Post_GBM086", "Pre_GBM088", "Post_GBM088")#"Pre_GBM046", 

missing <- samples[-which(samples %in% colnames(geneCN))]
dummycolumns <- setNames(data.frame(matrix(nrow = nrow(geneCN), ncol = length(missing),rep("Not available",each=nrow(geneCN)))),missing)
geneCN <- cbind(geneCN,dummycolumns)
samples <- factor(samples, levels = samples)
geneCN <- geneCN[,levels(samples)]
geneCN <- as.data.frame(geneCN)

# remove homozygous loss in PDGFRA in P076 pre-treatment -- it's called from a centromere region
geneCN["PDGFRA--4q12","Pre_GBM076"] <- NA

rownames(geneCN) <- unlist(lapply(strsplit(rownames(geneCN),split="--"), "[", 1))

toremove <- c()
for(r in 1:nrow(geneCN)){
  if(length(c(grep("Amplification",geneCN[r,]),grep("Homozygous deletion",geneCN[r,])))==0){
    toremove <- c(toremove,r)
  }
}
geneCN <- geneCN[-toremove,]

saveRDS(geneCN,"/Plots/GisticAllSamples-AmpHomDel-AllGenes-opiFormat.RDS")

genesetToPlot <- c("CDKN2A", "EGFR", "PDGFRA", "CDK4", "CDK6", "MDM2", "MET", "CCND2", "MYCN", "MDM4", "ATRX", "NF1", "RB1", "TP53")

geneCN <- geneCN[which(unlist(lapply(strsplit(as.character(rownames(geneCN)),split="--"), "[", 1)) %in% genesetToPlot),]

# reorder rows
newroworder <- genesetToPlot[which(genesetToPlot %in% rownames(geneCN))]
geneCN <- geneCN[newroworder,]

# add column split by patient
columnsplit <- gsub("Pre_GBM0|Post_GBM0|T","",colnames(geneCN))
colnames(geneCN) <- unlist(lapply(strsplit(colnames(geneCN),split="_GBM"), "[", 1))
colnames(geneCN)[which(colnames(geneCN)=="Pre")] <- "Initial"
colnames(geneCN)[which(colnames(geneCN)=="Post")] <- "Recurrent"
geneCN <- as.matrix(geneCN)

anno_width = unit(3, "cm")
fontsizeforall <- 9
colors_geneCN <- c("Amplification"="firebrick","Not available"="#D7D7D7")
geneCN_legend = Legend(labels = names(colors_geneCN), title = "", legend_gp = gpar(fill = colors_geneCN))
geneCN_plot <- Heatmap(geneCN, col = colors_geneCN,rect_gp = gpar(col = "grey"),
                       na_col = "white",name=" ",column_names_gp = gpar(fontsize = fontsizeforall),
                       show_heatmap_legend = FALSE,
                       row_names_gp = gpar(fontsize = fontsizeforall),column_split = columnsplit,column_title_gp = gpar(fontsize=fontsizeforall))
geneCN_plot
###########

#########################################################
#### make plot with mutations, geneCN and arm level
#########################################################
mutations_legend <- readRDS("/Plots/mutations_legend.RDS")
mutations_plot <- readRDS("/Plots/mutations_plot.RDS")

pd = packLegend(mutations_legend, geneCN_legend,arm_legend)#, groups_legend)
ht_list <- mutations_plot %v% geneCN_plot %v% arm_plot
draw(ht_list, padding = unit(c(2, 2, 20, 2), "mm"), heatmap_legend_list = pd, merge_legend = FALSE, heatmap_legend_side = "right",annotation_legend_side = "bottom")
save.image("/Plots/Figure1D.RData")
##############################################################################################################################################################################
##############################################################################################################################################################################

##############################################################################################################################################################################
##############################################################################################################################################################################
# Evolutionary analyses
# Fishplots -- Figure 2B
# Heatmaps, line plots and Fishplots -- Figure S7
#####################################
rm(list = ls())

setwd("/phyloHeatmaps/")

# which driver genes to annotate on the heatmap
gene_ann_file <- read.table("IntOGen-DriverGenes.tsv", header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=TRUE, quote = '""')


library(ComplexHeatmap)
library(ggplot2)
library(fishplot)
library(dplyr)

alltrees <- list.files("/phylogic/")
alltrees <- alltrees[grep(".timing.tsv", alltrees)]
allpatients <- gsub(".timing.tsv","",alltrees)

allpatients <- setNames(as.data.frame(allpatients),c("patient"))

thispatient <- allpatients$patient[1]

allmuts <- c()

for(thispatient in unique(allpatients$patient)){
  
  #############
  # The cluster CCF values in “<PatientID>_constrained_ccf.tsv” file are the ones used by phylogic to construct “Tree 1" and also to 
  # create the file with the clonal abundances (“<PatientID>_cell_population_abundances.tsv”). 
  # These don't always make sense (e.g. clusters with a mean CCF of 0.5 get forced down to a constrained CCF of ~0 for no obvious reason (e.g. no evidence that the mutations are all in CNA regions, are only INDELs or similar) )
  pathclusterccf <- paste0("/phylogic/",thispatient,".cluster_ccfs.txt")
  cluccf <- read.table(pathclusterccf, header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")
  cluccf <- setNames(cluccf[,c("Patient_ID","Sample_ID","Cluster_ID","postDP_ccf_mean")],c("Patient_ID","Sample_ID","Cluster","meanCCF"))
  cluccf$timepoint <- unlist(lapply(strsplit(cluccf$Sample_ID,split="_"), "[", 1))
  cluccf$timepoint <- gsub("Pre","Initial", cluccf$timepoint)
  cluccf$timepoint <- gsub("Post","Recurrent", cluccf$timepoint)
  cluccf$meanCCF <- as.numeric(as.character(cluccf$meanCCF))*100
  cluccf$Cluster <- as.character(cluccf$Cluster)
  
  pathmuts <- paste0("/phylogic/",thispatient,".mut_ccfs.txt")
  muts <- read.table(pathmuts, header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")
  muts$matchID <- paste(sep=":",muts$Chromosome,muts$Start_position,muts$Reference_Allele,muts$Tumor_Seq_Allele)
  
  # some clusters are <10% CCF in both samples and therefore not listed in “<PatientID>_meanCCF.tsv” --> remove them. We only use clusters with >0.1 CCF in at least one of the samples
  print("Removed clusters and sizes")
  print(table(muts[-which(muts$Cluster_Assignment %in% unique(cluccf$Cluster)),"Cluster_Assignment"]))
  muts <- muts[which(muts$Cluster_Assignment %in% unique(cluccf$Cluster)),]
  
  left <- muts[!duplicated(muts$matchID),c("Patient_ID","Hugo_Symbol","Protein_change","Variant_Classification","Variant_Type","Cluster_Assignment","matchID","Chromosome","Start_position","Reference_Allele","Tumor_Seq_Allele")]
  pre <- setNames(muts[which(muts$Sample_ID==paste0("Pre_",thispatient,"_T")),c("matchID","t_ref_count","t_alt_count","preDP_ccf_mean","preDP_ccf_CI_low","preDP_ccf_CI_high")],c("matchID","ref_count","alt_count","ccf_mean","ccf_CI_low","ccf_CI_high"))
  post <- setNames(muts[which(muts$Sample_ID==paste0("Post_",thispatient,"_T")),c("matchID","t_ref_count","t_alt_count","preDP_ccf_mean","preDP_ccf_CI_low","preDP_ccf_CI_high")],c("matchID","ref_count","alt_count","ccf_mean","ccf_CI_low","ccf_CI_high"))
  
  muts <- dplyr::left_join(left, dplyr::left_join(pre,post,by=c("matchID"="matchID"),suffix = c(".pre", ".post")))
  muts$clusterSize <- "NA"
  muts$clusterCCFpre <- 0
  muts$clusterCCFpost <- 0
  for(thiscluster in unique(cluccf$Cluster)){
    muts[which(muts$Cluster_Assignment == thiscluster),"clusterCCFpre"] <- cluccf[which(cluccf$Cluster == thiscluster & cluccf$Sample_ID == paste0("Pre_",thispatient,"_T")),"meanCCF"]
    muts[which(muts$Cluster_Assignment == thiscluster),"clusterCCFpost"] <- cluccf[which(cluccf$Cluster == thiscluster & cluccf$Sample_ID == paste0("Post_",thispatient,"_T")),"meanCCF"]
    muts[which(muts$Cluster_Assignment == thiscluster),"clusterSize"] <- nrow(muts[which(muts$Cluster_Assignment == thiscluster),])
  }
  muts$clonalitypre <- apply(as.data.frame(muts$clusterCCFpre), MARGIN = c(1), function(x)if(x>=80){x <- "clonal"}else if(x>=10 & x<80){x <- "subclonal"}else{x <- "missing"})
  muts$clonalitypost <- apply(as.data.frame(muts$clusterCCFpost), MARGIN = c(1), function(x)if(x>=80){x <- "clonal"}else if(x>=10 & x<80){x <- "subclonal"}else{x <- "missing"})
  
  clusterstoremove <- unique(muts[which(muts$clonalitypre == "missing" & muts$clonalitypost == "missing"),"Cluster_Assignment"])
  if(thispatient == "GBM064"){clusterstoremove <- 3}
  
  if(length(clusterstoremove)>0){muts <- muts[-which(muts$Cluster_Assignment %in% clusterstoremove),]}
  
  muts <- muts %>% dplyr::arrange(Cluster_Assignment)
  muts$clusterSize <- as.numeric(as.character(muts$clusterSize))
  allmuts <- rbind(muts,allmuts)
  ####################################
  
  
  ####################################
  # CCFs pre vs post - lineplot
  if(length(clusterstoremove)>0){cluccf <- cluccf[-which(cluccf$Cluster %in% clusterstoremove),]}
  clusterplot <- ggplot(cluccf,  aes(x = timepoint, y = meanCCF, group=Cluster)) +
    geom_line(aes(color=Cluster),size=1.5)+
    geom_point(aes(color=Cluster),size=2)+
    scale_y_continuous(breaks=seq(0, 100, by = 20))+
    scale_color_manual(values=c("#41ae76", "skyblue", "#E69F00","black","#8c6bb1","firebrick","#a8ddb5","salmon"), name = "Cluster")+
    theme_bw()+
    geom_hline(yintercept = c(10,80), color="darkgrey")+
    xlab("")+
    ylab("CCF")+
    theme(title=element_text(size=10),
          axis.text.x=element_text(size=10), axis.title.x=element_text(size=10),
          axis.text.y =element_text(size=10), axis.title.y=element_text(size=10),
          legend.position="right",legend.title = element_text(size=8), legend.text = element_text(size=8),plot.margin = unit(c(t = 0.5, r = 0.1, b = 0.0, l = 0.1), "cm"))
  clusterplot
  thisfilename <- paste0(thispatient,"-CCFs.pdf")
  ggsave(filename = thisfilename,plot=clusterplot,height = 2,width = 2.4,units=c("in"),device="pdf")
  ##################################
  
  
  #################################
  # Heatmap
  
  # matrix with CCFs
  opi <- setNames(muts[,c("ccf_mean.pre","ccf_mean.post")],c("Initial","Recurrent"))
  opi <- as.matrix(opi)*100
  
  fontsizeforall <- 9
  rowsplit <- muts$Cluster_Assignment
  
  # define colors
  ccf_col <- circlize::colorRamp2(c(0,10,20,40,60,80,100),c('white','#ffffb2','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026'))
  clusterID_col <- c("1"='#41ae76', "2"='skyblue', "3"='#E69F00', "4"='black', "5"='#8c6bb1', "6"='firebrick', "7"='#a8ddb5', "8"='salmon')
  size_col <- circlize::colorRamp2(c(0,max(muts$clusterSize)), c("white", "#dd1c77"))
  clonality_col <- c("missing"='#d9d9d9', "subclonal"='#e7d4e8', "clonal"='#762a83')
  
  inital_ann <- Heatmap(muts[,"clonalitypre"], name = "Initial", col = clonality_col,width = unit(5, "mm"),column_names_gp = grid::gpar(fontsize = fontsizeforall),show_row_names = FALSE, cluster_rows = FALSE,cluster_columns = FALSE) 
  recurrent_ann <- Heatmap(muts[,"clonalitypost"], name = "Recurrent", col = clonality_col,width = unit(5, "mm"),column_names_gp = grid::gpar(fontsize = fontsizeforall),show_row_names = FALSE, cluster_rows = FALSE,cluster_columns = FALSE) 
  size_ann <- Heatmap(muts[,"clusterSize"], name = "Size", col = size_col,width = unit(5, "mm"),column_names_gp = grid::gpar(fontsize = fontsizeforall), show_row_names = FALSE, cluster_rows = FALSE,cluster_columns = FALSE) 
  clust_ann <- Heatmap(muts[,"Cluster_Assignment"], name = "Cluster", col = clusterID_col,width = unit(5, "mm"),column_names_gp = grid::gpar(fontsize = fontsizeforall),row_split = rowsplit,show_row_names = FALSE, cluster_rows = FALSE,cluster_columns = FALSE) 
  
  ## Gene annotation on heatmap
  # only non-synoymous variants
  indels <- c("conservative_inframe_deletion", "conservative_inframe_insertion", "disruptive_inframe_deletion", "disruptive_inframe_deletion&splice_region_variant", 
              "disruptive_inframe_insertion", "frameshift_variant", "frameshift_variant&splice_region_variant", "frameshift_variant&start_lost", "frameshift_variant&stop_gained")
  nssnvs <- c("missense_variant", "missense_variant&splice_region_variant", "splice_region_variant", "splice_region_variant&stop_retained_variant", 
              "splice_region_variant&synonymous_variant", "start_lost", "start_lost&splice_region_variant", "stop_gained", "stop_gained&splice_region_variant", 
              "stop_lost", "stop_lost&conservative_inframe_deletion")
  allns <- c(indels, nssnvs)
  
  genelist <- c("IDH1","EGFR","NF1","RB1","PTEN","PIK3CA","PIK3R1","TP53","ATRX","MTOR")#gene_ann_file$Symbol
  
  # Generate mutationID (Gene + protein change) for the labeling of the heatmap.
  # There are a few duplicated mutationIDs in GBM006,GBM040, GBM059, GBM073, GBM076, GBM084 but they are all non-coding mutations (e.g. 2x'VPS35P1-NA' in GBM011 =two different 'non_coding_transcript_exon_variant' in 'VPS35P1') so this is not a problem.
  muts$mutationID <- paste0(muts$Hugo_Symbol,"-",muts$Protein_change)
  highlight.mutations <- which(muts$Hugo_Symbol %in% genelist & muts$Variant_Classification %in% allns)
  ha <- rowAnnotation(genes=anno_mark(at=highlight.mutations, labels = muts$mutationID[highlight.mutations],labels_gp = gpar(fontsize = fontsizeforall)))#fontsizeforall
  
  ht_list = clust_ann+
    size_ann+
    inital_ann+
    recurrent_ann+
    Heatmap(as.matrix(opi), border = FALSE, col = ccf_col, cluster_rows = FALSE,cluster_columns = FALSE,width = unit(20, "mm"),#rect_gp = gpar(col = "white", lwd = 0.00001), 
            right_annotation = ha, annotation_legend_size(fontsizeforall),row_names_gp = gpar(fontsize = fontsizeforall),column_names_gp = grid::gpar(fontsize = fontsizeforall),
            na_col = "white",name="CCF",show_row_names=FALSE, row_split = rowsplit,
            column_title = thispatient)
  
  ht_list
  phyloHeatmap <- paste0(thispatient,"_PhyloHeatmap.pdf")
  pdf(file=phyloHeatmap)
  draw(ht_list)
  dev.off()
  #################################
  
  #################################
  # Fishplot
  
  # use the '_cell_population_abundances.tsv' sheet only to extract the name of the parent clone for each subclone
  pathcloneccf <- paste0("/phylogic/",thispatient,"_cell_population_abundances.tsv")
  cloneccf <- read.table(pathcloneccf, header=TRUE,sep = "\t", row.names=NULL, na.string=".", stringsAsFactors=F, quote = "")
  cellpopulations <- paste0("CL0_",unique(cloneccf$Cell_population))
  parents <- as.numeric(as.character(gsub("CL","",unlist(sapply(strsplit(cellpopulations, split = "_"), function(x) x[length(x)-1])))))
  
  timepoints=c(0,1)
  
  frac.table  <- matrix(c(cluccf[grep("Pre",cluccf$Sample_ID),"meanCCF"],
                          cluccf[grep("Post",cluccf$Sample_ID),"meanCCF"]),
                        ncol = length(timepoints))
  
  # # Subclonal clusters are >10% and <=80% CCF
  # # Clonal clusters are >80% CCF
  # # Any clone <=10% CCF is treated as missing in that sample
  frac.table[which(frac.table<10)] <- 0
  if(thispatient=="GBM018"){
    # [,1] [,2]
    # [1,]   92   96 --> Daughter clone's CCF is 1% higher than parent; change 96 to 97
    # [2,]    0   97
    # [3,]    0   62
    # [4,]    0   32
    # [5,]   22    0
    # [6,]    0   14
    frac.table[,2]<- c(97,97,62,32,0,14)
  }
  if(thispatient=="GBM026"){
    # [,1] [,2]
    # [1,]  100   98
    # [2,]   77   88 --> Daughter clone's CCF is 4% higher than parent; change 88 to 92
    # [3,]    0   92
    # [4,]    0   11
    frac.table[,2] <- c(98, 92, 92, 11)
    # default tree is not plausible (1-2 & 1-3 can't be sibling clones),
    # tree nr.4 from the html was selected for this patient
    parents <- c(0,1,2,3)
  }
  if(thispatient=="GBM036"){
    # cluster 4 is missing in the phylogic '_cell_population_abundances.tsv'. It has a mean CCF of 10 in the pre-treatment and 0 post-treatment; phylogic excluded it because the constrained CCFs in both samples were below 10)
    # add it back to the clonal structure and add it to the tree on the html.
    # [,1] [,2]
    # [1,]  100  100
    # [2,]    0  100
    # [3,]    0   44
    # [4,]   10    0
    parents <- c(0,1,2,1)
    cellpopulations <- c(cellpopulations,"CL0_CL1_CL4")
  }
  if(thispatient=="GBM040"){
    # [,1] [,2]
    # [1,]   90   95 --> Daughter clone's CCF is 10% higher than parent; change 90 to 100 (it's the truncal cluster so it's likely 100)
    # [2,]  100    0
    # [3,]    0   72
    # [4,]   73    0
    # [5,]    0   28
    # [6,]    0   11
    frac.table[,1] <- c(100,100,0,73,0,0)
    # default tree is biologically not plausible (clone 1-3 can't be the parent of 5 and 6-2-4),
    # tree nr.9 from the html was selected for this patient (first biologically plausible solution)
    parents <- c(0,1,1,2,3,3)
  }
  if(thispatient=="GBM045"){
    # [,1] [,2]
    # [1,]   91  100 --> change 91 to 94 (it's the truncal cluster so it's likely even 100)
    # [2,]   53   95
    # [3,]    0   89
    # [4,]    0   59
    # [5,]    0   45
    # [6,]   41    0
    # [7,]    0   25
    # [8,]    0   10
    frac.table[,1] <- c(94,53,0,0,0,41,0,0)
    # default tree is biologically not plausible for many reasons,
    # tree nr.9 from the html was selected for this patient (first biologically plausible solution)
    parents <- c(0,1,2,3,4,1,5,5)
  }
  if(thispatient=="GBM052"){
    # [,1] [,2]
    # [1,]   92   95 --> Daughter clone's CCF is 8% resp. 3% higher than parent; change 92 to 100 and 95 to 98
    # [2,]  100    0
    # [3,]    0   98
    # [4,]   67    0
    # [5,]    0   28
    frac.table[,1] <- c(100,100,0,67,0)
    frac.table[,2] <- c(98,0,98,0,28)
    # default tree is biologically not plausible (cluster 2 can't be descendant of cluster 3: cluster 2 is clonal and cluster 3 is missing in the initial sample),
    # tree nr.6 from the html was selected for this patient (first biologically plausible solution)
    parents <- c(0,1,1,2,3)
  }
  if(thispatient=="GBM059"){
    # [,1] [,2]
    # [1,]   92   93 --> Daughter clone's CCF is 6% resp. 1% higher than parent; change 92 to 98 and 93 to 94 
    # [2,]   98    0
    # [3,]    0   94
    # [4,]    0   53
    # [5,]    0   10
    frac.table[,1] <- c(98,98,0,0,0)
    frac.table[,2] <- c(94,0,94,53,10)
    # default tree is biologically not plausible (cluster 2 can't be descendant of cluster 4: cluster 2 is clonal and cluster 4 is missing in the initial sample),
    # tree nr.7 from the html was selected for this patient (first biologically plausible solution)
    parents <- c(0,1,1,3,3)
  }
  if(thispatient=="GBM064"){
    # We removed the copy number driven cluster number 3 (all mutations from this cluster are in CNA regions).
    # remove it from the phylogeny
    parents <- c(0,1,2,3,2)# cluster 3 here is the cluster which was initially labelled as cluster 4
    cellpopulations <- cellpopulations[-3]
  }
  if(thispatient=="GBM073"){
    # [,1] [,2]
    # [1,]   88  100 --> Daughter clone's CCF is 7% higher than parent; change 88 to 95
    # [2,]   95    0
    # [3,]    0  100
    # [4,]    0   58
    # [5,]    0   16
    # [6,]   17    0
    frac.table[,1] <- c(95,95,0,0,0,17)
    # default tree is biologically not plausible (cluster 2 can't be descendant of cluster 4: cluster 2 is clonal and cluster 4 is missing in the initial sample),
    # tree nr.7 from the html was selected for this patient (first biologically plausible solution)
    parents <- c(0,1,1,3,3,2)
  }
  if(thispatient=="GBM084"){
    # [,1] [,2]
    # [1,]   98   92 --> Daughter clone's CCF is 2% higher than parent; change 98 to 100
    # [2,]  100    0
    # [3,]    0   85
    # [4,]   66    0
    # [5,]    0   40
    # [6,]   20    0
    frac.table[,1] <- c(100,100,0,66,0,20)
  }
  if(thispatient=="GBM086"){
    # default tree is biologically not plausible (cluster 2,4,5 can't all be siblings (combined CCF=114 in pre))
    # tree nr.2 from the html was selected for this patient (first biologically plausible solution)
    parents <- c(0,1,2,1,4)
  }
  if(thispatient=="GBM088"){
    # default tree is biologically not plausible (cluster 2 and 3 can't be siblings (CCFs in post: cl2=100;cl3=94)),
    # tree nr.13 from the html was selected for this patient (first biologically plausible solution)
    parents <- c(0,1,2,3,4)
  }
  print(thispatient)
  print(frac.table)
  print(parents)
  
  # does the case have a newly emerging subclone in the recurrence?
  # make this visible on the fishplot by adding a dummy timepoint to stretch the plots
  if(any(frac.table[,2]>0&frac.table[,1]==0)){
    timepointbeforerecurrence <- frac.table[,2]
    timepointbeforerecurrence[which(frac.table[,2]>0&frac.table[,2]<80&frac.table[,1]==0)] <- timepointbeforerecurrence[which(frac.table[,2]>0&frac.table[,2]<80&frac.table[,1]==0)]*0.01
    # does the case have a newly emerging subclone which appears clonal (illusion of clonality) in the recurrence?
    # make this also visible on the fishplot by adding a dummy timepoint to stretch the plots
    if(any(frac.table[,2]>=80&frac.table[,1]==0)){
      additionaltimepoint <- frac.table[,2]
      additionaltimepoint[which(frac.table[,2]>=80&frac.table[,1]==0)] <- additionaltimepoint[which(frac.table[,2]>=80&frac.table[,1]==0)]*0.01
      additionaltimepoint[which(frac.table[,2]>0&frac.table[,2]<80&frac.table[,1]==0)] <- 0
      frac.table <- cbind(frac.table[,1],additionaltimepoint,timepointbeforerecurrence,frac.table[,2])
      timepoints <- c(0,0.75,0.9,1)
    }else{
      frac.table <- cbind(frac.table[,1],timepointbeforerecurrence,frac.table[,2])
      timepoints <- c(0,0.75,1)
    }
  }
  
  # prepare plot
  fish <- createFishObject(frac.table,parents,timepoints=timepoints)
  fish <- layoutClones(fish)
  clustercolors <- c('#41ae76', 'skyblue', '#E69F00','black', '#8c6bb1', 'firebrick', '#a8ddb5', 'salmon')[1:length(cellpopulations)]
  fish <- setCol(fish,col=clustercolors)
  fishout <- paste0(thispatient,"_FishPlot.pdf")
  pdf(file=fishout, width = 4.25, height = 3.6)
  fishPlot(fish,shape="spline",title.btm=gsub("GBM","#",thispatient),cex.vlab = 1,
           cex.title=1, vlines=c(0,1),bg.col = c("white","white","white"),
           vlab=c("Initial","Recurrent"))
  dev.off()
  
}

allmuts$group <- paste0(allmuts$clonalitypre,"-",allmuts$clonalitypost)
saveRDS(allmuts,"/Tables/MutationsEvolutionaryAnalyses.RDS")
write.table(allmuts, "/Tables/MutationsEvolutionaryAnalyses.txt", sep = "\t", quote=FALSE, row.names=FALSE,col.names = TRUE)
##############################################################################################################################################################################
##############################################################################################################################################################################
