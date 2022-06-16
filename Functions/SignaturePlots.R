############################
### add aetiology to labels
############################
renameSigs <- function(tmp.keep) {
  if(length(grep("SBS",colnames(tmp.keep)))>0){
    tmp.keep <- t(tmp.keep)
  }
  if(length(which(rownames(tmp.keep)=="SBS1"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS1")] <- "SBS1-Deamination of 5-methylcytosine"}
  if(length(which(rownames(tmp.keep)=="SBS2"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS2")] <- "SBS2-APOBEC activity"}
  if(length(which(rownames(tmp.keep)=="SBS3"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS3")] <- "SBS3-Defective HR repair (BRCA1/2 mut)"}
  if(length(which(rownames(tmp.keep)=="SBS4"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS4")] <- "SBS4-Tobacco smoking"}
  if(length(which(rownames(tmp.keep)=="SBS5"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS5")] <- "SBS5-Unknown - clock-like"}
  if(length(which(rownames(tmp.keep)=="SBS6"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS6")] <- "SBS6-Defective MMR"}
  if(length(which(rownames(tmp.keep)=="SBS7a"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS7a")] <- "SBS7a-Ultraviolet light exposure"}
  if(length(which(rownames(tmp.keep)=="SBS7b"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS7b")] <- "SBS7b-Ultraviolet light exposure"}
  if(length(which(rownames(tmp.keep)=="SBS7c"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS7c")] <- "SBS7c-Ultraviolet light exposure"}
  if(length(which(rownames(tmp.keep)=="SBS7d"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS7d")] <- "SBS7d-Ultraviolet light exposure"}
  if(length(which(rownames(tmp.keep)=="SBS8"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS8")] <- "SBS8-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS9"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS9")] <- "SBS9-Pol eta activity"}
  if(length(which(rownames(tmp.keep)=="SBS10a"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS10a")] <- "SBS10a-Pol eps mut."}
  if(length(which(rownames(tmp.keep)=="SBS10b"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS10b")] <- "SBS10b-Pol eps mut."}
  if(length(which(rownames(tmp.keep)=="SBS11"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS11")] <- "SBS11-Temozolomide treatment"}
  if(length(which(rownames(tmp.keep)=="SBS12"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS12")] <- "SBS12-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS13"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS13")] <- "SBS13-APOBEC activity"}
  if(length(which(rownames(tmp.keep)=="SBS14"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS14")] <- "SBS14-Defective MMR + Pol eps mut."}
  if(length(which(rownames(tmp.keep)=="SBS15"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS15")] <- "SBS15-Defective MMR"}
  if(length(which(rownames(tmp.keep)=="SBS16"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS16")] <- "SBS16-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS17a"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS17a")] <- "SBS17a-Unknown (maybe 5FU; ROS damage)"}
  if(length(which(rownames(tmp.keep)=="SBS17b"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS17b")] <- "SBS17b-Unknown (maybe 5FU; ROS damage)"}
  if(length(which(rownames(tmp.keep)=="SBS18"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS18")] <- "SBS18-Reactive oxygen species"}
  if(length(which(rownames(tmp.keep)=="SBS19"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS19")] <- "SBS19-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS20"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS20")] <- "SBS20-Defective MMR + POLD1 mut."}
  if(length(which(rownames(tmp.keep)=="SBS21"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS21")] <- "SBS21-Defective MMR"}
  if(length(which(rownames(tmp.keep)=="SBS22"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS22")] <- "SBS22-Aristolochic acid exposure"}
  if(length(which(rownames(tmp.keep)=="SBS23"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS23")] <- "SBS23-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS24"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS24")] <- "SBS24-Aflatoxin exposure"}
  if(length(which(rownames(tmp.keep)=="SBS25"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS25")] <- "SBS25-Unknown (Chemotherapy?)"}
  if(length(which(rownames(tmp.keep)=="SBS26"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS26")] <- "SBS26-Defective MMR"}
  if(length(which(rownames(tmp.keep)=="SBS27"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS27")] <- "SBS27-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS28"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS28")] <- "SBS28-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS29"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS29")] <- "SBS29-Tobacco chewing"}
  if(length(which(rownames(tmp.keep)=="SBS30"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS30")] <- "SBS30-Defective BER (inactivating mut. in NTHL1)"}
  if(length(which(rownames(tmp.keep)=="SBS31"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS31")] <- "SBS31-Platinum treatment"}
  if(length(which(rownames(tmp.keep)=="SBS32"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS32")] <- "SBS32-Azathioprine treatment"}
  if(length(which(rownames(tmp.keep)=="SBS33"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS33")] <- "SBS33-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS34"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS34")] <- "SBS34-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS35"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS35")] <- "SBS35-Platinum treatment"}
  if(length(which(rownames(tmp.keep)=="SBS36"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS36")] <- "SBS36-Defective BER: MUTYH mut./ROS"}
  if(length(which(rownames(tmp.keep)=="SBS37"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS37")] <- "SBS37-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS38"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS38")] <- "SBS38-Unknown (Indir. effect of UV light?)"}
  if(length(which(rownames(tmp.keep)=="SBS39"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS39")] <- "SBS39-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS40"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS40")] <- "SBS40-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS41"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS41")] <- "SBS41-Unknown"}
  if(length(which(rownames(tmp.keep)=="SBS42"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS42")] <- "SBS42-Haloalkane exposure"}
  if(length(which(rownames(tmp.keep)=="SBS43"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS43")] <- "SBS43-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS44"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS44")] <- "SBS44-Defective MMR"}
  if(length(which(rownames(tmp.keep)=="SBS45"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS45")] <- "SBS45-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS46"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS46")] <- "SBS46-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS47"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS47")] <- "SBS47-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS48"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS48")] <- "SBS48-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS49"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS49")] <- "SBS49-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS50"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS50")] <- "SBS50-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS51"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS51")] <- "SBS51-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS52"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS52")] <- "SBS52-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS53"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS53")] <- "SBS53-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS54"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS54")] <- "SBS54-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS55"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS55")] <- "SBS55-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS56"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS56")] <- "SBS56-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS57"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS57")] <- "SBS57-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS58"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS58")] <- "SBS58-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS59"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS59")] <- "SBS59-Possible sequencing artefacts"}
  if(length(which(rownames(tmp.keep)=="SBS60"))>0){rownames(tmp.keep)[which(rownames(tmp.keep)=="SBS60")] <- "SBS60-Possible sequencing artefacts"}
  if(length(grep("SBS",rownames(tmp.keep)))>0){
    tmp.keep <- t(tmp.keep)
  }
  return(tmp.keep)
}


############################
### Define colors
############################
defineColorsGBMsubset <- function(tmp.keep) {
  if(length(grep("SBS",colnames(tmp.keep)))>0){
    tmp.keep <- t(tmp.keep)
  }
  color = rep(NA, length=length(factor(colnames(tmp.keep))));
  color[which(factor(colnames(tmp.keep))=="SBS1 Deamination of 5MeC")] <- "#882255";
  color[which(factor(colnames(tmp.keep))=="SBS3 Defective HR repair")] <- "#DDCC77";
  color[which(factor(colnames(tmp.keep))=="SBS5 Unknown (clock-like)")] <- "#332288";
  color[which(factor(colnames(tmp.keep))=="SBS8 Unknown")] <- "#117733";
  color[which(factor(colnames(tmp.keep))=="SBS11 Temozolomide treatment")] <- "#88CCEE";
  color[which(factor(colnames(tmp.keep))=="SBS15 Defective MMR")] <- "#CC6677";
  color[which(factor(colnames(tmp.keep))=="SBS16 Unknown")] <- "#44AA99";
  color[which(factor(colnames(tmp.keep))=="SBS30 Defective BER")] <- "#999933";
  color[which(factor(colnames(tmp.keep))=="SBS40 Unknown")] <- "#AA4499";
  if(length(grep("SBS",rownames(tmp.keep)))>0){
    tmp.keep <- t(tmp.keep)
  }
  return(color)
}

############################
### Barplots with optional mutations count on top
############################
makeSigBarplot <- function(sigdf,thislabel=NULL,sampleorder=NULL,mutationscounts=NULL,labelcolor = "grey20"){
  
  #######################################
  # note on mutationscounts label option:
  #######################################
  ### the dataframe 'mutationscounts' needs to have the following columns: c("SampleName","N mutations")
  ### it can e.g. be generated from the 'sigs.input.RDS' from deconstructSigs with the following code:
  # sigs.input <- readRDS("sigs.input.RDS")
  # rwsms <- rowSums(sigs.input)
  # rwsms <- setNames(cbind(names(rwsms),as.data.frame(rwsms)), c("SampleName","N mutations"))
  
  sigdf <- as.data.frame(sigdf)
  sigdf$SampleName <- rownames(sigdf)
  sigdf <- melt(sigdf, id.var = "SampleName")
  
  if(!is.null(sampleorder)){
    sigdf$SampleName <- as.factor(sigdf$SampleName)
    sigdf$SampleName <- factor(sigdf$SampleName,levels = sampleorder)
  }
  
  if(!is.null(mutationscounts)){
    sigdf <- dplyr::left_join(sigdf, mutationscounts)
    sigdf[which(duplicated(sigdf$SampleName)==TRUE),"N mutations"] <- NA
  }
  
  p <- ggplot(sigdf, aes(x=SampleName, y=value, fill=variable))+
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width=0.8)+
    ggtitle(thislabel)+
    scale_fill_manual(values=color)+
    xlab("")+
    ylab("% signature contribution")+
    theme_classic()+
    theme(legend.direction="vertical",
          legend.title=element_blank(),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(colour = labelcolor,size=9,margin = margin(t=0,r=5,b=0,l=0,unit = "pt")),
          legend.key.height=grid::unit(0.5,"cm"),
          legend.key.width=grid::unit(0.5,"cm"),
          legend.position = "bottom",
          axis.text.x=element_text(size=9,angle = 90, hjust = 1,vjust=0.5, color = labelcolor),
          axis.text.y=element_text(size=9, vjust=0.2, color = labelcolor),
          axis.title.y = element_text(size=10, vjust=-9, color = labelcolor),
          axis.ticks=element_blank(),
          plot.background=element_blank(),
          plot.title = element_text(size = 9),
          panel.border=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_blank(),
          strip.background = element_blank())+
    guides(fill = guide_legend(ncol = 6))
  
  # if given, add the mutation counts on top of the bars
  if(!is.null(mutationscounts)){
   p <- p + 
     
     geom_text(
      aes(label = `N mutations`),
      #y=1.1,
      y=110,
      angle=90,
      parse = TRUE,
      size = 3.2,
      color = labelcolor) +
      scale_y_continuous(limits = c(0,110),breaks=seq(0,110,10),
      labels=c("0","10","20","30","40","50","60","70","80","90","100","N mutations"))
  }
  return(p)
}






