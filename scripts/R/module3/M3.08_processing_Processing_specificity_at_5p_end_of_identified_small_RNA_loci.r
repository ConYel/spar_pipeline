processing_Processing_specificity_at_5p_end_of_identified_small_RNA_loci<- function(fprefix='input',wdir=".") {

start_time <- proc.time()
cat("Processing specificity at 5' end of all loci (M3.08) start", date(), "\n")
    
# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Segmentation characteristics  

# Module_3_Figure_8(Figure 3.08) 	
# Description: Processing specificity at the 5'end of all loci per sncRNA classes 
# input: input_annot.with_conservation.xls, input.unannot.final.with_conservation.xls
# output: Processing_specificity_at_5p_end_of_identified_small_RNA_loci.png
# output: Processing_specificity_at_5p_end_of_identified_small_RNA_loci.pdf  

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME="Processing_specificity_at_5p_end_of_identified_small_RNA_loci"
PLOTTITLE="Processing specificity at 5p end"
XTITLE="sncRNA_classes"
YTITLE="Normalized_processing_specificity"

# libraries needed 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))

#args<-commandArgs(TRUE)
#datafile1=args[1] # input file 1
#datafile2=args[2] # input file 2
#wdir=args[3] # output / working directory
#if (length(args)<1) { stop("ERROR: No input! USAGE: script inputfile <output-dir>")}
#if (length(args)<2) { stop("ERROR: No input! USAGE: script inputfile <output-dir>")}
#if (length(args)<3) { wdir="." } 
#use current dir if no 
#working dir has been specified

# output image file
pngfile= paste(wdir, "/../figures/", paste(BASENAME,".png",sep=""), sep="")
#pdffile= paste(wdir, "/", paste(BASENAME,".pdf",sep=""), sep="")

# Read in data
D = read.table(datafile1,sep='\t',header=T,comment.ch="")
E = read.table(datafile2,sep='\t',header=T,comment.ch="")
E$annotRNAclass = "unannot"
#DE = rbind(D[,c(1:20,25,29:31)],E) 
DE = rbind(D,E)

DE$Peak_length = DE$peakChrEnd-DE$peakChrStart
DX = DE[DE$Peak_length<=44,]
DXE = DE[,c("annotRNAclass","peakNormalizedEntropy5p")]
DXE_melt = melt(DXE,id=1) 
DXE_melt[,3]=1-DXE_melt[,3]
colnames(DXE_melt)= c(XTITLE,"Category",YTITLE)

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}


# make plot 
png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
options(scipen=10000)
print(	ggplot(DXE_melt, aes(x=DXE_melt[,XTITLE], y=DXE_melt[,YTITLE])) + 
		#geom_violin() + 
		stat_summary(fun.data=data_summary) + theme_classic()+
		# stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.5)+
		geom_boxplot(outlier.shape = NA,position=position_dodge(0.75))+
		ggtitle(PLOTTITLE) + xlab(XTITLE)+ylab(YTITLE)+
		ylim(0,1)+
		theme(text = element_text(size=20),axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)))
dev.off()

cat("Total time for processing specificity at 5' end of loci (M3.08) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")

}

# processing_Processing_specificity_at_5p_end_of_identified_small_RNA_loci(fprefix=fprefix,wdir=wdir)
