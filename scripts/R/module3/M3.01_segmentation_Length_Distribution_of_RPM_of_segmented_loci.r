segmentation_Length_Distribution_of_RPM_of_segmented_loci<- function(fprefix='input',wdir=".") {
start_time <- proc.time()
cat("Length distribution of loci RPM across lengths (M3.01) start", date(), "\n")
    
# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Segmentation characteristics  

# Module_3_Figure_1(Figure 3.01)
# Description: Distribution of loci RPM across different lengths after segmentation
# input: input_annot.with_conservation.xls, input.unannot.final.with_conservation.xls
# output: Length_Distribution_of_RPM_of_segmented_loci.png / Length_Distribution_of_RPM_of_segmented_loci.pdf

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME="Length_Distribution_of_RPM_of_segmented_loci"
PLOTTITLE="Distribution of loci RPM\n at different lengths after segmentation"
XTITLE="Peak_length(nt)"
YTITLE="Percentage"

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
pngfile = paste(wdir, "/../figures/", paste(BASENAME,".png",sep=""), sep="")
#pdffile= paste(wdir, "/", paste(BASENAME,".pdf",sep=""), sep="")

# Read in data
D = read.table(datafile1,sep='\t',header=T,comment.ch="")
E = read.table(datafile2,sep='\t',header=T,comment.ch="")
E$annotRNAclass = "unannot"
#DE = rbind(D[,c(1:20,25,29:31)],E) 
DE = rbind(D,E)

DE$Peak_length = DE$peakChrEnd-DE$peakChrStart
DX = DE[DE$Peak_length<=44,]
DXT = dcast(DX, DX$annotRNAclass ~ DX$Peak_length,sum)
DXT = colSums(DXT[,-1]) 
DXT = data.frame(DXT)
DXT$norm = DXT$DXT/sum(DXT$DXT)*100 
DXT$Var1 = rownames(DXT)

# make plot 
png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
options(scipen=10000)
print(	ggplot(data = DXT, aes(x = Var1, y= norm)) + 
		geom_bar(stat="identity") + theme_classic()+
		ggtitle(PLOTTITLE) + xlab(XTITLE)+ylab(YTITLE)+
		ylim(0,30)+
		theme(text = element_text(size=16),axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)))
dev.off()

cat("Total time for length distribution of loci RPM across lengths (M3.01) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")
}

# segmentation_Length_Distribution_of_RPM_of_segmented_loci(fprefix=fprefix,wdir=wdir)
