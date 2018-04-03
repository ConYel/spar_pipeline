proportion_of_mapped_reads_across_all_loci<- function(fprefix='input',wdir=".") {

## A pie chart summarizes the counts of loci per class at the front analysis page 

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME="proportion_of_mapped_reads_across_all_loci"
PLOTTITLE="Proportion of mapped reads all loci"
XTITLE="sncRNA_classes"
YTITLE="Count"

# libraries needed 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(RColorBrewer))

#args<-commandArgs(TRUE)
#datafile=args[1] # input file
#wdir=args[2] # output / working directory
#if (length(args)<1) { stop("ERROR: No input! USAGE: script inputfile <output-dir>")}
#if (length(args)<2) { wdir="." } 
#use current dir if no 
#working dir has been specified

# output image file
pngfile= paste(wdir, "/", paste(BASENAME,".png",sep=""), sep="")
pdffile= paste(wdir, "/", paste(BASENAME,".pdf",sep=""), sep="")

# Read in data
D = read.table(datafile1,sep='\t',header=T,comment.ch="")
E = read.table(datafile2,sep='\t',header=T,comment.ch="")
E$annotChr = "NA"
E$annotChrStart = 0
E$annotChrEnd = 0
E$annotID = "NA"
E$annotRNAclass = "unannot"
E$annotStrand = "+"
E$annotOverlap = 0
E$peakOverlap = 0
E = E[,c(1:20,24:31,21:23)]
DE = rbind(D,E) 


DE$Loci_length= DE$peakChrEnd-DE$peakChrStart
final = DE[,c("Loci_length","annotRNAclass","peakExpressionValue")]
colnames(final) = c("Loci_length",XTITLE,"Expression_log10_RAW")
final_Count = final[which(final$Loci_length<=100),]

# prepare to include the 'position' information for text 
zmelt = melt(final_Count[,-1])
zmelt = dcast(zmelt,sncRNA_class~., value.var="value",fun.aggregate=sum)
colnames(zmelt)=c(XTITLE,YTITLE)
combine_f = ddply(zmelt, .(), transform, weight=Count/sum(Count))[,-1]
m2 = ddply(zmelt, .(), transform, position = cumsum(Count) - 0.5*Count) 
combine_f =  cbind(combine_f,m2[,4])
colnames(combine_f)[4] = c("position")


## specifying the parameters for plot 
png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
colourCount = dim(table(final_Count[,XTITLE]))
getPalette = colorRampPalette(brewer.pal(colourCount, "Spectral"))
options(scipen=10000)
print(ggplot(combine_f, aes(x="", y=combine_f[,YTITLE], fill=combine_f[,XTITLE]))+
geom_bar(width = 1, stat = "identity")+ 
coord_polar("y", start=0)+
scale_fill_manual(values = getPalette(colourCount),name=XTITLE)+
ggtitle(PLOTTITLE) + xlab("")+ylab("")+
geom_text(aes(label = paste(sprintf("%1.1f%%", round(combine_f$weight,3)*100),sep="-"), y = combine_f$position))+
theme(axis.text = element_text(size = 12),axis.title = element_text(size =14),
legend.position = "bottom",plot.title = element_text(size = 16)))
dev.off()


}

# proportion_of_mapped_reads_across_all_loci(fprefix=fprefix,wdir=wdir)
