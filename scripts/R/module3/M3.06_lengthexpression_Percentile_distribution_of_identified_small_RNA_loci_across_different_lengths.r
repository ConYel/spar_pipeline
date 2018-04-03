lengthexpression_Percentile_distribution_of_identified_small_RNA_loci_across_different_lengths<-function(fprefix='input',wdir=".") {

start_time <- proc.time()
cat("Percentile distribution across all characterized and uncharacterized loci (M3.06) start", date(), "\n")
    
# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Expression characteristics  

# Module_3_Figure_6(Figure 3.06) 	
# Description: Percentile distribution across all characterized and uncharacterized loci
# input: input_annot.with_conservation.xls, input.unannot.final.with_conservation.xls
# output: Percentile_distribution_of_identified_small_RNA_loci_across_different_lengths.png
# output: Percentile_distribution_of_identified_small_RNA_loci_across_different_lengths.pdf 

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME="Percentile_distribution_of_identified_small_RNA_loci_across_different_lengths"
PLOTTITLE="Percentile distribution of all loci \n(normalized by each sncRNA class)"
XTITLE="Peak_length(nt)"
YTITLE="sncRNA_classes"

# libraries needed 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))

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

# Read in data (both characterized and uncharacterized) 
D = read.table(datafile1,sep='\t',header=T,comment.ch="")
E = read.table(datafile2,sep='\t',header=T,comment.ch="")
E$annotRNAclass = "unannot"
#DE = rbind(D[,c(1:20,25,29:31)],E) 
DE = rbind(D,E)
DE$Peak_length_nt= DE$peakChrEnd-DE$peakChrStart
DEX = DE[,c("Peak_length_nt","peakExprPercentile","annotRNAclass")]

Z = ddply(DEX,Peak_length_nt~annotRNAclass,summarise,sum=sum(peakExprPercentile))

# Read in data (characterized only) 
#D = read.table(datafile1,sep='\t',header=T,comment.ch="")
#DEX = D[,c("peakExprPercentile","annotRNAclass")]
#Z = ddply(DEX,Peak_length_nt~annotRNAclass,summarise,sum=sum(peakExprPercentile))

colnames(Z) = c(XTITLE,YTITLE,"sum_of_ranks")
Z = Z[which(Z[,XTITLE]<=44),]

zmelt = dcast(Z,Z[,2]~Z[,1])
zmelt[is.na(zmelt)] <- 0 

# col normlization 
#cs = colSums(zmelt[,-1])
#u_z = (t((t(zmelt[,-1])/cs)))*100

# row normalization
rs = rowSums(zmelt[,-1])
u_z = zmelt[,-1]/rs*100
u_z = data.frame(zmelt[,1],u_z)
colnames(u_z) = c(YTITLE,as.numeric(colnames(zmelt)[-1]))
final_count = melt(u_z,ID=1)
colnames(final_count)=c(YTITLE,XTITLE,"Norm_sum_of_rank")
final_count[,XTITLE]=as.factor(final_count[,XTITLE])
final_count$Norm_sum_of_rank=as.numeric(final_count$Norm_sum_of_rank)

levels(final_count[,1])[1] = paste(levels(final_count[,1])[1],"(",as.character(table(DEX$annotRNAclass)[1]),")")
levels(final_count[,1])[2] = paste(levels(final_count[,1])[2],"(",as.character(table(DEX$annotRNAclass)[2]),")")
levels(final_count[,1])[3] = paste(levels(final_count[,1])[3],"(",as.character(table(DEX$annotRNAclass)[3]),")")
levels(final_count[,1])[4] = paste(levels(final_count[,1])[4],"(",as.character(table(DEX$annotRNAclass)[4]),")")
levels(final_count[,1])[5] = paste(levels(final_count[,1])[5],"(",as.character(table(DEX$annotRNAclass)[5]),")")
levels(final_count[,1])[6] = paste(levels(final_count[,1])[6],"(",as.character(table(DEX$annotRNAclass)[6]),")")
levels(final_count[,1])[7] = paste(levels(final_count[,1])[7],"(",as.character(table(DEX$annotRNAclass)[7]),")")
levels(final_count[,1])[8] = paste(levels(final_count[,1])[8],"(",as.character(table(DEX$annotRNAclass)[8]),")")
levels(final_count[,1])[9] = paste(levels(final_count[,1])[9],"(",as.character(table(DEX$annotRNAclass)[9]),")")
levels(final_count[,1])[10] = paste(levels(final_count[,1])[10],"(",as.character(table(DEX$annotRNAclass)[10]),")")
levels(final_count[,1])[11] = paste(levels(final_count[,1])[11],"(",as.character(table(DEX$annotRNAclass)[11]),")")
levels(final_count[,1])[12] = paste(levels(final_count[,1])[12],"(",as.character(table(DEX$annotRNAclass)[12]),")")
levels(final_count[,1])[13] = paste(levels(final_count[,1])[13],"(",as.character(table(DEX$annotRNAclass)[13]),")")

# make plot 
png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
print(	ggplot(final_count,aes(x = final_count[,XTITLE], y = final_count[,YTITLE], fill = Norm_sum_of_rank))+
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ 
		geom_tile()+ 
		#scale_fill_gradientn(colours=brewer.pal(10,"Paired"),limits=c(0,100), breaks=seq(0,100,by=10),name = "Normalized RPM (%)")+
		scale_fill_gradientn(colours=c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"),limits=c(0,100), breaks=seq(0,100,by=10),name = "Normalized RPM (%)")+
		ggtitle(PLOTTITLE) + xlab("Peak Length")+ylab(YTITLE) +
		theme(legend.position="bottom",legend.direction="horizontal",legend.text=element_text(size=6),
		axis.text = element_text(size = 14),axis.text.x = element_text(angle = 90, hjust = 1),
		axis.title = element_text(size =18),plot.title = element_text(size = 18)))
dev.off()

cat("Total time for percentile distribution across all characterized and uncharacterized loci analysis (M3.06):", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")
}

# lengthexpression_Percentile_distribution_of_identified_small_RNA_loci_across_different_lengths(fprefix=fprefix,wdir=wdir)

