wdir="test_out"
fprefix="peaks"
partition_ref_path="./annot/partition_files/hg19"
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
segmentation_Genomic_length_distribution_of_identified_small_RNA_loci<- function(fprefix='input',wdir=".") {
start_time <- proc.time()
cat("Length distribution of loci after segmentation (M3.02) start", date(), "\n")

# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Segmentation characteristics  

# Module_3_Figure_2(Figure 3.02) 	
# Description: Length distribution of loci after segmentation
# input: input_annot.with_conservation.xls, input.unannot.final.with_conservation.xls
# output: Genomic_length_distribution_of_identified_small_RNA_loci.png
# output: Genomic_length_distribution_of_identified_small_RNA_loci.pdf

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME="Genomic_length_distribution_of_identified_small_RNA_loci"
PLOTTITLE="Length distribution of all loci \n after segmentation"
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
DXS = table(DX$Peak_length)
DXS = data.frame(DXS)
DXS$norm = DXS$Freq/sum(DXS$Freq)*100 

# make plot 
png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
options(scipen=10000)
print(	ggplot(data = DXS, aes(x = Var1, y = norm)) + 
		geom_bar(stat="identity") + theme_classic()+
		ggtitle(PLOTTITLE) + xlab(XTITLE)+ylab(YTITLE)+
		ylim(0,30)+
		theme(text = element_text(size=16),axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)))
dev.off()

cat("Total time for length distribution of loci after segmentation (M3.02) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")
}

# segmentation_Genomic_length_distribution_of_identified_small_RNA_loci(fprefix=fprefix,wdir=wdir)
lengthexpression_RPM_distribution_of_identified_small_RNA_loci<- function(fprefix='input',wdir=".") {

start_time <- proc.time()
cat("Distribution of RPM values per sncRNA class (M3.05) started", date(), "\n")
    
# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Expression characteristics  

# Module_3_Figure_5(Figure 3.05) 	
# Description: Distribution of log 10 RPM values of all loci per sncRNA class
# input: input_annot.with_conservation.xls, input.unannot.final.with_conservation.xls
# output: RPM_distribution_of_identified_small_RNA_loci.png / RPM_distribution_of_identified_small_RNA_loci.pdf

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME = "RPM_distribution_of_identified_small_RNA_loci"
PLOTTITLE="RPM distribution of all loci \nacross all RNA class"
XTITLE="Log10(RPM)"
YTITLE="Density"

# libraries needed 
suppressPackageStartupMessages(library(ggplot2))
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

# Read in data
D = read.table(datafile1,sep='\t',header=T,comment.ch="")
E = read.table(datafile2,sep='\t',header=T,comment.ch="")
E$annotRNAclass = "unannot"
#DE = rbind(D[,c(1:20,25,29:31)],E) 
DE = rbind(D,E)
DE$annotRNAclass2 = DE$annotRNAclass 
levels(DE$annotRNAclass2)[1:4]<-"miRNA"
levels(DE$annotRNAclass2)[5:6]<-"snoRNA"

## making necessary quantiles 
u = log10(DE$peakRPM) 
quantiles1 = quantile(u, probs=seq(0,1,.1)); quantiles1 = quantiles1[c(2,10)]
quantiles2 = quantile(u, probs=seq(0,1,.25)); quantiles2 = quantiles2[c(2:4)]
quantile_line = c(quantiles1, quantiles2)[c(2,4,5)]; quantile_line = sort(quantile_line); quantile_line = data.frame(t(quantile_line))


p = ggplot(DE, aes(x=log10(DE$peakRPM))) + geom_density(aes(group=DE$annotRNAclass2))
py_plot = max(ggplot_build(p)$data[[1]]$y)


## change yaxis axis + title names 
## p$labels$colour="Categories"		

#scale_fill_manual(values = getPalette(colourCount))
#scale_fill_manual 

#ggplot(D, aes(log10(peakRPM), ..count.., fill = annotRNAclass2))+
#geom_density(position = "fill")+
#scale_fill_brewer(palette="Set2")+
#geom_vline(data = quantile_line,aes(xintercept = quantile_line$X10.),color="black",linetype = "dashed",size=1.25)+
#geom_vline(data = quantile_line,aes(xintercept = quantile_line$X90.),color="black",linetype = "dashed",size=1.25)

# make plot 
int = sort.int(ggplot_build(p)$data[[1]]$x, index.return=TRUE) 
front <- data.frame(position = quantile_line) 

front = rbind(quantile_line,rep(py_plot*0.75,5))
front = t(front)
colnames(front)=c("x","y")
front = data.frame(front)

#rownames(front)[1] = paste("10%, log10(RPM)=",round(quantile_line$X10,2))
#rownames(front)[2] = paste("25%, log10(RPM)=",round(quantile_line$X25,2))
#rownames(front)[3] = paste("50%, log10(RPM)=",round(quantile_line$X50,2))
#rownames(front)[4] = paste("75%, log10(RPM)=",round(quantile_line$X75,2))
#rownames(front)[5] = paste("90%, log10(RPM)=",round(quantile_line$X90,2))

rownames(front)[1] = paste("50%, log10(RPM)=",round(quantile_line$X50,2))
rownames(front)[2] = paste("75%, log10(RPM)=",round(quantile_line$X75,2))
rownames(front)[3] = paste("90%, log10(RPM)=",round(quantile_line$X90,2))


png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
print(ggplot(DE, aes(x=log10(DE$peakRPM),colour=as.factor(DE$annotRNAclass2))) + 
		geom_density(aes(group=DE$annotRNAclass2),size=1.5)+
		scale_color_brewer(palette = "Spectral",name = "sncRNA classes")+ 
		coord_cartesian(ylim=c(0,py_plot+0.5))+ theme_classic()+
		ggtitle(PLOTTITLE) + xlab(XTITLE)+ylab(YTITLE)+
		geom_vline(xintercept = as.numeric(quantile_line),color="black",linetype = "dashed",size=1.25)+
		geom_text(data=front,aes(x=x-0.1,y=y+0.5,label=rownames(front)), colour="blue", angle=90,size=4)+
		theme(legend.position="bottom")
		)
dev.off()

cat("Total time for RPM per sncRNA class (M3.05) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")
}
# RPM_distribution_of_identified_small_RNA_loci(fprefix=fprefix,wdir=wdir)
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

processing_5p_end_positional_offset_between_identified_loci_and_small_RNA_gene<- function(fprefix='input',wdir=".") {

start_time <- proc.time()
cat("Offset between loci and genes at 5' end (M3.07) start", date(), "\n")
    
# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Processing characteristics  

# Module_3_Figure_7 (Figure 3.07) 	
# Description: Normalized offset between loci and genes at the 5' end (displaying relative positions of loci on the genes) 
# input: input_annot.with_conservation.xls
# output: 5p_end_positional_offset_between_identified_loci_and_small_RNA_gene.png
# output: 5p_end_positional_offset_between_identified_loci_and_small_RNA_gene.pdf 

# parameters for the plot 
datafile=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
BASENAME="5p_end_positional_offset_between_identified_loci_and_small_RNA_gene"
PLOTTITLE="Normalized offset at 5' end, between annotated loci and sncRNA gene"
XTITLE="Normalized_offset"
YTITLE="Number_of_loci"

# libraries needed 
suppressPackageStartupMessages(library(ggplot2))

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
D = read.table(datafile,sep='\t',header=T,comment.ch="")
D_pos = subset(D,peakStrand=="+")
D_neg = subset(D,peakStrand=="-")

# for positive strand (peak start - annotation start) 
D_pos$Offset_5p = D_pos$peakChrStart-D_pos$annotChrStart + D_pos$peakMostCommon5pPosition - 1
D_pos$Loci_length = D_pos$peakChrEnd-D_pos$peakChrStart 
D_pos$Annot_length = D_pos$annotChrEnd-D_pos$annotChrStart 

# for negative strand (peak end - annotation end)
D_neg$Offset_5p = D_neg$annotChrEnd-D_neg$peakChrEnd + D_neg$peakMostCommon5pPosition - 1
D_neg$Loci_length = D_neg$peakChrEnd-D_neg$peakChrStart 
D_neg$Annot_length = D_neg$annotChrEnd-D_neg$annotChrStart 

D_comb = rbind(D_pos,D_neg)
A_final = D_comb[,c("annotRNAclass","annotStrand","Offset_5p","Loci_length","Annot_length")]
colnames(A_final) = c("sncRNA_class", "Strand", "Offset_5p","Loci_length","Annot_length")

A_final = A_final[which(A_final$Loci_length<=44),]
A_final[,6] = A_final[,3]/A_final[,5]
colnames(A_final)[6] = XTITLE

# hist(mir_3p[,6],breaks=dim(mir_3p)[1],xlim=c(-0.2,0.2))

# A histogram
png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
options(scipen=10000)
print(ggplot(A_final, aes(x=Normalized_offset))+ 
geom_histogram(colour="white", breaks=seq(round(min(A_final[,6])-0.05,1), 1, by= 0.05)+0.025,right = TRUE)+
ggtitle(PLOTTITLE)+theme_classic()+
theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 8))+
facet_grid(sncRNA_class~.,scales="free_y"))
dev.off()

cat("Total time for offset between loci and genes at 5' end (M3.07) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")
}

# processing_5p_end_positional_offset_between_identified_loci_and_small_RNA_gene(fprefix=fprefix,wdir=wdir)
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
processing_Read_pipeup_at_5p_end_of_identified_small_RNA_loci<- function(fprefix='input',wdir=".") {

start_time <- proc.time()
cat("Read pileup at 5' end of segmented loci (M3.09) start", date(), "\n")    
    
# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Segmentation characteristics  

# Module_3_Figure_9(Figure 3.09) 	
# Description: Read pipeup at the 5p end point of segmented loci across all sncRNA classes
# input: input_annot.with_conservation.xls, input.unannot.final.with_conservation.xls
# output: Read_pipeup_at_5p_end_of_identified_small_RNA_loci.png / Read_pipeup_at_5p_end_of_identified_small_RNA_loci.pdf

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME="Read_pipeup_at_5p_end_of_identified_small_RNA_loci"
PLOTTITLE="Normalized entropy at 5p end of loci"
XTITLE="sncRNA_classes"
YTITLE="Normalized_peak_entropy_5p_read_end"

# libraries needed 
#suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(plyr))
library(reshape2)
library(ggplot2)
library(plyr)
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
print(pngfile)
flush.console()

# Read in data
D = read.table(datafile1,sep='\t',header=T,comment.ch="")
E = read.table(datafile2,sep='\t',header=T,comment.ch="")
E$annotRNAclass = "unannot"
#DE = rbind(D[,c(1:20,25,29:31)],E) 
DE = rbind(D,E)

DE$Peak_length = DE$peakChrEnd-DE$peakChrStart
DX = DE[DE$Peak_length<=44,]
DXE = DE[,c("annotRNAclass","peakProportionOfReadsAtMostCommon5pPosition")]
DXE_melt = melt(DXE,id=1) 
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
print(ggplot(DXE_melt, aes(x=DXE_melt[,XTITLE], y=DXE_melt[,YTITLE])) + 
		stat_summary(fun.data=data_summary) + theme_classic()+ 
		# stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.5)+
		# geom_boxplot(width=.1)
		geom_boxplot(outlier.shape = NA,position=position_dodge(0.75))+
		ggtitle(PLOTTITLE) + xlab(XTITLE)+ylab(YTITLE)+
		ylim(0,1)+
		theme(text = element_text(size=20),axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)))
#ylim1=boxplot.stats(DXE_melt[,YTITLE])$stats[c(1, 5)]

#print(p0+coord_cartesian(ylim = ylim1*1.05))
dev.off()

cat("Total time for read pileup at 5' end of loci (M3.09) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")

}

# processing_Read_pipeup_at_5p_end_of_identified_small_RNA_loci(fprefix=fprefix,wdir=wdir)
genomewide_Proportion_of_expressed_annotated_small_RNA_genes<- function(fprefix='input',wdir=".") {

start_time <- proc.time()
cat("Proportion of characterized RNA loci per class (M3.10) start", date(), "\n")    
    
# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Genome wide characteristics  

# Module_3_Figure_10 (Figure 3.10) 	
# Description: Proportion of characterized RNA loci per sncRNA class 
# input: input_annot.with_conservation.xls
# output: Proportion_of_expressed_annotated_small_RNA_genes.png / Proportion_of_expressed_annotated_small_RNA_genes.pdf

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
BASENAME="Proportion_of_expressed_annotated_small_RNA_genes"
PLOTTITLE="Percentage of characterized RNA loci in dataset \n as compared to annotations"
XTITLE="sncRNA_classes"
YTITLE="Percentage"

# libraries needed 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(plyr))
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
pngfile= paste(wdir, "/../figures/", paste(BASENAME,".png",sep=""), sep="")
#pdffile= paste(wdir, "/", paste(BASENAME,".pdf",sep=""), sep="")

# Read in data
D = read.table(datafile1,sep='\t',header=T,comment.ch="")
df <- data.frame(matrix(ncol = 13, nrow = 1))
df[1,] = c("645","653","926","1582","9205","3","5","398","718","12","625","625","625")
colnames(df)=c("mir-3p","mir-5p","mir-5p3pno","miRNAprimary","piRNA","rRNA","scRNA","snoRNA","snoRNAnar","snRNA","tRF3","tRF5","tRNA")

c = D[,c("annotID","annotRNAclass")]
z = transform(c, annotID = colsplit(annotID, pattern  = "\\:.:", names = c('chr', 'name')))
z = z[!duplicated(z$annotID[,2]),]

count = table(z$annotRNAclass)
tb = rbind(count,df)
tbp = as.numeric(tb[1,])/as.numeric(tb[2,])
#tbp = tbp-0.01
tbp_comp = 1-tbp
tbp_plot = rbind(tbp,tbp_comp)*99.8
rownames(tbp_plot)=c("Present in data","Overall annotations")
colnames(tbp_plot)=c("mir-3p","mir-5p","mir-5p3pno","miRNAprimary","piRNA","rRNA","scRNA","snoRNA","snoRNAnar","snRNA","tRF3","tRF5","tRNA")
colSums(tbp_plot)


melt_tbp_plot = melt(tbp_plot,ID=1)
colnames(melt_tbp_plot) = c("Category",XTITLE,YTITLE)

# make plot 
png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
print(	ggplot(melt_tbp_plot,aes(x = melt_tbp_plot[,XTITLE], y = melt_tbp_plot[,YTITLE], fill = Category))+
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ 
		geom_bar(stat = "identity")+scale_fill_grey()+scale_y_continuous(limits = c(0, 100.0))+theme_classic()+
		ggtitle(PLOTTITLE) + xlab(XTITLE)+ylab(YTITLE) +
		theme(legend.position="bottom",legend.text=element_text(size=10),
		axis.text = element_text(size = 14),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
		axis.title = element_text(size =18),plot.title = element_text(size = 18)))
dev.off()

cat("Total time for proportion of characterized loci per class (M3.10) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")

}

# genomewide_Proportion_of_expressed_annotated_small_RNA_genes(fprefix=fprefix,wdir=wdir)
## Genomewide_distribution_patterns_of_small_RNA_loci.r
## alex amlie-wolf
## computes genomic partition of all SPAR-called loci (annot and unannot)

genomewide_Genomewide_distribution_patterns_of_small_RNA_loci <- function(fprefix='input',wdir=".",partition_ref_path="./partition_files/") {
    start_time <- proc.time()
    cat("Genomic partition analysis (M3.11) start", date(), '\n')
    
    BASENAME="Genomewide_distribution_patterns_of_small_RNA_loci"
    PLOTTITLE="Genomic partitioning of all loci\nacross sncRNA classes"
    XTITLE="Genomic_Element"
    YTITLE="Count"

    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(RColorBrewer))
    suppressPackageStartupMessages(library(parallel))
    
    ## output image file
    pngfile= paste(wdir, "/../figures/", paste(BASENAME,".png",sep=""), sep="")
#    pdffile= paste(wdir, "/", paste(BASENAME,".pdf",sep=""), sep="")

    ## first make the mRNA partition proportion plot   
    load(paste0(wdir, "/", BASENAME, "_mRNA_partition_props_byclass.Rdata"))

    partition_order <- c("5' UTR Exon", "5' UTR Intron", "3' UTR Exon", "3' UTR Intron",
                         "Promoter", "mRNA Exon", "mRNA Intron", "Intergenic")    
    
    png(pngfile, width = 7, height = 7, units = 'in', res = 300, type="cairo")
    print(ggplot(full_mRNA_class_prop_tab, aes(x=class, y=Proportion, fill=factor(Genomic_Element, levels=partition_order))) +
          xlab("sncRNA Class") + ylab("Genomic Partition Proportions") +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          scale_fill_manual(values=setNames(brewer.pal(n=9, name="Paired"), partition_order), name=XTITLE) +
          theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1)) +
          guides(fill=guide_legend(nrow=3, byrow=TRUE)) + 
          ggtitle(paste("Genomic Proportions of annotated and unannotated loci")) #, BASENAME))
          )
    dev.off()
    
    # ggsave(paste0(wdir, "/", BASENAME, "_partition_byclass.png"), class_barplots, dpi=600, type="cairo-png")
    # ggsave(paste0(wdir, "/", BASENAME, "_partition_byclass.pdf"), class_barplots, dpi=600)

    cat("Total time for genomic partition (M3.11) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")
}

# genomewide_Genomewide_distribution_patterns_of_small_RNA_loci(fprefix=fprefix,wdir=wdir, partition_ref_path=partition_ref_path)




genomewide_Genomewide_distribution_of_expressed_small_RNA_loci<- function(fprefix='input',wdir=".") {

start_time <- proc.time()
cat("Genomewide scatter (Philadelphia) plot (M3.13) start", date(), "\n")    
    
# Plots for Module 3 Analysis and visualization of SPAR output 
# Session: Genomewide characteristics  

# Module_3_Figure_13(Figure 3.13) 	
# Description: Genomewide scatter plot on RPM values (log10) for all loci 
# input: input_annot.with_conservation.xls, input.unannot.final.with_conservation.xls
# output: Genomewide_distribution_of_expressed_small_RNA_loci.png / Genomewide_distribution_of_expressed_small_RNA_loci.pdf

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME="Genomewide_distribution_of_expressed_small_RNA_loci"
PLOTTITLE="Genome wide scatter plot \nof RPM of all loci"
XTITLE="Chromosome Number"
YTITLE="Log10(RPM)"

# libraries needed 
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(scales))

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
DE$annotRNAclass2 = DE$annotRNAclass 
DE$annotRNAclass2 = revalue(DE$annotRNAclass2, c("mir-3p"="miRNA", "mir-5p"="miRNA", "mir-5p3pno"="miRNA", "miRNAprimary"="miRNA", "snoRNAnar"="snoRNA", "tRF3"="tRF","tRF5"="tRF"))

DE$plotchr = DE$X.peakChr
DE$plotchr = revalue(DE$plotchr,c("chr1"=1,"chr2"=2,"chr3"=3,"chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chr21"=21,"chr22"=22,"chrX"="X","chrY"="Y","chrM"="M"))
DE = DE[nchar(as.character(DE$plotchr))<=2,]

DE$plotchr = droplevels(DE$plotchr)
         
DE$plotpos = ifelse(DE$plotchr==1, DE$peakChrStart,
ifelse(DE$plotchr==2, 249250621 +DE$peakChrStart,
ifelse(DE$plotchr==3, 492449994 +DE$peakChrStart,
ifelse(DE$plotchr==4, 690472424 +DE$peakChrStart,
ifelse(DE$plotchr==5, 881626700 +DE$peakChrStart,
ifelse(DE$plotchr==6, 1062541960 +DE$peakChrStart,
ifelse(DE$plotchr==7, 1233657027 +DE$peakChrStart,
ifelse(DE$plotchr==8, 1392795690 +DE$peakChrStart,
ifelse(DE$plotchr==9, 1548066250 +DE$peakChrStart,
ifelse(DE$plotchr==10, 1694430272 +DE$peakChrStart,
ifelse(DE$plotchr==11, 1835643703 +DE$peakChrStart,
ifelse(DE$plotchr==12, 1971178450 +DE$peakChrStart,
ifelse(DE$plotchr==13, 2106184966 +DE$peakChrStart,
ifelse(DE$plotchr==14, 2240036861 +DE$peakChrStart,
ifelse(DE$plotchr==15, 2355206739 +DE$peakChrStart,
ifelse(DE$plotchr==16, 2462556279 +DE$peakChrStart,
ifelse(DE$plotchr==17, 2565087671 +DE$peakChrStart,
ifelse(DE$plotchr==18, 2655442424 +DE$peakChrStart,
ifelse(DE$plotchr==19, 2736637634 +DE$peakChrStart,
ifelse(DE$plotchr==20, 2814714882 +DE$peakChrStart,
ifelse(DE$plotchr==21, 2877740402 +DE$peakChrStart,
ifelse(DE$plotchr==22, 2937113968 +DE$peakChrStart,
ifelse(DE$plotchr=="X", 2996242951 +DE$peakChrStart,
ifelse(DE$plotchr=="Y", 3047547517 +DE$peakChrStart,
ifelse(DE$plotchr=="M", 3095677412 +DE$peakChrStart,
DE$plotchr)))))))))))))))))))))))))


DE$plotchr[DE$plotchr=="23"]<-"X"
DE$plotchr[DE$plotchr=="24"]<-"Y"
DE$plotchr[DE$plotchr=="25"]<-"M"

ChrLoc = c(1,249250621,492449994,690472424,881626700,1062541960,1233657027,1392795690,1548066250,1694430272,1835643703,1971178450,2106184966,2240036861,2355206739,2462556279,2565087671,2655442424,2736637634,2814714882,2877740402,2937113968,2996242951,3047547517,3095677412)
DE2=DE[order(DE$plotpos),]

# DE$plotpos = revalue(DE$plotchr,c("2"=249250622,"3"=492449995,"4"=690472425,"5"=881626701,"6"=1062541961,"7"=1233657028,"8"=1392795691,"9"=1548066251,"10"=1694430273,"11"=1835643704,"12"=1971178451,"13"=2106184967,"14"=2240036862,"15"=2355206740,"16"=2462556280,"17"=2565087672,"18"=2655442425,"19"=2736637635,"20"=2814714883,"21"=2877740403,"22"=2937113969,"X"=2996242952,"Y"=3047547518,"M"=3095677413))

png(pngfile,width = 14, height = 7, units = 'in', res = 300, type="cairo")

print(ggplot(data = DE2, aes(x = plotpos, y = log10(peakRPM))) + 
		geom_point(aes(colour=annotRNAclass2)) + facet_wrap( ~ annotRNAclass2, shrink=FALSE)+
		scale_x_continuous(breaks=ChrLoc, labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","M"), expand=c(0.02,0))+
		ggtitle(PLOTTITLE) + xlab(XTITLE)+ylab(YTITLE)+theme_classic()+
		scale_color_discrete(name = "sncRNA classes")+
		theme(legend.position="bottom",legend.text=element_text(size=10),
                      axis.text = element_text(size = 7),
                      axis.text.x = element_text(size=6, angle = 90, hjust=1, vjust=0.5),
                      axis.title = element_text(size =14),plot.title = element_text(size = 14)))
dev.off()

cat("Total time for RPM scatter/Phliadelphia plot (M3.13) analysis:", (proc.time() - start_time)[['elapsed']], "seconds (", date(), ")\n")

}
# genomewide_Genomewide_distribution_of_expressed_small_RNA_loci(fprefix=fprefix,wdir=wdir)



proportion_of_mapped_reads_across_all_loci<- function(fprefix='input',wdir=".") {

## A pie chart summarizes the counts of loci per class at the front analysis page 

# parameters for the plot 
datafile1=paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
datafile2=paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")
BASENAME="Proportion_of_mapped_reads_across_all_loci"
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
pngfile= paste(wdir, "/../figures/", paste(BASENAME,".png",sep=""), sep="")
pdffile= paste(wdir, "/../figures/", paste(BASENAME,".pdf",sep=""), sep="")

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
#E = E[,c(1:20,24:31,21:23)]
DE = rbind(D,E) 


DE$Loci_length= DE$peakChrEnd-DE$peakChrStart
final = DE[,c("Loci_length","annotRNAclass","peakExpressionValue")]
colnames(final) = c("Loci_length",XTITLE,"Expression_log10_RAW")
final_Count = final[which(final$Loci_length<=100),]

# prepare to include the 'position' information for text 
zmelt = melt(final_Count[,-1])
zmelt = dcast(zmelt,sncRNA_classes~., value.var="value",fun.aggregate=sum)
colnames(zmelt)=c(XTITLE,YTITLE)
combine_f = ddply(zmelt, .(), transform, weight=Count/sum(Count))[,-1]
m2 = ddply(zmelt, .(), transform, position = cumsum(Count) - 0.5*Count) 
combine_f =  cbind(combine_f,m2[,4])
colnames(combine_f)[4] = c("position")
combine_f[,XTITLE] = paste(combine_f[,XTITLE], sprintf("(%1.1f%%)", round(combine_f$weight,3)*100))

## specifying the parameters for plot 
png(pngfile,width = 7, height = 7, units = 'in', res = 300, type="cairo")
colourCount = dim(table(final_Count[,XTITLE]))
getPalette = colorRampPalette(brewer.pal(colourCount, "Paired"))
options(scipen=100000)
print(ggplot(combine_f, aes(x="", y=combine_f[,YTITLE], fill=combine_f[,XTITLE]))+
geom_bar(width = 1, stat = "identity")+ 
coord_polar("y", start=0)+
scale_fill_manual(values = getPalette(colourCount),name=XTITLE)+
ggtitle(PLOTTITLE) + xlab("")+ylab("")+
## geom_text(aes(label = paste(sprintf("%1.1f%%", round(combine_f$weight,3)*100),sep="-"), y = combine_f$position))+
theme(axis.text = element_text(size = 10),axis.title = element_text(size =12),
legend.position = "bottom",plot.title = element_text(size = 16), panel.background = element_blank()))
dev.off()


}

# proportion_of_mapped_reads_across_all_loci(fprefix=fprefix,wdir=wdir)

## these are our arguments:
#fprepix <- "fname"
#wdir <- "."
#partition_ref_path <- "./partition_files/"

library(parallel)

## these are our functions we want to run
func_list_list <- list(
    ## circos plot first since it's the slowest
    ##list("genomewide_Circular_genome_data_visualization",list(fprefix=fprefix,wdir=wdir)),
    ## then genomic partition
    list("genomewide_Genomewide_distribution_patterns_of_small_RNA_loci",list(fprefix=fprefix,wdir=wdir,partition_ref_path=partition_ref_path)),
    ## then the philadelphia plot
    list("genomewide_Genomewide_distribution_of_expressed_small_RNA_loci",list(fprefix=fprefix,wdir=wdir)),
    ## Segmentation characteristics 
    list("segmentation_Length_Distribution_of_RPM_of_segmented_loci",list(fprefix=fprefix,wdir=wdir)),
    list("segmentation_Genomic_length_distribution_of_identified_small_RNA_loci",list(fprefix=fprefix,wdir=wdir)),
    ## Length and Expression characteristics
    ## list("lengthexpression_RPM_and_length_of_all_loci_across_RNA_classes",list(fprefix=fprefix,wdir=wdir)),
    ## list("lengthexpression_proportion_of_mapped_reads_all_loci",list(fprefix=fprefix,wdir=wdir)),
    list("lengthexpression_RPM_distribution_of_identified_small_RNA_loci",list(fprefix=fprefix,wdir=wdir)),
    list("lengthexpression_Percentile_distribution_of_identified_small_RNA_loci_across_different_lengths",list(fprefix=fprefix,wdir=wdir)),
    ## Processing characteristics
    list("processing_5p_end_positional_offset_between_identified_loci_and_small_RNA_gene",list(fprefix=fprefix,wdir=wdir)),
    list("processing_Processing_specificity_at_5p_end_of_identified_small_RNA_loci",list(fprefix=fprefix,wdir=wdir)),
    list("processing_Read_pipeup_at_5p_end_of_identified_small_RNA_loci",list(fprefix=fprefix,wdir=wdir)),
    ## Genomewide characteristics 
    list("proportion_of_mapped_reads_across_all_loci",list(fprefix=fprefix,wdir=wdir)),
    list("genomewide_Proportion_of_expressed_annotated_small_RNA_genes",list(fprefix=fprefix,wdir=wdir))
    )

## source everything (these cant have function calls in them)
r_files <- list.files(path="SPAR-master/scripts/R/module3",pattern="^M.+\\.r$",full.names=TRUE)
lapply(r_files, source)

## make our cluster
clust <- makeCluster(detectCores(), outfile=paste0(wdir, "/../logs/parallel_r_output.log"))
clusterExport(cl=clust, varlist=unlist(lapply(func_list_list, "[[", 1)))

## parLapply(clust, func_list_list,
##           function(func_list) { func <- func_list[[1]]
##                                 func_args <- func_list[[2]]
##                                 do.call(func, func_args) })

clusterApply(clust, func_list_list,
             function(func_list) { func <- func_list[[1]]
                                   func_args <- func_list[[2]]
                                   do.call(func, func_args) })

stopCluster(clust)
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
