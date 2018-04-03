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
