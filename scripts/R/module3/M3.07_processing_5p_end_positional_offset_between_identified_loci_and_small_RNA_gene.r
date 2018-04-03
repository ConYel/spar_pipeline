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
