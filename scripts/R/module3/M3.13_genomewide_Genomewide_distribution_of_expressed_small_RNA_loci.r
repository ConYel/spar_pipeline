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



