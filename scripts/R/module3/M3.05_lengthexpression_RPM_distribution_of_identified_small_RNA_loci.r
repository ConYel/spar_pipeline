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
