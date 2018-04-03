## combine_partition_results.R
## alex amlie-wolf

args <- commandArgs(trailingOnly=TRUE) 
fprefix <- args[1]
wdir <- args[2]

suppressPackageStartupMessages(library(plyr))

BASENAME <- "Genomewide_distribution_patterns_of_small_RNA_loci"

## read in all the partition results and merge them together
annot_mRNA_genomic_partition <- read.table(paste0(wdir, "/", BASENAME, "_annot_mRNA_genomic_partition.txt"), header=T, sep="\t", quote="", as.is=T)
unannot_mRNA_genomic_partition <- read.table(paste0(wdir, "/", BASENAME, "_unannot_mRNA_genomic_partition.txt"), header=T, sep="\t", quote="", as.is=T)

annot_lncRNA_genomic_partition <- read.table(paste0(wdir, "/", BASENAME, "_annot_lncRNA_genomic_partition.txt"), header=T, sep="\t", quote="", as.is=T)
unannot_lncRNA_genomic_partition <- read.table(paste0(wdir, "/", BASENAME, "_unannot_lncRNA_genomic_partition.txt"), header=T, sep="\t", quote="", as.is=T)

annot_repeat_genomic_partition <- read.table(paste0(wdir, "/", BASENAME, "_annot_repeat_genomic_partition.txt"), header=T, sep="\t", quote="", as.is=T)
unannot_repeat_genomic_partition <- read.table(paste0(wdir, "/", BASENAME, "_unannot_repeat_genomic_partition.txt"), header=T, sep="\t", quote="", as.is=T)

## first merge the mRNA and lncRNA results; they should have the same number of rows
## the first 31 columns are from the input locus file
annot_mRNA_lncRNA_partition <- merge(annot_mRNA_genomic_partition, annot_lncRNA_genomic_partition, by=1:31, suffixes=c(".mRNA", ".lncRNA"))
unannot_mRNA_lncRNA_partition <- merge(unannot_mRNA_genomic_partition, unannot_lncRNA_genomic_partition, by=1:31, suffixes=c(".mRNA", ".lncRNA"))

## next, parse the repeat results so we have a summary for each input locus
annot_repeat_locus_summary <- ddply(annot_repeat_genomic_partition, colnames(annot_repeat_genomic_partition)[1:31], function(x) {
    if(nrow(x)==1) {
        return(cbind(numRepeatOverlaps=ifelse(x$repeatBpOverlap==0, 0, 1),
                     x[,c("repeat_family", "repeat_id", "repeatChr", "repeatStart", "repeatEnd", "repeatStrand", "repeatBpOverlap")])) }
    else { 
        return(data.frame(numRepeatOverlaps=nrow(x),
                          repeat_family=paste0(x$repeat_family, collapse=";"),
                          repeat_id=paste0(x$repeat_id, collapse=";"),
                          repeatChr=paste0(x$repeatChr, collapse=";"),
                          repeatStart=paste0(x$repeatStart, collapse=";"),
                          repeatEnd=paste0(x$repeatEnd, collapse=";"),
                          repeatStrand=paste0(x$repeatStrand, collapse=";"),
                          repeatBpOverlap=paste0(x$repeatBpOverlap, collapse=";"),
                          stringsAsFactors = FALSE))
    }})

unannot_repeat_locus_summary <- ddply(unannot_repeat_genomic_partition, colnames(unannot_repeat_genomic_partition)[1:31], function(x) { 
    if(nrow(x)==1) {
        return(cbind(numRepeatOverlaps=ifelse(x$repeatBpOverlap==0, 0, 1),
                     x[,c("repeat_family", "repeat_id", "repeatChr", "repeatStart", "repeatEnd", "repeatStrand", "repeatBpOverlap")])) }
    else {
        return(data.frame(numRepeatOverlaps=nrow(x),
                          repeat_family=paste0(x$repeat_family, collapse=";"),
                          repeat_id=paste0(x$repeat_id, collapse=";"),
                          repeatChr=paste0(x$repeatChr, collapse=";"),
                          repeatStart=paste0(x$repeatStart, collapse=";"),
                          repeatEnd=paste0(x$repeatEnd, collapse=";"),
                          repeatStrand=paste0(x$repeatStrand, collapse=";"),
                          repeatBpOverlap=paste0(x$repeatBpOverlap, collapse=";"),
                          stringsAsFactors = FALSE))
    }})

## merge the repeat information with the mRNA and lncRNA information
annot_mRNA_lncRNA_repeat_partition <- merge(annot_mRNA_lncRNA_partition, annot_repeat_locus_summary, by=1:31)
unannot_mRNA_lncRNA_repeat_partition <- merge(unannot_mRNA_lncRNA_partition, unannot_repeat_locus_summary, by=1:31)

## now write them out
write.table(annot_mRNA_lncRNA_repeat_partition,
            paste0(wdir, "/", BASENAME, "_all_partition_summary_per_annot_locus.txt"),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(unannot_mRNA_lncRNA_repeat_partition,
            paste0(wdir, "/", BASENAME, "_all_partition_summary_per_unannot_locus.txt"),
            quote=F, sep="\t", row.names=F, col.names=T)

