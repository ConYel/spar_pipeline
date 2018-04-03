args <- commandArgs(trailingOnly=TRUE) 
fprefix <- args[1]
wdir <- args[2]

suppressPackageStartupMessages(library(plyr))

BASENAME <- "Genomewide_distribution_patterns_of_small_RNA_loci"

## read in the original data
## in this script, we process both the annotated and unannotated loci at the same time
annot_datafile <- paste(wdir, "/", fprefix, "_annot.with_conservation.xls", sep="")
unannot_datafile <- paste(wdir, "/", fprefix, "_unannot.with_conservation.xls", sep="")

## annotated data
annot_locus_data <- read.table(annot_datafile, header=T, sep="\t", quote="", as.is=T,
                               comment.char="", strip.white=T)
## unannotated data
unannot_locus_data <- read.table(unannot_datafile, header=T, sep="\t", quote="", as.is=T,
                               comment.char="", strip.white=T)

## read in the full repeat overlap results
annot_partition_data <- read.table(paste0(wdir, "/", BASENAME, "_repeat_annot_overlaps.txt"), header=F, sep="\t", quote="", as.is=T)
unannot_partition_data <- read.table(paste0(wdir, "/", BASENAME, "_repeat_unannot_overlaps.txt"), header=F, sep="\t", quote="", as.is=T)

## we want to eventually merge these tables, using unannotated as a 'class'   
## add a class column (for annotated loci)
annot_partition_data$class <- unlist(lapply(strsplit(annot_partition_data[,4], ","), '[', 1))

## add more specific annotations for specific repeat type (in addition to V11 which has the broader families)
annot_partition_data$repeat_family <- unlist(lapply(strsplit(gsub("repFamily=", "", annot_partition_data$V10), ";"), "[", 2))
annot_partition_data$repeat_id <- unlist(lapply(strsplit(gsub("ID=", "", annot_partition_data$V10), ";"), "[", 1))
intergenic_loci <- annot_partition_data$repeat_id=="."
annot_partition_data$repeat_id[intergenic_loci] <- "Non-repeat"
annot_partition_data$repeat_family[intergenic_loci] <- "Non-repeat"

## also add a peak ID column so that we can merge on that
annot_partition_data$peakID <- unlist(lapply(strsplit(annot_partition_data[,4], ","), '[', 2))

## now write this information out as extra columns on the original file, and also include the
## overlapping repeat and the number of base pairs of overlap
## need to duplicate rows for multiple repeat hits
merged_annot_locus_data <- merge(annot_locus_data, annot_partition_data[,c("V1", "V2", "V3", "class", "repeat_family", "repeat_id", "V7", "V8", "V9", "V12", "V13", "peakID")], by.x=1:4, by.y=c("V1", "V2", "V3", "peakID"))
merged_annot_locus_data <- rename(merged_annot_locus_data, c(V7 = "repeatChr",  V8 = "repeatStart", V9 = "repeatEnd", V12 = "repeatStrand", V13 = "repeatBpOverlap"))

write.table(merged_annot_locus_data, paste0(wdir, "/", BASENAME, "_annot_repeat_genomic_partition.txt"), quote=F, sep="\t", row.names=F, col.names=T)

## also do this for the unannotated loci
unannot_partition_data$class <- "unannotated"

## add more specific annotations for specific repeat type (in addition to V11 which has the broader families)
unannot_partition_data$repeat_family <- unlist(lapply(strsplit(gsub("repFamily=", "", unannot_partition_data$V10), ";"), "[", 2))
unannot_partition_data$repeat_id <- unlist(lapply(strsplit(gsub("ID=", "", unannot_partition_data$V10), ";"), "[", 1))
intergenic_loci <- unannot_partition_data$repeat_id=="."
unannot_partition_data$repeat_id[intergenic_loci] <- "Non-repeat"
unannot_partition_data$repeat_family[intergenic_loci] <- "Non-repeat"

## also add a peak ID column so that we can merge on that
unannot_partition_data$peakID <- unlist(lapply(strsplit(unannot_partition_data[,4], ","), '[', 2))

## now write this information out as extra columns on the original file, and also include the
## overlapping repeat and the number of base pairs of overlap
merged_unannot_locus_data <- merge(unannot_locus_data, unannot_partition_data[,c("V1", "V2", "V3", "class", "repeat_family", "repeat_id", "V7", "V8", "V9", "V12", "V13", "peakID")], by.x=1:4, by.y=c("V1", "V2", "V3", "peakID"))
merged_unannot_locus_data <- rename(merged_unannot_locus_data, c(V7 = "repeatChr",  V8 = "repeatStart", V9 = "repeatEnd", V12 = "repeatStrand", V13 = "repeatBpOverlap"))

write.table(merged_unannot_locus_data, paste0(wdir, "/", BASENAME, "_unannot_repeat_genomic_partition.txt"), quote=F, sep="\t", row.names=F, col.names=T)

## merge the tables to save as Rdata for the next stage of analysis
full_repeat_partition_data <- rbind(annot_partition_data, unannot_partition_data) 
save(full_repeat_partition_data, file=paste0(wdir, "/", BASENAME, "_repeat_partition_results.Rdata"))

