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

## assign to partitions
annot_partition_data <- read.table(paste0(wdir, "/", BASENAME, "_annot_overlaps.txt"), header=F, sep="\t", quote="", as.is=T)
unannot_partition_data <- read.table(paste0(wdir, "/", BASENAME, "_unannot_overlaps.txt"), header=F, sep="\t", quote="", as.is=T)

## we want to eventually merge these tables, using unannotated as a 'class'   
## add a class column (for annotated loci)
annot_partition_data$class <- unlist(lapply(strsplit(annot_partition_data[,4], ","), '[', 1))

## do partition assignments (more nicely) across all the loci
## because we used a 51% filter, there is only one possible annotation for each locus
annot_partition_data$partition <- unlist(lapply(strsplit(annot_partition_data$V10, ";"), "[", 1))
annot_partition_data$overlap_id <- unlist(lapply(strsplit(annot_partition_data$V10, ";"), "[", 2))
intergenic_loci <- annot_partition_data$partition=="."
annot_partition_data$partition[intergenic_loci] <- "Intergenic"
annot_partition_data$overlap_id[intergenic_loci] <- "Intergenic"

## rename the partition
annot_partition_data$partition <- ifelse(annot_partition_data$partition=="5utr_exon", "5' UTR Exon",
                                 ifelse(annot_partition_data$partition=="5utr_intron", "5' UTR Intron",
                                 ifelse(annot_partition_data$partition=="3utr_exon", "3' UTR Exon",
                                 ifelse(annot_partition_data$partition=="3utr_intron", "3' UTR Intron",
                                 ifelse(annot_partition_data$partition=="promoter", "Promoter",
                                 ifelse(annot_partition_data$partition=="exon", "mRNA Exon",
                                 ifelse(annot_partition_data$partition=="intron", "mRNA Intron",
                                 ifelse(annot_partition_data$partition=="repeat", "Repeat", "Intergenic"))))))))

## now write this out as extra columns on the original file, and also include the
## overlapping partition and the number of base pairs of overlap
write.table(cbind(annot_locus_data, genomicPartition=annot_partition_data$partition, overlapIds=annot_partition_data$overlap_id, partitionChr=annot_partition_data[,7], partitionStart=annot_partition_data[,8], partitionEnd=annot_partition_data[,9], partitionStrand=annot_partition_data[,12], partitionBpOverlap=annot_partition_data[,13]), paste0(wdir, "/", BASENAME, "_annot_genomic_partition.txt"), quote=F, sep="\t", row.names=F, col.names=T)

## also do this for the unannotated loci
unannot_partition_data$class <- "unannotated"
unannot_partition_data$partition <- unlist(lapply(strsplit(unannot_partition_data$V10, ";"), "[", 1))
unannot_partition_data$overlap_id <- unlist(lapply(strsplit(unannot_partition_data$V10, ";"), "[", 2))
intergenic_loci <- unannot_partition_data$partition=="."
unannot_partition_data$partition[intergenic_loci] <- "Intergenic"
unannot_partition_data$overlap_id[intergenic_loci] <- "Intergenic"

## rename the partition
unannot_partition_data$partition <- ifelse(unannot_partition_data$partition=="5utr_exon", "5' UTR Exon",
                                 ifelse(unannot_partition_data$partition=="5utr_intron", "5' UTR Intron",
                                 ifelse(unannot_partition_data$partition=="3utr_exon", "3' UTR Exon",
                                 ifelse(unannot_partition_data$partition=="3utr_intron", "3' UTR Intron",
                                 ifelse(unannot_partition_data$partition=="promoter", "Promoter",
                                 ifelse(unannot_partition_data$partition=="exon", "mRNA Exon",
                                 ifelse(unannot_partition_data$partition=="intron", "mRNA Intron",
                                 ifelse(unannot_partition_data$partition=="repeat", "Repeat", "Intergenic"))))))))                                     

## now write this out as extra columns on the original file, and also include the
## overlapping partition and the number of base pairs of overlap
write.table(cbind(unannot_locus_data, genomicPartition=unannot_partition_data$partition, overlapIds=unannot_partition_data$overlap_id, partitionChr=unannot_partition_data[,7], partitionStart=unannot_partition_data[,8], partitionEnd=unannot_partition_data[,9], partitionStrand=unannot_partition_data[,12], partitionBpOverlap=unannot_partition_data[,13]), paste0(wdir, "/", BASENAME, "_unannot_genomic_partition.txt"), quote=F, sep="\t", row.names=F, col.names=T)

## merge the tables for visualization
full_partition_data <- rbind(annot_partition_data, unannot_partition_data) 

class_col <- which(colnames(full_partition_data)=="partition")
id_col <- which(colnames(full_partition_data)=="overlap_id")

## now get the proportions and counts
## define the desired order of (known) partition classes
partition_order <- c("5' UTR Exon", "5' UTR Intron", "3' UTR Exon", "3' UTR Intron",
                     "Promoter", "mRNA Exon", "mRNA Intron", "Repeat", "Intergenic")

full_prop_tab <- as.data.frame(table(full_partition_data$partition), stringsAsFactors = F)
full_prop_tab$prop <- full_prop_tab$Freq / sum(full_prop_tab$Freq)

## re-order    
full_prop_tab <- full_prop_tab[match(c(partition_order, as.character(full_prop_tab$Var1[!(full_prop_tab$Var1 %in% partition_order)])), full_prop_tab$Var1),]

colnames(full_prop_tab) <- c("Genomic_Element", "Count", "Proportion")

write.table(full_prop_tab, paste0(wdir, "/", BASENAME, "_allclass_partition_props.txt"), quote=F, sep="\t", row.names=F, col.names=T)

## then look at assignment per sncRNA class
class_prop_tab <- ddply(full_partition_data, .(class), function(x) {
    count_df <- as.data.frame(table(x$partition), stringsAsFactors = F)
    count_df$prop <- count_df$Freq / sum(count_df$Freq)
    ## add entries for those classes that we didn't observe
    if (sum(!(partition_order %in% count_df$Var1)) > 0) {
        count_df <- rbind(count_df,
                          data.frame(Var1=partition_order[!(partition_order %in% count_df$Var1)], Freq=0, prop=0))
    }

    ## re-order    
    count_df <- count_df[match(c(partition_order, as.character(count_df$Var1[!(count_df$Var1 %in% partition_order)])), count_df$Var1),]

    colnames(count_df) <- c("Genomic_Element", "Count", "Proportion")
    return(count_df)
})

## save the classes without total
unique_classes <- sort(unique(class_prop_tab$class))

## plot the proportion of each partition by class as stacked barplots
## add a column for the total proportions
tot_counts <- ddply(class_prop_tab, .(Genomic_Element), summarize, Count=sum(Count))
## get the correct ordering of the classes
tot_ordering <- match(partition_order, tot_counts$Genomic_Element)

full_class_prop_tab <- rbind(class_prop_tab,
                        data.frame(class="Total", Genomic_Element=tot_counts$Genomic_Element[tot_ordering],
                                   Count=tot_counts$Count[tot_ordering],
                                   Proportion=tot_counts$Count[tot_ordering]/sum(tot_counts$Count)))
## order it so 'total' is last
full_class_prop_tab$class <- factor(full_class_prop_tab$class,
                               levels=c(unique_classes, "Total"))

#save(full_class_prop_tab, file=paste0(wdir, "/", BASENAME, "_partition_props_byclass.Rdata"))
save.image(file=paste0(wdir, "/", BASENAME, "_all_data.Rdata")

write.table(full_class_prop_tab, paste0(wdir, "/", BASENAME, "_partition_props_byclass.txt"), quote=F, sep="\t", row.names=F, col.names=T)

