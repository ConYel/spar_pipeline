wdir="test_out/results"
fprefix="peaks"
partition_ref_path="./annot/partition_files/hg19"
SPAR_path="."

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
