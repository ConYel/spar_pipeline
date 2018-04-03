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




