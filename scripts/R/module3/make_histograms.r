library(tidyverse)
setwd("/home")
files_hist <- dir(full.names = TRUE,
                  recursive = TRUE,
                  pattern = "histogram_table.txt")

sapply(files_hist, function(x)
      read_tsv(x, col_names = c("count", "length"),
              col_types = cols_only(
              count = col_integer(),
              length = col_integer()
              )) %>%
    ggplot(aes(x = length )) +
    geom_histogram(
      binwidth = 0.5,
      color = "red",
      fill = "red",
      alpha = .2) +
    theme_minimal() +
    labs(title="filtered softclipped reads Length Distribution",
         y = "count")+
    scale_x_continuous(limits = c(13,51), breaks = seq(13,51,by = 1)) + 
    ggsave(filename = "histogram_length_distribution.pdf",
           device = "pdf",
           dpi = "retina",
           path = str_replace(x,"/results/histogram_table.txt","/figures"))
)
