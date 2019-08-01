library(tidyverse)
setwd("/home")
files_hist <- dir(full.names = TRUE,
                  recursive = TRUE,
                  pattern = "histogram_table.txt")

sapply(files_hist, function(x)read_delim(x, col_names = c("count", "length"),
           delim = " ") %>% 
         mutate(count = as.numeric(str_trim(count))) %>% 
    ggplot(aes(x = length, y =  count)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    labs(title="Filtered softclipped reads Length Distribution",
         y = "count")+
    scale_x_continuous(breaks = seq(13,52,by = 2)) + 
    ggsave(filename = "histogram_length_distribution.pdf",
           device = "pdf",
           dpi = "retina",
           path = str_replace(x,"/results/histogram_table.txt","/figures"))
)


files_hist[1] %>% 
    read_delim(col_names = c("count", "length"),
                             delim = " ") %>% 
    mutate(count = as.numeric(str_trim(count))) %>% 
    mutate(countm = count / 1e6)%>% 
    ggplot(aes(x = length, y =  countm)) +
    geom_bar(stat="identity") +
    theme_minimal() 
