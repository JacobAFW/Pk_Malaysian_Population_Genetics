# depth plot
library(tidyverse)
library(janitor)

plot_depth <- function(SAMPLE){
    depth_plot <- read_tsv("depth_with_header.tsv", col_names = TRUE) %>% 
        mutate_if(is.logical, as.numeric) %>% 
        mutate(across(2:ncol(.), ~ ifelse(is.na(.x), 0, .x))) %>%
        select(Contig, Bases, SAMPLE) %>%
        rename(Depth = 3) %>%
        ggplot(aes(x = Bases, y = Depth, fill = Contig)) +
        geom_col() +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
        scale_fill_viridis_d()
    ggsave(paste0(SAMPLE,"_depth.pdf"), depth_plot, dpi = 300, width = 50)
}

plot_depth("M-23-4013_NKM4-15_C8-C8_L002")
plot_depth("M-23-4007_SK23_E7-E7_L002")
plot_depth("M-23-4017_NKM4-84_H8-H8_L002")
plot_depth("M-23-4010_NKM1-18_H7-H7_L002")
plot_depth("Merged_MK29_L002")