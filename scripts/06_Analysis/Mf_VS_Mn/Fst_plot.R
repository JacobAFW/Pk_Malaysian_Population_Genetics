library(tidyverse)

Fst_matrix <- read_table("Pk.fst", col_names=T) %>%
    mutate(CHR = str_remove(SNP, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2.*"))

# Sliding window
window_size <- 10000

Fst_plot_window <- Fst_matrix %>%
    filter(FST != "nan") %>%
    mutate(FST = as.numeric(FST)) %>%
    mutate(Window = (floor(POS/window_size) * window_size)+ (window_size/2)) %>%
    group_by(CHR, Window) %>%
    summarise(Window_Fst = mean(FST), count = n()) %>% # caculate the mean Fst for each window in each chr
    filter(count > 10) %>% # remove windows with low counts
    ggplot(aes(x = Window, y = Window_Fst, colour = CHR)) +
    geom_point() +
    facet_wrap(~CHR) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
    ylab("Fst") +
    xlab("Windows (10000)")  +
    scale_color_viridis_d("Chr") 

ggsave("Pk_Mn_vs_Mf/Fst_sliding_window_plot.png", dpi = 600, width = 14, Fst_plot_window)


Fst_plot <- Fst_matrix %>%
    filter(FST != "nan") %>%
    mutate(FST = as.numeric(FST)) %>%
    filter(FST > 0) %>%
    ggplot(aes(x = POS, y = FST, colour = CHR)) +
    geom_point() +
    facet_wrap(~CHR) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
    ylab("Fst") +
    xlab("Windows (10000)")  +
    scale_color_viridis_d("Chr") 

ggsave("Pk_Mn_vs_Mf/Fst_manhattan.png", dpi = 600, width = 20, Fst_plot)

Fst_plot <- Fst_matrix %>%
    filter(FST != "nan") %>%
    mutate(FST = as.numeric(FST)) %>%
    arrange(CHR, POS) %>%
    add_column(ROW = 1:nrow(.)) %>%
    filter(FST > 0) %>%
    ggplot(aes(x = ROW, y = FST, colour = CHR)) +
    geom_point() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
    ylab("Fst") +
    xlab("Genome")  +
    scale_color_viridis_d("Chr") +
    guides(color=guide_legend(nrow=1, byrow=TRUE)) 

ggsave("Pk_Mn_vs_Mf/Fst_manhattan_whole_genome.png", dpi = 600, width = 20, Fst_plot)

Fst_plot <- Fst_matrix %>%
    filter(FST != "nan") %>%
    filter(CHR == "14") %>%
    mutate(FST = as.numeric(FST)) %>%
    ggplot(aes(x = POS, y = FST)) +
    geom_point() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
    ylab("Fst") +
    xlab("SNPs")  +
    scale_color_viridis_d("Chr") 

ggsave("Pk_Mn_vs_Mf/Fst_manhattan_14.png", dpi = 600, width = 14, Fst_plot)

