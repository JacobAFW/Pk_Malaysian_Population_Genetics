library(tidyverse)

Fst_matrix <- read_table("Pk.fst", col_names=T) %>%
    mutate(CHR = str_remove(SNP, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2.*"))

# Sliding window
window_size <- 500

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
    #geom_hline(yintercept = 0.4, colour = "#1F968BFF")

ggsave("Mn_plots/Fst_sliding_window_plot.png", dpi = 600, width = 14, Fst_plot_window)

Fst_plot_window <- Fst_matrix %>%
    filter(FST != "nan") %>%
    mutate(FST = as.numeric(FST)) %>%
    mutate(Window = (floor(POS/window_size) * window_size)+ (window_size/2)) %>%
    group_by(CHR, Window) %>%
    summarise(Window_Fst = mean(FST), count = n()) %>% # caculate the mean Fst for each window in each chr
    filter(count > 2) %>% # remove windows with low counts
    add_column(ROW = 1:nrow(.)) %>%
    ggplot(aes(x = ROW, y = Window_Fst, colour = CHR)) +
    geom_point() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
    ylab("Fst") +
    xlab("Genome")  +
    scale_color_viridis_d("Chr") +
    ylim(-.01, 1.2) +
    guides(color=guide_legend(nrow=1, byrow=TRUE)) 

ggsave("Mn_plots/Fst_sliding_window_plot_full_genome.png", dpi = 600, width = 14, Fst_plot_window)

Fst_plot <- Fst_matrix %>%
    filter(FST > 0) %>%
    filter(FST != "nan") %>%
    mutate(FST = as.numeric(FST)) %>%
    ggplot(aes(x = POS/1000000, y = FST, colour = CHR)) +
    geom_point() +
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylab("Fst") +
    xlab("Chromsome position (Mb)")  +
    scale_color_viridis_d("Chr") +
    facet_wrap(~CHR)

ggsave("Mn_plots/Fst_manhattan.png", dpi = 600, width = 20, Fst_plot)

Fst_plot <- Fst_matrix %>%
    filter(FST > 0) %>%
    filter(FST != "nan") %>%
    mutate(FST = as.numeric(FST)) %>%
    arrange(CHR, POS) %>%
    add_column(ROW = 1:nrow(.)) %>%
    ggplot(aes(x = ROW, y = FST, colour = CHR)) +
    geom_point() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
    ylab("Fst") +
    xlab("Genome")  +
    scale_color_viridis_d("Chr") +
    guides(color=guide_legend(nrow=1, byrow=TRUE)) 

ggsave("Mn_plots/Fst_manhattan_whole_genome.png", dpi = 600, width = 20, Fst_plot)


Fst_plot <- Fst_matrix %>%
    filter(FST > 0) %>%
    filter(FST != "nan") %>%
    filter(CHR == "01") %>%
    mutate(FST = as.numeric(FST)) %>%
    ggplot(aes(x = POS/1000000, y = FST, colour = CHR)) +
    geom_point() +
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylab("Fst") +
    xlab("Chromsome position (Mb)")  +
    scale_color_viridis_d("Chr") 

ggsave("Mn_plots/Fst_manhattan_01.png", dpi = 600, width = 14, Fst_plot)

# Identify SNPs in outlier regions

Mn_Fst_towers <- Fst_matrix %>%
    mutate(FST = as.numeric(FST)) %>%
    filter(FST != "nan" & CHR == "06" & FST > 0.75 & POS < 700000 & POS >  500000) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "08" & FST > 0.875)
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "10" & FST > 0.875 & POS > 1000000) 
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "11" & FST > 0.875 & POS < 1500000) 
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "13" & FST > 0.875 & POS < 1500000)
    )

write_tsv(Mn_Fst_towers, "Mn_plots/Mn_Fst_towers.tsv")
