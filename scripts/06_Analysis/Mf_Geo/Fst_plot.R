library(tidyverse)

Fst_matrix <- read_table("Pk.fst", col_names=T) %>%
    mutate(CHR = str_remove(SNP, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2.*"))

# Sliding window
window_size <- 1000

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
    xlab("Windows (1000)")  +
    scale_color_viridis_d("Chr") 
    #geom_hline(yintercept = 0.4, colour = "#1F968BFF")

ggsave("Mf_plots/Fst_sliding_window_plot.png", dpi = 600, width = 14, Fst_plot_window)


Fst_plot_window <- Fst_matrix %>%
    filter(FST != "nan") %>%
    mutate(FST = as.numeric(FST)) %>%
    mutate(Window = (floor(POS/window_size) * window_size)+ (window_size/2)) %>%
    group_by(CHR, Window) %>%
    summarise(Window_Fst = mean(FST), count = n()) %>% # caculate the mean Fst for each window in each chr
    filter(count > 10) %>% # remove windows with low counts
    add_column(ROW = 1:nrow(.)) %>%
    ggplot(aes(x = ROW, y = Window_Fst, colour = CHR)) +
    geom_point() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
    ylab("Fst") +
    xlab("Genome")  +
    scale_color_viridis_d("Chr") +
    guides(color=guide_legend(nrow=1, byrow=TRUE)) 

ggsave("Mf_plots/Fst_sliding_window_plot_full_genome.png", dpi = 600, width = 14, Fst_plot_window)



Fst_plot <- Fst_matrix %>%
    filter(FST != "nan") %>%
    mutate(FST = as.numeric(FST)) %>%
    ggplot(aes(x = POS/1000000, y = FST, colour = CHR)) +
    geom_point() +
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylab("Fst") +
    xlab("Chromsome position (Mb)")  +
    scale_color_viridis_d("Chr") +
    facet_wrap(~CHR)

ggsave("Mf_plots/Fst_manhattan.png", dpi = 600, width = 20, Fst_plot)

Fst_plot <- Fst_matrix %>%
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

ggsave("Mf_plots/Fst_manhattan_whole_genome.png", dpi = 600, width = 20, Fst_plot)

Fst_plot <- Fst_matrix %>%
    filter(FST != "nan") %>%
    filter(CHR == "08" | CHR == "11" | CHR == "12") %>%
    mutate(CHR = paste0("Chromosome ", .$CHR)) %>%
    mutate(FST = as.numeric(FST)) %>%
    ggplot(aes(x = POS/1000000, y = FST, colour = CHR)) +
    geom_point() +
        theme(legend.position = "none", 
        legend.title = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    ylab("Fst") +
    xlab("Chromsome position (Mb)")  +
    scale_color_manual(values = c("#404788FF", "#1F968BFF", "#55C667FF")) +
    facet_grid(~CHR, scales = "free", space = "free")

ggsave("Mf_plots/Fst_manhattan_facet.png", dpi = 600, Fst_plot)

# Identify SNPs in outlier regions
Fst_plot <- Fst_matrix %>%
    filter(FST != "nan") %>%
    filter(CHR == "12") %>%
    mutate(FST = as.numeric(FST)) %>%
    ggplot(aes(x = POS/1000000, y = FST, colour = CHR)) +
    geom_point() +
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylab("Fst") +
    xlab("Chromsome position (Mb)")  +
    scale_color_viridis_d("Chr") +
    scale_x_continuous(breaks = c(0.2, .4, .6, .8, 2, 2.05, 2.1, 2.15, 2.2, 2.3, 2.4, 2.6, 2.8, 3))

ggsave("Mf_plots/Fst_manhattan_12.png", dpi = 600, width = 20, Fst_plot)

Fst_matrix %>%
    mutate(FST = as.numeric(FST)) %>%
    filter(FST != "nan" & CHR == "12" & FST > 0.25 & POS > 2075000 & POS < 2125000) %>% 
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "08" & FST > 0.25 & POS > 940000 & POS < 950000) 
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "08" & FST > 0.25 & POS > 500000 & POS < 510000) 
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "11" & FST > 0.5 & POS < 1800000) 
    )   %>%
    write_tsv("Mf_plots/Mf_Fst_towers_2.tsv") 


Mf_Fst_towers <- Fst_matrix %>%
    mutate(FST = as.numeric(FST)) %>%
    filter(FST != "nan" & CHR == "08" & FST > 0.5 ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "11" & FST > 0.75)
    ) %>%
        rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "11" & FST > 0.5 & POS < 1500000 ) 
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "12" & FST > 0.5 & FST < 0.65 ) 
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "12" & FST > 0.5 & POS > 1500000 & POS < 2500000)
    )

write_tsv(Mf_Fst_towers, "Mf_plots/Mf_Fst_towers.tsv")



