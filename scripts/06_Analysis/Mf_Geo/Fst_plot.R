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
    #geom_hline(yintercept = 0.4, colour = "#1F968BFF")

ggsave("Mf_plots/Fst_sliding_window_plot.png", dpi = 600, width = 14, Fst_plot_window)


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

Mf_Fst_towers <- Fst_matrix %>%
    mutate(FST = as.numeric(FST)) %>%
    filter(FST != "nan" & CHR == "04" & FST > 0.5 & POS > 1000000 ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "07" & FST > 0.75)
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "08" & FST > 0.5 & FST < 0.65 ) 
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "11" & FST > 0.5 & POS < 1500000 ) 
    ) %>%
    rbind(
        Fst_matrix %>%
            mutate(FST = as.numeric(FST)) %>%
            filter(FST != "nan" & CHR == "12" & FST > 0.5 & POS > 1500000 & POS < 2500000)
    )

write_tsv(Mf_Fst_towers, "Mf_plots/Mf_Fst_towers.tsv")



