# Load packages
library(tidyverse)

# Read in data
snp_data <- read_table("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/cleaned.bim", col_names = F) %>%
    select(X2) %>%
    separate(X2, into = c("Contig", "Pos", "Major", "Minor", "Nothing"), sep = ":") %>%
    mutate(Pos = as.numeric(Pos)) %>%
    select(1:2) %>%
    arrange(Contig, Pos)

# Create window 
window_size <- 1000

snp_plot <- snp_data %>%
    mutate(Window = (floor(Pos/window_size) * window_size)+ (window_size/2)) %>%
    group_by(Contig, Window) %>%
    summarise(SNPs = n()) %>%
    ggplot(aes(x = Window, y = SNPs, colour = Contig)) +
    geom_point() +
    facet_wrap(~Contig) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
    xlab("Chromosome position")

ggsave("SNP_density.png", dpi = 600, width = 14, snp_plot)


# Histogram
snp_plot <- snp_data %>%
    mutate(Window = (floor(Pos/window_size) * window_size)+ (window_size/2)) %>%
    group_by(Contig, Window) %>%
    summarise(SNPs = n()) %>%
    arrange(SNPs) %>%
    add_column(window = 1:nrow(.)) %>%
    ggplot(aes(x = window, y = SNPs)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
    geom_col()

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/SNP_col.png", dpi = 600, width = 14, snp_plot)


snp_plot <- snp_data %>%
    mutate(Window = (floor(Pos/window_size) * window_size)+ (window_size/2)) %>%
    group_by(Contig, Window) %>%
    summarise(SNPs = n()) %>%
    arrange(SNPs) %>%
    ggplot(aes(x = SNPs)) +
    geom_histogram(binwidth = 5) +
    ylab("Window count")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/SNP_hist.png", dpi = 600, width = 14, snp_plot)
