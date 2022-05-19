# Load packages
library(tidyverse)

# Read in data
pk_maf <- read_table("Pk.frq", col_names=T) %>%
    mutate(CHR = str_remove(SNP, ":.*"))

# Plot MAF by SNP
pk_maf_SNP <- pk_maf  %>%
    ggplot(aes(x = SNP, y = MAF, fill = CHR)) +
    geom_col() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_hline(yintercept = 0.01, colour = "#1F968BFF") +
    theme(axis.text.x = element_text(angle = 45), legend.position = "none") +
    scale_fill_viridis_d() 

ggsave("Pk_maf_plot.png", dpi = 600, pk_maf_SNP)

# Density of MAF by SNP
pk_maf_density <- pk_maf  %>%
    ggplot(aes(x = MAF)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept = 0.01, colour = "#1F968BFF") 

ggsave("Pk_maf_density.png", dpi = 600, pk_maf_density)

# Plot MAF by chr
pk_maf_chr <- pk_maf %>% 
    mutate(SNP = str_remove(SNP, ":.*"),
        SNP = str_remove(SNP, "ordered_"),
        SNP = str_remove(SNP, "new_"),
        SNP = str_remove(SNP, "_v2"),
        SNP = str_remove(SNP, "A1.*")) %>%
    group_by(SNP) %>%
    summarise(MAF = mean(MAF)) %>%
    ggplot(aes(x = SNP, y = MAF)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45), legend.position = "none") +
    scale_fill_viridis_d() 

ggsave("Pk_maf_plot_chr.png", dpi = 600, pk_maf_chr)

# Summary stats

pk_maf %>% 
    mutate(SNP = str_remove(SNP, ":.*")) %>%
    group_by(SNP) %>%
    summarise(MAF_mean = mean(MAF),
            MAF_SD = sd(MAF),
            MAF_SE = (sd(MAF))/sqrt(n())) %>%
rbind(
    pk_maf %>%
        summarise(MAF_mean = mean(MAF),
                MAF_SD = sd(MAF),
                MAF_SE = (sd(MAF))/sqrt(n())) %>%
        add_column(SNP = "Total")
    ) %>%
    write_csv(col_names = T, "pk_maf_summary.csv")
