# Load packages
library(tidyverse)

# Read in data
pk_maf <- read_table("Pk.frq.counts", col_names = T) %>%
    mutate(CHR = str_remove(SNP, ":.*")) %>%
    rename(minor_count = C1, major_count = C2)

# Plot MAF by chr
pk_maf_plot <- pk_maf %>% 
    mutate(CHR = str_remove(CHR, ":.*"),
        CHR = str_remove(CHR, "ordered_"),
        CHR = str_remove(CHR, "new_"),
        CHR = str_remove(CHR, "_v2"),
        CHR = str_remove(CHR, "A1.*")) %>%
    group_by(CHR) %>%
    summarise(minor_count = mean(minor_count)) %>%
    ggplot(aes(x = CHR, y = minor_count, fill = CHR)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45), legend.position = "none") +
    scale_fill_viridis_d() 

ggsave("Pk_maf_counts_plot_chr.png", dpi = 600, pk_maf_plot, width = 10)

# Summary stats

pk_maf %>% 
    mutate(SNP = str_remove(SNP, ":.*")) %>%
    group_by(SNP) %>%
    summarise(MAC_mean = mean(minor_count),
            MAC_SD = sd(minor_count),
            MAC_SE = (sd(minor_count))/sqrt(n()),
            MAC_min = min(minor_count),
            MAC_max = max(minor_count)) %>%
rbind(
    pk_maf %>%
        summarise(MAC_mean = mean(minor_count),
                MAC_SD = sd(minor_count),
                MAC_SE = (sd(minor_count))/sqrt(n()),
                MAC_min = min(minor_count),
            MAC_max = max(minor_count)) %>%
        add_column(SNP = "Total")
    ) %>%
    write_csv(col_names = T, "pk_maf_counts_summary.csv")
