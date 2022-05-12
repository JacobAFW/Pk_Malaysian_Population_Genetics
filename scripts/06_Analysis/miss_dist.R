# Load packages
library(tidyverse)
library(data.table)


# Sample-based missing data
## Read in data
sample_mis <- read_table("Pk.imiss", col_names=T) 

## Plot dist by sample
sample_miss_plot <- sample_mis %>% 
    rename(Sample = IID,
        Geno_Miss = F_MISS) %>%
    mutate(Sample = str_remove(Sample, "DKD.*"),
            Sample = str_remove(Sample, "PK_SB_DNA_")) %>%
    select(2:ncol(.)) %>%
    ggplot(aes(x = Sample, y = Geno_Miss, fill = ifelse(Geno_Miss > 0.75, ">0.75", "<0.75"))) +
    geom_col() +
    scale_fill_manual(values = c("#440154FF", "#1F968BFF"))  +
    theme(axis.text.x = element_text(angle = 90))  +
    labs(fill = "Threshold")
    

ggsave("sample_miss_plot.png", width = 12, dpi = 600, sample_miss_plot)

# Variant-based missing data
## Read in data
var_miss <- read_table("Pk.lmiss", col_names=T) %>%
    mutate(CHR = str_remove(SNP, ":.*")) %>%
    rename(Geno_Miss = F_MISS,
            Chr = CHR)

## Plot dist by variant and coloured by chromosome
var_miss_plot <- var_miss %>% 
    ggplot(aes(x = SNP, y = Geno_Miss, fill =  Chr)) +
    geom_col() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_viridis_d(option = "D") 

ggsave("var_miss_plot.png", width = 12, dpi = 600, var_miss_plot)


## Plot dist by chromomsome averges for variants
chr_miss_plot <- var_miss %>% 
    group_by(Chr) %>%
    summarise(Geno_mean = mean(Geno_Miss)) %>%
    ggplot(aes(x = Chr, y = Geno_mean, fill =  Chr)) +
    geom_col() +
    theme(axis.text.x = element_text (angle = 90), legend.position="none") +
    scale_fill_viridis_d(option = "D") +
    ylab("Mean variant geno-miss per chromosome ")

ggsave("chr_miss_plot.png", width = 12, dpi = 600, chr_miss_plot)

##########################################################################################################################
# Combined summary stats - MAF and GM

maf_miss_summary <- read_csv("pk_maf_summary.csv") %>% 
    rename(Chr = SNP) %>%
    left_join(
        var_miss %>% 
            group_by(Chr) %>%
            summarise(Miss_mean = mean(Geno_Miss),
            Miss_SD = sd(Geno_Miss),
            Miss_SE = (sd(Geno_Miss))/sqrt(n())) %>%
        rbind(
            var_miss %>%
                summarise(Miss_mean = mean(Geno_Miss),
                Miss_SD = sd(Geno_Miss),
                Miss_SE = (sd(Geno_Miss))/sqrt(n())) %>%
            add_column(Chr = "Total")) 
    )

write_csv(maf_miss_summary, "pk_maf_miss_summary.csv")