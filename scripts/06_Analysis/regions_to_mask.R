library(tidyverse)
library(readr)
library(stringr)
library(data.table)

# Read in and wrangle VCF - same approach as hmmIBD genotype file creation - working in /g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/

genotype_matrix <- read_tsv("hmmIBD/grep_patterns.tsv", col_names = c("CHROM", "POS")) %>%
        left_join(
                read_table("PK_consensus_no_MOI.vcf", skip = 87) %>%
                rename("CHROM" = `#CHROM`) # this first section selects only variants that passed the PLINK filters we applied
        ) %>%
        mutate(CHROM = str_remove(CHROM, "ordered_PKNH_"),
        CHROM = str_remove(CHROM, "_v2")) %>%
        mutate_at(c(10:ncol(.)), ~str_remove(., ":.*")) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "1/1" | . %like% "1/0" | . %like% "0/1" , "1", .)) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "2/2" | . %like% "2/0" | . %like% "0/2" , "2", .)) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "3/3" | . %like% "3/0" | . %like% "0/3" , "3", .)) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "4/4" | . %like% "4/0" | . %like% "0/4" , "4", .)) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "0/0" , "0", .)) %>% # 0 for reference allele
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "./." , "-1", .)) %>% # for missing 
        select(-c(3:9)) %>%
        arrange(CHROM, POS) 

# Create sliding windows with SNP counts
genotype_matrix %>%
    pivot_longer(3:ncol(.), names_to = "Sample", values_to = "SNP") %>%
    mutate(SNP = as.integer(SNP)) %>%
    group_by(CHROM, POS) %>%
    summarise(SNP_count = sum(SNP > 0), SNP_Fre = sum(SNP > 0)/length(unique(.$Sample))*100) %>%
    write_tsv("SICAvar_KIRR_exclusion_table.tsv")