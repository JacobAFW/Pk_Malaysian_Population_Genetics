# Packages
library(dplyr)
library(readr)
library(stringr)
library(data.table)
library(tidyverse)

# Read in and wrangle VCF

SNP_matrix <- read_tsv("variant_positions_combined.tsv", col_names = c("CHROM", "POS")) %>%
        left_join(
                read_table("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/PK_consensus_no_MOI.vcf", skip = 87, n_max = 20000) %>%
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



metadata <- read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.csv", col_names = FALSE) %>% 
        select(1, 6, 7) %>%
        rename(Sample = X1) %>%
        rename(Location = X6) %>%
        rename(Cluster = X7) %>% 
        as.data.frame()

metadata <- metadata %>% 
    mutate(Cluster = ifelse(grepl("PK_SB_DNA_011", Sample) | # does the ID match any of these
                                grepl("PK_SB_DNA_091", Sample) | 
                                grepl("PK_SB_DNA_043", Sample) |
                                grepl("PK_SB_DNA_016", Sample) |  
                                grepl("PK_SB_DNA_092", Sample) | 
                                grepl("PK_SB_DNA_030", Sample) | 
                                grepl("PK_SB_DNA_093", Sample) | 
                                grepl("PK_SB_DNA_042", Sample) | 
                                grepl("PK_SB_DNA_063", Sample), "Mn", .$Cluster)) %>% # if not, just use values from X7 - clusters and Sabah
    mutate(Cluster = ifelse(Cluster == "Sabah", "Mf", .$Cluster)) # if its the remaining Sabah samples, make them Mn, else keep them the same


SNP_matrix %>% 
    pivot_longer(3:ncol(.), names_to = "Sample", values_to = "SNP") %>% 
    filter(SNP != "-1") %>% # remove missing data
    mutate(SNP = ifelse(SNP > 0, "1", .$SNP)) %>% # convert to presence-absence of alt allele
    left_join(
        metadata %>% select(Sample, Cluster)
    ) %>% 
    group_by(CHROM, POS, Cluster) %>% 
    count(SNP) %>%
    pivot_wider(names_from = Cluster, values_from = n) %>% 
    filter(SNP != "0") %>% # filter for variants only 
    filter(Mn == nrow(metadata %>% filter(Cluster == "Mn")) & is.na(Mf) & is.na(Peninsular)) %>% 
    select(CHROM, POS) %>%
    write_tsv("SNPs_segregated_Mn.txt")

