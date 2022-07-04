# Packages
library(dplyr)
library(readr)
library(stringr)
library(data.table)

# Read in and wrangle VCF
read_table("PK_consensus_unzipped.bcf", skip = 93) %>%
    rename("CHROM" = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, "ordered_PKNH_"),
            CHROM = str_remove(CHROM, "_v2")) %>%
    mutate_at(c(10:ncol(.)), ~str_remove(., ":.*")) %>% 
    mutate_at(c(10:ncol(.)), ~ifelse(. %like% "1/1", 1, .)) %>%
    mutate_at(c(10:ncol(.)), ~ifelse(. %like% "1/0" | . %like% "0/1" , 1, .)) %>%  # heterozygous - 1 for alt allele
    mutate_at(c(10:ncol(.)), ~ifelse(. %like% "0/0" , 0, .)) %>% # 0 for reference allele
    mutate_at(c(10:ncol(.)), ~ifelse(. %like% "./." , -1, .)) %>% # for missing 
    mutate_at(c(10:ncol(.)), ~as.numeric(.)) %>%
    mutate_at(c(10:ncol(.)), ~ifelse(POS == lag(POS) & . != 0 & . != -1, lag(.) + 1, .)) %>%  # if an alt allele is present, but the row above is the same position, assign the genotype code to +1 of the row above - this will give us a unique genotype number for each allele
    select(-c(3:9)) %>%
    na.omit()# %>%
    write_tsv("hmmIBD_BCF.tsv")   



