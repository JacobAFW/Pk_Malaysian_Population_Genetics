# Packages
library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(data.table)


# Read in and wrangle VCF

pkdhfr_filtered <- read_table("pkdhfr_filtered.vcf", skip = 90) %>%
        rename("CHROM" = `#CHROM`) %>%
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

pkdhfr_unfiltered <- read_table("pkdhfr_unfiltered.vcf", skip = 89) %>%
        rename("CHROM" = `#CHROM`) %>%
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

annotation <- read_table("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/PK_consensus_filtered_annotation.txt", skip = 82) %>%
  rename(Variation = `#Uploaded_variation`) %>%
  separate(Variation, sep = "_", c("ordered", "PKNH", "CHROM", "v2", "POS", "Variant")) %>%
  select(-c(ordered, PKNH, v2, Allele, Location, Variant, Feature, Feature_type, cDNA_position, CDS_position, Codons, Existing_variation)) %>%
  mutate(POS = as.numeric(POS))

pkdhfr_filtered %>%
  left_join(annotation) %>%
  filter(grepl("missense", Consequence)) %>%
  pivot_longer(3:204, names_to = "Sample", values_to = "SNP") %>% 
  mutate(SNP = as.numeric(SNP)) %>%
  filter(SNP > 0) %>% # filter for alt alleles 
  group_by(CHROM, POS, SNP) %>%
  summarise(n = n()) %>% 
  filter(n > 1) %>% write_tsv("SNP_positions.tsv") # filter for MAC > 2

pkdhfr_filtered %>%
  left_join(annotation) %>%
  filter(grepl("missense", Consequence)) %>%
  pivot_longer(3:204, names_to = "Sample", values_to = "SNP") %>% 
  mutate(SNP = as.numeric(SNP)) %>%
  filter(SNP > 0) %>%
  filter(POS == "446152") #%>% write_tsv("POS_446152.tsv")

Metadata <- readxl::read_xlsx("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/metadata/PK_Sabah_Sample_naming_indexes.xlsx") %>% 
      select(1) %>%
      add_column(Location = "Sabah") %>%
      rename(Sample = sampleid) %>%
      add_column(Cluster = "Sabah") %>% 
      rbind(
        read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/data/metadata/Pk_clusters_metadata.csv") %>%
        select(1, 3, 5) %>% 
        rename(Location = area) %>%
        mutate(Cluster = str_remove(Group, "-Pk")) %>% 
        select(-Group) 
      ) %>%
      rbind(
        read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/data/metadata/Pk_clusters_peninsular_metadata.csv") %>%
        select(2, 8) %>%
        rename(Sample = ENA_accession_no_ES) %>%
        add_column(Cluster = "Peninsular") 
      ) %>%
    mutate(Location = str_replace(Location, " ", "_")) %>%
    mutate(Cluster = ifelse(grepl("PK_SB_DNA_011", Sample) | # does the ID match any of these
                                grepl("PK_SB_DNA_091", Sample) | 
                                grepl("PK_SB_DNA_043", Sample) |
                                grepl("PK_SB_DNA_016", Sample) |  
                                grepl("PK_SB_DNA_092", Sample) | 
                                grepl("PK_SB_DNA_030", Sample) | 
                                grepl("PK_SB_DNA_093", Sample) | 
                                grepl("PK_SB_DNA_042", Sample) | 
                                grepl("PK_SB_DNA_063", Sample), "Mn", .$Cluster)) %>% # if not, just use values from Cluster
    mutate(Cluster = ifelse(Cluster == "Sabah", "Mf", .$Cluster))



pkdhfr_meta <- pkdhfr_filtered %>%
  left_join(annotation) %>%
  filter(grepl("missense", Consequence)) %>%
  pivot_longer(3:204, names_to = "Sample", values_to = "SNP") %>% 
  mutate(SNP = as.numeric(SNP)) %>%
  #filter(SNP > 0) %>% # filter for alt alleles 
  mutate(Sample = str_remove(Sample, "_DK.*")) %>%
  left_join(Metadata)

# AC filter
pkdhfr_filter <- pkdhfr_meta %>%
        filter(SNP >= 1) %>% 
        group_by(CHROM, POS) %>%
        summarise(n = n()) %>%
        filter(n > 1)

pkdhfr_meta <- pkdhfr_meta %>%
  filter(POS %in% pkdhfr_filter$POS)

pkdhfr_meta %>%
        filter(SNP != -1) %>%
        group_by(Cluster, CHROM, POS) %>%
        summarise(SNP_freq = sum(SNP >= 1)/n()*100) %>%
        na.omit()  %>%
        arrange(desc(SNP_freq)) %>%
        write_tsv("pkdhfr_cluster_stats.tsv")

pkdhfr_meta %>%
        filter(SNP != -1) %>%
        group_by(CHROM, POS) %>%
        summarise(SNP_freq = sum(SNP >= 1)/n()*100) %>%
        na.omit()  %>%
        arrange(desc(SNP_freq)) %>%
        write_tsv("pkdhfr_stats.tsv")

pkdhfr_meta %>% 
  filter(SNP >= 1) %>%
  select(CHROM, POS, Gene, Consequence, Protein_position, Amino_acids) %>% 
  unique() %>%
  write_tsv("orthologue_variants_in_Pk_pkdhfr.tsv")