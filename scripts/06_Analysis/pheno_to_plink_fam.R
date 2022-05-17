# Packages
library(tidyverse)
library(janitor)
library(data.table)

# Meatdata
Location_data <- readxl::read_xlsx("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/metadata/PK_Sabah_Sample_naming_indexes.xlsx") %>% 
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
  mutate(Location = str_replace(Location, " ", "_"))


Plink_fam <- read_csv("Pk.csv", col_names=FALSE) %>%
  select(-ncol(.)) %>%
  mutate(Sample = str_remove(X2, "_DKD.*")) %>% 
  left_join(Location_data) %>% 
  select(-Sample) 

write_csv(Plink_fam, "Pk.csv", col_names = FALSE)
