library(tidyverse)
library(janitor)

# Read in depth and bases data
depth_data <- read_tsv("batch1/read_depth.tsv") %>% 
    rbind(
        read_tsv("batch2/read_depth.tsv")
    ) %>% 
    rbind(
        read_tsv("batch3/read_depth.tsv")
    ) %>% 
    rbind(
        read_tsv("batch4/read_depth.tsv")
    ) %>% 
    rbind(
        read_tsv("batch5/read_depth.tsv")
    ) %>% 
    select(-c("new_API_strain_A1_H.1", "PKNH_MIT_v2", "Data", "Alignment")) %>% 
    pivot_longer(!Sample, names_to = "Contig", values_to = "Depth") %>% 
    mutate(Threshold = ifelse(Depth > 5, "Yes", "No"))
    #rename_with( ~ paste0("Depth_", .x), starts_with("ordered"))


missing_bases <- read_tsv("batch1/bases.tsv") %>% 
    rbind(
        read_tsv("batch2/bases.tsv")
    ) %>% 
    rbind(
        read_tsv("batch3/bases.tsv")
    ) %>% 
    rbind(
        read_tsv("batch4/bases.tsv")
    ) %>% 
    rbind(
        read_tsv("batch5/bases.tsv")
    ) %>% 
    select(-c("new_API_strain_A1_H.1", "PKNH_MIT_v2", "Data", "Alignment")) %>% 
    pivot_longer(!Sample, names_to = "Contig", values_to = "Bases") %>% 
    mutate(Threshold = ifelse(Bases < 25, "Yes", "No"))
    #rename_with( ~ paste0("Bases_", .x), starts_with("ordered")) 

combined_data <- depth_data %>% 
    left_join(missing_bases) 


# Add metadata
sequencing_with_meta <- combined_data %>% 
    mutate(Sample = str_remove(Sample, "_DK.*")) %>% 
    mutate(Sample = str_replace(Sample, "_", "-")) %>% 
    left_join(
        read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/metadata/full_metadata.tsv") %>% 
            rename(Sample = studycode) %>% 
            select(Sample, sevemal, parasitemia) 
    ) %>% 
    mutate(sevemal = ifelse(grepl("UM", sevemal), "UM", .$sevemal)) %>% 
    mutate(sevemal = ifelse(sevemal == "1", "SM", ifelse(sevemal == "0", "UM", .$sevemal))) %>% 
    filter(!grepl("pos-ctrl", Sample)) %>% 
    filter(!grepl("...99", Sample))

# How many samples < 5000 pass?
sequencing_with_meta %>% 
     group_by(Sample) %>% 
     summarise(parasitemia = mean(parasitemia), Depth = mean(Depth), Bases = mean(Bases)) %>% 
     filter(parasitemia < 5000) %>% 
     mutate(Bases = ifelse(is.na(Bases), 0, .$Bases)) %>% 
     mutate(Threshold = ifelse(Depth > 5 & Bases < 25, "Yes", "No")) %>%
     count(Threshold)

sequencing_with_meta %>% 
     group_by(Sample) %>% 
     summarise(parasitemia = mean(parasitemia), Depth = mean(Depth), Bases = mean(Bases)) %>% 
     filter(parasitemia < 5000) %>% 
     mutate(Bases = ifelse(is.na(Bases), 0, .$Bases)) %>% 
     mutate(Threshold = ifelse(Depth > 30 & Bases < 10, "Yes", "No")) %>% #write_csv("low_para_samples.csv")
     count(Threshold)
     
# How many severe samples fail our thresholds?
sequencing_with_meta %>% 
    filter(sevemal == "SM") %>%
     group_by(Sample) %>% 
     summarise(parasitemia = mean(parasitemia), Depth = mean(Depth), Bases = mean(Bases)) %>% 
     mutate(Bases = ifelse(is.na(Bases), 0, .$Bases)) %>% 
     mutate(Threshold = ifelse(Depth > 5 & Bases < 25, "Yes", "No")) %>% 
     count(Threshold)

sequencing_with_meta %>% 
    filter(sevemal == "SM") %>%
     group_by(Sample) %>% 
     summarise(parasitemia = mean(parasitemia), Depth = mean(Depth), Bases = mean(Bases)) %>% 
     mutate(Bases = ifelse(is.na(Bases), 0, .$Bases)) %>% 
     mutate(Threshold = ifelse(Depth > 5 & Bases < 25, "Yes", "No")) %>% 
     filter(Threshold == "No") %>% 
     write_tsv("Failed_severe_samples.tsv")

sequencing_with_meta %>% 
    filter(sevemal == "SM") %>%
     group_by(Sample) %>% 
     summarise(parasitemia = mean(parasitemia), Depth = mean(Depth), Bases = mean(Bases)) %>% 
     mutate(Bases = ifelse(is.na(Bases), 0, .$Bases)) %>% 
     mutate(Threshold = ifelse(Depth > 30 & Bases < 10, "Yes", "No")) %>% 
     count(Threshold)