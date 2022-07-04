library(tidyverse)

read_csv("Pk.csv", col_names = F) %>% 
    mutate(Fst_cluster = ifelse(grepl("PK_SB_DNA_011", X1) | # does the ID match any of these
                                grepl("PK_SB_DNA_091", X1) | 
                                grepl("PK_SB_DNA_043", X1) |
                                grepl("PK_SB_DNA_016", X1) |  
                                grepl("PK_SB_DNA_092", X1) | 
                                grepl("PK_SB_DNA_030", X1) | 
                                grepl("PK_SB_DNA_093", X1) | 
                                grepl("PK_SB_DNA_042", X1) | 
                                grepl("PK_SB_DNA_063", X1) |
                                grepl("PK_SB_DNA_028", X1) , "Mn", .$X7)) %>% # if not, just use values from X7 - clusters and Sabah
    mutate(Fst_cluster = ifelse(Fst_cluster == "Sabah", "Mf", .$Fst_cluster)) %>% # if its the remaining Sabah samples, make them Mn, else keep them the same
    select(1, 2, 8) %>%
    filter(Fst_cluster != "Peninsular") %>%
    write_csv("sample_clusters.csv", col_names = F)