library(tidyverse)

read_csv("Pk.csv", col_names = F) %>% 
    mutate(Fst_cluster = ifelse(grepl("PK_SB_DNA_021", X1) | grepl("PK_SB_DNA_088", X1) | grepl("PK_SB_DNA_049", X1) | grepl("PK_SB_DNA_094", X1) | grepl("PK_SB_DNA_046", X1) | # does the ID match any of these
                    grepl("PK_SB_DNA_054", X1) | grepl("PK_SB_DNA_040", X1) | grepl("PK_SB_DNA_005", X1) | grepl("PK_SB_DNA_090", X1) | grepl("PK_SB_DNA_014", X1) |
                    grepl("PK_SB_DNA_036", X1) | grepl("PK_SB_DNA_029", X1) | grepl("PK_SB_DNA_044", X1) | grepl("PK_SB_DNA_072", X1) | grepl("PK_SB_DNA_067", X1) |
                    grepl("PK_SB_DNA_082", X1) | grepl("PK_SB_DNA_070", X1) | grepl("PK_SB_DNA_052", X1) | grepl("PK_SB_DNA_002", X1) | grepl("PK_SB_DNA_009", X1) |
                    grepl("PK_SB_DNA_035", X1) | grepl("PK_SB_DNA_001", X1) | grepl("PK_SB_DNA_031", X1) | grepl("PK_SB_DNA_084", X1) | grepl("PK_SB_DNA_077", X1) |
                    grepl("PK_SB_DNA_074", X1) | grepl("PK_SB_DNA_015", X1) | grepl("PK_SB_DNA_023", X1) | grepl("PK_SB_DNA_013", X1) | grepl("PK_SB_DNA_047", X1) |
                    grepl("PK_SB_DNA_038", X1) | grepl("PK_SB_DNA_081", X1) | grepl("PK_SB_DNA_003", X1) | grepl("PK_SB_DNA_086", X1) | grepl("PK_SB_DNA_089", X1) |
                    grepl("PK_SB_DNA_034", X1) | grepl("PK_SB_DNA_057", X1) | grepl("PK_SB_DNA_045", X1), "Mf", .$X7)) %>% # if not, just use values from X7 - clusters and Sabah
    mutate(Fst_cluster = ifelse(Fst_cluster == "Sabah", "Mn", .$Fst_cluster)) %>% # if its the remaining Sabah samples, make them Mn, else keep them the same
    select(1, 2, 8) %>%
    write_csv("sample_clusters.csv", col_names = F)