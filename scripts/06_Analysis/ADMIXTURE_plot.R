# Packages
library(tidyverse)

# Read in metadata
metadata <- read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.csv", col_names = FALSE) %>% 
        select(1, 6, 7) %>%
        rename(Sample = X1) %>%
        rename(Location = X6) %>%
        rename(Cluster = X7) %>% 
        as.data.frame()

### add a priori classfication of Sabah clusters based on trees

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

# read in admix and bind to metadata
admixture_table <- read.table("cleaned.3.Q") %>% 
  cbind(metadata)

# Plot
admix_plot <- admixture_table %>% 
  pivot_longer(1:3, names_to = "Ancestry", values_to = "Proportion") %>% 
  mutate(Cluster = ifelse(Ancestry == "V1", "Mf",
                          ifelse(Ancestry == "V2", "Peninsular", "Mn"))) %>% 
  mutate(Ancestry = Proportion) %>% 
  ggplot(aes(x = Sample, y = Ancestry, fill = Cluster)) +
  geom_col() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("#440154FF", "#39568CFF", "#73D055FF"))

ggsave("admix_plot.png", width = 16, dpi = 600)

# Plot mixed samples (> 0.05)
admix_plot <- admixture_table %>% 
  filter(V1 > 0.005 & V2 > 0.005 |
           V1 > 0.005 & V3 > 0.005 |
           V2 > 0.005 & V3 > 0.005) %>% 
  pivot_longer(1:3, names_to = "Ancestry", values_to = "Proportion") %>% 
  mutate(Sample = str_remove(Sample, "_DKD.*")) %>%
  mutate(Cluster = ifelse(Ancestry == "V1", "Mf",
                          ifelse(Ancestry == "V2", "Peninsular", "Mn"))) %>% 
  mutate(Ancestry = Proportion) %>% 
  ggplot(aes(x = Sample, y = Ancestry, fill = Cluster)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("#440154FF", "#39568CFF", "#73D055FF"))

ggsave("admix_plot_mixed_samples.png", dpi = 600)
