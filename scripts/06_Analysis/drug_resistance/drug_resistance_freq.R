# Packages
library(tidyverse)
library(janitor)
library(dplyr)
library(readr)
library(stringr)
library(data.table)

# Read in genotype file previously created and subset to orthologue genes
genotype_file <- read_tsv("hmmIBD.tsv")
orthologue_genes <- genotype_file %>% 
  filter(CHROM == "01" & POS > 381348 & POS < 384824) %>%
  rbind(genotype_file %>% filter(CHROM == "05" & POS > 445338 & POS < 447218)) %>%
  rbind(genotype_file %>% filter(CHROM == "08" & POS > 1408170 & POS < 1411360)) %>%
  rbind(genotype_file %>% filter(CHROM == "10" & POS > 448974 & POS < 453368)) %>%
  rbind(genotype_file %>% filter(CHROM == "12" & POS > 2627516 & POS < 2629654)) %>%
  rbind(genotype_file %>% filter(CHROM == "13" & POS > 2302686 & POS < 2304038)) %>%
  rbind(genotype_file %>% filter(CHROM == "14" & POS > 1308259 & POS < 1310954)) %>%
  rbind(genotype_file %>% filter(CHROM == "14" & POS > 2116423 & POS < 2122473))


# Read in annotation text file, left join to subset genotype file & subset to non-synomynous
annotation <- read_table("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/PK_consensus_filtered_annotation.txt", skip = 82) %>%
  rename(Variation = `#Uploaded_variation`) %>%
  separate(Variation, sep = "_", c("ordered", "PKNH", "CHROM", "v2", "POS", "Variant")) %>%
  select(-c(ordered, PKNH, v2, Allele, Location, Variant, Feature, Feature_type, cDNA_position, CDS_position, Codons, Existing_variation)) %>%
  mutate(POS = as.numeric(POS))

# Left join annotation & genotype file & filter out synomynous variants 
annotated_orthologue_genes <- orthologue_genes %>%
  left_join(annotation) %>%
  filter(!grepl("synonymous_variant", Consequence)) %>%
  filter(!grepl("downstream_gene_variant", Consequence)) %>%
  filter(!grepl("upstream_gene_variant", Consequence)) %>%
  filter(!grepl("intron_variant", Consequence)) %>%
  filter(!grepl("stop_retained_variant", Consequence))


# plot

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

SNP_Frequency <- annotated_orthologue_genes %>%
  pivot_longer(3:187, names_to = "Sample", values_to = "Allele")  %>%
  mutate(Sample = str_remove(Sample, "_DK.*")) %>%
  #filter(grepl("PK", Sample)) %>%
  left_join(Metadata) %>%
  filter(Allele != -1) %>% # filter out samples with missing calles
  group_by(Cluster, CHROM, POS) %>%
  summarise(SNP_freq = sum(Allele >= 1)/n()*100) %>%
  na.omit()

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "01") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Cluster)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_01.png", dpi = 600, width = 12, SNP_Frequency_plot)


SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "05") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Cluster)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_05.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "08") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Cluster)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_08.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "10") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Cluster)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_10.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "13") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Cluster)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_13.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "14" & POS < 2000000) %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Cluster)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_14_1.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "14" & POS > 2000000) %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Cluster)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_14_2.png", dpi = 600, width = 12, SNP_Frequency_plot)

# Sub-populations within Mf
SNP_Frequency <- annotated_orthologue_genes %>%
  pivot_longer(3:187, names_to = "Sample", values_to = "Allele")  %>%
  mutate(Sample = str_remove(Sample, "_DK.*")) %>%
  #filter(grepl("PK", Sample)) %>%
  left_join(Metadata) %>%
  filter(Allele != -1) %>% # filter out samples with missing calles
  filter(Cluster == "Mf") %>%
  mutate(Location = ifelse(Location != "Sabah", "Kapit-Betong-Sarikei", "Sabah")) %>% # mutate location to be based on two sub-populations identified within Mf
  group_by(Location, CHROM, POS) %>%
  summarise(SNP_freq = sum(Allele >= 1)/n()*100)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "01") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mf_location_01.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "05") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mf_location_05.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "08") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mf_location_08.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "10") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mf_location_10.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "14" & POS < 2000000) %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mf_location_14.1.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "14" & POS > 2000000) %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mf_location_14.2.png", dpi = 600, width = 12, SNP_Frequency_plot)

# Sub-populations within Mn
SNP_Frequency <- annotated_orthologue_genes %>%
  pivot_longer(3:187, names_to = "Sample", values_to = "Allele")  %>%
  mutate(Sample = str_remove(Sample, "_DK.*")) %>%
  #filter(grepl("PK", Sample)) %>%
  left_join(Metadata) %>%
  filter(Allele != -1) %>% # filter out samples with missing calles
  filter(Cluster == "Mn") %>%
  mutate(Location = ifelse(Location != "Sabah", "Kapit-Betong-Sarikei", "Sabah")) %>% # mutate location to be based on two sub-populations identified within Mf
  group_by(Location, CHROM, POS) %>%
  summarise(SNP_freq = sum(Allele >= 1)/n()*100)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "01") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mn_location_01.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "05") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mn_location_05.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "08") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mn_location_08.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "10") %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mn_location_10.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "14" & POS < 2000000) %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mn_location_14.1.png", dpi = 600, width = 12, SNP_Frequency_plot)

SNP_Frequency_plot <- SNP_Frequency %>%
  filter(CHROM == "14" & POS > 2000000) %>%
  ggplot(aes(x = POS, y = SNP_freq, colour = Location)) +
  geom_point(size = 3) +
  geom_line(aes(group = POS), colour = "grey") +
  scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
  ylab("SNP Frequency") +
  xlab("Chromosome Postion")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/drug_resistance/SNP_Frequency_plot_Mn_location_14.2.png", dpi = 600, width = 12, SNP_Frequency_plot)

