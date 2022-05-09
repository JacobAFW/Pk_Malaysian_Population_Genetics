# Wrangle variant calling data - version 2 - using merged VCF file

# Load Packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(tibble)

Names <- c("X1") # X1 is the first columns
for (i in 1:100){
  Names <- append(Names, print(paste0("Sample_", i)))
}

read_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/04_Variant_calling/consensus/TEST/CHROMOSOME.gatk.tsv", col_names =  Names) %>% 
  select(X1) %>% 
  inner_join(read_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/04_Variant_calling/consensus/TEST/CHROMOSOME.bcftools.tsv", col_names =  Names) %>% 
               select(X1)) %>% 
  separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt")) %>% 
  write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/04_Variant_calling/consensus/TEST/CHROMOSOME_variant_names.tsv")