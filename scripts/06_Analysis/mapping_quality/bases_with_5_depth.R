# define read depth function to get read depth per contig
library(tidyverse)
library(janitor)
library(data.table)

depth_chr1 <- read_tsv("ordered_PKNH_01_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr1_depth = sum(Depth > 5), chr1_total = n(), chr1_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr2 <- read_tsv("ordered_PKNH_02_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr2_depth = sum(Depth > 5), chr2_total = n(), chr2_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr3 <- read_tsv("ordered_PKNH_03_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(ch3_depth = sum(Depth > 5), chr3_total = n(), chr3_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr4 <- read_tsv("ordered_PKNH_04_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr4_depth = sum(Depth > 5), chr4_total = n(), chr4_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr5 <- read_tsv("ordered_PKNH_05_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr5_depth = sum(Depth > 5), chr5_total = n(), chr5_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr6 <- read_tsv("ordered_PKNH_06_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr6_depth = sum(Depth > 5), chr6_total = n(), chr6_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr7 <- read_tsv("ordered_PKNH_07_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr7_depth = sum(Depth > 5), chr7_total = n(), chr7_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr8 <- read_tsv("ordered_PKNH_08_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr8_depth = sum(Depth > 5), chr8_total = n(), chr8_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr9 <- read_tsv("ordered_PKNH_09_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr9_depth = sum(Depth > 5), chr9_total = n(), chr9_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr10 <- read_tsv("ordered_PKNH_10_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr10_depth = sum(Depth > 5), chr10_total = n(), chr10_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr11 <- read_tsv("ordered_PKNH_11_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr11_depth = sum(Depth > 5), chr11_total = n(), chr11_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr12 <- read_tsv("ordered_PKNH_12_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr12_depth = sum(Depth > 5), chr12_total = n(), chr12_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr13 <- read_tsv("ordered_PKNH_13_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr13_depth = sum(Depth > 5), chr13_total = n(), chr13_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chr14 <- read_tsv("ordered_PKNH_14_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr14_depth = sum(Depth > 5), chr14_total = n(), chr14_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chrAPI <- read_tsv("new_API_strain_A1_H.1_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr_API_depth = sum(Depth > 5), chr_API_total = n(), chr_API_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_chrMIT <- read_tsv("PKNH_MIT_v2_depth_with_header.tsv", col_names = TRUE) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "Sample", values_to = "Depth") %>%
    group_by(Sample) %>% 
    summarise(chr_MIT_depth = sum(Depth > 5), chr_MIT_total = n(), , chr_MIT_prop = sum(Depth > 5)/n()*100) %>%
    na.omit() 

depth_summary <- depth_chr1 %>%
    left_join(depth_chr2) %>%
    left_join(depth_chr3) %>%
    left_join(depth_chr4) %>%
    left_join(depth_chr5) %>%
    left_join(depth_chr6) %>%
    left_join(depth_chr7) %>%
    left_join(depth_chr8) %>%
    left_join(depth_chr9) %>%
    left_join(depth_chr10) %>%
    left_join(depth_chr11) %>%
    left_join(depth_chr12) %>%
    left_join(depth_chr13) %>%
    left_join(depth_chr14) %>%
    left_join(depth_chrAPI) %>%
    left_join(depth_chrMIT)
 
write_tsv(depth_summary, "Bases_above_depth5.tsv")

depth_data <- list.files(pattern = "*read_depth.tsv") %>% 
  map_dfc(~fread(.)) %>%
  rename("Sample" = "Sample...1") %>%
  rename("Data" = "Data...3") %>%
  rename("Alignment" = "Alignment...4") %>%
  select(Sample, Data, Alignment, ends_with("v2"), ends_with("H.1"))  %>%
  rename_with(function(x){gsub("ordered","depth",x)}) %>%
  rename("depth_PKNH_MIT" = "PKNH_MIT_v2") %>%
  rename("depth_new_API" = "new_API_strain_A1_H.1") %>%
  na.omit()

porportion_above_5_and_depth <- depth_data %>%
    pivot_longer(cols = 4:ncol(.), names_to = "Contig", values_to = "Depth") %>%
    select(1, 4, 5) %>%
    left_join(
        depth_summary %>%
            select(Sample, ends_with("prop")) %>%
            pivot_longer(cols = 2:ncol(.), names_to = "Contig", values_to = "Prop_above_5") %>%
            mutate(
                Contig = ifelse(Contig == "chr1_prop", "depth_PKNH_01_v2", 
                    ifelse(Contig == "chr2_prop", "depth_PKNH_02_v2",
                        ifelse(Contig == "chr3_prop", "depth_PKNH_03_v2", 
                            ifelse(Contig == "chr4_prop", "depth_PKNH_04_v2",
                                ifelse(Contig == "chr5_prop", "depth_PKNH_05_v2", 
                                    ifelse(Contig == "chr6_prop", "depth_PKNH_06_v2", 
                                        ifelse(Contig == "chr7_prop", "depth_PKNH_07_v2",
                                            ifelse(Contig == "chr8_prop", "depth_PKNH_08_v2", 
                                                ifelse(Contig == "chr9_prop", "depth_PKNH_09_v2", 
                                                    ifelse(Contig == "chr10_prop", "depth_PKNH_10_v2", 
                                                        ifelse(Contig == "chr11_prop", "depth_PKNH_11_v2", 
                                                            ifelse(Contig == "chr12_prop", "depth_PKNH_12_v2",
                                                                ifelse(Contig == "chr13_prop", "depth_PKNH_13_v2", 
                                                                    ifelse(Contig == "chr14_prop", "depth_PKNH_14_v2", 
                                                                        ifelse(Contig == "chr_API_prop", "depth_new_API", 
                                                                            ifelse(Contig == "chr_MIT_prop", "depth_PKNH_MIT", .$Contig)))))))))))))))))
    )

porportion_above_5_and_depth %>%
    group_by(Sample) %>%
    summarise(Depth = mean(Depth), Prop_above_5 = mean(Prop_above_5)) %>%
    arrange(Depth) %>%
    write_tsv("Proportion_of_bases_above_5.tsv")
    
 
