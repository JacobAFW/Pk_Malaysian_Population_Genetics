# define read depth function to get read depth per contig
library(tidyverse)
library(janitor)

read_depth_data <- function(file_path, dataset, alignment){
read_tsv(file_path, col_names = TRUE) %>% 
  select(!Bases) %>% 
  mutate_if(is.logical, as.numeric) %>% 
  mutate(across(2:ncol(.), ~ ifelse(is.na(.x), 0, .x))) %>%
  group_by(Contig) %>% 
  summarise_all(mean) %>%
  t() %>% 
  as.data.frame() %>% 
  row_to_names(1) %>% 
  add_column(Data = dataset, Alignment = alignment) %>% 
  rownames_to_column("Sample") %>%
  na.omit()
}

# define base_pairs function to change names and get percentage of bases that are NA
base_pairs <- function(file_path, dataset, alignment){
read_tsv(file_path, col_names = TRUE) %>% 
  select(!Bases) %>% 
  mutate_if(is.logical, as.numeric) %>% 
  mutate(across(2:ncol(.), ~ ifelse(is.na(.x), 0, .x))) %>%
  mutate(across(2:ncol(.), ~ na_if(.x, 0))) %>%
  group_by(Contig) %>% 
  summarise_all(funs(sum(is.na(.))/length(.) * 100)) %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(1) %>% 
  add_column(Data = dataset, Alignment = alignment) %>% 
  rownames_to_column("Sample") 
}

depth_data <- read_depth_data("depth_with_header.tsv", "Pk", "Direct")
base_data <- base_pairs("depth_with_header.tsv", "Pk", "Direct")

write_tsv(depth_data, "read_depth.tsv")
write_tsv(base_data, "bases.tsv")
