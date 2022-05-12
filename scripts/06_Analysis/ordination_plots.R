# Load packages
library(tidyverse)
library(ape)
library(ggtree)

# PCA

## Read in data
pca <- read_table("Pk.eigenvec", col_names=F) %>%
    select(-1) %>%
    mutate(X2 = str_remove(X2, "_DK.*")) %>%
    rename("sampleid" = "X2") %>%
    rename_at(vars(starts_with("X")), ~str_replace(., "X.*", paste0("PC", seq_along(.))))

eigenval <- scan("Pk.eigenval")

## Add metadata 
pca <- pca %>%
    left_join(read_csv("Pk.csv", col_names = FALSE) %>% 
        select(1,6) %>%
        rename(sampleid = X1) %>%
        rename(Location = X6) %>%
        mutate(sampleid = str_remove(sampleid, "_DK.*")))

## Percentage variance explained by each PC
pve_plot <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100) %>%
    ggplot(aes(PC, pve)) + 
    geom_bar(stat = "identity") + 
    ylab("Percentage variance explained") + 
    theme_light()

ggsave("Percentage_variance_explained.png", dpi=600, pve_plot)

## Plot PCA
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Location)) + 
    geom_point(size = 3) +
    scale_color_viridis_d() +
    coord_equal() + 
    theme_light() + 
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("PCA.png", dpi = 600, pca_plot)


# MDS

## Read in data
mds <- read_table("Pk.mds", col_names=T) %>%
    select(-c("FID", "X14")) %>%
    mutate(IID = str_remove(IID, "_DK.*")) %>%
    rename("sampleid" = "IID") %>%
    rename_at(vars(starts_with("C")), ~str_replace(., "C", "MDS"))

## Add metadata 
mds <- mds %>%
    left_join(read_csv("Pk.csv", col_names = FALSE) %>% 
        select(1,6) %>%
        rename(sampleid = X1) %>%
        rename(Location = X6) %>%
        mutate(sampleid = str_remove(sampleid, "_DK.*")))

## Plot MDS
mds_plot <- ggplot(mds, aes(MDS1, MDS2, colour = Location)) + 
    geom_point(size = 3) +
    scale_color_viridis_d() +
    coord_equal() + 
    theme_light()

ggsave("MDS.png", dpi=600, mds_plot)


# Neighbour-joining tree

## Create distance matrix from PLINK data and built NJT
NJT_ID <- read_table("Pk.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("Pk.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

## Plot tree
# options(ignore.negative.edge=TRUE)
NJT_metadata <- read_csv("Pk.csv", col_names = FALSE) %>%
    select(1, 6 ) %>%
    rename(Sample = X1) %>%
    rename(Location = X6) %>%
    mutate(Sample = str_remove(Sample, "_DKD.*")) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(Sabah = ifelse(Location == "Sabah", "Sabah", "Other"))

options(ignore.negative.edge=TRUE)

NJT_tree_plot <- ggtree(NJT_tree, layout = "circular", size = 0.5, aes(colour = Location)) %<+% NJT_metadata +
    geom_tippoint(aes(colour = Location, shape = Sabah)) + 
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) 
    
ggsave("NJT_tree.png", dpi = 600, height = 8, width = 16, NJT_tree_plot)