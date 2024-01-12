# Load packages
library(tidyverse)

# Read in metadata
metadata <- read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.csv", col_names = FALSE) %>% 
        dplyr::select(1, 6, 7) %>%
        rename(Sample = X1) %>%
        rename(Location = X6) %>%
        rename(Cluster = X7) %>% 
        as.data.frame() 

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


############################################################################## WORK IN PROGRESS - using VCF to identify dominant alleles
# Get allele call for each sample
# Add metadata/cluster information
# Create column that has dominant allele call for each cluster at each position in the chromosome 
# use ifelse statement - for each sample, does it have the Mn or Mf allele

genotype_table <- read_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD/hmmIBD.tsv")


# create long version of genotype data
genotype_table_long <- genotype_table %>%
    pivot_longer(3:ncol(.), names_to = "SAMPLE", values_to = "SNP") %>%
    left_join(
        metadata %>%
            dplyr::select(1,3) %>%
            rename(SAMPLE = Sample)
    ) #%>%
    #filter(CHROM == "08") # removed when finish


dominant_allele <- genotype_table_long %>% 
    filter(SNP >= 0) %>% # remove missing data
    mutate(SNP = as.factor(SNP)) %>% # convert to factor so we can create a summary of the observations
    group_by(CHROM, POS, Cluster, SNP) %>% 
    summarise(Allele_count = n()) %>% # get the counts for each unique combination you've grouped by
    #filter(POS==56510) %>%
    filter(Allele_count == max(Allele_count)) %>% # filter to get the allele with max count for each combination (including the SNP/allele)
    select(-5) %>% # remove the allele count column
    ungroup() %>% 
    group_by(CHROM, POS) %>% 
    mutate(row = cur_group_id()) %>% # create an id specific to each combo of CHROM & POS
    na.omit() %>% # remove rows with missing data
    group_by(row) %>% # group by the new id
    summarise(CHROM, POS, Cluster, SNP, n = n()) %>% # get the number of rows for this new id - any with 4 or more have duplicate clusters for the same position, ie they have a 50-50 split in allele frequency 
    filter(n < 4) %>% # filter out positions that have a cluster that has 50-50 allele split
    pivot_wider(names_from = Cluster, values_from = SNP) %>% 
    ungroup() %>%
    dplyr::select(-c("row", "n"))

# ASSUMPTION = if it is not equal to the dominant allele in any of the clusters, just make it the allele of the cluster it belongs to - this only occurs in ~0.4% of SNPs
introgression_table <- genotype_table_long %>%
    left_join(dominant_allele) %>%
    add_column(Intro_clust = .$SNP) %>%
    mutate(Intro_clust = ifelse(Cluster == "Mn" & Intro_clust == Mn, "Mn", # Mn
        ifelse(Cluster == "Mn" & Intro_clust == Mf & Intro_clust != Peninsular, "Mf", 
        ifelse(Cluster == "Mn" & Intro_clust == Peninsular & Intro_clust != Mf, "Peninsular", 
        ifelse(Cluster == "Mn" & Intro_clust == Mf & Intro_clust == Peninsular, "Mf_Peninsular",
        ifelse(Cluster == "Mn" & (Intro_clust != Mn | Intro_clust != Mf | Intro_clust != Peninsular), "Mn", .$Intro_clust)))))) %>% 
    mutate(Intro_clust = ifelse(Cluster == "Mf" & Intro_clust == Mf, "Mf", # Mf
        ifelse(Cluster == "Mf" & Intro_clust == Mn & Intro_clust != Peninsular, "Mn", 
        ifelse(Cluster == "Mf" & Intro_clust == Peninsular & Intro_clust != Mn, "Peninsular", 
        ifelse(Cluster == "Mf" & Intro_clust == Mn & Intro_clust == Peninsular, "Mn_Peninsular",
        ifelse(Cluster == "Mf" & (Intro_clust != Mn | Intro_clust != Mf | Intro_clust != Peninsular), "Mf", .$Intro_clust)))))) %>%
    mutate(Intro_clust = ifelse(Cluster == "Peninsular" & Intro_clust == Peninsular, "Peninsular", # Peninsular
        ifelse(Cluster == "Peninsular" & Intro_clust == Mf & Intro_clust != Mn, "Mf",
        ifelse(Cluster == "Peninsular" & Intro_clust == Mn & Intro_clust != Mf, "Mn", 
        ifelse(Cluster == "Peninsular" & Intro_clust == Mf & Intro_clust == Mn, "Mf_Mn",
        ifelse(Cluster == "Peninsular" & (Intro_clust != Mn | Intro_clust != Mf | Intro_clust != Peninsular), "Peninsular", .$Intro_clust)))))) 

#write_tsv(introgression_table, "Introgression/introgression_table.tsv")
introgression_table <- read_tsv("Introgression/introgression_table.tsv")

# genetic distance per sliding window of each sample to each cluster (dominant allele) - the number of times the SNP is not the dominant allele of that cluster

window_size <- 10000 # should give us ~180 windows per sample per chrom

#WINDOWS Table for downstream
PKA1H1_windows <- read.table("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/ref_genomes/PKA1H1/strain_A1_H.1.Icor.fasta.bed") %>% 
    mutate(CHROM = str_remove(V1, "ordered_PKNH_")) %>% 
    mutate(CHROM = str_remove(CHROM, "_v2")) %>% 
    rename(Start = V2, End = V3) %>% 
    select(-V1) %>% 
    filter(CHROM != "PKNH_MIT", CHROM != "new_API_strain_A1_H.1") 

library(plyr)
PKA1H1_windows <- ddply(PKA1H1_windows, "CHROM", summarise, POS = seq(Start, End)) %>% 
    mutate(TMP_WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% 
    mutate(WINDOW = group_indices(., CHROM, TMP_WINDOW)) %>% 
    select(-TMP_WINDOW)
detach("package:plyr", unload=TRUE)

#write_tsv(PKA1H1_windows , "Introgression/PKA1H1_windows.tsv")
#PKA1H1_windows <- read_tsv("Introgression/PKA1H1_windows.tsv")


introgression_table_window <- introgression_table %>% 
    left_join(PKA1H1_windows) %>% 
    filter(SNP >= 0) %>% # filter out missing calls - messes with the distance calculation
    group_by(SAMPLE, WINDOW) %>%
    summarise(Mf_distance = sum(SNP != Mf, na.rm = TRUE)/n()*100, # the number of times the sample allele doesn't match the dominant allele = genetic distance
        Mn_distance = sum(SNP != Mn, na.rm = TRUE)/n()*100,
        Pen_distance = sum(SNP != Peninsular, na.rm = TRUE)/n()*100)  %>%
    left_join(
        metadata %>% rename(SAMPLE = Sample)
    ) %>%
    filter(SAMPLE %in% metadata$Sample)

#write_tsv(introgression_table_window, "Introgression/introgression_table_window.tsv")
#introgression_table_window <- read_tsv("Introgression/introgression_table_window.tsv")



######################################## CREATE POLYGON PLOTS TO REPRESENT INTROGRESSION PER SAMPLE ########################################






# Plot function
## for unfiltered data - entire genome (doesn't contain point labels - too messy)
genetic_distance_plot <- function(data, sample){
subset <- data %>%  # for title
    filter(SAMPLE == sample)

data %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = introgression_table_window %>%
        filter(SAMPLE == sample),
        size = 0.75, alpha = 0.75) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    ggtitle(paste0(sample, ", ", subset[1,"Cluster"], ", ", subset[1,"Location"])) 
}

## for data filtered to a particular CHROM (contain point labels with window number) - filtering must be done prior - line 97
genetic_distance_plot_chrom <- function(data, sample){

subset <- data %>%  # subset to deal with smaller dataset
    filter(SAMPLE == sample)

max_value <- subset %>% # get max values to use for x and y limits
    summarise(Mf = max(Mf_distance), Mn = max(Mn_distance)) 

subset %>% # create plot
    ggplot(aes(x = Mf_distance, y = Mn_distance, colour = Cluster)) +
        geom_point() +
        ggtitle(paste0(sample, ", ", subset[1,"Cluster"], ", ", subset[1,"Location"])) + # extract cluster and location information to include in title
        ylim(0, ifelse(max_value$Mf > max_value$Mn, max_value$Mf + 5, max_value$Mn + 5)) + # use the max values calculated above to inform the axis values
        xlim(0, ifelse(max_value$Mf > max_value$Mn, max_value$Mf + 5, max_value$Mn + 5)) +
        theme(legend.position = "none") +
        scale_colour_manual(values = ifelse(subset$Cluster == "Mf", "#440154FF", ifelse(subset$Cluster == "Mn", "#39568CFF", "#73D055FF") )) + # colour based on cluster
        geom_vline(xintercept = 100, linetype = "dashed", size = 1) +
        geom_hline(yintercept = 100, linetype = "dashed", size = 1) +
        geom_text(aes(label = ifelse(Mn_distance > 100, as.character(WINDOW), "")), hjust = 0, vjust = 2) + # add 'window' text to values in outlier quandrants
        geom_text(aes(label = ifelse(Mf_distance > 100, as.character(WINDOW), "")), hjust = 0, vjust = 2)
}


# Create a loop to plot all samples that are Mn or Mf
sample_names <- introgression_table_window %>% 
    filter(Cluster == "Mn" | Cluster == "Mf") %>%
    dplyr::select(SAMPLE) %>%
    unique() %>%
    mutate(SAMPLE = as.factor(SAMPLE)) 
     

for(i in levels(sample_names$SAMPLE)){
    distance_plot <- genetic_distance_plot(introgression_table_window, paste0(i))
    ggsave(paste0("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/",i,"MfvsMn_dot.png"), dpi = 300, distance_plot)
}      





######################################## USE CONTOURS TO IDENTIFY INTROGRESSED REGIONS ########################################






# Function for identifying regions on introgression from a contours and point plot 
find_introgressed_regions <- function(SAMPLENAME){
# libaries
library(sp) 
library(MASS)

# Create plot with points and polygons
raster_plot <- introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = introgression_table_window %>%
        filter(SAMPLE == SAMPLENAME),
        size = 0.75, alpha = 0.75) +
    geom_density_2d(mapping = aes(x = Mf_distance, y = Mn_distance, colour = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    scale_colour_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 

# Extract info on points and polygons from ggplot
polygon_data <- ggplot_build(raster_plot) # create matrix for contours based on plot above - includes all datapoints
polygon_points <- polygon_data$data[[1]]
polygon_contours <- polygon_data$data[[2]] # extract contour data

# Filter to cluster-specific polygon/contours and identify whether or not points fall within that cluster
Mf_contours <- polygon_contours %>% 
    filter(colour == "#440154FF") %>% 
    filter(group == "1-002-001")

Mf_contours <- point.in.polygon(pol.x = as.numeric(unlist(Mf_contours[3])), 
    pol.y = as.numeric(unlist(Mf_contours[4])),
    point.x = as.numeric(unlist(polygon_points[1])), 
    point.y = as.numeric(unlist(polygon_points[2]))) %>% 
    as.data.frame() %>% 
    rename("Mf" = ".")

Mn_contours <- polygon_contours %>% 
    filter(colour == "#73D055FF") %>% 
    filter(group == "2-002-001")

Mn_contours <- point.in.polygon(pol.x = as.numeric(unlist(Mn_contours[3])), 
    pol.y = as.numeric(unlist(Mn_contours[4])),
    point.x = as.numeric(unlist(polygon_points[1])), 
    point.y = as.numeric(unlist(polygon_points[2]))) %>% 
    as.data.frame() %>% 
    rename("Mn" = ".")

Pen_contours <- polygon_contours %>% 
    filter(colour == "#39568CFF") %>% 
    filter(group == "3-002-001") 

Pen_contours <- point.in.polygon(pol.x = as.numeric(unlist(Pen_contours[3])), 
    pol.y = as.numeric(unlist(Pen_contours[4])),
    point.x = as.numeric(unlist(polygon_points[1])), 
    point.y = as.numeric(unlist(polygon_points[2]))) %>% 
    as.data.frame() %>% 
    rename("Pen" = ".")

# Combine newly generated point presence-abscence info together
points_in_contours_data <- Mf_contours %>% 
    cbind(Mn_contours) %>% 
    cbind(Pen_contours) 

# Cmbine with original introgression table to filter to "introgressed" points - those that fall within another clusters polygon but not its own
subset <- introgression_table_window %>%
    filter(SAMPLE == SAMPLENAME) %>% 
    cbind(points_in_contours_data) %>% 
    filter(ifelse(Cluster == "Mn", Mf == 1 & Mn == 0 & Pen == 0, # Mf to Mn introgression
        ifelse(Cluster == "Mf", Mf == 0 & Mn == 1 & Pen == 0, # Mn to Mf introgression
        ifelse(Cluster == "Peninsular", Mf == 1 & Mn == 1 & Pen == 0)))) # Mn or Mf to Peninsular introgression

return(subset)
}


# Create a loop to identify introgressed windows
## extract sample names
sample_names <- introgression_table_window %>% 
    dplyr::select(SAMPLE) %>%
    unique() %>%
    mutate(SAMPLE = as.factor(SAMPLE)) 
     
## create empty df
introgressed_windows <-  introgression_table_window %>% 
    add_column(Mf = 1, Mn = 1, Pen = 1) %>% 
    slice(-(1:nrow(.)))

## loop into df
for(i in levels(sample_names$SAMPLE)){
    introgressed_windows <- introgressed_windows %>%
        rbind(find_introgressed_regions(paste0(i)))
}      


#write_tsv(introgressed_windows, "/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/introgressed_windows.tsv")
#introgressed_windows <- read_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/introgressed_windows.tsv")



######################################## APPLY FILTERS ########################################



# Filter out low n windows - full datatset

## shoulder plot to come up with n threshold 
shoulder_plot <- introgressed_windows %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>%
    arrange(n) %>% 
    add_column(WINDOWS = 1:nrow(.)) %>%
    ggplot(aes(x = WINDOWS, y = n)) +
    geom_col() + 
    scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40))

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/shoulder_plot.png", dpi = 300,  shoulder_plot)

# filter out windows that only appear in a single sample - removes ~500 windows
introgressed_windows_filter <- introgressed_windows %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()/152*100) %>% 
    filter(n > 5) 

introgressed_windows <- introgressed_windows %>% 
    filter(WINDOW %in% introgressed_windows_filter$WINDOW) 

# Filter out low n windows - per cluster

Mf_filter <- introgressed_windows %>%
    filter(Cluster == "Mf") %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()/81*100) %>% 
    filter(n > 5) 


Mn_filter <- introgressed_windows %>%
    filter(Cluster == "Mn") %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% # total numbers are too low to use proportions - 5% is one sample
    filter(n > 5) 

Pen_filter <- introgressed_windows %>%
    filter(Cluster == "Peninsular") %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% # total numbers are too low to use proportions - 5% is 1.5 samples
    filter(n > 5) 

introgressed_windows <- introgressed_windows %>% 
    filter(Cluster == "Mf") %>% 
    filter(WINDOW %in% Mf_filter$WINDOW) %>% 
    rbind(
        introgressed_windows %>% 
            filter(Cluster == "Mn") %>%
            filter(WINDOW %in% Mn_filter$WINDOW)
    ) %>%
    rbind(
        introgressed_windows %>% 
            filter(Cluster == "Peninsular") %>%
            filter(WINDOW %in% Pen_filter$WINDOW)
    )


# Filter out regions with SICAvar or Kir genes

## Read in GFF file
GFF_annotation <- read.table("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/ref_genomes/PKA1H1/gff/updated_version/PlasmoDB-63_PknowlesiA1H1.gff_names_updated", skip = 18, sep = "\t") %>% 
    rename(CHROM = V1,
    source = V2,
    feature= V3,
    start = V4,
    end = V5,
    score = V6,
    strand = V7,
    frame = V8,
    attribute = V9) %>% 
    filter(grepl("ordered", CHROM)) %>% 
    mutate(CHROM = str_remove(CHROM, "ordered_PKNH_")) %>% 
    mutate(CHROM = str_remove(CHROM, "_v2"))

SICAvar_KIR <- GFF_annotation %>% 
    filter(grepl("SICA", attribute)) %>% 
    rbind(
        GFF_annotation %>% 
            filter(grepl("KIR", attribute))
    ) %>% 
    as_tibble() %>% 
    pivot_longer(c("start", "end"), names_to = "range", values_to = "POS") %>% 
    unite(CHROM_POS, c("CHROM", "POS"), sep = "_") %>% 
    select(CHROM_POS) %>%
    unique()

SICAvar_KIR <- PKA1H1_windows %>% 
    unite(CHROM_POS, c("CHROM", "POS"), sep = "_") %>% 
    filter(CHROM_POS %in% SICAvar_KIR$CHROM_POS)

introgressed_windows <- introgressed_windows %>% 
    filter(!(WINDOW %in% SICAvar_KIR$WINDOW)) # filter out known hypervariable regions - removes another ~1K windows

# Filter windows that appear in multiple clusters as potential introgression events - suggests hypervariability

hypervariable_filter <- introgressed_windows %>%
    select(WINDOW, Cluster) %>% 
    unique() %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    filter(n > 1)

introgressed_windows <- introgressed_windows %>% 
    filter(!(WINDOW %in% hypervariable_filter$WINDOW)) # filter out unknown hypervariable regions - those windows that appear in multiple clusters


#write_tsv(introgressed_windows, "/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/introgressed_windows_filtered.tsv")


######################################## Explore windows of introgression ########################################





# Compare degree of introgression between different clusters
introgressed_windows %>%
    dplyr::select(SAMPLE, Cluster) %>%
    group_by(Cluster, SAMPLE) %>% 
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(Cluster) %>% 
    summarise(mean = mean(n), median = median(n))

# Add state and district info
metadata <- read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.csv", col_names = FALSE) %>% 
        dplyr::select(1, 6, 7) %>%
        rename(Sample = X1) %>%
        rename(Location = X6) %>%
        rename(Cluster = X7) %>% 
        as.data.frame() 

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

# add state and district info
metadata <- metadata %>%
  mutate(State = ifelse(Cluster == "Peninsular", "Peninsular", Location)) %>%
  mutate(State = ifelse(State != "Sabah" & State != "Peninsular", "Sarawak", State)) %>%
  mutate(Sample2 = Sample) %>%
  mutate(Sample = str_remove(Sample2, "_DK.*")) %>%
  left_join(
    read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD_geo_clusters/districts_metadata.csv") %>%
    dplyr::select(Sample, district)
  ) %>%
  mutate(district = ifelse(is.na(district), Location, district)) %>%
  dplyr::select(-Sample) %>%
  rename(SAMPLE = Sample2) 

introgressed_windows_meta <- introgressed_windows %>%
    left_join(metadata)

# Summary info for TOTAL introgression across states and districts
introgressed_windows_meta %>% 
    dplyr::select(SAMPLE, district, State) %>%
    group_by(district, SAMPLE, State) %>% 
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(district, State) %>% 
    summarise(mean = mean(n), median = median(n), sd = sd(n), n = n()) %>% # number of windows
    arrange(desc(median)) #%>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/average_windows_for_districts.tsv")

# Summary info for TOTAL introgression across clusters
introgressed_windows_meta %>% 
    dplyr::select(SAMPLE, Cluster) %>%
    group_by(SAMPLE, Cluster) %>% 
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(Cluster) %>% 
    summarise(mean = mean(n), median = median(n), sd = sd(n), n = n()) %>% # number of windows
    arrange(desc(median)) #%>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/average_windows_for_clusters.tsv")

# Number of introgression events per sample
introgressed_windows_meta %>% 
    group_by(SAMPLE) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    mutate(quantile = ntile(n, 3)) 

# Summary by groupings - proprtion of introgression levels per group
intro_per_sample_summary <- introgressed_windows_meta %>% 
    group_by(SAMPLE) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    mutate(intro_level = ntile(n, 3)) %>% 
    mutate(intro_level = ifelse(intro_level == 3, "High", ifelse(intro_level == 2, "Medium", "Low"))) %>% 
    left_join(metadata) %>%
    dplyr::select(-c(Location)) #%>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/intro_per_sample_summary.tsv")

intro_per_sample_summary %>% 
    dplyr::select(3:ncol(.)) %>% 
    group_by(Cluster, intro_level) %>% 
    summarise(n = n()) %>% 
    rename(Group = Cluster) %>% 
    add_column(Variable = "Cluster") %>%
    rbind(
        intro_per_sample_summary %>% 
            dplyr::select(3:ncol(.)) %>% 
            group_by(State, intro_level) %>% 
            summarise(n = n()) %>% 
            rename(Group = State) %>% 
            add_column(Variable = "State") 
    ) %>% 
    rbind(
        intro_per_sample_summary %>% 
            dplyr::select(3:ncol(.)) %>% 
            group_by(district, intro_level) %>% 
            summarise(n = n()) %>% 
            rename(Group = district) %>% 
            add_column(Variable = "District")
    ) %>% 
    ungroup() %>% 
    pivot_wider(names_from = intro_level, values_from = n, values_fn = mean) %>% 
    relocate(Group, Variable, High, Medium, Low) #%>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/levels_of_introgression_summary.tsv")






######################## Introgressed windows for each cluster





######################## Mf samples
introgression_of_Mn_into_Mf <- introgressed_windows %>% # position of windows
    filter(Cluster == "Mf" & Mf == 0 & Mn == 1) %>% 
    dplyr::select(WINDOW) %>% 
    unique() %>% 
    dplyr::left_join(PKA1H1_windows)

### introgressed windows per chrom
introgression_of_Mn_into_Mf %>% 
    dplyr::select(WINDOW, CHROM) %>% 
    unique() %>% 
    group_by(CHROM) %>% 
    summarise(n_windows = n()) %>% 
    arrange(desc(n_windows)) #%>% write_tsv("Introgression/Mf_windows_across_chrom.tsv")

# Summary table of samples and windows - presence/absence of introgression
Mf_introgression_matrix <- introgressed_windows %>% 
    filter(Cluster == "Mf" & Mf == 0 & Mn == 1)  %>%
    dplyr::select(SAMPLE, WINDOW) %>% 
    add_column(PRESENCE = 1) %>% 
    mutate(WINDOW = as.factor(as.character(WINDOW))) %>%
    pivot_wider(names_from = SAMPLE, values_from = PRESENCE) %>% 
    replace(is.na(.), 0 )

# looking into CHROM 08
## counts of samples per window, with CHROM and POS
chrom_08 <- introgression_of_Mn_into_Mf %>% 
    filter(CHROM == "08") %>% 
    arrange() %>%
    group_by(WINDOW) %>% 
    summarise(min = min(POS), max = max(POS)) %>% 
    left_join(
        introgressed_windows %>% 
            filter(Cluster == "Mf" & Mf == 0 & Mn == 1) %>% 
            dplyr::select(WINDOW) %>% 
            group_by(WINDOW) %>% 
            summarise(n = n())
    ) #%>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/CHROM_08_Introgression_Mn_into_Mf.tsv")

# looking into CHROM 11
chrom_11 <- introgression_of_Mn_into_Mf %>% 
    filter(CHROM == "11") %>% 
    arrange() %>%
    group_by(WINDOW) %>% 
    summarise(min = min(POS), max = max(POS)) %>% 
    left_join(
        introgressed_windows %>% 
            filter(Cluster == "Mf" & Mf == 0 & Mn == 1) %>% 
            dplyr::select(WINDOW) %>% 
            group_by(WINDOW) %>% 
            summarise(n = n())
    ) #%>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/CHROM_11_Introgression_Mn_into_Mf.tsv")

# Plot across chromosome 
CHROM_plot <- chrom_08 %>% 
    add_column(Cluster = "Mf") %>%
    ggplot(aes(xmin = min/1000000, xmax = max/1000000, ymin = 0.25, ymax = 0.75, fill = Cluster)) + 
        geom_hline(yintercept = 0.5, size = 5, colour = "#440154FF") + 
        geom_rect(stat = "identity") +
        ylim(0.2, 0.8) + 
        scale_fill_manual(values= c ("#73D055FF")) + 
        theme(axis.title.y = element_blank(), 
            axis.text.y=element_blank(), 
            axis.ticks.y = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = "none") +
        xlab("Chromosome 8 (Mb)") 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/chromosome_8.png", dpi = 300, width = 10, height = 3, CHROM_plot)

CHROM_plot <- chrom_11 %>% 
    add_column(Cluster = "Mf") %>%
    ggplot(aes(xmin = min/1000000, xmax = max/1000000, ymin = 0.25, ymax = 0.75, fill = Cluster)) + 
        geom_hline(yintercept = 0.5, size = 5, colour = "#440154FF") + 
        geom_rect(stat = "identity") +
        ylim(0.2, 0.8) + 
        scale_fill_manual(values= c ("#73D055FF")) + 
        theme(axis.title.y = element_blank(), 
            axis.text.y=element_blank(), 
            axis.ticks.y = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = "none") +
        xlab("Chromosome 11 (Mb)") 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/chromosome_11.png", dpi = 300, width = 11.4, height = 3, CHROM_plot)





# Subset metadata
introgressed_windows_meta_Mn_to_Mf <- introgressed_windows_meta %>% # metadata of Mn to Mf introgression
    filter(Cluster == "Mf" & Mf == 0 & Mn == 1)


# How many samples have introgressed windows on chrom 8 and 11?
CHROM_filter_08 <- introgression_of_Mn_into_Mf %>% 
    filter(CHROM == "08") %>% 
    dplyr::select(WINDOW) %>%
    unique()

introgressed_windows_meta_Mn_to_Mf %>% 
    filter(WINDOW %in% CHROM_filter_08$WINDOW) %>% 
    dplyr::select(SAMPLE) %>% 
    unique()

### min number of windows in a sample
introgressed_windows_meta_Mn_to_Mf %>% 
    filter(WINDOW %in% CHROM_filter_08$WINDOW) %>% 
    group_by(SAMPLE) %>% 
    summarise(n = n()) %>% 
    arrange(n) %>% 
    filter(n > 1)

CHROM_filter_11 <- introgression_of_Mn_into_Mf %>% 
    filter(CHROM == "11") %>% 
    dplyr::select(WINDOW) %>%
    unique()

introgressed_windows_meta_Mn_to_Mf %>% 
    filter(WINDOW %in% CHROM_filter_11$WINDOW) %>% 
    dplyr::select(SAMPLE) %>% 
    unique()

# Sabah windows
introgressed_windows_meta_Mn_to_Mf  %>% 
    mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
    group_by(State, WINDOW) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = State, values_from = n) %>% 
    arrange(desc(Sabah)) %>% 
    slice(1:10) %>% 
    left_join(
        introgression_of_Mn_into_Mf %>% 
            group_by(WINDOW, CHROM) %>% 
            summarise(start = min(POS), end = max(POS))
    ) #%>% write_tsv("Introgression/Mf_major_windows_Sabah_samples.tsv")

# Sarawak samples
introgressed_windows_meta_Mn_to_Mf  %>% 
    mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
    group_by(State, WINDOW) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = State, values_from = n) %>% 
    arrange(desc(Sarawak)) %>% 
    slice(1:10) %>% 
    left_join(
        introgression_of_Mn_into_Mf %>% 
            group_by(WINDOW, CHROM) %>% 
            summarise(start = min(POS), end = max(POS))
    ) #%>% write_tsv("Introgression/Mf_major_windows_Sarawak_samples.tsv")

# Windows found in both states
introgressed_windows_meta_Mn_to_Mf  %>% 
    mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
    group_by(State, WINDOW) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = State, values_from = n) %>% 
    filter(Sabah >= 1 & Sarawak >= 1) 

# district breakdown of interesting windows 
introgressed_windows_meta_Mn_to_Mf  %>% 
    mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
    group_by(State, district, WINDOW) %>% 
    summarise(n = n()) %>% 
    filter(WINDOW == 854)


######################## Mn Samples
introgression_of_Mf_into_Mn <- introgressed_windows %>% # position of windows
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>%
    dplyr::select(WINDOW) %>% 
    unique() %>% 
    dplyr::left_join(PKA1H1_windows)

### introgressed windows per chrom
introgression_of_Mf_into_Mn  %>% 
    dplyr::select(WINDOW, CHROM) %>% 
    unique() %>% 
    group_by(CHROM) %>% 
    summarise(n_windows = n()) %>% 
    arrange(desc(n_windows))

# Summary table of samples and windows - presence/absence of introgression
Mn_introgression_matrix <- introgressed_windows %>% 
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>%
    dplyr::select(SAMPLE, WINDOW) %>% 
    add_column(PRESENCE = 1) %>% 
    mutate(WINDOW = as.factor(as.character(WINDOW))) %>%
    pivot_wider(names_from = SAMPLE, values_from = PRESENCE) %>% 
    replace(is.na(.), 0 )

# Subset metadata
introgressed_windows_meta_Mf_to_Mn <- introgressed_windows_meta %>% # metadata of Mn to Mf introgression
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) 

# Windows with greatest n
introgressed_windows_meta_Mf_to_Mn %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n))

# Sabah windows
introgressed_windows_meta_Mf_to_Mn  %>% 
    mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
    group_by(State, WINDOW) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = State, values_from = n) %>% 
    arrange(desc(Sabah)) %>% 
    slice(1:10) %>% 
    left_join(
        introgression_of_Mf_into_Mn %>% 
            group_by(WINDOW, CHROM) %>% 
            summarise(start = min(POS), end = max(POS))
    ) #%>% write_tsv("Introgression/Mn_major_windows_Sabah_samples.tsv")

# Sarawak samples
introgressed_windows_meta_Mf_to_Mn  %>% 
    mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
    group_by(State, WINDOW) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = State, values_from = n) %>% 
    arrange(desc(Sarawak)) %>% 
    slice(1:10) %>% 
    left_join(
        introgression_of_Mf_into_Mn %>% 
            group_by(WINDOW, CHROM) %>% 
            summarise(start = min(POS), end = max(POS))
    ) #%>% write_tsv("Introgression/Mn_major_windows_Sarawak_samples.tsv")

# district breakdown of interesting windows 
introgressed_windows_meta_Mf_to_Mn  %>% 
    mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
    group_by(State, district, WINDOW) %>% 
    summarise(n = n()) %>% 
    filter(WINDOW == 1347)


# Windows with greatest n for both clades      
window_freq <- introgressed_windows_meta_Mn_to_Mf  %>% 
    mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
    group_by(State, WINDOW) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = State, values_from = n) %>% 
    mutate(Sabah = replace_na(Sabah, 0)) %>%
    mutate(Sarawak = replace_na(Sarawak, 0)) %>%
    mutate(Total = Sabah + Sarawak) %>% 
    arrange(desc(Total)) %>% 
    left_join(introgression_of_Mn_into_Mf %>% 
        group_by(WINDOW, CHROM) %>% 
        summarise(start_POS = min(POS), end_POS = max(POS))) %>% 
    add_column(Cluster = "Mf") %>%
    rbind(
        introgressed_windows_meta_Mf_to_Mn  %>% 
            mutate(Sabah = ifelse(State == "Sabah", "Sabah", "Other")) %>% 
            group_by(State, WINDOW) %>% 
            summarise(n = n()) %>% 
            pivot_wider(names_from = State, values_from = n) %>% 
            mutate(Sabah = replace_na(Sabah, 0)) %>%
            mutate(Total = Sabah + Sarawak) %>% 
            arrange(desc(Total)) %>% 
            left_join(introgression_of_Mf_into_Mn %>% 
                group_by(WINDOW, CHROM) %>% 
                summarise(start_POS = min(POS), end_POS = max(POS))) %>% 
            add_column(Cluster = "Mn") 
        ) %>% 
    relocate(WINDOW, CHROM, start_POS, end_POS, Cluster) %>% 
    arrange(Cluster, WINDOW) #%>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Windows_with_greatest_freq.tsv")


window_freq_plot <- window_freq %>% 
    mutate(WINDOW = as.factor(WINDOW)) %>%
    filter(Cluster == "Mf") %>%
    pivot_longer(c(Sabah, Sarawak), names_to = "State", values_to = "Count") %>% 
    ggplot(aes(x = WINDOW, y = Count, fill = State)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
        theme_dark() +
        xlab("Window") +
        ylab("Sample count")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_frequency_across_states_Mf.png", dpi = 300, window_freq_plot)


window_freq_plot <- window_freq %>% 
    mutate(WINDOW = as.factor(WINDOW)) %>%
    filter(Cluster == "Mn") %>%
    pivot_longer(c(Sabah, Sarawak), names_to = "State", values_to = "Count") %>% 
    ggplot(aes(x = WINDOW, y = Count, fill = State)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
        theme_dark() +
        xlab("Window") +
        ylab("Sample count")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_frequency_across_states_Mn.png", dpi = 300, window_freq_plot)

# Does frequency of samples increase over time for any of the windows?
## Mf
introgressed_windows_meta_Mn_to_Mf  %>%
    mutate(SAMPLE = str_remove(SAMPLE, "_DK.*")) %>%
    left_join(
        read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD_geo_clusters/districts_metadata.csv")  %>% 
            dplyr::select(Sample, date) %>%
            rename(SAMPLE = Sample)
    ) %>% 
    filter(grepl("PK", SAMPLE)) %>% 
    mutate(year = lubridate::year(date)) %>% 
    group_by(WINDOW, year) %>% 
    summarise(intro_n = n()) %>% 
    left_join(
        introgressed_windows_meta %>%
        filter(Cluster == "Mf") %>%
        mutate(SAMPLE = str_remove(SAMPLE, "_DK.*")) %>%
        left_join(
            read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD_geo_clusters/districts_metadata.csv")  %>% 
                dplyr::select(Sample, date) %>%
                rename(SAMPLE = Sample)
        ) %>% 
        filter(grepl("PK", SAMPLE)) %>% 
        mutate(year = lubridate::year(date)) %>% 
        dplyr::select(SAMPLE, year) %>% 
        unique() %>%
        group_by(year) %>% 
        summarise(total_n = n())
    ) %>%
    mutate(prop_intro = intro_n/total_n*100) %>% 
    filter(total_n > 1) %>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mf_freq_over_years.tsv")

## Mn

  

# Plot major windows
## Plot function
genetic_distance_window_plot <- function(data, window){
data %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    geom_point(data = introgression_table_window %>%
        filter(WINDOW == window),
        size = 1, alpha = 0.75, aes(shape = Cluster)) +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100))
}

# Mf 
distance_plot <- genetic_distance_window_plot(introgression_table_window, 1504)
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_plot_1504_Mf.png", dpi = 300, distance_plot)

distance_plot <- genetic_distance_window_plot(introgression_table_window, 852)
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_plot_852_Mf.png", dpi = 300, distance_plot)

distance_plot <- introgression_table_window %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    geom_point(data = introgression_table_window %>%
        filter(WINDOW == 1440 | WINDOW == 1441),
        size = 1, alpha = 0.75, aes(shape = Cluster)) +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_plot_1440_1441_Mf.png", dpi = 300, distance_plot)

# Mn
distance_plot <- genetic_distance_window_plot(introgression_table_window, 1347)
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_plot_1347_Mn.png", dpi = 300, distance_plot)

distance_plot <- genetic_distance_window_plot(introgression_table_window, 557)
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_plot_557_Mn.png", dpi = 300, distance_plot)

distance_plot <- genetic_distance_window_plot(introgression_table_window, 1869)
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_plot_1869_Mn.png", dpi = 300, distance_plot)



## counts of windows per sample - high counts should be outliers in NJT
### chrom 08
Mf_introgression_matrix %>% 
    filter(WINDOW %in% CHROM_filter_08$WINDOW) %>%
    pivot_longer(2:ncol(.), names_to = "SAMPLE", values_to = "INTROGRESSION") %>% 
    group_by(SAMPLE) %>% 
    summarise(n_windows = sum(INTROGRESSION == 1)) %>% 
    arrange(desc(n_windows))

### chrom 11
Mf_introgression_matrix %>% 
    filter(WINDOW %in% CHROM_filter_11$WINDOW) %>%
    pivot_longer(2:ncol(.), names_to = "SAMPLE", values_to = "INTROGRESSION") %>% 
    group_by(SAMPLE) %>% 
    summarise(n_windows = sum(INTROGRESSION == 1)) %>% 
    arrange(desc(n_windows))


# NJT of regions of interest 
library(ape)
library(ggtree)
library(tidyverse)

# Neighbour-joining tree

## Create distance matrix from PLINK data and build NJT - subset to this cluster of windows: 880027 - 1339957
NJT_ID <- read_table("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

## Plot tree - large region on chrom 08
options(ignore.negative.edge=TRUE)

NJT_metadata <- metadata %>% 
    relocate(SAMPLE) %>% 
    mutate(SAMPLE = str_remove(SAMPLE, "_DK.*"))

NJT_tree_plot <- ggtree(NJT_tree, layout = "daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#73D055FF", "#39568CFF"))
    
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/tree_plots/NJT_tree_chrom_08_unrooted.png", dpi = 300, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
    legend.title = element_blank(), 
    legend.key = element_blank()) +
    geom_tiplab(size = 2) +
    scale_color_manual(values = c("#440154FF", "#73D055FF", "#39568CFF"))

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/tree_plots/NJT_tree_chrom_08_unrooted_labelled.png", dpi = 300, height = 20, width = 20, limitsize = FALSE, NJT_tree_plot)

# Do "jumping" samples align with introgression results?
tree_windows <- introgression_of_Mn_into_Mf %>% 
    filter(CHROM == "08") %>% 
    filter(POS > 880027	& POS < 1339957) %>% 
    dplyr::select(WINDOW) %>%
    unique()

introgressed_windows_meta_Mn_to_Mf %>% 
    filter(WINDOW %in% tree_windows$WINDOW) %>% 
    group_by(SAMPLE) %>% 
    summarise(n_windows = n()) %>% 
    arrange(desc(n_windows)) %>% 
    head(n = 14)



## Plot tree - individual windows
## get inputs for filters to be used in trees
### data needs to be generated in PLINK for the specific windows
tree_windows <- read_tsv("Introgression/Mf_major_windows_Sabah_samples.tsv") %>% 
    rbind(
        read_tsv("Introgression/Mf_major_windows_Sarawak_samples.tsv")
    ) %>% 
    rbind(
        read_tsv("Introgression/Mn_major_windows_Sarawak_samples.tsv")
    ) %>% 
    rbind(
        read_tsv("Introgression/Mn_major_windows_Sabah_samples.tsv")
    ) %>% 
    mutate(start = start - 20000) %>% 
    mutate(end = end + 20000) %>% 
    select(1,4:6) %>% 
    mutate(CHROM = paste0("ordered_PKNH_",CHROM,"_v2")) %>% 
    mutate(WINDOW = as.factor(WINDOW)) %>% # needed to paste into loop downstream
    write_tsv("Introgression/tree_windows.tsv")

generate_tree_plots <- function(WINDOW){
library(ape)
library(ggtree)
library(tidyverse)

NJT_metadata <- metadata %>% 
    relocate(SAMPLE) %>% 
    mutate(SAMPLE = str_remove(SAMPLE, "_DK.*"))

options(ignore.negative.edge=TRUE)

NJT_ID <- read_table(paste0("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/NJT_data/",WINDOW,".NJT.dist.id"), col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table(paste0("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/NJT_data/",WINDOW,".NJT.dist"), col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

NJT_tree_plot <- ggtree(NJT_tree, layout = "daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#73D055FF", "#39568CFF"))
    
ggsave(paste0("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/tree_plots/NJT_tree_window_",WINDOW,"_unrooted.png"), dpi = 300, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
    legend.title = element_blank(), 
    legend.key = element_blank()) +
    geom_tiplab(size = 2) +
    scale_color_manual(values = c("#440154FF", "#73D055FF", "#39568CFF"))

ggsave(paste0("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/tree_plots/NJT_tree_window_",WINDOW,"_unrooted_labelled.png"), dpi = 300, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)
}


for(i in levels(tree_windows$WINDOW)){
    generate_tree_plots(paste0(i))
}      



# SNP barcode plots
SNP_barcode_table <- genotype_table_long %>% 
    select(-Cluster) %>% 
    left_join(metadata 
        ) %>% 
    filter(!is.na(Cluster)) %>% # those samples that have been filtered out for various reasons wont have this data
  mutate(SNP = ifelse(SNP < 0, NA, .$SNP)) # don't want missing data to be factored into clustering, so, convert -1 to NA

############################## FUNCTION
clustered_snp_barcode <- function(WINDOW_NUM, CHROM_NUM, START, END, CLUSTER){
  SNP_barcode_table <- SNP_barcode_table %>%
    filter(CHROM == paste0(CHROM_NUM) & POS > START	& POS < END) %>% 
    mutate(SNP = ifelse(SNP < 0, NA, .$SNP))

  clust_data <- SNP_barcode_table %>% 
    select(CHROM, POS, SAMPLE, SNP) %>% 
    pivot_wider(names_from = SAMPLE, values_from = SNP) 
  
  clust <- clust_data %>%
    select(-c(1:2)) %>% 
    t() %>% # each row is a sample and each column a SNP positon 
    dist() %>%
    hclust() #euclidean distance
 
  sample_order <- as.data.frame(clust$labels[c(clust$order)]) %>% 
    rename(SAMPLE = "clust$labels[c(clust$order)]") %>% 
    filter(SAMPLE != "CHROM" & SAMPLE != "POS") %>% 
    add_column(SAMPLE_2 = 1:nrow(.)) %>% 
    mutate(SAMPLE_2 = as.factor(SAMPLE_2))

  SNP_barcode_table <- SNP_barcode_table %>% 
    left_join(sample_order)
  
  introgressed_samples <- introgressed_windows  %>% 
    filter(WINDOW == WINDOW_NUM & Cluster == paste0(CLUSTER)) 
    
  SNP_barcode_table %>% 
    ggplot(aes(x = POS, y = SAMPLE_2, fill = Cluster, alpha = SNP, group = Cluster)) +
      geom_col() +
      scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      scale_alpha(guide = 'none') +
      ylab("Sample") +
      xlab("Chromosome Position") 
}

snp_barcode_plot <- clustered_snp_barcode(826, "08", 950000, 959999, "Mf")
ggsave("Introgression/snp_barcodes/oocyte_window_snp_barcode_clustered.png", dpi = 300, width = 10, snp_barcode_plot)

snp_barcode_plot <- clustered_snp_barcode(1504, "11", 2080000, 2089999, "Mf")
ggsave("Introgression/snp_barcodes/1504_snp_barcode_clustered.png", dpi = 300, width = 10, snp_barcode_plot)

snp_barcode_plot <- clustered_snp_barcode(1440, "11", 1440000, 1449999, "Mf")
ggsave("Introgression/snp_barcodes/1440_snp_barcode_clustered.png", dpi = 300, width = 10, snp_barcode_plot)









## Mf to Mn 
######################## Mf to Mn
introgressed_windows_filter <- introgressed_windows %>% # filter for n samples > 1 for each window
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    filter(n > 1) 

introgression_of_Mf_into_Mn <- introgressed_windows %>% # position of windows
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>% 
    filter(WINDOW %in% introgressed_windows_filter$WINDOW) %>%
    dplyr::select(WINDOW) %>% 
    unique() %>% 
    dplyr::left_join(PKA1H1_windows)

### introgressed windows per chrom
introgression_of_Mf_into_Mn %>% 
    dplyr::select(WINDOW, CHROM) %>% 
    unique() %>% 
    group_by(CHROM) %>% 
    summarise(n_windows = n()) %>% 
    arrange(desc(n_windows))

# looking into CHROM 08
introgression_of_Mf_into_Mn %>% 
    filter(CHROM == "08") %>% 
    arrange() %>%
    group_by(WINDOW) %>% 
    summarise(min = min(POS), max = max(POS)) %>% 
    left_join(
        introgressed_windows %>% 
            filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>% 
            filter(WINDOW %in% introgressed_windows_filter$WINDOW) %>%
            dplyr::select(WINDOW) %>% 
            group_by(WINDOW) %>% 
            summarise(n = n())
    ) #%>% write_tsv("CHROM_08_Introgression_Mf_into_Mn.tsv")

# looking into CHROM 11
introgression_of_Mf_into_Mn %>% 
    filter(CHROM == "11") %>% 
    arrange() %>%
    group_by(WINDOW) %>% 
    summarise(min = min(POS), max = max(POS)) %>% 
    left_join(
        introgressed_windows %>% 
            filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>% 
            filter(WINDOW %in% introgressed_windows_filter$WINDOW) %>%
            dplyr::select(WINDOW) %>% 
            group_by(WINDOW) %>% 
            summarise(n = n())
    ) #%>% write_tsv("CHROM_11_Introgression_Mf_into_Mn.tsv")


# Add metadata
introgressed_windows_meta_Mf_to_Mn <- introgressed_windows_meta %>% # meatdata of Mn to Mf introgression
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0)




# Create table with overlapping windows for 028 and 053 - outlier samples that are highly similar
introgressed_windows %>% 
    filter(grepl("PK_SB_DNA_028", SAMPLE) | grepl("PK_SB_DNA_053", SAMPLE)) %>% 
    arrange(WINDOW) %>% 
    select(SAMPLE, WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    filter(n > 1) %>%
    left_join(
        PKA1H1_windows %>% 
            group_by(CHROM, WINDOW) %>% 
            summarise(start = min(POS), end = max(POS))
    ) %>%
    rename(n_1 = n) %>%
    left_join(
        introgressed_windows_meta %>% 
            filter(State == "Sarawak") %>% 
            group_by(WINDOW) %>% 
            summarise(n = n()) %>% 
            rename(n_2 = n)
    )












