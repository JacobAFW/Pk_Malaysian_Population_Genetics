
introgression_table_window  <- introgression_table %>% 
    mutate(TMP_WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% # creates sliding windows based on window_size - POS is specific to the location within the CHROMOSOME, and so, windows with the same POS on different CHROM will end up with the same value 
    mutate(WINDOW = group_indices(., CHROM, TMP_WINDOW)) %>%  # this uses the window generated above to give us a unique ID for every window across all chromosomes 
    group_by(SAMPLE, WINDOW) %>%
    summarise(Mf_distance = sum(SNP != Mf, na.rm = TRUE), # the number of times the sample allele doesn't match the dominant allele = genetic distance
        Mn_distance = sum(SNP != Mn, na.rm = TRUE),
        Pen_distance = sum(SNP != Peninsular, na.rm = TRUE))  %>%
    left_join(
        metadata %>% rename(SAMPLE = Sample)
    ) %>%
    filter(SAMPLE %in% metadata$Sample) # only include samples in metadata - these are samples that have pass all filters



# Raster type plots
raster_plot <- introgression_table_window %>%
    filter(Cluster == "Mn") %>% 
    ggplot(aes(x = Mf_distance, y = Mn_distance)) +
    #geom_point() +
    geom_density_2d_filled(contour_var = "density", alpha = 0.5) +
    theme(legend.position = "none")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/raster_tests.png", dpi = 300, raster_plot)


raster_plot <- introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = introgression_table_window %>%
        filter(SAMPLE == "ERR366425"),
        size = 0.75, alpha = 0.75) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window %>% 
            filter(Cluster == "Mn" | Cluster == "Mf" | Cluster == "Peninsular"),
        contour_var = "density") +
    geom_smooth(mapping = aes(x = Mf_distance, y = Mn_distance),
        data = introgression_table_window %>% 
            filter(Cluster == "Mn"),
        method = "lm",
        colour = "black") +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/raster_tests_2.png", dpi = 300,  raster_plot)

raster_plot <- introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = introgression_table_window %>%
        filter(SAMPLE == "ERR366425"),
        size = 0.75, alpha = 0.75) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/raster_tests_3.png", dpi = 300,  raster_plot)

raster_plot <- introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = introgression_table_window %>%
        filter(SAMPLE == "ERR366425"),
        size = 0.75, alpha = 0.75) +
    geom_density_2d(mapping = aes(x = Mf_distance, y = Mn_distance, colour = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    scale_colour_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/raster_tests_4.png", dpi = 300,  raster_plot)


# Plot initial Sarikei samples (loop below to plot all samples)
distance_plot <- genetic_distance_plot(introgression_table_window, "ERR274221")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274221_MfvsMn_dot.png", dpi = 300, distance_plot)

distance_plot <- genetic_distance_plot(introgression_table_window, "ERR274222")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274222_MfvsMn_dot.png", dpi = 300, distance_plot)

distance_plot <- genetic_distance_plot(introgression_table_window, "ERR274224")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274224_MfvsMn_dot.png", dpi = 300, distance_plot)

distance_plot <- genetic_distance_plot(introgression_table_window, "ERR274225")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274225_MfvsMn_dot.png", dpi = 300, distance_plot)

distance_plot <- genetic_distance_plot(introgression_table_window, "ERR366425")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR366425_MfvsMn_dot.png", dpi = 300, distance_plot)

distance_plot <- genetic_distance_plot(introgression_table_window, "ERR366426")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR366426_MfvsMn_dot.png", dpi = 300, distance_plot)




# Interogate regions of interest 

## Outliers in Mf plots that have a higher degree of windows with "Mn-distances" to Mf
Outlier_windows <- introgression_table_window %>%
    filter(SAMPLE == "ERR274221" |
        SAMPLE == "ERR274222" | 
        SAMPLE == "ERR985372" |
        SAMPLE == "ERR985373" |
        SAMPLE == "ERR985374" |
        SAMPLE == "ERR985375" |
        SAMPLE == "ERR985377" |
        SAMPLE == "ERR985378" |
        SAMPLE == "ERR985379" |
        SAMPLE == "ERR985380" |
        SAMPLE == "ERR985381" |
        SAMPLE == "ERR985384" |
        SAMPLE == "PK_SB_DNA_028_DKDL210002157-1a_HWHGKDSXY_L4") %>% 
    filter(Mf_distance > 12.5 & Mn_distance < 12.5) %>%
    group_by(WINDOW) %>%
    summarise(n = n()) %>% 
    filter(n >= 6) # filter for appearance in at least 50% of samples

### Are these interesting windows in other Mf samples?
introgression_table_window %>%
    filter(SAMPLE != "ERR274221" &
        SAMPLE != "ERR274222" &
        SAMPLE != "ERR985372" &
        SAMPLE != "ERR985373" &
        SAMPLE != "ERR985374" &
        SAMPLE != "ERR985375" &
        SAMPLE != "ERR985377" &
        SAMPLE != "ERR985378" &
        SAMPLE != "ERR985379" &
        SAMPLE != "ERR985380" &
        SAMPLE != "ERR985381" &
        SAMPLE != "ERR985384" &
        SAMPLE != "PK_SB_DNA_028_DKDL210002157-1a_HWHGKDSXY_L4") %>% 
    filter(Cluster == "Mf") %>%
    filter(Mf_distance > 12.5 & Mn_distance < 12.5) %>%
    group_by(WINDOW) %>%
    summarise(n = n()) %>% 
    filter(n >= 6) 

introgression_table_window %>% 
    filter(SAMPLE == "PK_SB_DNA_028_DKDL210002157-1a_HWHGKDSXY_L4") %>% 
    filter(Mf_distance > 12.5 & Mn_distance < 12.5) %>% 
    filter(WINDOW > 1260)

### create table with original CHROM and POS for exploring outliers
introgression_table_window_with_POS <- introgression_table %>%
    mutate(TMP_WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% # creates sliding windows based on window_size - POS is specific to the location within the CHROMOSOME, and so, windows with the same POS on different CHROM will end up with the same value 
    mutate(WINDOW = group_indices(., CHROM, TMP_WINDOW)) %>%
    select(CHROM, POS, SAMPLE, SNP, Cluster, WINDOW)

introgression_table_window_with_POS %>%
    filter(WINDOW %in%  Outlier_windows$WINDOW) %>% 
    select(CHROM, POS, WINDOW) %>%
    distinct() %>%
    group_by(CHROM) %>% 
    summarise(min = min(POS), max = max(POS))


###################### Exploring
introgression_table_window_with_POS %>%
    filter(WINDOW %in%  Outlier_windows$WINDOW) %>% 
    select(CHROM, POS, WINDOW) %>%
    distinct()  %>%
    filter(CHROM == "08") %>%
    filter(POS > 924000)

introgression_table_window_with_POS %>%
    filter(WINDOW %in%  Outlier_windows$WINDOW) %>% 
    select(CHROM, POS, WINDOW) %>%
    distinct()  %>%
    filter(CHROM == "08") %>%
    group_by(WINDOW) %>%
    summarise(min = min(POS), max = max(POS)) %>% 
    filter(min > 1290062)


## Look at outliers across all samples 

Outlier_windows <- introgression_table_window %>% 
    filter(Cluster == "Mf") %>%
    filter(Mf_distance > 15 & Mn_distance < 15) %>% # we need to divy this up slightly differently 
    group_by(WINDOW) %>%
    summarise(n = n()) %>% 
    arrange(desc(n)) 





## Outliers in Mn plots that have a higher degree of windows with "Mf-distances" to Mn
Outlier_windows <- introgression_table_window %>%
    filter(grepl("DNA_030", SAMPLE) |
        grepl("DNA_063", SAMPLE) |
        grepl("DNA_093", SAMPLE)) %>% 
    filter(Mn_distance > 24) %>%
    group_by(WINDOW) %>%
    summarise(n = n()) %>% 
    filter(n == 3)

introgression_table_window_with_POS %>%
    filter(WINDOW %in%  Outlier_windows$WINDOW) %>% 
    select(CHROM, POS, WINDOW) %>%
    distinct() 

# Use linear model and residuals to find outlier points
library(broom)
Mf_glm <- introgression_table_window %>%
    filter(Cluster == "Mf") %>%
    glm(Mf_distance ~ WINDOW, data = .) 

resid_plot <- augment(Mf_glm) %>%
    ggplot(aes(x = .fitted, y = .resid)) + 
    geom_point()

augment(Mf_glm) %>%
    filter(Mf_distance < 100 & .resid > (mean(.resid) + (3*sd(.resid)))) %>% # filter to 3*SD above the mean and below distance of 100 to its own cluster seems like an error, so, filter out (maybe because they're proportions) 
    group_by(Mf_distance, WINDOW) %>% 
    summarise(n = n()) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mf_resid.png", dpi = 300, resid_plot)

augment(Mf_glm) %>%
    filter(Mf_distance < 100 & .resid > (mean(.resid) + (3*sd(.resid)))) %>% # filter to 3*SD above the mean and below distance of 100 to its own cluster seems like an error, so, filter out (maybe because they're proportions) 
    group_by(Mf_distance, WINDOW) %>% 
    summarise(n = n()) %>% 
    filter(n > 5) 


Cluster_glm <- introgression_table_window %>%
    glm(Mf_distance ~ Cluster, data = .) 














# Identifying contours that points are present in

library(sp) # for point.in.polygon
library(MASS)

raster_plot <- introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = introgression_table_window %>%
        filter(SAMPLE == "PK_SB_DNA_053_DKDL210002182-1a_HWHGKDSXY_L4"),
        size = 0.75, alpha = 0.75) +
    geom_density_2d(mapping = aes(x = Mf_distance, y = Mn_distance, colour = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    scale_colour_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 

polygon_data <- ggplot_build(raster_plot) # create matrix for contours based on plot above - includes all datapoints

polygon_points <- polygon_data$data[[1]]
polygon_contours <- polygon_data$data[[2]] # extract contour data

# Do points fall within lower polygons of its own cluster AND within major polygons of another cluster 
## Do points fall within Major Mf
Mf_contours <- polygon_contours %>% 
    filter(colour == "#440154FF") %>% 
    #filter(level == 0.0045) # exclude bottom two polygon layers
    filter(group == "1-002-001")

Mf_contours <- point.in.polygon(pol.x = as.numeric(unlist(Mf_contours[3])), 
    pol.y = as.numeric(unlist(Mf_contours[4])),
    point.x = as.numeric(unlist(polygon_points[1])), 
    point.y = as.numeric(unlist(polygon_points[2]))) %>% 
    as.data.frame() %>% 
    rename("Mf" = ".")

## Do points fall OUTSIDE Major Mn - its own cluster
Mn_contours <- polygon_contours %>% 
    filter(colour == "#73D055FF") %>% 
    filter(group == "2-002-001")

Mn_contours <- point.in.polygon(pol.x = as.numeric(unlist(Mn_contours[3])), 
    pol.y = as.numeric(unlist(Mn_contours[4])),
    point.x = as.numeric(unlist(polygon_points[1])), 
    point.y = as.numeric(unlist(polygon_points[2]))) %>% 
    as.data.frame() %>% 
    rename("Mn" = ".")

## Do points faill within Major Peninsular
Pen_contours <- polygon_contours %>% 
    filter(colour == "#39568CFF") %>% #dplyr::select(level, nlevel, group) %>% unique() 
    filter(group == "3-002-001") # exclude bottom two polygon layers

Pen_contours <- point.in.polygon(pol.x = as.numeric(unlist(Pen_contours[3])), 
    pol.y = as.numeric(unlist(Pen_contours[4])),
    point.x = as.numeric(unlist(polygon_points[1])), 
    point.y = as.numeric(unlist(polygon_points[2]))) %>% 
    as.data.frame() %>% 
    rename("Pen" = ".")

points_in_contours_data <- Mf_contours %>% 
    cbind(Mn_contours) %>% 
    cbind(Pen_contours) 

# Identify and plot Mf to Mn exchange events
subset <- introgression_table_window %>%
    filter(SAMPLE == "PK_SB_DNA_053_DKDL210002182-1a_HWHGKDSXY_L4") %>% 
    cbind(points_in_contours_data) %>% 
    filter(Mf == 0 & Mn == 1 & Pen == 0)


raster_plot <- introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = subset, size = 0.75, alpha = 0.75) +
    geom_density_2d(mapping = aes(x = Mf_distance, y = Mn_distance, colour = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    scale_colour_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/PK_SB_DNA_053_DKDL210002182-1a_HWHGKDSXY_L4.png", dpi = 300,  raster_plot)





# Use polygon data to define outliers
polygon_plot <- introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) 


ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/polygon_plot.png", dpi = 300, polygon_plot)

polygon_data <- layer_data(polygon_plot)

polygon_data %>% 
    select(fill) %>%
    unique()

introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    stat_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window)



polygon_plot <- introgression_table_window  %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = introgression_table_window %>%
        filter(SAMPLE == "ERR366425"),
        size = 0.75, alpha = 0.75) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window,
        contour_var = "density") 

polygon_data <- ggplot_build(polygon_plot)

head(polygon_data$data[[1]]) # points
head(polygon_data$data[[2]]) # contours


library(sp) # for point.in.polygon
library(MASS)

dens <- with(introgression_table_window, MASS::kde2d(Mf_distance, Mn_distance, n = 50))

dens <- kde2d(polygon_data$data[[1]]$x, polygon_data$data[[1]]$y) 

Mf <- introgression_table_window %>%
    filter(Cluster == "Mf")
dens <- kde2d(Mf$Mf_distance, Mf$Mn_distance, n = 200)
ls <- contourLines(dens, level = levels)
inner <- point.in.polygon(Mf$Mf_distance, Mf$Mn_distance, ls[[2]]$x, ls[[2]y])





# This tells uf if the points in this sample fall within a polygon
point.in.polygon(pol.x = as.numeric(unlist(polygon_contours[3])), 
    pol.y = as.numeric(unlist(polygon_contours[4])),
    point.x = as.numeric(unlist(polygon_points[1])), 
    point.y = as.numeric(unlist(polygon_points[2])))







points_in_polygon_data <- introgression_table_window %>% # use full dataset to identify contours that points are present in
    split(.$SAMPLE) %>% 
    sapply(function(x) 
        point.in.polygon(as.numeric(unlist(polygon_contours[3])), as.numeric(unlist(polygon_contours[4])), # contours
        as.numeric(unlist(x$Mf_distance)), as.numeric(unlist(x$Mn_distance)))) %>% # points
    as.data.frame() 




points_in_polygon_data <- polygon_contours %>%
    cbind(points_in_polygon_data) %>% 
    select(1,2,13:ncol(.))
    write_tsv("points_in_polygon.tsv")





# Mf
Mf_contours <- polygon_contours %>% 
    filter(colour == "#440154FF")

# Mn
Mn_contours <- polygon_contours %>% 
    filter(colour == "#73D055FF")

# Peninsular
Pen_contours <- polygon_contours %>% 
    filter(colour == "#39568CFF")







# Get the genome position for the windows of interest
introgression_table %>%
    mutate(WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>%
    filter(WINDOW == 35000 | WINDOW ==  25000 | WINDOW == 415000) %>% 
    group_by(WINDOW) %>%
    summarise(min(POS), max(POS))


# distance across windows
distance_plot <- introgression_table_window %>% 
    filter(SAMPLE == "ERR274221") %>% 
    ggplot(aes(x = WINDOW, y = Mf_distance, colour = "#404788FF", shape = Mn_distance > 100)) +
        geom_point() +
        scale_colour_viridis_d() +
        ggtitle("ERR274221, Mf Cluster, Sarikei")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274221_Mf_dot.png", dpi = 300, distance_plot)

distance_plot <- introgression_table_window %>% 
    filter(SAMPLE == "ERR274221") %>% 
    ggplot(aes(x = WINDOW, y = Mn_distance)) +
        geom_point() +
        scale_colour_viridis_d() +
        ggtitle("ERR274221, Mf Cluster, Sarikei")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274221_Mn_dot.png", dpi = 300, distance_plot)



































################### TESTING - using original data to create sliding windows
## create sliding windows using original genotype data table (not long format)

sliding_windows <- genotype_table %>% 
     mutate(WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>%
     select(CHROM, POS, WINDOW) 

genotype_table %>% 
    select(CHROM, POS) %>%
    mutate()
    mutate(WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% 

    

sliding_windows %>% filter(CHROM == "08") %>% ungroup() %>% select(WINDOW) %>% unique()

introgression_table %>%













# Delete? this method was incorrect but might still be needed? Probably not...
introgression_table_window <- introgression_table %>%
    mutate(WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>%
    group_by(SAMPLE, WINDOW) %>%
    summarise(Mf_true = sum(SNP == Mf, na.rm=TRUE), Mf_false = sum(SNP != Mf, na.rm=TRUE),
        Mn_true = sum(SNP == Mn, na.rm=TRUE), Mn_false = sum(SNP != Mn, na.rm=TRUE),
        Pen_true = sum(SNP == Peninsular, na.rm=TRUE), Pen_false = sum(SNP != Peninsular, na.rm=TRUE)) %>%
    mutate(Mf_distance = ifelse(Mf_true > Mf_false, (Mf_true - Mf_false), (Mf_false - Mf_true)),
        Mn_distance = ifelse(Mn_true > Mn_false, (Mn_true - Mn_false), (Mn_false - Mn_true)),
        Pen_distance = ifelse(Pen_true > Pen_false, (Pen_true - Pen_false), (Pen_false - Pen_true))) %>%
    select(-c(3:8)) %>%
    left_join(
        metadata %>% rename(SAMPLE = Sample)
    ) 











# Original Method - Probably not 
########################################################################### TESTING - sample vs sample - nucleotide diveristy
window_size <- 50000
genotype_density <- read_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD/hmmIBD.tsv") %>%
    mutate(WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% # create window column
    relocate(CHROM, POS, WINDOW) %>% 
    pivot_longer(4:ncol(.), names_to = "SAMPLE", values_to = "SNP") %>%
    filter(SNP >= 1) %>% # filter for SNPs only - so that we can count the number of SNPs in each window (and sample)
    group_by(CHROM, WINDOW, SAMPLE) %>%
    summarise(Density = n()) # count the number of rows/SNPs per grouping

# Calculate pairwise diveristy - difference in SNP density between samples divided by the window size
pairwise_diversity <- genotype_density %>%
    inner_join(
        genotype_density, by = c("CHROM", "WINDOW")
    ) %>%
    filter(Density.x != Density.y) %>%
    mutate(Diversity = ifelse(Density.x > Density.y, (Density.x - Density.y)/window_size , (Density.y - Density.x)/window_size)) %>% # calculates paired SNP diveristy perwindow
    select(-c(Density.x, Density.y)) 

# Add metadata
pairwise_diversity_meta <- pairwise_diversity %>%
    left_join(
        metadata %>%
        rename(CLUSTER.x = Cluster) %>%
        rename(SAMPLE.x = Sample) %>%
        rename(LOCATION.x = Location)
    ) %>%
    left_join(
        metadata %>%
        rename(CLUSTER.y = Cluster) %>%
        rename(SAMPLE.y = Sample) %>%
        rename(LOCATION.y = Location)
    )

# Dot plot for cluster 
## outliers for a given window should be visible & we will be able to see the average
diversity_plot <- pairwise_diversity_meta %>%
    filter(CLUSTER.x == "Mn" & CLUSTER.y == "Mn") %>% # diversity for Mn to Mn pairwise comparison
    group_by(CHROM, WINDOW) %>%
    add_column(x_axis = 1:nrow(.)) %>%
    ungroup() %>%
    ggplot(aes(x = x_axis, y = Diversity, colour = CHROM)) +
        geom_point() +
        scale_colour_viridis_d() +
        theme(legend.position = "bottom", axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Diversity") + 
        guides(colour = guide_legend(nrow = 1))

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mn_dot.png", dpi = 300, width = 20, diversity_plot)

diversity_plot <- pairwise_diversity_meta %>%
    filter(CLUSTER.x == "Mf" & CLUSTER.y == "Mf") %>%
    group_by(CHROM, WINDOW) %>%
    add_column(x_axis = 1:nrow(.)) %>%
    ungroup() %>%
    ggplot(aes(x = x_axis, y = Diversity, colour = CHROM)) +
        geom_point() +
        scale_colour_viridis_d() +
        theme(legend.position = "bottom", axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Diversity") + 
        guides(colour = guide_legend(nrow = 1))

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mf_dot.png", dpi = 300, width = 20, diversity_plot)



# Distribution 
diversity_plot <- pairwise_diversity_meta %>%
    filter(CLUSTER.x == "Mn" & CLUSTER.y == "Mn") %>%
    ggplot(aes(x = Diversity)) +
        geom_histogram() 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mn_hist.png", dpi = 300, width = 20, diversity_plot)

diversity_plot <- pairwise_diversity_meta %>%
    filter(CLUSTER.x == "Mf" & CLUSTER.y == "Mf") %>%
    ggplot(aes(x = Diversity)) +
        geom_histogram() 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mf_hist.png", dpi = 300, width = 20, diversity_plot)

# Stats
pairwise_diversity_meta %>%
    filter(CLUSTER.x == "Mf" & CLUSTER.y == "Mf") %>%
    group_by(CLUSTER.x) %>%
    summarise(mean = mean(Diversity), median = median(Diversity), sd = sd(Diversity)) %>%
    rbind(
        pairwise_diversity_meta %>%
            filter(CLUSTER.x == "Mn" & CLUSTER.y == "Mn") %>%
            group_by(CLUSTER.x) %>%
            summarise(mean = mean(Diversity), median = median(Diversity), sd = sd(Diversity)) 
    ) %>%
    mutate(sd_3 = 3*sd)

# Filter based on 3*SD and see what comes out
diversity_plot <- pairwise_diversity_meta %>%
    filter(CLUSTER.x == "Mf" & CLUSTER.y == "Mf") %>%
    filter(Diversity > 0.00868) %>% # mean + 3*sd
    group_by(CHROM, WINDOW) %>%
    add_column(x_axis = 1:nrow(.)) %>%
    ungroup() %>%
    ggplot(aes(x = x_axis, y = Diversity, colour = CHROM)) +
        geom_point() +
        scale_colour_viridis_d() +
        theme(legend.position = "bottom", axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Diversity") + 
        guides(colour = guide_legend(nrow = 1))

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mf_dot_filtered.png", dpi = 300, width = 20, diversity_plot)

diversity_plot <- pairwise_diversity_meta %>%
    filter(CLUSTER.x == "Mf" & CLUSTER.y == "Mf") %>%
    filter(Diversity > 0.00868) %>% # mean + 3*sd
    filter(CHROM == "08") %>%
    group_by(CHROM, WINDOW) %>%
    add_column(x_axis = 1:nrow(.)) %>%
    ungroup() %>%
    ggplot(aes(x = x_axis, y = Diversity, colour = CHROM)) +
        geom_point() +
        scale_colour_viridis_d() +
        theme(legend.position = "bottom", axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Diversity") + 
        guides(colour = guide_legend(nrow = 1))

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mf_dot_filtered_chr8.png", dpi = 300, width = 20, diversity_plot)


diversity_plot <- pairwise_diversity_meta %>%
    filter(CLUSTER.x == "Mn" & CLUSTER.y == "Mn") %>%
    filter(Diversity > 0.00247) %>% # mean + 3*sd
    group_by(CHROM, WINDOW) %>%
    add_column(x_axis = 1:nrow(.)) %>%
    ungroup() %>%
    ggplot(aes(x = x_axis, y = Diversity, colour = CHROM)) +
        geom_point() +
        scale_colour_viridis_d() +
        theme(legend.position = "bottom", axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Diversity") + 
        guides(colour = guide_legend(nrow = 1))

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Mn_dot_filtered.png", dpi = 300, width = 20, diversity_plot)










########################################## Delete below?
introgression_table <- genotype_table %>% 
    filter(CHROM == "08") %>% # remove when finished
    left_join(dominant_allele) %>% 
    mutate_at(3:187, ~ifelse(. == Mn, "Mn", ifelse(. == "Mf", "Mf", "Peninsular"))) %>% # specific to number of samples - doesn't work, because if multiple clusters have the same value then the sample will just equal the first cluster
    pivot_longer(3:187, names_to = "SAMPLE", values_to = "Allele")
########################################## 








###############################################################################


# Subset for comparisons
## Mf
## ERR274221 - other cluster
diversity_plot <- pairwise_diversity_meta %>%
    filter(SAMPLE.x == "ERR274222" & CLUSTER.y == "Mn") %>%
    group_by(CHROM, WINDOW) %>%
    summarise(Diversity = mean(Diversity, na.rm = T)) %>%
    filter(CHROM == "01") %>%
    ggplot(aes(x = WINDOW, y = Diversity)) +
        geom_col() +
        #scale_fill_manual(values = c("#440154FF", "#3CBB75FF")) +
        theme(legend.position = "none", axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Diversity") 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274221_Mn.png", dpi = 600, width = 12, diversity_plot)

## ERR274221 - same cluster
diversity_plot <- pairwise_diversity_meta %>%
    filter(SAMPLE.x == "ERR274222" & CLUSTER.y == "Mf" & SAMPLE.y != "ERR274222") %>%
    group_by(CHROM, WINDOW) %>%
    summarise(Diversity = mean(Diversity, na.rm = T)) %>%
    filter(CHROM == "01") %>%
    ggplot(aes(x = WINDOW, y = Diversity)) +
        geom_col() +
        #scale_fill_manual(values = c("#440154FF", "#3CBB75FF")) +
        theme(legend.position = "none", axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Diversity") 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274221_Mf.png", dpi = 600, width = 12, diversity_plot)

## PK_SB_DNA_028 - other cluster

diversity_plot <- pairwise_diversity_meta %>%
    filter(grepl("PK_SB_DNA_028", SAMPLE.x) | CLUSTER.y == "Mn") %>%
    group_by(CHROM, WINDOW) %>%
    summarise(Diversity = mean(Diversity, na.rm = T)) %>%
    filter(CHROM == "08") %>%
    ggplot(aes(x = WINDOW, y = Diversity, fill = ifelse(Diversity > 0.002, "Above", "Below"))) +
        geom_col() +
        scale_fill_manual(values = c("#440154FF", "#3CBB75FF")) +
        theme(legend.position = "none", axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Diversity") 





















########################################################################### TESTING - works, but may not be needed 

# Read in genotype file created for IBD, convert to sliding window summary and join metadata
window_size <- 50000
genotype_density <- read_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD/hmmIBD.tsv") %>%
    mutate(WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% # create window column
    relocate(CHROM, POS, WINDOW) %>% 
    pivot_longer(4:ncol(.), names_to = "SAMPLE", values_to = "SNP") %>%
    filter(SNP >= 1) %>% # filter for SNPs only - so that we can count the number of SNPs in each window (and sample)
    group_by(CHROM, WINDOW, SAMPLE) %>%
    summarise(Density = n()) %>% # count the number of rows/SNPs per grouping
    left_join(
        metadata %>%
        rename(SAMPLE = Sample)
        )






########################################################################### TESTING - sample vs cluster - density 
 

## Mf
## ERR274221
geno_plot <- genotype_density %>%
    filter(Cluster == "Mn" | SAMPLE == "ERR274222") %>%
    group_by(Cluster, CHROM, WINDOW) %>%
    summarise(Density = mean(Density, na.rm = T)) %>%
    filter(CHROM == "10") %>%
    ggplot(aes(x = WINDOW, y = Density, fill = Cluster, group = Cluster)) +
       # geom_point() +
       # geom_line() + 
        geom_col(position = "dodge") +
        theme(axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Density") +
        scale_fill_manual(values = c("#440154FF", "#3CBB75FF")) 
        #facet_grid(~CHROM, scales = "free", space = "free")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274221_Mn.png", dpi = 600, width = 12, geno_plot)

geno_plot <- genotype_density %>%
    filter(Cluster == "Mf") %>%
    mutate(Cluster = ifelse(SAMPLE == "ERR274222", "Mf_Mixed", .$Cluster)) %>%
    group_by(Cluster, CHROM, WINDOW) %>%
    summarise(Density = mean(Density, na.rm = T)) %>%
    filter(CHROM == "10") %>%
    ggplot(aes(x = WINDOW, y = Density, fill = Cluster, group = Cluster)) +
       # geom_point() +
       # geom_line() + 
        geom_col(position = "dodge") +
        theme(axis.text.x = element_blank()) +
        xlab("Sliding Window") +
        ylab("SNP Density") +
        scale_fill_manual(values = c( "#3CBB75FF", "#440154FF")) 
        #facet_grid(~CHROM, scales = "free", space = "free")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/ERR274221_Mf.png", dpi = 600, width = 12, geno_plot)

########################################################################### TESTING - sample vs cluster for nucleotide diversity (in same object)

## Differences between sample and cluster

## MF

geno_plot <- genotype_density %>%
    filter(Cluster == "Mn" | SAMPLE == "PK_SB_DNA_028_DKDL210002157-1a_HWHGKDSXY_L4") %>%
    filter(SAMPLE != "ERR274221" & SAMPLE != "ERR274222" & SAMPLE != "ERR985375" & SAMPLE != "ERR985378" & SAMPLE != "ERR985381" & SAMPLE != "ERR985384" & SAMPLE != "PK_SB_DNA_042_DKDL210002171-1a_HWHGKDSXY_L4" &  SAMPLE != "PK_SB_DNA_063_DKDL210002192-1a_HWHGKDSXY_L4" & SAMPLE != "PK_SB_DNA_093_DKDL210002222-1a_HWHGKDSXY_L4") %>% # filter out other mixed samples
    group_by(Cluster, CHROM, WINDOW) %>%
    summarise(Density = mean(Density, na.rm = T)) %>%
    pivot_wider(names_from = Cluster, values_from = Density) %>%
    mutate(Other_cluster = ifelse(Mf > Mn, (Mf - Mn)/window_size , (Mn - Mf)/window_size)) %>%
    select(-c(3:4)) %>%
    left_join(
        genotype_density %>%
            filter(Cluster == "Mf") %>%
            filter(SAMPLE != "ERR274221" & SAMPLE != "ERR274222" & SAMPLE != "ERR985375" & SAMPLE != "ERR985378" & SAMPLE != "ERR985381" & SAMPLE != "ERR985384" & SAMPLE != "PK_SB_DNA_042_DKDL210002171-1a_HWHGKDSXY_L4" &  SAMPLE != "PK_SB_DNA_063_DKDL210002192-1a_HWHGKDSXY_L4" & SAMPLE != "PK_SB_DNA_093_DKDL210002222-1a_HWHGKDSXY_L4") %>% # filter out other mixed samples
            mutate(Cluster = ifelse(SAMPLE == "PK_SB_DNA_028_DKDL210002157-1a_HWHGKDSXY_L4", "Mf_Mixed", .$Cluster)) %>%
            group_by(Cluster, CHROM, WINDOW) %>%
            summarise(Density = mean(Density, na.rm = T)) %>%
            pivot_wider(names_from = Cluster, values_from = Density) %>%
            mutate(Same_cluster = ifelse(Mf > Mf_Mixed, (Mf - Mf_Mixed)/window_size, (Mf_Mixed - Mf)/window_size)) %>% # window size to diversity within window
            select(-c(3:4)) 
    ) %>%
    pivot_longer(3:4, names_to = "Comparison", values_to = "Density") %>%
    filter(CHROM == "10") %>%
    ggplot(aes(x = WINDOW, y = Density, fill = Comparison, group = Comparison)) +
       # geom_point() +
       # geom_line() + 
        geom_col(position = "dodge") +
        theme(axis.text.x = element_blank()) +
        xlab("Sliding Window Position") +
        ylab("Difference in SNP diversity") +
        scale_fill_manual(values = c("#440154FF", "#3CBB75FF")) 
        #facet_grid(~CHROM, scales = "free", space = "free")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/PK_SB_DNA_028_DKDL210002157-1a_HWHGKDSXY_L4_difference.png", dpi = 600, width = 12, geno_plot)

## Mn

geno_plot <- genotype_density %>%
    filter(Cluster == "Mf" | SAMPLE == "PK_SB_DNA_042_DKDL210002171-1a_HWHGKDSXY_L4") %>%
    filter(SAMPLE != "ERR274221" & SAMPLE != "ERR274222" & SAMPLE != "ERR985375" & SAMPLE != "ERR985378" & SAMPLE != "ERR985381" & SAMPLE != "ERR985384" & SAMPLE !=  "PK_SB_DNA_028_DKDL210002157-1a_HWHGKDSXY_L4" & SAMPLE != "PK_SB_DNA_063_DKDL210002192-1a_HWHGKDSXY_L4" & SAMPLE != "PK_SB_DNA_093_DKDL210002222-1a_HWHGKDSXY_L4") %>% # filter out other mixed samples
    group_by(Cluster, CHROM, WINDOW) %>%
    summarise(Density = mean(Density, na.rm = T)) %>%
    pivot_wider(names_from = Cluster, values_from = Density) %>%
    mutate(Other_cluster = ifelse(Mf > Mn, (Mf - Mn)/window_size , (Mn - Mf)/window_size)) %>%
    select(-c(3:4)) %>%
    left_join(
        genotype_density %>%
            filter(Cluster == "Mn") %>%
            filter(SAMPLE != "ERR274221" & SAMPLE != "ERR274222" & SAMPLE != "ERR985375" & SAMPLE != "ERR985378" & SAMPLE != "ERR985381" & SAMPLE != "ERR985384" & SAMPLE != "PK_SB_DNA_028_DKDL210002157-1a_HWHGKDSXY_L4" &  SAMPLE != "PK_SB_DNA_063_DKDL210002192-1a_HWHGKDSXY_L4" & SAMPLE != "PK_SB_DNA_093_DKDL210002222-1a_HWHGKDSXY_L4") %>% # filter out other mixed samples
            mutate(Cluster = ifelse(SAMPLE == "PK_SB_DNA_042_DKDL210002171-1a_HWHGKDSXY_L4", "Mn_Mixed", .$Cluster)) %>%
            group_by(Cluster, CHROM, WINDOW) %>%
            summarise(Density = mean(Density, na.rm = T)) %>%
            pivot_wider(names_from = Cluster, values_from = Density) %>%
            mutate(Same_cluster = ifelse(Mn > Mn_Mixed, (Mn - Mn_Mixed)/window_size, (Mn_Mixed - Mn)/window_size)) %>% # window size to diversity within window
            select(-c(3:4)) 
    ) %>%
    pivot_longer(3:4, names_to = "Comparison", values_to = "Density") %>%
    filter(CHROM == "08") %>%
    ggplot(aes(x = WINDOW, y = Density, fill = Comparison, group = Comparison)) +
       # geom_point() +
       # geom_line() + 
        geom_col(position = "dodge") +
        theme(axis.text.x = element_blank()) +
        xlab("Sliding Window Position") +
        ylab("Difference in SNP diversity") +
        scale_fill_manual(values = c("#440154FF", "#3CBB75FF")) 
        #facet_grid(~CHROM, scales = "free", space = "free")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/PK_SB_DNA_042_DKDL210002171-1a_HWHGKDSXY_L4_difference.png", dpi = 600, width = 12, geno_plot)

################################################################################################





## ERR274222 
## ERR985375
## ERR985378
## ERR985381
## ERR985384
## PK_SB_DNA_028

## Mn
## PK_SB_DNA_042
## PK_SB_DNA_063
## PK_SB_DNA_093





# Ploting windows with contours designed just for that contour 

distance_plot <- introgression_table_window %>% 
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window%>%
        filter(WINDOW == 1269 | WINDOW == 1268),
        contour_var = "density") +
    geom_point(data = introgression_table_window %>%
        filter(WINDOW == 1269 | WINDOW == 1268),
        size = 1, alpha = 0.75, aes(shape = Cluster)) +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") 
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/window_plot_1268_1269_Mf.png", dpi = 300, distance_plot)





# OLD window method

# WINDOWS with introgression data
## Using proportions - removes bias introduced through number of SNPs in a window
introgression_table_window  <- introgression_table %>% 
    mutate(TMP_WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% # creates sliding windows based on window_size - POS is specific to the location within the CHROMOSOME, and so, windows with the same POS on different CHROM will end up with the same value 
    mutate(WINDOW = group_indices(., CHROM, TMP_WINDOW)) %>%  # this uses the window generated above to give us a unique ID for every window across all chromosomes 
    filter(SNP >= 0) %>% # filter out missing calls - messes with the distance calculation
    group_by(SAMPLE, WINDOW) %>%
    summarise(Mf_distance = sum(SNP != Mf, na.rm = TRUE)/n()*100, # the number of times the sample allele doesn't match the dominant allele = genetic distance
        Mn_distance = sum(SNP != Mn, na.rm = TRUE)/n()*100,
        Pen_distance = sum(SNP != Peninsular, na.rm = TRUE)/n()*100)  %>%
    left_join(
        metadata %>% rename(SAMPLE = Sample)
    ) %>%
    filter(SAMPLE %in% metadata$Sample) # only include samples in metadata - these are samples that have pass all filters



##################### REMOVE? PKA1H1_windows should work for this 
# combine with POS and CHROM info and export
introgression_table_window_with_POS <- introgression_table %>%
    mutate(TMP_WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% # creates sliding windows based on window_size - POS is specific to the location within the CHROMOSOME, and so, windows with the same POS on different CHROM will end up with the same value 
    mutate(WINDOW = group_indices(., CHROM, TMP_WINDOW)) %>%
    dplyr::select(CHROM, POS, WINDOW) 

introgression_table_window_with_POS <- introgression_table_window_with_POS %>% unique()

introgressed_windows_with_pos <- introgressed_windows %>%
    dplyr::select(WINDOW) %>% 
    unique() %>% 
    dplyr::left_join(introgression_table_window_with_POS)

#write_tsv(introgressed_windows_with_pos, "/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/introgressed_windows_with_pos.tsv")
#introgressed_windows_with_pos <- read_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/introgressed_windows_with_pos.tsv")
##################### REMOVE?




# Summary info for TOTAL introgression across clusters
introgressed_windows_meta %>% 
    dplyr::select(SAMPLE, Cluster) %>%
    group_by(SAMPLE, Cluster) %>% 
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(Cluster) %>% 
    summarise(mean = mean(n), median = median(n), sd = sd(n), n = n()) %>% # number of windows
    arrange(desc(median)) #%>% write_tsv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/average_windows_for_clusters.tsv")