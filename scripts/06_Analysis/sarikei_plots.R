raster_plot <- introgression_table_window %>%
    filter(SAMPLE == "ERR274225") %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(size = 0.75, alpha = 0.75) +
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
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Sarikei_ERR274225.png", dpi = 300,  raster_plot)

raster_plot <- introgression_table_window %>%
    filter(SAMPLE == "ERR366425") %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(size = 0.75, alpha = 0.75) +
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
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Sarikei_ERR366425.png", dpi = 300,  raster_plot)

raster_plot <- introgression_table_window %>%
    filter(SAMPLE == "ERR366426") %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(size = 0.75, alpha = 0.75) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window %>% 
            filter(Cluster == "Mn" | Cluster == "Mf" | Cluster == "Peninsular"),
        contour_var = "density") +
    geom_smooth(mapping = aes(x = Mf_distance, y = Mn_distance),
        data = introgression_table_window %>% 
            filter(Cluster == "Mf"),
        method = "lm",
        colour = "black") +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Sarikei_ERR366426.png", dpi = 300,  raster_plot)

raster_plot <- introgression_table_window %>%
    filter(SAMPLE == "ERR274221") %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(size = 0.75, alpha = 0.75) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window %>% 
            filter(Cluster == "Mn" | Cluster == "Mf" | Cluster == "Peninsular"),
        contour_var = "density") +
    geom_smooth(mapping = aes(x = Mf_distance, y = Mn_distance),
        data = introgression_table_window %>% 
            filter(Cluster == "Mf"),
        method = "lm",
        colour = "black") +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Sarikei_ERR274221.png", dpi = 300,  raster_plot)

raster_plot <- introgression_table_window %>%
    filter(SAMPLE == "ERR274222") %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance)) +
    geom_point(size = 0.75, alpha = 0.75) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window %>% 
            filter(Cluster == "Mn" | Cluster == "Mf" | Cluster == "Peninsular"),
        contour_var = "density") +
    geom_smooth(mapping = aes(x = Mf_distance, y = Mn_distance),
        data = introgression_table_window %>% 
            filter(Cluster == "Mf"),
        method = "lm",
        colour = "black") +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Sarikei_ERR274222.png", dpi = 300,  raster_plot)

raster_plot <- introgression_table_window %>%
    filter(SAMPLE == "ERR274224") %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(size = 0.75, alpha = 0.75) +
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
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200)) 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Introgression/Sarikei_ERR274224.png", dpi = 300,  raster_plot)