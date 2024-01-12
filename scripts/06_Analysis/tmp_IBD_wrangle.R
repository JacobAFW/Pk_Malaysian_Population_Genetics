Pk.IB.txt  <- read_table("Pk.hmm.txt")

# samples of interest in sample1
subset_IBD <- Pk.IB.txt %>% 
    filter(grepl("PK_SB_DNA_028", sample1) | grepl("PK_SB_DNA_053", sample1)) %>% 
    left_join(
      metadata %>%
        rename(sample2 = Sample),
      by = "sample2"
    ) %>%
    filter(Cluster == "Mn") %>% 
    add_column(data = "sample1") %>%
    select(sample1, sample2, chr, start, end, different, data) %>% 
    rbind(
        # samples of interest in sample2
        Pk.IB.txt %>% 
            filter(grepl("PK_SB_DNA_028", sample2) | grepl("PK_SB_DNA_053", sample2)) %>% 
            left_join(
            metadata %>%
                rename(sample1 = Sample),
            by = "sample1"
            ) %>%
            filter(Cluster == "Mn") %>% 
            add_column(data = "sample2") %>%
            select(sample1, sample2, chr, start, end, different, data)
    ) %>%
    rbind(
        # Mn-to-Mn
        Pk.IB.txt %>% 
            left_join(
            metadata %>%
                select(Sample, Cluster) %>% 
                rename(sample1 = Sample) %>%
                rename(cluster1 = Cluster),
            by = "sample1"
            ) %>% 
            left_join(
            metadata %>%
                select(Sample, Cluster) %>% 
                rename(sample2 = Sample) %>%
                rename(cluster2 = Cluster),
            by = "sample2" 
            ) %>%
            filter(cluster1 == "Mn" & cluster2 == "Mn") %>% 
            add_column(data = "Mn-Mn") %>%
            select(sample1, sample2, chr, start, end, different, data)
    )

windows_data <- subset_IBD %>% 
    group_by(data, chr) %>% 
    summarise(start = mean(start), end = mean(end)) %>%
    arrange(chr, start) 

windows_plot <- windows_data %>% 
    pivot_longer(cols = start:end, names_to = "pos", values_to = "pos_value") %>%
    ggplot(aes(x = pos_value, y = data, group = data)) +
    geom_point() +
    geom_line() +
    theme(legend.position = "none") +
    facet_grid(~chr, scales = "free", space = "free") +
    xlab("Position in chromosome") 

ggsave("regions_of_shared_IBD.png", dpi = 300, width = 20, windows_plot)










# Plot Fraction of sites that are IBD
windows_data <- subset_IBD %>%
    left_join(
      metadata %>%
        rename(sample2 = Sample),
      by = "sample2"
    ) %>%
    filter(Cluster == "Mn") %>%
    filter(different == 1) %>%
    arrange(chr, start) %>% 
    mutate(sample1 = str_remove(sample1, "_DK.*")) %>% 
    mutate(sample1 = str_remove(sample1, "PK_SB_DNA_")) %>% 
    group_by(sample1, Cluster, chr) %>% 
    summarise(start = mean(start), end = mean(end)) %>%
    arrange(chr, start) 

windows_plot <- windows_data %>% 
    pivot_longer(cols = start:end, names_to = "pos", values_to = "pos_value") %>%
    ggplot(aes(x = pos_value, y = sample1, group = sample1)) +
    geom_point() +
    geom_line() +
    theme(legend.position = "none") +
    facet_grid(~chr, scales = "free", space = "free") +
    xlab("Position in chromosome") 

ggsave("regions_of_shared_IBD.png", dpi = 300, windows_plot)