library(igraph)

.reorder.metadata <-
  function(IBD, metadata, sample.col = "Sample") {
    samples <- unique(c(IBD[, "sample1"], IBD[, "sample2"]))
    sample.order <- match(metadata[, sample.col], samples)
    actual.length <- length(na.omit(sample.order))
    
    metadata <- metadata[order(sample.order), ]
    metadata[1:actual.length, ]
  }


.get.labels <-
  function(metadata,
           label.col,
           label.palette = NULL) {
    labels <- metadata[, label.col]
    labels[is.na(labels)] <- NaN
    
    labels <- factor(labels)
    
    if (is.null(label.palette)) {
      labels
    } else {
      factor(labels, levels = names(label.palette))
    }
  }


.generate.label.palette <-
  function(labels) {
    label.names <- levels(labels)
    label.palette <- # c("#440154FF", "#73D055FF", "#39568CFF") # for clusters
      viridis(length(label.names), option = ifelse(length(label.names) < 5, "D", "H")) # for districts
    names(label.palette) <- label.names
    
    label.palette
  }


.generate.label.colours <-
  function(labels, label.palette)
    label.palette[labels]


.create.edgelist <-
  function(IBD)
    IBD[, c("sample1", "sample2", "fract_sites_IBD")]


.create.vertices <-
  function(metadata, sample.col = "Sample") {
    vertices <- data.frame(metadata[, sample.col])
    names(vertices) <- sample.col
    
    vertices
  }


.internal.plot.IBD <- function(IBD.graph,
                               unlabelled = TRUE,
                               coords = layout_nicely,
                               percent.cutoff = NA) {
  if (unlabelled) {
    plot(
      IBD.graph,
      layout = coords,
      vertex.size = 4,
      vertex.label = NA,
      main = paste0("IBD >=", percent.cutoff, "%")
    )
  } else {
    plot(
      IBD.graph,
      layout = coords,
      vertex.size = 4,
      vertex.label.cex = 0.3,
      main = paste0("IBD >=", percent.cutoff, "%")
    )
  }
}


.plot.IBD <-
  function(file,
           IBD.graph,
           unlabelled = TRUE,
           label.palette = NULL,
           labels = NULL,
           legend.title = NA,
           coords = layout_nicely,
           percent.cutoff = NA) {
    pdf(file)
    par(mar = rep.int(1, 4) + 0.1)
    
    .internal.plot.IBD(
      IBD.graph,
      unlabelled = unlabelled,
      coords = layout_nicely,
      percent.cutoff = percent.cutoff
    )
    
    if (!is.null(label.palette))
      legend(
        "bottomleft",
        legend = levels(labels),
        fill = label.palette,
        title = legend.title,
        cex = 0.6
      )
    
    dev.off()
  }


.plot.IBDs <- function(prefix,
                       IBD.graph,
                       label.palette,
                       labels,
                       legend.title,
                       coords,
                       percent.cutoff) {
  labelled.file <-
    paste0(prefix, "_labelled_IBD", percent.cutoff, ".pdf")
  .plot.IBD(
    labelled.file,
    IBD.graph,
    unlabelled = FALSE,
    label.palette = label.palette,
    labels = labels,
    legend.title = legend.title,
    coords = coords,
    percent.cutoff = percent.cutoff
  )
  
  unlabelled.file <-
    paste0(prefix, "_unlabelled_IBD", percent.cutoff, ".pdf")
  .plot.IBD(
    unlabelled.file,
    IBD.graph,
    unlabelled = TRUE,
    label.palette = label.palette,
    labels = labels,
    legend.title = legend.title,
    coords = coords,
    percent.cutoff = percent.cutoff
  )
}


plot.IBD <- function(file,
                     IBD.cutoffs,
                     IBD,
                     metadata,
                     label.col,
                     unlabelled = TRUE,
                     label.palette = NULL,
                     legend.title = NA,
                     sample.col = "Sample") {
  edgelist <- .create.edgelist(IBD)
  
  metadata <- .reorder.metadata(IBD, metadata)
  
  vertices <- .create.vertices(metadata)
  
  sqrt.n.plots <- ceiling(sqrt(length(IBD.cutoffs)))
  size <- 7 * sqrt.n.plots
  
  legend.position <- sqrt.n.plots ^ 2 - sqrt.n.plots + 1
  legend.plot <- FALSE
  
  pdf(file, width = size, height = size)
  par(mar = rep.int(1, 4) + 0.1,
      mfrow = c(sqrt.n.plots, sqrt.n.plots))
  
  if (is.null(label.palette)) {
    labels <- .get.labels(metadata, label.col)
    label.palette <- .generate.label.palette(labels)
    
  } else {
    labels <-
      .get.labels(metadata, label.col, label.palette = label.palette)
  }
  
  label.colours <-
    .generate.label.colours(labels, label.palette)
  
  for (i in seq_along(IBD.cutoffs)) {
    d <- edgelist[edgelist[, "fract_sites_IBD"] >= IBD.cutoffs[i], ]
    
    IBD.graph <-
      graph_from_data_frame(d, directed = FALSE, vertices = vertices)
    
    coords <- layout_(IBD.graph, nicely())
    percent.cutoff <- round(IBD.cutoffs[i] * 100, digits = 2)
    
    IBD.graph <-
      set_vertex_attr(IBD.graph, "color", value = label.colours)
    
    .internal.plot.IBD(IBD.graph,
                       unlabelled,
                       coords,
                       percent.cutoff)
    
    if (i == legend.position) {
      legend(
        "bottomleft",
        legend = levels(labels),
        fill = label.palette,
        title = legend.title,
        cex = 0.6
      )
      
      legend.plot <- TRUE
    }
  }
  
  if (!legend.plot) {
    while (i != legend.position) {
      plot.new()
      i <- i + 1
    }
    
    legend(
      "bottomleft",
      legend = levels(labels),
      fill = label.palette,
      title = legend.title,
      cex = 0.6
    )
  }
  
  dev.off()
}


plot.IBDs <- function(prefixes,
                      IBD.cutoffs,
                      IBD,
                      metadata,
                      label.cols,
                      label.palettes = NULL,
                      legend.titles = NA,
                      sample.col = "Sample") {
  edgelist <- .create.edgelist(IBD)
  
  metadata <- .reorder.metadata(IBD, metadata)
  
  vertices <- .create.vertices(metadata)
  
  for (IBD.cutoff in IBD.cutoffs) {
    d <- edgelist[edgelist[, "fract_sites_IBD"] >= IBD.cutoff, ]
    
    IBD.graph <-
      graph_from_data_frame(d, directed = FALSE, vertices = vertices)
    
    coords <- layout_(IBD.graph, nicely())
    percent.cutoff <- round(IBD.cutoff * 100, digits = 2)
    
    for (label.col in label.cols) {
      label.palette <- label.palettes[[label.col]]
      
      if (is.null(label.palette)) {
        labels <- .get.labels(metadata, label.col)
        label.palette <- .generate.label.palette(labels)
        
      } else {
        labels <-
          .get.labels(metadata, label.col, label.palette = label.palette)
      }
      
      label.colours <-
        .generate.label.colours(labels, label.palette)
      
      IBD.graph <-
        set_vertex_attr(IBD.graph, "color", value = label.colours)
      
      legend.title <- legend.titles[[label.col]]
      
      prefix <- prefixes[[label.col]]
      
      .plot.IBDs(prefix,
                 IBD.graph,
                 label.palette,
                 labels,
                 legend.title,
                 coords,
                 percent.cutoff)
    }
  }
}


read.pairwise.matrix <- function(file) {
  as.matrix(read.delim(file, check.names = FALSE))
}


IBD.cutoffs <- c(0.125, 0.25, 0.5)

# include +-5% IBD
relaxed.IBD.cutoffs <- c(IBD.cutoffs, 1)
relaxed.IBD.cutoffs <- relaxed.IBD.cutoffs * 0.95


## Pk
library(igraph)
library(tidyverse)
library(MetamapsDB)
library(viridis)

Pk.IBD.file <- "Pk_subset.hmm_fract.txt"
Pk.IBD <- read.delim(Pk.IBD.file)

metadata <- read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.csv", col_names = FALSE) %>% 
        select(1, 6, 7) %>%
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
    select(Sample, district)
  ) %>%
  mutate(district = ifelse(is.na(district), Location, district)) %>%
  select(-Sample) %>%
  rename(Sample = Sample2)

# Plot Fraction of sites that are IBD
IBD_meta_combined <- Pk.IBD %>%
    as.data.frame() %>%
    select(sample1, sample2, fract_sites_IBD) %>%
    left_join(
      metadata %>%
        rename(sample1 = Sample) %>%
        rename(cluster1 = Cluster) %>%
        rename(district1 = district) %>% 
        rename(state1 = State), 
      by = "sample1"
    ) %>% 
    left_join(
      metadata %>%
        rename(sample2 = Sample) %>%
        rename(cluster2 = Cluster) %>%
        rename(district2 = district) %>% 
        rename(state2 = State), 
      by = "sample2" 
    )
  
## Fraction IBD within and between clusters
IBD_fract_plot <- IBD_meta_combined %>%
  unite("Clusters", c("cluster1", "cluster2"), sep = "-") %>% 
  ggplot(aes(x = Clusters, y = fract_sites_IBD, colour = Clusters)) +
  geom_boxplot() +
  theme(legend.position = "none", 
    legend.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  coord_flip() +
  scale_color_viridis_d() +
  ylab("Fraction of IBD sites")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD/fract_IBD_clusters.png", dpi = 300, height = 7, width = 7, IBD_fract_plot)

## Fraction IBD within and between districts
IBD_fract_plot <- IBD_meta_combined %>%
  unite("Districts", c("district1", "district2"), sep = "-") %>% 
  ggplot(aes(x = Districts, y = fract_sites_IBD, colour = Districts)) +
  geom_boxplot() +
  theme(legend.position = "none", 
    legend.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  coord_flip() +
  scale_color_viridis_d() +
  ylab("Fraction of IBD sites")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD/fract_IBD_locations.png", height = 30, dpi = 300, IBD_fract_plot)


# Base plots
Pk.label.cols <- c("Cluster")
Pk.legend.titles <- list("Cluster" = "Cluster")
Pk.prefixes <- metadata %>% 
  select(Cluster) %>%
  unique() %>% 
  as.list()

plot.IBDs(
  Pk.prefixes,
  IBD.cutoffs,
  Pk.IBD,
  metadata,
  Pk.label.cols,
  legend.titles = Pk.legend.titles
)

Pk.median.file <- "Pk_pairwise_IBD_median.txt"
Pk.median <- read.pairwise.matrix(Pk.median.file)

Pk.IBD.cutoffs <- c(fivenum(Pk.median), relaxed.IBD.cutoffs)

IBDs.plot.file <- "Pk_major_IBDs.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  metadata,
  "Cluster",
  legend.title = "Cluster"
)

# Base plots - districts
Pk.label.cols <- c("district")
Pk.legend.titles <- list("district" = "district")
Pk.prefixes <- metadata %>% 
  select(district) %>%
  unique() %>% 
  as.list()

plot.IBDs(
  Pk.prefixes,
  IBD.cutoffs,
  Pk.IBD,
  metadata,
  Pk.label.cols,
  legend.titles = Pk.legend.titles
)

IBDs.plot.file <- "Major_IBDs_district.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  metadata,
  "district",
  legend.title = "Districts"
)


# Base plots - states - might want to change the palette for fewer levels
Pk.label.cols <- c("States")
Pk.legend.titles <- list("States" = "States")
Pk.prefixes <- metadata %>% 
  select(State) %>%
  unique() %>% 
  as.list()

plot.IBDs(
  Pk.prefixes,
  IBD.cutoffs,
  Pk.IBD,
  metadata,
  Pk.label.cols,
  legend.titles = Pk.legend.titles
)


IBDs.plot.file <- "Major_IBDs_States.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  metadata,
  "State",
  legend.title = "States"
)

# Identify the two outlier at 3.42% - that connect the two clusters
## Filter for state to state, calc mean and arrange in desc order
IBD_meta_combined %>%
  filter(cluster1 == "Mf" & cluster2 == "Mn") %>%
  group_by(sample1) %>%
  summarise(mean = mean(fract_sites_IBD)) %>%
  arrange(desc(mean)) %>%
  head(n=12)

IBD_meta_combined %>%
  filter(cluster1 == "Mf" & cluster2 == "Mn") %>%
  filter(fract_sites_IBD >= 0.034) %>%
  arrange(desc(fract_sites_IBD)) # the top two samples should be two in the plot that are closest to Mn


############################################# WORK IN PORGRESS
## Plotting with ggplot
.create.vertices <- function(metadata, sample.col = "Sample") {
    vertices <- data.frame(metadata[, sample.col])
    names(vertices) <- sample.col
    
    vertices
  }

vertices <- .create.vertices(metadata)

IBD.graph <- graph_from_data_frame(Pk.IBD, directed = FALSE, vertices = vertices)

mylayout <- igraph::layout_as_tree(IBD.graph, circular = T) 
IBD.graph$layout = mylayout #store layout as a graph attribute


#another gotcha is that ig2ggplot needs both vertex names and vertex labels. 
#as of now you just have vertex names. 
V(IBD.graph)$label = V(IBD.graph)$name #store label as a vertex attrbute

IBD_plot <- MetamapsDB::ig2ggplot(IBD.graph, 
                      dfOnly = FALSE, 
                      labels = FALSE, 
                      metab = TRUE ) + 
    theme(legend.position = 'none')

ggsave("IBD_plot.png", dpi = 600, IBD_plot)


######################################################## TROUBLESHOOTING

Pk.IBD %>%
  select(1) %>%
  rename(sample2 = sample1) %>%
  rbind(
    Pk.IBD %>%
    select(2)
  ) %>%
  unique()


graph_from_data_frame(Pk.IBD, directed = FALSE, vertices)


Pk.IBD.file <-
  "Pk.hmm_fract.txt"
Pk.IBD <- read.delim(Pk.IBD.file)
IBD.graph <- graph_from_data_frame(Americas.IBD, directed = FALSE, vertices)


######################################################## 