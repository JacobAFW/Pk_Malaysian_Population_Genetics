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
  function(labels, palette = "Tableau 10") {
    label.names <- levels(labels)
    label.palette <-
      palette.colors(length(label.names), palette = "Tableau 10")
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

metadata1 <- read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.csv", col_names = FALSE) %>% 
        select(1, 6, 7) %>%
        rename(Sample = X1) %>%
        rename(Location = X6) %>%
        rename(Cluster = X7) %>% 
        as.data.frame()

### add a priori classfication of Sabah clusters based on trees

metadata2 <- metadata1 %>% 
    mutate(Cluster = ifelse(grepl("PK_SB_DNA_011", Sample) | # does the ID match any of these
                                grepl("PK_SB_DNA_091", Sample) | 
                                grepl("PK_SB_DNA_043", Sample) |
                                grepl("PK_SB_DNA_016", Sample) |  
                                grepl("PK_SB_DNA_092", Sample) | 
                                grepl("PK_SB_DNA_030", Sample) | 
                                grepl("PK_SB_DNA_093", Sample) | 
                                grepl("PK_SB_DNA_042", Sample) | 
                                grepl("PK_SB_DNA_063", Sample), "Mn", .$Cluster)) %>% # if not, just use values from X7 - clusters and Sabah
    mutate(Cluster = ifelse(Cluster == "Sabah", "Mf", .$Cluster)) %>% # if its the remaining Sabah samples, make them Mn, else keep them the same
    mutate(Sample = str_remove(Sample, "_DK.*")) %>%
    left_join(
      read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/hmmIBD_geo_clusters/geo_clusters.csv") %>%
        select(Sample, Geo_cluster)
    ) 

# wrangle the metadata to re-introduce the full sample names (removed in the last step to align with geo clusters) so that they align with the IBD data
metadata <- metadata2 %>%
  left_join(
    metadata1 %>%
      select(1) %>%
      rename(Sample2 = Sample) %>%
      mutate(Sample = str_remove(Sample2, "_DK.*"))
  ) %>% 
  mutate(Sample = Sample2) %>%
  select(-Sample2) 

#Mf
Pk.IBD.file <- "Pk_Mf.hmm_fract.txt"
Pk.IBD <- read.delim(Pk.IBD.file)

# Base plots
Pk.label.cols <- c("Geo_cluster")
Pk.legend.titles <- list("Geo_cluster" = "Geo_cluster")
Pk.prefixes <- metadata %>% 
  select(Geo_cluster) %>%
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

Pk.median.file <- "Pk_Mf_pairwise_IBD_median.txt"
Pk.median <- read.pairwise.matrix(Pk.median.file)

Pk.IBD.cutoffs <- c(fivenum(Pk.median), relaxed.IBD.cutoffs)

IBDs.plot.file <- "Mf_major_IBDs.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  metadata,
  "Geo_cluster",
  legend.title = "Geographic clusters within Mf"
)


#Mn
Pk.IBD.file <- "Pk_Mn.hmm_fract.txt"
Pk.IBD <- read.delim(Pk.IBD.file)

# Base plots
Pk.label.cols <- c("Geo_cluster")
Pk.legend.titles <- list("Geo_cluster" = "Geo_cluster")
Pk.prefixes <- metadata %>% 
  select(Geo_cluster) %>%
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

Pk.median.file <- "Pk_Mn_pairwise_IBD_median.txt"
Pk.median <- read.pairwise.matrix(Pk.median.file)

Pk.IBD.cutoffs <- c(fivenum(Pk.median), relaxed.IBD.cutoffs)

IBDs.plot.file <- "Mn_major_IBDs.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  metadata,
  "Geo_cluster",
  legend.title = "Geographic clusters within Mn"
)