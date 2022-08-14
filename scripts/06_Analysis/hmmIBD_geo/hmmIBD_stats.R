# Packages
library(dplyr)
library(readr)
library(stringr)
library(data.table)

# Define Edwins Function
.add.label.to.IBD <-
  function(IBD, metadata, label.col, sample.col = "Sample") {
    metadata <- metadata[, c(sample.col, label.col)]
    
    IBD <- merge(
      IBD,
      metadata,
      by.x = "sample1",
      by.y = sample.col,
      all.x = TRUE,
      sort = FALSE
    )
    colnames(IBD)[length(IBD)] <- paste0(label.col, 1)
    
    IBD <- merge(
      IBD,
      metadata,
      by.x = "sample2",
      by.y = sample.col,
      all.x = TRUE,
      sort = FALSE
    )
    colnames(IBD)[length(IBD)] <- paste0(label.col, 2)
    
    IBD
  }


.create.pairwise.matrix <- function(metadata, label.col) {
  unique.labels <- sort(unique(metadata[, label.col]))
  labels.length <- length(unique.labels)
  matrix(
    nrow = labels.length,
    ncol = labels.length,
    dimnames = list(unique.labels, unique.labels)
  )
}


.populate.pairwise.matrix.with.statistic <-
  function(pairwise.matrix,
           IBD,
           label.col,
           statistic) {
    label.cols1 <- paste0(label.col, 1)
    label.cols2 <- paste0(label.col, 2)
    labels <- rownames(pairwise.matrix)
    
    for (i in 1:nrow(pairwise.matrix)) {
      label1 <- labels[i]
      
      for (j in 1:ncol(pairwise.matrix)) {
        if (i < j)
          break
        label2 <- labels[j]
        
        e11 <- IBD[, label.cols1] == label1
        e12 <- IBD[, label.cols1] == label2
        e21 <- IBD[, label.cols2] == label1
        e22 <- IBD[, label.cols2] == label2
        
        pairwise.matrix[i, j] <-
          statistic(IBD[(e11 &
                           e22) | (e12 & e21), "fract_sites_IBD"]) 
      }
    }
    
    pairwise.matrix
  }


calculate.pairwise.matrix <-
  function(IBD,
           metadata,
           label.col,
           statistic,
           sample.col = "Sample") {
    IBD <-
      .add.label.to.IBD(IBD, metadata, label.col, sample.col = sample.col)
    
    pairwise.matrix <-
      .create.pairwise.matrix(IBD, paste0(label.col, 1))
    
    .populate.pairwise.matrix.with.statistic(pairwise.matrix, IBD, label.col, statistic)
  }


save.pairwise.matrix <- function(pairwise.matrix, file) {
  write.table(pairwise.matrix,
              file,
              quote = FALSE,
              sep = "\t",
              na = "")
}


.add.sentinel.to.upper.triangle <- function(pairwise.matrix) {
  random <- rnorm(1)
  while (any(pairwise.matrix == random, na.rm = TRUE))
    random <- rnorm(1)
  
  pairwise.matrix[upper.tri(pairwise.matrix)] <- random
  
  pairwise.matrix
}


.transform.pairwise.matrix.to.table <- function(pairwise.matrix) {
  pairwise.matrix <- .add.sentinel.to.upper.triangle(pairwise.matrix)
  
  randoms <- pairwise.matrix[upper.tri(pairwise.matrix)]
  random <- randoms[1]
  if (any(is.na(randoms)) ||
      !all(randoms == random))
    warning("Matrix upper triangle does not contain sentinels.")
  
  pairwise.table <- as.data.frame(as.table(pairwise.matrix))
  pairwise.table[!(pairwise.table[, "Freq"] %in% random), ]
}


.concatenate.pairwise.labels <-
  function(pairwise.table, sep = "-") {
    labels <-
      do.call(paste, c(pairwise.table[, c("Var1", "Var2")], sep = sep))
    as.data.frame(cbind(labels, pairwise.table[, "Freq"]))
  }


calculate.pairwise.table <- function(IBD,
                                     metadata,
                                     label.col,
                                     statistics,
                                     sample.col = "Sample") {
  pairwise.tables <- NULL
  # safeguard against non-vector
  statistics <- c(statistics)
  stats <- names(statistics)
  
  
  for (i in seq_along(statistics)) {
    pairwise.matrix <- calculate.pairwise.matrix(IBD,
                                                 metadata,
                                                 label.col,
                                                 statistics[[i]],
                                                 sample.col = sample.col)
    pairwise.table <-
      .transform.pairwise.matrix.to.table(pairwise.matrix)
    
    pairwise.table <- .concatenate.pairwise.labels(pairwise.table)
    if (is.null(stats)) {
      names(pairwise.table) <- c(label.col, LETTERS[i])
    } else {
      names(pairwise.table) <- c(label.col, stats[i])
    }
    
    if (is.null(pairwise.tables)) {
      pairwise.tables <- pairwise.table
      next
    }
    
    pairwise.tables <-
      merge(pairwise.tables,
            pairwise.table,
            by = label.col,
            all = TRUE)
  }
  
  pairwise.tables
}


save.pairwise.table <- function(pairwise.table, file) {
  write.table(
    pairwise.table,
    file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}


statistics <- c(min, function(x)
  quantile(x, probs = 0.25),
  median, mean, function(x)
    quantile(x, probs = 0.75), max)

names(statistics) <- c("Min.", "1st Qu.", "Median", "Mean",
                       "3rd Qu.", "Max.")


# Pk Clusters in Malaysia

## stat medians
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
      read_csv("geo_clusters.csv") %>%
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
  


# Mf
Pk.IBD <- read_table("Pk_Mf.hmm_fract.txt") 
Pk.pairwise.median <- calculate.pairwise.matrix(Pk.IBD, metadata, "Geo_cluster", median)
save.pairwise.matrix(Pk.pairwise.median, "Pk_Mf_pairwise_IBD_median.txt")

## stat summary
Pk.pairwise.summary <- calculate.pairwise.table(Pk.IBD, metadata, "Geo_cluster", statistics)
save.pairwise.table(Pk.pairwise.summary, "Pk_Mf_pairwise_summary.txt")


# Mn
Pk.IBD <- read_table("Pk_Mn.hmm_fract.txt") 
Pk.pairwise.median <- calculate.pairwise.matrix(Pk.IBD, metadata, "Geo_cluster", median)
save.pairwise.matrix(Pk.pairwise.median, "Pk_Mn_pairwise_IBD_median.txt")

## stat summary
Pk.pairwise.summary <- calculate.pairwise.table(Pk.IBD, metadata, "Geo_cluster", statistics)
save.pairwise.table(Pk.pairwise.summary, "Pk_Mn_pairwise_summary.txt")

