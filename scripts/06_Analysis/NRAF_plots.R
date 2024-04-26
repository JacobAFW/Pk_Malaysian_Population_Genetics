# Converting VCF to GDS for estMOI
library(tidyverse)
library(SeqArray)
library(gridExtra)
library(ggpubr)
library(moimix)

seqVCF2GDS("moimix/PK_consensus_masked.vcf.gz", "moimix/PK_consensus_filtered_pass.gds", storage.option = "LZ4_RA")
isolates <- seqOpen("moimix/PK_consensus_filtered_pass.gds") # read into R
seqSummary(isolates)
sample.id <- seqGetData(isolates, "sample.id")
seqCheck(isolates)
str(seqGetData(isolates, "annotation/info/AC"))

head(seqGetData(isolates, "sample.id"))

# moimix analysis
coords <- getCoordinates(isolates) # get genomic coordinates of all variants

## estimating BAF matrix 
isolate_baf <- bafMatrix(isolates)

baf_df <- isolate_baf$coords %>%
  cbind(
    as.data.frame(isolate_baf$baf_site) %>% 
      rename(baf = "isolate_baf$baf_site")
    ) %>%
  left_join(
    isolate_baf$baf_matrix %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("variant.id") %>%
      mutate(variant.id = as.numeric(variant.id)) 
    ) %>%
  mutate(chromosome = str_remove(chromosome, ".*PKNH_")) %>%
  mutate(chromosome = str_remove(chromosome, "_v2"))

## Plot BAF for different Fws samples with ggplot
plot_baf <- function(DATA, SAMPLE){
DATA %>%
    ggplot(aes_string(x = "variant.id", y = SAMPLE, colour = "chromosome")) + # use aes string so that we can paste the str in from the function arguments - passing the sample in not as a string does not work
    geom_point() +
    scale_colour_manual(values = c(rep_len(c("#404788FF", "#1F968BFF"), length(unique(baf_df$chromosome))))) +
    xlab("Genome Position") +
    ylab("NRAF") + 
    theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 25)) +
    guides(colour = guide_legend(nrow = 1))
}

baf_plot <- plot_baf(baf_df, "ERR2214850")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR2214850_72.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR2214856")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR2214856_80.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR3374041")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR3374041_78.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985386")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985386_67.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985395")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985395_54.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985396")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985396_77.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985397")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985397_82.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985405")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985405_62.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985410")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985410_62.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985417")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985417_84.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_006 = "PK_SB_DNA_006_DKDL210002135-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_006")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_006_65.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_032 = "PK_SB_DNA_032_DKDL210002161-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_032")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_032_78.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_039 = "PK_SB_DNA_039_DKDL210002168-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_039")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_039_81.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_048 = "PK_SB_DNA_048_DKDL210002177-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_048")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_048_82.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_051 = "PK_SB_DNA_051_DKDL210002180-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_051")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_051_81.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_055 = "PK_SB_DNA_055_DKDL210002184-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_055")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_055_77.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_083 = "PK_SB_DNA_083_DKDL210002212-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_083")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_083_59.png", width = 14, dpi = 300, baf_plot)



# Low MOI
baf_plot <- baf_df %>%
    rename(PK_SB_DNA_083 = "PK_SB_DNA_083_DKDL210002212-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_083")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_083_59.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_069 = "PK_SB_DNA_069_DKDL210002198-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_069")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_069_86.png", width = 14, dpi = 300, baf_plot)

baf_plot <- baf_df %>%
    rename(PK_SB_DNA_068 = "PK_SB_DNA_068_DKDL210002197-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_068")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_068_92.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR2214842")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR2214842_87.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985375")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985375_87.png", width = 14, dpi = 300, baf_plot)

baf_plot <- plot_baf(baf_df, "ERR985376")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/ERR985376_93.png", width = 14, dpi = 300, baf_plot)


# No MOI
baf_plot <- baf_df %>%
    rename(PK_SB_DNA_004 = "PK_SB_DNA_004_DKDL210002133-1a_HWHGKDSXY_L4") %>%
    plot_baf(., "PK_SB_DNA_004")
ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/PK_SB_DNA_004_99.png", width = 14, dpi = 300, baf_plot)



# Combined 

plot_baf_grid <- function(DATA, SAMPLE){
DATA %>%
    ggplot(aes_string(x = "variant.id", y = SAMPLE, colour = "chromosome")) + # use aes string so that we can paste the str in from the function arguments - passing the sample in not as a string does not work
    geom_point(size = 1) +
    scale_colour_manual(values = c(rep_len(c("#404788FF", "#1F968BFF"), length(unique(baf_df$chromosome))))) +
    ylab("NRAF") + 
    theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none", 
        legend.title = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15)) +
    guides(colour = guide_legend(nrow = 1)) 
}

plot_baf_grid2 <- function(DATA, SAMPLE){
x_axis <- DATA %>%
    group_by(chromosome) %>%
    summarise(variant.id = median(variant.id)) %>%
    mutate(chromosome = str_remove(chromosome, "^0"))

DATA %>%
    ggplot(aes_string(x = "variant.id", y = SAMPLE, colour = "chromosome")) + # use aes string so that we can paste the str in from the function arguments - passing the sample in not as a string does not work
    geom_point(size = 1) +
    scale_colour_manual(values = c(rep_len(c("#404788FF", "#1F968BFF"), length(unique(baf_df$chromosome))))) +
    ylab("NRAF") + 
    xlab("Chromosome") +
    theme(legend.position = "none", 
        legend.title = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(nrow = 1)) +
    scale_x_continuous(breaks = x_axis$variant.id, labels = x_axis$chromosome)
}



# Figure in manuscript
ERR985395 <- plot_baf_grid(baf_df, "ERR985395") + scale_x_continuous(name = "ERR985395 (Fws 0.54)", position = "top") 
ERR985405 <- plot_baf_grid(baf_df, "ERR985405") + scale_x_continuous(name = "ERR985405 (Fws 0.62)", position = "top") 
ERR2214850 <- plot_baf_grid(baf_df, "ERR2214850") + scale_x_continuous(name = "ERR2214850 (Fws 0.72)", position = "top") 
ERR985417 <- plot_baf_grid(baf_df, "ERR985417") + scale_x_continuous(name = "ERR985417 (Fws 0.84)", position = "top") 
ERR985376 <- plot_baf_grid(baf_df, "ERR985376") + scale_x_continuous(name = "ERR985376 (Fws 0.93)", position = "top") 
PK_SB_DNA_004 <- baf_df %>% 
    rename(PK_SB_DNA_004 = "PK_SB_DNA_004_DKDL210002133-1a_HWHGKDSXY_L4") %>% 
    plot_baf_grid(., "PK_SB_DNA_004") + scale_x_continuous(name = "PK_SB_DNA_004 (Fws 0.99)", position = "top") 

grid_plot <- grid.arrange(nrow = 3,
    ERR985395,
    ERR985405,
    ERR2214850,
    ERR985417,
    ERR985376,
    PK_SB_DNA_004
    )

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/grid_plot.png", width = 12, height = 10, dpi = 600, grid_plot)


# Figure in Sup

fws_high <- read_tsv("moimix/fws_MOI.tsv") %>% filter(Proportion < 0.95)


baf_df <- baf_df %>% 
    rename_at(vars(starts_with("PK_SB_DNA")), ~str_remove(., "_DKD.*")) 

nraf_plots <- list()

for (i in pull(fws_high, Sample)){
    name <- str_remove(paste0(i), "_DK.*")

    Fws <- fws_high %>% filter(Sample == paste0(i)) %>% pull(Proportion) %>% round(2)
    
    nraf_plots[[name]]<- plot_baf_grid2(baf_df, paste0(name)) + 
        ggtitle(paste0(name, " (Fws ", Fws, ")"))
}

grid_plot <- arrangeGrob(
    grobs = nraf_plots,
    ncol = 3
)

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/MAF-Fws_plots/grid_plot_sup2.png", width = 16, height = 30, dpi = 300, grid_plot)
