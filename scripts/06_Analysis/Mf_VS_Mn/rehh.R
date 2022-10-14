# Packages
library(rehh)
library(tidyverse)
library(data.table)
library(R.utils)

# WHOLE GENOME iHS

## need to read in chr files independently and build up an object resulting from the scan_hh outputs
### first create an empty df to be appended in the loop and a list of chr names
### the loop then creates the vcf file name (hap_file), and reads in the file to hh, which is then parsed into the scan_hh function to produce an object i that is appended to the growing wgscan object
### the resulting wgscan produces the same results as the scan_hh object above
### this can then be parsed into ihh2ihs

chromsomes_files <- c("ordered_PKNH_01_v2", "ordered_PKNH_02_v2", "ordered_PKNH_03_v2", "ordered_PKNH_04_v2", "ordered_PKNH_05_v2",
    "ordered_PKNH_06_v2", "ordered_PKNH_07_v2", "ordered_PKNH_08_v2", "ordered_PKNH_09_v2", "ordered_PKNH_10_v2", "ordered_PKNH_11_v2",
    "ordered_PKNH_12_v2", "ordered_PKNH_13_v2", "ordered_PKNH_14_v2")

wgscan <- tibble("CHR" = NA, "POSITION" = NA, "FREQ_A" = NA, "FREQ_D" = NA, "NHAPLO_A" = NA, "NHAPLO_D" = NA, "IHH_A" = NA, "IHH_D" = NA, "IES" = NA, "INES" = NA)

for(i in chromsomes_files) {
    # haplotype file name for each chromosome
    hap_file = paste(i, ".vcf.gz", sep = "") # filename pattern
    # create internal representation
    hh <- data2haplohh(hap_file = hap_file, 
                    chr.name = i, 
                    polarize_vcf = FALSE, # if the AA key is absent
                    min_perc_geno.mrk = 50, # discard markers genotyped on < 50% of haplotypes
                    min_maf = 0.05, # discard markers with a minor allele frequency of < 0.05)
                    vcf_reader = "data.table") # use this package
    scan <- scan_hh(hh) # perform scan on a single chromosome (calculate iHH/iES values)
    wgscan <- wgscan %>% # append the wgscan object every round with additional data
    rbind(scan)
}

wgs_ihs <- wgscan %>%
    na.omit() %>%
    ihh2ihs(min_maf = 0.05, # default
            freqbin = 0.025 # default
            )

### Output is a list with two elements:
### ihs: a data frame with iHS and corresponding p-value piHS (p-value in a −log10 scale assuming the iHS are normally distributed under the neutral hypothesis)
### frequency.class: a data frame summarizing the derived allele frequency bins used for standardization and mean and standard deviation of the un-standardized values
wgs_ihs$ihs
wgs_ihs$frequency.class


# Scan for candidate regions 
## First update iHS table so that its compatible
wgs_ihs_2 <- wgs_ihs
wgs_ihs_2$ihs <- wgs_ihs_2$ihs %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    mutate(CHR = as.numeric(CHR)) %>%
    as.data.frame()

candidate_regions <- calc_candidate_regions(wgs_ihs_2,
                                 threshold = 4,
                                 pval = TRUE,
                                 window_size = 1000,
                                 overlap = 100,
                                 min_n_extr_mrk = 3) %>%
                                 add_column(Stat = "iHS") 


## Plot iHS and iHS p-values

iHS <- wgs_ihs$ihs %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    ggplot(aes(x = CHR, y = IHS, colour = CHR)) +
    geom_jitter() + 
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    scale_x_discrete() +
    xlab("Chromsome") +
    ylab("iHS") 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/iHS.png", dpi = 600, width = 10, iHS)

iHS_pvalue <- wgs_ihs$ihs %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    ggplot(aes(x = CHR, y = LOGPVALUE, colour = CHR)) +
    geom_jitter() + 
    geom_hline(yintercept = 4, linetype = 2, colour = "#1F968BFF") + # adjusted with bonferonni
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    xlab("Chromosomes") +
    ylab("iHS p-value")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/iHS_pvalue.png", dpi = 600, width = 10, iHS_pvalue)


## Plotting using rehh

png('/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/iHS_rehh.png')
manhattanplot(wgs_ihs)
dev.off()

png('/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/iHS_p_rehh.png')
manhattanplot(wgs_ihs, pval = TRUE, threshold = 4)
dev.off()


## QQ-Plot
png('/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/iHS_QQ.png')
distribplot(wgs_ihs$ihs$IHS, 
            xlab = "iHS", 
            qqplot = TRUE)
dev.off()

#######################################################################################


# Comparing EHH between clusters

# CALCULATE Rsb  and XP-EHH

## pairwise population statistic 
## ratio of EHHS between populations

chromsomes_files <- c("ordered_PKNH_01_v2", "ordered_PKNH_02_v2", "ordered_PKNH_03_v2", "ordered_PKNH_04_v2", "ordered_PKNH_05_v2",
    "ordered_PKNH_06_v2", "ordered_PKNH_07_v2", "ordered_PKNH_08_v2", "ordered_PKNH_09_v2", "ordered_PKNH_10_v2", "ordered_PKNH_11_v2",
    "ordered_PKNH_12_v2", "ordered_PKNH_13_v2", "ordered_PKNH_14_v2")

### Mn
wgscan_Mn <- tibble("CHR" = NA, "POSITION" = NA, "FREQ_A" = NA, "FREQ_D" = NA, "NHAPLO_A" = NA, "NHAPLO_D" = NA, "IHH_A" = NA, "IHH_D" = NA, "IES" = NA, "INES" = NA)

for(i in chromsomes_files) {
    # haplotype file name for each chromosome
    hap_file = paste(i, "_Mn.vcf.gz", sep = "") # filename pattern
    # create internal representation
    hh <- data2haplohh(hap_file = hap_file, 
                    chr.name = i, 
                    polarize_vcf = FALSE, # if the AA key is absent
                    min_perc_geno.mrk = 25, # discard markers genotyped on < 25% of haplotypes
                    min_maf = 0.05, # discard markers with a minor allele frequency of < 0.05)
                    vcf_reader = "data.table") # use this package
    scan <- scan_hh(hh) # perform scan on a single chromosome (calculate iHH/iES values)
    wgscan_Mn <- wgscan_Mn %>% # append the wgscan object every round with additional data
    rbind(scan)
}

### Mf
wgscan_Mf <- tibble("CHR" = NA, "POSITION" = NA, "FREQ_A" = NA, "FREQ_D" = NA, "NHAPLO_A" = NA, "NHAPLO_D" = NA, "IHH_A" = NA, "IHH_D" = NA, "IES" = NA, "INES" = NA)

for(i in chromsomes_files) {
    # haplotype file name for each chromosome
    hap_file = paste(i, "_Mf.vcf.gz", sep = "") # filename pattern
    # create internal representation
    hh <- data2haplohh(hap_file = hap_file, 
                    chr.name = i, 
                    polarize_vcf = FALSE, # if the AA key is absent
                    min_perc_geno.mrk = 50, # discard markers genotyped on < 50% of haplotypes
                    min_maf = 0.05, # discard markers with a minor allele frequency of < 0.05)
                    vcf_reader = "data.table") # use this package
    scan <- scan_hh(hh) # perform scan on a single chromosome (calculate iHH/iES values)
    wgscan_Mf <- wgscan_Mf %>% # append the wgscan object every round with additional data
    rbind(scan)
}

# Rsb
Rsb <- ines2rsb(
        scan_pop1 = wgscan_Mn,
        scan_pop2 = wgscan_Mf,
        popname1 = "Mn",
        popname = "Mf"
    )

# XP-EHH
XP_EHH <- ies2xpehh(
        scan_pop1 = wgscan_Mn,
        scan_pop2 = wgscan_Mf,
        popname1 = "Mn",
        popname = "Mf"
    )


# Scan for candidate regions 
candidate_regions_Rsb <- calc_candidate_regions(Rsb,
                                 threshold = 4,
                                 pval = TRUE,
                                 window_size = 1000,
                                 overlap = 100,
                                 min_n_extr_mrk = 3) %>%
                                 add_column(Stat = "Rsb")

candidate_regions_XP_EHH <- calc_candidate_regions(XP_EHH,
                                 threshold = 4,
                                 pval = TRUE,
                                 window_size = 1000,
                                 overlap = 100,
                                 min_n_extr_mrk = 3) %>%
                                 add_column(Stat = "XP_EHH")

## Save candidate regions for all stats
candidate_regions %>%
    rbind(candidate_regions_Rsb) %>%
    rbind(candidate_regions_XP_EHH) %>%
    write_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/candidate_regions.csv")


# Plot Rsb and p-value
Rsb_plot <- Rsb %>%
    na.omit() %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    rename("Rsb" = "RSB_Mn_Mf") %>%
    ggplot(aes(x = CHR, y = Rsb, colour = CHR)) +
    geom_jitter() + 
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    xlab("Chromsome") +
    ylab("Rsb") 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/Rsb.png", dpi = 600, width = 10, Rsb_plot)

Rsb_pvalue <- Rsb %>%
    na.omit() %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    ggplot(aes(x = CHR, y = LOGPVALUE, colour = CHR)) +
    geom_jitter() + 
    geom_hline(yintercept = 4, linetype = 2, colour = "#1F968BFF") + # adjusted with bonferonni
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    xlab("Chromosomes") +
    ylab("Rsb p-value")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/Rsb_pvalue.png", dpi = 600, width = 10, Rsb_pvalue)

#Plot XP-EHH and p-value

XP_EHH_plot <- XP_EHH %>%
    na.omit() %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    ggplot(aes(x = CHR, y = XPEHH_Mn_Mf, colour = CHR)) +
    geom_jitter() + 
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    xlab("Chromsome") +
    ylab("XP_EHH") 

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/XP_EHH.png", dpi = 600, width = 10, XP_EHH_plot)

XP_EHH_pvalue <- XP_EHH %>%
    na.omit() %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    ggplot(aes(x = CHR, y = LOGPVALUE, colour = CHR)) +
    geom_jitter() + 
    geom_hline(yintercept = 4, linetype = 2, colour = "#1F968BFF") + # adjusted with bonferonni
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    xlab("Chromosomes") +
    ylab("XP_EHH p-value")

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/XP_EHH_pvalue.png", dpi = 600, width = 10, XP_EHH_pvalue)

# Rsb vs XP-EHH

Rsb_XP_comaprison <- Rsb %>%
    left_join(
        XP_EHH %>%
            select(1:3)
        ) %>%
    na.omit() %>%
    ggplot(aes(x = RSB_Mn_Mf, y = XPEHH_Mn_Mf, colour = CHR)) +
    geom_point() + 
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    xlab("Rsb") +
    ylab("XP_EHH") +
    geom_abline()

ggsave("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk_Mn_vs_Mf/XP_Rsb_comparison.png", dpi = 600, width = 10, Rsb_XP_comaprison)

