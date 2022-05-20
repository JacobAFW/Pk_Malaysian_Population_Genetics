library(tidyverse)

Fst_matrix <- read_table("Pk.fst", col_names=T) %>%
    mutate(CHR = str_remove(SNP, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2.*"))

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

# Dot plots
Fst_plot <- Fst_matrix %>% 
    filter(CHR == "08") %>%
    ggplot(aes(x = POS, y = FST, colour = CHR)) +
    geom_point() 
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 

ggsave("Fst_dot_plot.png", dpi = 600, Fst_plot)

# Sliding window
window_size <- 10000

Fst_plot_window <- Fst_matrix %>%
    mutate(Window = (floor(POS/window_size) * window_size)+ (window_size/2)) %>%
    group_by(CHR, Window) %>%
    summarise(Window_Fst = mean(FST), count = n()) %>% # caculate the mean Fst for each window in each chr
    filter(count > 10) %>% # remove windows with low counts
    ggplot(aes(x = Window, y = Window_Fst, colour = CHR)) +
    geom_point() +
    facet_wrap(~CHR) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
    ylab("Fst") +
    xlab("Windows (10000)") +
    scale_color_discrete("Chr") +
    geom_hline(yintercept = 0.3, colour = "#1F968BFF")

ggsave("Fst_sliding_window_plot.png", dpi = 600, width = 14, Fst_plot_window)