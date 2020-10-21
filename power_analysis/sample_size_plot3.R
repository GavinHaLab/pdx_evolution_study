#Sample Size Plot 3
#August 6th, 2020
#Anna Hoge
#Ha Lab


library(tidyverse)
library(ggplot2)
library(scales)
library(data.table)


#make sample size by known frequency dataframe
calculate_sample_size <- function(p0, difference) {
    alpha = 0.05
    beta = 0.2
    p = p0 + difference
    n=p0*(1-p0)*((qnorm(1-alpha)+qnorm(1-beta)*sqrt(p*(1-p)/p0/(1-p0)))/(p-p0))^2
    return(ceiling(n))
}
p0 <- seq(0.01, 0.99, 0.01)
n_diff_neg_0.1 <- calculate_sample_size(p0, -0.1)
power_analysis_df <- data.frame(known_frequency = p0,
                                n_diff_neg_0.1 = n_diff_neg_0.1)
                    
#count sample size of PT-PDX pairings in our study for each cohort
pairings_df <- fread("pairing_discordance_values_1mb_thresholds.txt")
names(pairings_df) <- c("cohort", "sample1", "sample2", "discordance")
cohort_num_pairings <- pairings_df %>%
                          filter(str_detect(sample1, "PT")) %>%
                          group_by(cohort) %>%
                          count()

#count sample size of PT-PDX pairings in our study for each cancer type
breast_counts <- cohort_num_pairings %>%
                filter(cohort == "bcm_snp_breast" | cohort == "hci_snp_breast" | cohort == "europdx_wgs_brca" | cohort == "hci_wes_breast" | cohort == "jax_snp_breast" | cohort == "wust_wes_breast")
breast_n <- sum(breast_counts$n) #BRCA
colorectal_counts <- cohort_num_pairings %>%
                    filter(cohort == "europdx_wgs_crc" | cohort == "jax_snp_colorectal" | cohort == "pdmr_wes_colorectal")
colorectal_n <- sum(colorectal_counts$n) #COAD
bladder_counts <- cohort_num_pairings %>%
                      filter(cohort == "jax_snp_bladder" | cohort == "pdmr_wes_bladder")
bladder_n <- sum(bladder_counts$n) #BLCA
gbm_counts <- cohort_num_pairings %>%
                  filter(cohort == "jax_snp_gbm")
gbm_n <- sum(gbm_counts$n) #GBM
melanoma_counts <- cohort_num_pairings %>%
                       filter(cohort == "jax_snp_melanoma" | cohort == "wistar_wes_melanoma")
melanoma_n <- sum(melanoma_counts$n) #SKCM
skin_squamous_counts <- cohort_num_pairings %>%
                       filter(cohort == "pdmr_wes_skin")
skin_squamous_n <- sum(skin_squamous_counts$n) #N/A
sarcoma_counts <- cohort_num_pairings %>%
                      filter(cohort == "jax_snp_sarcoma" | cohort == "pdmr_wes_sarcoma")
sarcoma_n <- sum(sarcoma_counts$n) #SARC
gastric_counts <- cohort_num_pairings %>%
                      filter(cohort == "jax_snu_wes_gastric")
gastric_n <- sum(gastric_counts$n) #STAD
luad_counts <- cohort_num_pairings %>%
                   filter(cohort == "mdacc_wes_luad" | cohort == "pdmr_wes_lung")
luad_n <- sum(luad_counts$n) #LUAD
lusc_counts <- cohort_num_pairings %>%
                    filter(cohort == "mdacc_wes_lusc")
lusc_n <- sum(lusc_counts$n) #LUSC
headneck_counts <- cohort_num_pairings %>%
                       filter(cohort == "pdmr_wes_headneck")
headneck_n <- sum(headneck_counts$n) #HNSC
renal_counts <- cohort_num_pairings %>%
                    filter(cohort == "pdmr_wes_renal")
renal_n <- sum(renal_counts$n) #KIRC ###check this
hcc_counts <- cohort_num_pairings %>%
                  filter(cohort == "sibs_snp_hcc")
hcc_n <- sum(hcc_counts$n) #LIHC
pancreatic_counts <- cohort_num_pairings %>%
                         filter(cohort == "wustl_wes_pancreatic")
pancreatic_n <- sum(pancreatic_counts$n) #PAAD
sample_size_df <- data.table(cancer_type = c("breast", "colorectal", "bladder", "gbm", "melanoma", "skin_squamous", "sarcoma", "gastric", "luad", "lusc", "head_neck", "renal", "hcc", "pancreatic"),
                             sample_size = c(breast_n, colorectal_n, bladder_n, gbm_n, melanoma_n, skin_squamous_n, sarcoma_n, gastric_n, luad_n, lusc_n, headneck_n, renal_n, hcc_n, pancreatic_n))

#convert arm counts table to frequencies
arm_counts_df <- fread("tumors_with_arm_calls.txt")
arm_freq_df <- arm_counts_df %>%
    gather(arm, call, 14:52) %>%
    select(Type, arm, call) %>%
    group_by(Type, arm) %>%
    count(call) %>%
    mutate(arm_status = paste0(arm, ".", call)) %>%
    mutate(frequency = n / sum(n))

#get recurrent arms with frequencies per cancer type
recurrent_arms_df <- arm_freq_df %>%
                         filter(call == 1 | call == -1) %>%
                         filter(frequency >= 0.25)     
recurrent_arms_frequencies_df <- recurrent_arms_df %>%
                                     mutate(sample_size = ifelse(Type == "BRCA", breast_n, ifelse(Type == "COAD", colorectal_n, ifelse(Type == "BLCA", bladder_n, ifelse(Type == "GBM", gbm_n, ifelse(Type == "SKCM", melanoma_n, ifelse(Type == "SARC", sarcoma_n, ifelse(Type == "STAD", gastric_n, ifelse(Type == "LUAD", luad_n, ifelse(Type == "LUSC", lusc_n, ifelse(Type == "HNSC", headneck_n, ifelse(Type == "KIRC", renal_n, ifelse(Type == "LIHC", hcc_n, ifelse(Type == "PAAD", pancreatic_n, 0)))))))))))))) %>%
                                     filter(sample_size > 0) %>%
                                     arrange(-sample_size)

#plot with power analysis
ggplot() + 
    geom_point(data = recurrent_arms_frequencies_df, aes(x = frequency, y = sample_size, color = Type)) +
    scale_color_manual("Cancer Types",
                       breaks = c("colorectal", "breast", "gastric", "luad", "bladder", "head_neck", "sarcoma", "lusc", "gbm", "melanoma", "hcc", "pancreatic", "renal"),
                       values = c("#0b43bd", "#e364d6", "#a6b7e3", "#D3D3D3", "#e0c02f", "#7d0e0a", "#e4e80c", "#778899", "black", "#808080", "#82d7b4", "#881db3", "#c97712")) +
    geom_line(data = power_analysis_df, aes(x = p0, y = n_diff_neg_0.1, size = "")) +
    scale_size_manual("Minimum Sample Size to Detect\n10% Absolute Decrease in\nFrequency of Recurrent Arm",
                       breaks = c(""),
                       values = c(0.5)) +
    guides(size = guide_legend(order = 1),
           color = guide_legend(order = 2)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          panel.grid.minor = element_blank()) +
    xlab("Frequency of Recurrently Altered Arm") +
    ylab("Sample Size") +
    scale_x_continuous(limit = c(0.1, 1),
                       breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                       expand = c(0, 0.005)) +
    scale_y_continuous(limit = c(0, 160),
                       breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160),
                       expand = c(0, 1)) +
    ggsave(file = "power_analysis_plot.pdf", width = 15, height = 8)



