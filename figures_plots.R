##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##            CREATE FIGURES AND TABLES FOR SIMULATION EXPERIMENT           ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
library(tidyverse)
library(patchwork)
library(knitr)
library(readr)
library(ggtext)
library(latex2exp)

# Date and working directory
cur_date <- "20230314"
wd <- "C:/r_proj/tsln_sralc_sims/"

# Load data
pr_all <- readRDS(paste0(wd, "results/", cur_date, "/pr_all.rds"))
sim_list <- readRDS(paste0(wd, "results/", cur_date, "/sim_list.rds"))

# Create plot folders
loc_plots <- paste0(wd, "results/", cur_date, "/plots")
if(!file.exists(loc_plots))dir.create(loc_plots)

## ---- Summarise objects ---- ## ----------------------------------------------

# SR - median across 100 reps ####
SR <- bind_rows(lapply(pr_all$SR, 
                       FUN = function(x){cbind(x, rep_counter = 1:50)}), 
                .id = "QaS") %>% 
  group_by(QaS) %>% 
  summarise(SR = median(SR))

# ACL - median across 100 reps ####
foo <- function(x){cbind(ACL = x, rep_counter = 1:50)}
ACL <- bind_rows(lapply(pr_all$WOLSB, 
                       FUN = function(x){data.frame(ACL = x, rep_counter = 1:50)}), 
                .id = "QaS") %>% 
  group_by(QaS) %>% 
  summarise(ACL = median(ACL))
rm(foo)

# Coverage - QaS by missing ####
coverage <- bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  mutate(in_ci = ifelse(prop > s2LN.lower & prop < s2LN.upper, 1, 0)) %>% 
  group_by(QaS, rep_counter, missing) %>% 
  summarise(in_ci = mean(in_ci), .groups = "drop") %>% 
  group_by(QaS, missing) %>% 
  summarise(coverage = median(in_ci), .groups = "drop")

## ---- Create Plots ---- ## ---------------------------------------------------

# Plot data
plot_df <- bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  group_by(QaS, rep_counter, missing) %>% 
  summarise(RRMSE = mean(s2LN_RRMSE),
            ARB = mean(s2LN_ARB), 
            .groups = "drop") %>% 
  group_by(QaS, missing) %>% 
  summarise(MRRMSE = median(RRMSE),
            MRRMSE_l = quantile(RRMSE, probs = 0.25),
            MRRMSE_u = quantile(RRMSE, probs = 0.75),
            MARB = median(ARB),
            MARB_l = quantile(ARB, probs = 0.25),
            MARB_u = quantile(ARB, probs = 0.75),
            .groups = "drop") %>% 
  left_join(.,SR, by = "QaS") %>% 
  left_join(.,ACL, by = "QaS") %>% 
  left_join(.,coverage, by = c("QaS", "missing"))

# SR vs ACL ----
plot_df %>% 
  ggplot(aes(x = SR, y = ACL, col = missing, group = missing))+
  theme_bw()+
  geom_line()+
  geom_point()
ggsave(paste0(loc_plots, "/SR vs ACL.png"), width = 10, height = 8.35)

# MRRMSE and MARB vs ACL ----
(plot_df %>% 
  ggplot(aes(y = MRRMSE, ymin = MRRMSE_l, ymax = MRRMSE_u, 
             x = ACL, col = missing, group = missing))+
  theme_bw()+
  geom_errorbar()+
  geom_line()+
  geom_point()+
    theme(legend.position = "bottom"))+
(plot_df %>% 
   ggplot(aes(y = MARB, ymin = MARB_l, ymax = MARB_u, 
              x = ACL, col = missing, group = missing))+
   theme_bw()+
   geom_errorbar()+
   geom_line()+
   geom_point()+
   theme(legend.position = "none"))
ggsave(paste0(loc_plots, "/MRRMSE and MARB_ACL.png"), width = 10, height = 8.35)

# MRRMSE and MARB vs SR ----
(plot_df %>% 
    ggplot(aes(y = MRRMSE, ymin = MRRMSE_l, ymax = MRRMSE_u, 
               x = SR, col = missing, group = missing))+
    theme_bw()+
    geom_errorbar()+
    geom_line()+
    geom_point()+
    theme(legend.position = "bottom"))+
(plot_df %>% 
   ggplot(aes(y = MARB, ymin = MARB_l, ymax = MARB_u, 
              x = SR, col = missing, group = missing))+
   theme_bw()+
   geom_errorbar()+
   geom_line()+
   geom_point()+
   theme(legend.position = "none"))
ggsave(paste0(loc_plots, "/MRRMSE and MARB_SR.png"), width = 10, height = 8.35)

# Coverage vs SR and ALC ----
(plot_df %>% 
  ggplot(aes(y = coverage, 
             x = SR, col = missing, group = missing))+
  theme_bw()+
  geom_hline(yintercept = 0.95)+
  geom_line()+
  geom_point()+
  theme(legend.position = "bottom"))+
(plot_df %>% 
   ggplot(aes(y = coverage, 
              x = ACL, col = missing, group = missing))+
   theme_bw()+
   geom_hline(yintercept = 0.95)+
   geom_line()+
   geom_point()+
   theme(legend.position = "none"))
ggsave(paste0(loc_plots, "/Coverage.png"), width = 10, height = 8.35)

# Bias reduction ----
foo <- function(x){cbind(x, JaS = 1:100)}
bind_rows(lapply(pr_all$bias_red, foo), .id = "QaS") %>% 
  group_by(QaS, JaS) %>% 
  summarise(ratio = median(ratio, na.rm = T),.groups = "drop") %>%
  group_by(QaS) %>% 
  summarise(Ratio_median = median(ratio),
            Ratio_l = quantile(ratio, probs = 0.025),
            Ratio_u = quantile(ratio, probs = 0.975)) %>% 
  left_join(.,SR, by = "QaS") %>% 
  ggplot(aes(y = Ratio_median, ymin = Ratio_l, ymax = Ratio_u, 
             x = SR, group = 1))+
  theme_bw()+
    geom_line()+
    geom_errorbar()+
    geom_point()+
  labs(y = "Reduction in MAB between\ndirect and s1 estimates")
ggsave(paste0(loc_plots, "/Reduction_MAB.png"), width = 10, height = 8.35)

# Increase in Stage 1 sampling variance ----
bind_rows(pr_all$samp_var_sm, .id = "QaS") %>% 
  group_by(QaS, JaS) %>% 
  summarise(ratio = median(ratio, na.rm = T),.groups = "drop") %>% 
  group_by(QaS) %>%
  summarise(Ratio_median = median(ratio, na.rm = T),
            Ratio_l = quantile(ratio, probs = 0.025, na.rm = T),
            Ratio_u = quantile(ratio, probs = 0.975, na.rm = T)) %>% 
  left_join(.,SR, by = "QaS") %>% 
  left_join(.,ACL, by = "QaS") %>%
  pivot_longer(-c(QaS, Ratio_median, Ratio_l, Ratio_u)) %>% 
  ggplot(aes(y = Ratio_median, ymin = Ratio_l, ymax = Ratio_u, 
             x = value, group = 1))+
  geom_line()+
  geom_errorbar()+
  geom_point()+
  facet_grid(name~.)+
  labs(y = "Increase in sampling variance\nS1 vs direct")
ggsave(paste0(loc_plots, "/Sampling_var.png"), width = 10, height = 8.35)

# HDI size ----
bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  group_by(QaS, rep_counter) %>% 
  summarise(cisize = median(s2LN_ci_size, na.rm = T),.groups = "drop") %>% 
  group_by(QaS) %>%
  summarise(cisize_median = median(cisize, na.rm = T),
            cisize_l = quantile(cisize, probs = 0.25, na.rm = T),
            cisize_u = quantile(cisize, probs = 0.75, na.rm = T)) %>% 
  left_join(.,SR, by = "QaS") %>% 
  left_join(.,ACL, by = "QaS") %>%
  pivot_longer(-c(QaS, cisize_median, cisize_l, cisize_u)) %>% 
  ggplot(aes(y = cisize_median, ymin = cisize_l, ymax = cisize_u, 
             x = value, group = 1))+theme_bw()+
  geom_line()+
  geom_errorbar()+
  geom_point()+
    facet_wrap(name~., scales = "free")+
  labs(y = "Size of HDI")

## ---- Frequentist Plots ---- ## -------------------------------------------
plot_df_freq <- bind_rows(sim_list$fp_global, .id = "QaS") %>%
  filter(model == "s2LN") %>% 
  left_join(.,SR, by = "QaS") %>% 
  left_join(.,ACL, by = "QaS") %>% 
  dplyr::select(SR, ACL, QaS, missing, MSE, Bias, Variance) %>% 
  pivot_longer(-c(SR, ACL, QaS, missing))

# MSE, Variance and bias
(plot_df_freq %>% 
  ggplot(aes(y = value, x = SR, group = missing, col = missing))+
  theme_bw()+
  geom_line()+
  geom_point()+
  facet_grid(name~., scales = "free_y"))+
(plot_df_freq %>% 
   ggplot(aes(y = value, x = ACL, group = missing, col = missing))+
   theme_bw()+
   geom_line()+
   geom_point()+
   facet_grid(name~., scales = "free_y"))



