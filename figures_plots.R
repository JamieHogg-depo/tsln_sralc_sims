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
library(quantreg)
export = TRUE

# Date and working directory
cur_date <- "20230411_0710"
wd <- "C:/r_proj/tsln_sralc_sims/"

# Load data
pr_all <- readRDS(paste0(wd, "results/", cur_date, "/pr_all.rds"))
sim_list <- readRDS(paste0(wd, "results/", cur_date, "/sim_list.rds"))

# Create plot folders
loc_plots <- paste0(wd, "results/", cur_date, "/plots")
if(!file.exists(loc_plots))dir.create(loc_plots)

## ---- Functions ---- ## ------------------------------------------------------

getQuantilePlot <- function(outcome, dep, df){
  
  f <- as.formula(paste0(outcome, " ~ ", dep, " + I(", dep, "^2)"))
  qfit <- rq(f, data=df, tau=c(0.05, 0.2, 0.5, 0.8, 0.95))
  newdata <- data.frame(sm = seq(0,1,0.001)) %>% setNames(c(dep))
  pr <- predict(qfit, newdata = newdata) %>% as.data.frame()
  nd <- cbind(sm = seq(0,1,0.001), pr) %>% 
    setNames(c("sm", "bottom", "lower", "median", "upper", "top"))
  
  # return the plot
  in_df2 <- df
  names(in_df2)[names(in_df2) == outcome] <- "met"
  names(in_df2)[names(in_df2) == dep] <- "sm"
  ggplot()+theme_bw()+
    geom_point(data = in_df2, aes(y = met, x = sm), col = "grey")+
    geom_line(data = nd, aes(y = lower, x = sm), col = "tomato2")+
    geom_line(data = nd, aes(y = upper, x = sm), col = "tomato2")+
    geom_line(data = nd, aes(y = top, x = sm), col = "sienna1")+
    geom_line(data = nd, aes(y = bottom, x = sm), col = "sienna1")+
    geom_line(data = nd, aes(y = median, x = sm), col = "red")+
    labs(y = outcome, x = dep)
}

jsave <- function(filename, square = T, square_size = 5000, ratio = c(6,9)){
  if(square){
    ggsave(filename = filename,
           path = loc_plots,
           dpi = 1000,
           width = square_size,
           height = square_size,
           scale = 1,
           units = "px")
  }else{
    total = square_size^2
    a <- sqrt((total*ratio[1])/ratio[2])
    b <- (ratio[2]*a)/ratio[1]
    ggsave(filename = filename,
           path = loc_plots,
           dpi = 1000,
           width = round(b),
           height = round(a),
           scale = 1,
           units = "px")
  }
}

make_numeric_decimal <- function(.data){
  df <- .data
  cols_to_format <- unlist(lapply(df, is.numeric))
  df[,cols_to_format] <- bind_cols(lapply(df[,cols_to_format], sprintf, fmt = '%#.2f'))
  return(df)
}

## ---- Summarise objects ---- ## ----------------------------------------------

options(scipen = 999)

# Concordance
concor <- data.frame(QaS = as.character(1:13),
                     sigma_res = c(0.01, 0.1, 0.25, 0.5, 0.75, 1,
                                   1.25, 1.5, 1.75, 2, 2.5, 3, 3.5))

# SR - median across 100 reps ####
SR_i <- bind_rows(lapply(pr_all$SR, 
                       FUN = function(x){cbind(x, rep_counter = 1:100)}), 
                .id = "QaS")
SR <- SR_i %>% 
  group_by(QaS) %>% 
  summarise(SR = median(SR))

# ALC - median across 100 reps ####
ALC_i <- bind_rows(lapply(pr_all$ALC, 
                       FUN = function(x){data.frame(ALC = x, rep_counter = 1:100)}), 
                .id = "QaS") 
ALC <- ALC_i %>% 
  group_by(QaS) %>% 
  summarise(ALC = median(ALC))

# Coverage - QaS by missing ####
coverage_i <- bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  mutate(in_ci = ifelse(prop > s2LN.lower & prop < s2LN.upper, 1, 0)) %>% 
  group_by(QaS, rep_counter, missing) %>% 
  summarise(in_ci = mean(in_ci), .groups = "drop") 
coverage <- coverage_i %>% 
  group_by(QaS, missing) %>% 
  summarise(coverage = median(in_ci), .groups = "drop")

# All reps and scenarios
spm_global <- bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  mutate(in_ci = ifelse(prop > s2LN.lower & prop < s2LN.upper, 1, 0)) %>% 
  group_by(QaS, rep_counter) %>%
  summarise(MRRMSE = mean(s2LN_RRMSE), 
            MARB = mean(s2LN_ARB), 
            Coverage = mean(in_ci), 
            cisize = mean(s2LN_ci_size),
            .groups = "drop")
in_df <- spm_global %>% 
  left_join(.,dplyr::select(SR_i, QaS, rep_counter, SR), by = c("QaS", "rep_counter")) %>% 
  left_join(.,dplyr::select(ALC_i, QaS, rep_counter, ALC), by = c("QaS", "rep_counter")) %>% 
  left_join(.,concor, by = "QaS")

## ---- Quantile regression to individual data ---- ##

qfit <- rq(Coverage ~ SR + I(SR^2), data=in_df, tau=c(0.5))
predict(qfit, newdata = data.frame(SR = c(0.1,0.5)))

# Specific table
qfit <- rq(MRRMSE ~ ALC + I(ALC^2), data=in_df, tau=c(0.5,0.80,0.95))
newdata <- data.frame(ALC = c(0,0.4,0.5, 0.6))
pp <- t(predict(qfit, newdata = newdata)) %>% as.data.frame() %>% 
  setNames(c("ALC_0", "ALC_0.4", "ALC_0.5", "ALC_0.6")) %>% 
  mutate(Percentile = paste0(c(50, 80, 95), "th")) %>% relocate(Percentile)
rownames(pp) = NULL
pp$ratio0 <- pp$ALC_0/pp$ALC_0
pp$ratio1 <- pp$ALC_0/pp$ALC_0.4
pp$ratio2 <- pp$ALC_0/pp$ALC_0.5
pp$ratio3 <- pp$ALC_0/pp$ALC_0.6
pp %>% make_numeric_decimal() %>% 
  mutate(ALC_0 = paste0(ALC_0, " (", ratio0, ")"),
         ALC_0.4 = paste0(ALC_0.4, " (", ratio1, ")"),
         ALC_0.5 = paste0(ALC_0.5, " (", ratio2, ")"),
         ALC_0.6 = paste0(ALC_0.6, " (", ratio3, ")")) %>% 
  dplyr::select(Percentile, contains("ALC_")) %>% 
  setNames(c("Percentile", str_replace(names(.)[-1], "_", " = "))) %>% 
  kable(format = "latex")

# RRMSE
getQuantilePlot("MRRMSE", "SR", in_df) / getQuantilePlot("MRRMSE", "ALC", in_df)
if(export) jsave("MRRMSE_quantreg.png")

# ARB
(getQuantilePlot("MARB", "SR", in_df)+
    ylim(0,0.8))/
(getQuantilePlot("MARB", "ALC", in_df)+
   ylim(0,0.8))
if(export) jsave("MARB_quantreg.png")

# Coverage and SR
(getQuantilePlot("Coverage", "SR", in_df)+
    ylim(0.25,1)+
    geom_hline(yintercept = 0.95))/
(getQuantilePlot("Coverage", "ALC", in_df)+
   ylim(0.25,1)+
   geom_hline(yintercept = 0.95))
if(export) jsave("Coverage_quantreg.png")

# HDI size
(getQuantilePlot("cisize", "SR", in_df)+labs(y = "Mean width of 95% HDI"))/
  (getQuantilePlot("cisize", "ALC", in_df)+labs(y = "Mean width of 95% HDI"))
if(export) jsave("HDIsize_quantreg.png")

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
  left_join(.,ALC, by = "QaS") %>% 
  left_join(.,coverage, by = c("QaS", "missing")) %>% 
  left_join(.,concor, by = "QaS")

# sigma_res vs SR, ALC ----
(plot_df %>% 
  ggplot(aes(y = SR, x = sigma_res, group = 1))+
  theme_bw()+
  geom_line()+
  geom_point())+
(plot_df %>% 
  ggplot(aes(y = ALC, x = sigma_res, group = 1))+
  theme_bw()+
  geom_line()+
  geom_point())

# SR vs ALC ----
in_df %>% 
  ggplot(aes(x = SR, y = ALC))+
  theme_bw()+
  geom_point()+
  geom_abline()
if(export) ggsave(paste0(loc_plots, "/SR vs ALC.png"), width = 10, height = 8.35)

# MRRMSE and MARB vs ALC ----
(plot_df %>% 
  ggplot(aes(y = MRRMSE, ymin = MRRMSE_l, ymax = MRRMSE_u, 
             x = ALC, col = missing, group = missing))+
  theme_bw()+
  geom_errorbar()+
  geom_line()+
  geom_point()+
    theme(legend.position = "bottom"))+
(plot_df %>% 
   ggplot(aes(y = MARB, ymin = MARB_l, ymax = MARB_u, 
              x = ALC, col = missing, group = missing))+
   theme_bw()+
   geom_errorbar()+
   geom_line()+
   geom_point()+
   theme(legend.position = "none"))
if(export) ggsave(paste0(loc_plots, "/MRRMSE and MARB_ALC.png"), width = 10, height = 8.35)

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
if(export) ggsave(paste0(loc_plots, "/MRRMSE and MARB_SR.png"), width = 10, height = 8.35)

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
              x = ALC, col = missing, group = missing))+
   theme_bw()+
   geom_hline(yintercept = 0.95)+
   geom_line()+
   geom_point()+
   theme(legend.position = "none"))
if(export) ggsave(paste0(loc_plots, "/Coverage.png"), width = 10, height = 8.35)

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
if(export) ggsave(paste0(loc_plots, "/Reduction_MAB.png"), width = 10, height = 8.35)

# Increase in Stage 1 sampling variance ----
bind_rows(pr_all$samp_var_sm, .id = "QaS") %>% 
  group_by(QaS, JaS) %>% 
  summarise(ratio = median(ratio, na.rm = T),.groups = "drop") %>% 
  group_by(QaS) %>%
  summarise(Ratio_median = median(ratio, na.rm = T),
            Ratio_l = quantile(ratio, probs = 0.025, na.rm = T),
            Ratio_u = quantile(ratio, probs = 0.975, na.rm = T)) %>% 
  left_join(.,SR, by = "QaS") %>% 
  left_join(.,ALC, by = "QaS") %>%
  pivot_longer(-c(QaS, Ratio_median, Ratio_l, Ratio_u)) %>% 
  ggplot(aes(y = Ratio_median, ymin = Ratio_l, ymax = Ratio_u, 
             x = value, group = 1))+
  geom_line()+
  geom_errorbar()+
  geom_point()+
  facet_grid(name~.)+
  labs(y = "Increase in sampling variance\nS1 vs direct")
if(export) ggsave(paste0(loc_plots, "/Sampling_var.png"), width = 10, height = 8.35)

# HDI size ----
bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  group_by(QaS, rep_counter) %>% 
  summarise(cisize = median(s2LN_ci_size, na.rm = T),.groups = "drop") %>% 
  group_by(QaS) %>%
  summarise(cisize_median = median(cisize, na.rm = T),
            cisize_l = quantile(cisize, probs = 0.25, na.rm = T),
            cisize_u = quantile(cisize, probs = 0.75, na.rm = T)) %>% 
  left_join(.,SR, by = "QaS") %>% 
  left_join(.,ALC, by = "QaS") %>%
  pivot_longer(-c(QaS, cisize_median, cisize_l, cisize_u)) %>% 
  ggplot(aes(y = cisize_median, ymin = cisize_l, ymax = cisize_u, 
             x = value, group = 1))+theme_bw()+
  geom_line()+
  geom_errorbar()+
  geom_point()+
    facet_wrap(name~., scales = "free")+
  labs(y = "Size of HDI")
if(export) ggsave(paste0(loc_plots, "/HDI_size.png"), width = 10, height = 8.35)

## ---- Tables ---- ## -------------------------------------------

in_df %>% 
  group_by(sigma_res) %>% 
  summarise(MRRMSE = median(MRRMSE),
            MARB = median(MARB),
            SR = median(SR),
            ALC = median(ALC),
            coverage = median(Coverage)) %>% 
  arrange(SR) %>% 
  ggplot(aes(y = MRRMSE, x = ALC))+
  geom_line()



