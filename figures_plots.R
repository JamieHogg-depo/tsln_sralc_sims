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
export = T

# Date and working directory
#cur_date <- "20230411_0710"
cur_date <- "20230412_0627"
wd <- "C:/r_proj/tsln_sralc_sims/"

# Load data
pr_all <- readRDS(paste0(wd, "results/", cur_date, "/pr_all.rds"))
sim_list <- readRDS(paste0(wd, "results/", cur_date, "/sim_list.rds"))
pr_all_nare <- readRDS(paste0(wd, "results/20230411_0710/pr_all.rds"))

# Create plot folders
loc_plots <- paste0(wd, "results/", cur_date, "/plots")
if(!file.exists(loc_plots))dir.create(loc_plots)

## ---- Functions ---- ## ------------------------------------------------------

getQuantilePlot <- function(outcome, dep, df){
  
  f <- as.formula(paste0(outcome, " ~ ", dep, " + I(", dep, "^2)"))
  #f <- as.formula(paste0(outcome, " ~ ", dep))
  qfit <- rq(f, data=df, tau=c(0.05, 0.2, 0.5, 0.8, 0.95))
  newdata <- data.frame(sm = seq(0,1,0.001)) %>% setNames(c(dep))
  pr <- predict(qfit, newdata = newdata) %>% as.data.frame()
  nd <- cbind(sm = seq(0,1,0.001), pr) %>% 
    setNames(c("sm", "bottom", "lower", "median", "upper", "top"))
  
  getIndexMin <- function(x){which(x == min(x))}
  
  low_value1 <- nd[getIndexMin(nd$upper),1]
  low_value2 <- nd[getIndexMin(nd$median),1]
  low_value3 <- nd[getIndexMin(nd$top),1]
  
  # return the plot
  in_df2 <- df
  names(in_df2)[names(in_df2) == outcome] <- "met"
  names(in_df2)[names(in_df2) == dep] <- "sm"
  
  if(outcome == "MRRMSE"){
  ggplot()+theme_bw()+
    geom_point(data = in_df2, aes(y = met, x = sm), col = "grey", shape = 1)+
    geom_line(data = nd, aes(y = lower, x = sm), col = "tomato2")+
    geom_line(data = nd, aes(y = upper, x = sm), col = "tomato2")+
    geom_line(data = nd, aes(y = top, x = sm), col = "sienna1")+
    geom_line(data = nd, aes(y = bottom, x = sm), col = "sienna1")+
    geom_line(data = nd, aes(y = median, x = sm), col = "red")+
    geom_vline(xintercept = low_value1, linetype = "dotted")+
    geom_vline(xintercept = low_value2, linetype = "dotted")+
    geom_vline(xintercept = low_value3, linetype = "dotted")+
    labs(y = outcome, x = dep)+
    scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  }else{
    ggplot()+theme_bw()+
      geom_point(data = in_df2, aes(y = met, x = sm), col = "grey", shape = 1)+
      geom_line(data = nd, aes(y = lower, x = sm), col = "tomato2")+
      geom_line(data = nd, aes(y = upper, x = sm), col = "tomato2")+
      geom_line(data = nd, aes(y = top, x = sm), col = "sienna1")+
      geom_line(data = nd, aes(y = bottom, x = sm), col = "sienna1")+
      geom_line(data = nd, aes(y = median, x = sm), col = "red")+
      #geom_vline(xintercept = low_value1, linetype = "dotted")+
      #geom_vline(xintercept = low_value2, linetype = "dotted")+
      #geom_vline(xintercept = low_value3, linetype = "dotted")+
      labs(y = outcome, x = dep)+
      scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  }
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
                .id = "QaS") %>% mutate(areaRE = 1)
SR_i_nare <- bind_rows(lapply(pr_all_nare$SR, 
                         FUN = function(x){cbind(x, rep_counter = 1:100)}), 
                  .id = "QaS") %>% mutate(areaRE = 0)
SR <- bind_rows(SR_i, SR_i_nare)

# ALC - median across 100 reps ####
ALC_i <- bind_rows(lapply(pr_all$ALC, 
                       FUN = function(x){data.frame(ALC = x, rep_counter = 1:100)}), 
                .id = "QaS") %>% mutate(areaRE = 1)
ALC_i_nare <- bind_rows(lapply(pr_all_nare$ALC, 
                          FUN = function(x){data.frame(ALC = x, rep_counter = 1:100)}), 
                   .id = "QaS") %>% mutate(areaRE = 0)
ALC <- bind_rows(ALC_i, ALC_i_nare)

# Coverage - QaS by missing ####
coverage_i <- bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  mutate(in_ci = ifelse(prop > s2LN.lower & prop < s2LN.upper, 1, 0)) %>% 
  group_by(QaS, rep_counter, missing) %>% 
  summarise(in_ci = mean(in_ci), .groups = "drop") %>% mutate(areaRE = 1)
coverage_i_nare <- bind_rows(pr_all_nare$spm_pa, .id = "QaS") %>% 
  mutate(in_ci = ifelse(prop > s2LN.lower & prop < s2LN.upper, 1, 0)) %>% 
  group_by(QaS, rep_counter, missing) %>% 
  summarise(in_ci = mean(in_ci), .groups = "drop") %>% mutate(areaRE = 0)
coverage <- bind_rows(coverage_i, coverage_i_nare)

# All reps and scenarios
spm_global_i <- bind_rows(pr_all$spm_pa, .id = "QaS") %>% 
  mutate(in_ci = ifelse(prop > s2LN.lower & prop < s2LN.upper, 1, 0)) %>% 
  group_by(QaS, rep_counter) %>%
  summarise(MRRMSE = mean(s2LN_RRMSE), 
            MARB = mean(s2LN_ARB), 
            Coverage = mean(in_ci), 
            cisize = mean(s2LN_ci_size),
            .groups = "drop") %>% mutate(areaRE = 1)
spm_global_nare <- bind_rows(pr_all_nare$spm_pa, .id = "QaS") %>% 
  mutate(in_ci = ifelse(prop > s2LN.lower & prop < s2LN.upper, 1, 0)) %>% 
  group_by(QaS, rep_counter) %>%
  summarise(MRRMSE = mean(s2LN_RRMSE), 
            MARB = mean(s2LN_ARB), 
            Coverage = mean(in_ci), 
            cisize = mean(s2LN_ci_size),
            .groups = "drop") %>% mutate(areaRE = 0)
in_df <- bind_rows(spm_global_i, spm_global_nare) %>% 
  left_join(.,dplyr::select(SR, QaS, rep_counter, SR, areaRE), by = c("QaS", "rep_counter", "areaRE")) %>% 
  left_join(.,dplyr::select(ALC, QaS, rep_counter, ALC, areaRE), by = c("QaS", "rep_counter", "areaRE")) %>% 
  left_join(.,concor, by = "QaS")

## ---- Quantile regression to individual data ---- ##

qfit <- rq(Coverage ~ SR + I(SR^2), data=in_df, tau=c(0.5))
predict(qfit, newdata = data.frame(SR = c(0.1,0.4)))

# Specific table
qfit <- rq(MRRMSE ~ ALC + I(ALC^2), data=in_df, tau=c(0.5,0.80,0.95))
newdata <- data.frame(ALC = c(0,0.4,0.5, 0.6, 0.7))
pp <- t(predict(qfit, newdata = newdata)) %>% as.data.frame() %>% 
  setNames(c("ALC_0", "ALC_0.4", "ALC_0.5", "ALC_0.6", "ALC_0.7")) %>% 
  mutate(Percentile = paste0(c(50, 80, 95), "th")) %>% relocate(Percentile)
rownames(pp) = NULL
pp$ratio0 <- pp$ALC_0/pp$ALC_0
pp$ratio1 <- pp$ALC_0/pp$ALC_0.4
pp$ratio2 <- pp$ALC_0/pp$ALC_0.5
pp$ratio3 <- pp$ALC_0/pp$ALC_0.6
pp$ratio4 <- pp$ALC_0/pp$ALC_0.7
pp %>% make_numeric_decimal() %>% 
  mutate(ALC_0 = paste0(ALC_0, " (", ratio0, ")"),
         ALC_0.4 = paste0(ALC_0.4, " (", ratio1, ")"),
         ALC_0.5 = paste0(ALC_0.5, " (", ratio2, ")"),
         ALC_0.6 = paste0(ALC_0.6, " (", ratio3, ")"),
         ALC_0.7 = paste0(ALC_0.7, " (", ratio4, ")")) %>% 
  dplyr::select(Percentile, contains("ALC_")) %>% 
  setNames(c("Percentile", str_replace(names(.)[-1], "_", " = "))) %>% 
  kable(format = "latex")

# Find global minimum
qfit <- rq(MRRMSE ~ SR + I(SR^2), data=in_df, tau=c(0.05, 0.2, 0.5, 0.8, 0.95))
newdata <- data.frame(SR = seq(0,1,0.001))
pr <- predict(qfit, newdata = newdata) %>% as.data.frame()
nd <- cbind(SR = seq(0,1,0.001), pr) %>% 
  setNames(c("SR", "bottom", "lower", "median", "upper", "top"))

# RRMSE
getQuantilePlot("MRRMSE", "SR", in_df) + getQuantilePlot("MRRMSE", "ALC", in_df)
if(export) jsave("MRRMSE_quantreg.png", square = F)

# ARB
(getQuantilePlot("MARB", "SR", in_df)+
    ylim(0,0.8))+
(getQuantilePlot("MARB", "ALC", in_df)+
   ylim(0,0.8))
if(export) jsave("MARB_quantreg.png", square = F)

# Coverage and SR
(getQuantilePlot("Coverage", "SR", in_df)+
    ylim(0.25,1)+
    geom_hline(yintercept = 0.95)+ 
    scale_x_continuous(position="top"))+
(getQuantilePlot("Coverage", "ALC", in_df)+
   ylim(0.25,1)+
   geom_hline(yintercept = 0.95)+ 
   scale_x_continuous(position="top"))
if(export) jsave("Coverage_quantreg.png", square = F)

# HDI size
(getQuantilePlot("cisize", "SR", in_df)+labs(y = "Mean width of 95% HDI"))+
(getQuantilePlot("cisize", "ALC", in_df)+labs(y = "Mean width of 95% HDI"))
if(export) jsave("HDIsize_quantreg.png", square = F)

# Compare SR and ALC
summary(lm(ALC ~ SR, in_df[in_df$areaRE == 0,]))
in_df %>% 
  filter(areaRE == 0) %>% 
  ggplot(aes(y=ALC, x = SR))+
  theme_bw()+
  geom_point(col = "grey")+
  geom_abline()+
  geom_smooth(method = "lm", col = "red")+
  ylim(0,1)+xlim(0,1)
if(export) jsave("ALCvsSR.png", square = F)


