##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                        SCRIPT CONTROL: GET RUNNING                       ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start_time_full <- Sys.time()
library(tidyverse)
library(rstan)
library(tidybayes)
library(MASS)
library(MCMCvis)
library(loo)

wd <- "tsln_sralc_sims/r_src/"
source(paste0(wd, "functions.R"))
source(paste0(wd, "generateCENSUS.R"))
source(paste0(wd, "takeSAMPLE.R"))

## GRAND PARAMETERS ##
n.iter <- 2000
n.warm <- 1000 
n.ch <- 2

## ---- Define the sequence of residual errors ---- ##
sim_grid <- c(0.01, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5)
sigma_res <- sim_grid[QaS]

## ---- Compile the stan models ---- ## 
NORM_comp <- stan_model(file=paste0(wd, "modelcode/NORM_stan.stan"))
s1LN_comp <- stan_model(file=paste0(wd, "modelcode/s1LN_stan.stan"))

# start the timer
start_time_full <- Sys.time()

# Create lists
source(paste0(wd, "create_list.R"))

## STEP 1: Setup census ## 
GC <- generateCENSUS(seed = 47, 
                     alpha_survey = 2, #x1 - categorical, 
                     alpha_nonsurvey = 2, #x2, 
                     gamma = 0.01)

## ---- Start taking repeated sampled from census ---- ##
JaS <- 1
repeat{
  
message("Now starting rep ", JaS, " of ", no_reps)
  
## STEP 2: Take stratified sample ## 
ss_list <- takeSAMPLE(GC, input_seed = JaS)
pr$true_prop[[JaS]] <- ss_list$true_prop

# rm(GC, ss_list)
# GC <- generateCENSUS(seed = 47, 
#                      alpha_survey = 2, #8, #x1 - categorical, 
#                      alpha_nonsurvey = 2, #x2, 
#                      gamma = 0.01)
# ss_list <- takeSAMPLE(GC, input_seed = JaS)
# 
# fit <- glm(y ~ x1 + x2, family = binomial, ss_list$sample, 
#            weights = (w/sum(w)) * nrow(ss_list$sample))
# summary(fit)
# exp(coef(fit))
# with(ss_list$sample, getSR(y, predict(fit, type = "response"), w, area))
# with(ss_list$sample, getWOLS(y, predict(fit, type = "response"), w, area))

## STEP 3: Fit the models ##
source(paste0(wd, "s1LN_modelrun.R"))
message("Finished model: ", MC)

source(paste0(wd, "s2LN_modelrun.R"))
message("Finished model: ", MC)

## STEP 4: Add results to PER REP LIST ##
source(paste0(wd, "step4.R"))

## Stopping rule
if(JaS == no_reps){
  break
}

# Next JaS value
JaS <- JaS + 1

}

## STEP 6: Collapse repetitions ## 
source(paste0(wd, "step6.R"))

total_run_time <- as.numeric(Sys.time() - start_time_full, units = "mins")
message("Total run time was ", 
        round(total_run_time, 2), 
        " mins. Each repetition took ", round(total_run_time/no_reps, 2), 
        " mins on average.")

# Save output
QaS_f <- sprintf("%02d", QaS)
saveRDS(sim_list, file = paste0("tsln_sralc_sims/outputs/", cur_date, "/r/QaS", QaS_f, "_sim_list.rds"))
#saveRDS(ar, file = paste0("tsln_sralc_sims/outputs/", cur_date, "/r/QaS", QaS_f, "_ar.rds"))
saveRDS(pr, file = paste0("tsln_sralc_sims/outputs/", cur_date, "/r/QaS", QaS_f, "_pr.rds"))
saveRDS(pr_all, file = paste0("tsln_sralc_sims/outputs/", cur_date, "/r/QaS", QaS_f, "_pr_all.rds"))








