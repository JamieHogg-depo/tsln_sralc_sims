##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                        STEP 6: COLLAPSE REPETITIONS                      ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Collect c matrixes for each model
pr_all$s1LN_c[[QaS]] <- bind_rows(pr$s1LN_c, .id = "rep_counter")
pr_all$s2LN_c[[QaS]] <- bind_rows(pr$s2LN_c, .id = "rep_counter")

# Add median of posterior median of model coefficients
sim_list$mpl[[QaS]] <- suppressWarnings(bind_rows(mpl) %>% 
  group_by(i, .variable, model) %>% 
  summarise(median = median(.value),
            .groups = "drop") %>% 
  mutate(parameter = ifelse(is.na(i), .variable, paste0(.variable, "[", i, "]"))) %>% 
  dplyr::select(-c(i, .variable)) %>% 
  relocate(model, parameter) %>% 
  arrange(model, parameter))

# SPM: global measures
pr_all$spm_global[[QaS]] <- bind_rows(pr$spm_global)
sim_list$spm_global[[QaS]] <- bind_rows(pr$spm_global) %>% 
  dplyr::select(-rep_counter) %>% 
  group_by(model, missing) %>% 
  summarise_all(median, na.rm = TRUE) %>% 
  arrange(missing, model)

# SPM: variance of SPM RRMSE AND ARB across repetitions
pr_all$spm_pa[[QaS]] <- bind_rows(pr$spm_pa)
variances <- bind_rows(pr$spm_pa) %>%
  dplyr::select(-rep_counter) %>% 
  dplyr::select(area, contains("RRMSE"), contains("ARB")) %>% 
  pivot_longer(-area,
               names_to = c("model", "metric"),
               names_pattern = "(.*)_(.*)") %>% 
  group_by(area, model, metric) %>% 
  summarise(variance_reps = var(value), .groups = "drop") %>% 
  pivot_wider(values_from = variance_reps, 
              names_from = c(model, metric),
              names_glue = "{model}_{metric}_repvar")

# SPM: per area
sim_list$spm_pa[[QaS]] <- bind_rows(pr$spm_pa) %>%
  dplyr::select(-rep_counter) %>% 
  group_by(area) %>% 
  summarise_at(vars(-c(missing, ps_area)), median, na.rm = TRUE) %>% 
  left_join(.,(bind_rows(pr$spm_pa) %>%
                 dplyr::select(area, missing) %>% 
                 group_by(area) %>% 
                 summarise(missing = mean(as.numeric(missing)))), 
            by = "area") %>% 
  left_join(.,variances, by = "area") %>% 
  relocate(area, missing)

# CD: per area
 # Posterior median and HDI
	tt_s1LN <- bind_rows(ar$s1LN) %>% 
      dplyr::select(-prop) %>% 
      median_hdci(s1LN, .simple_names = F) %>% 
      dplyr::select(-.width, -.point, -.interval) %>% 
      rename(s1LN.median = s1LN) %>% 
      mutate(
        # Confidence interval size
        s1LN_ci_size = s1LN.upper - s1LN.lower,
      )
	tt_s2LN <- bind_rows(ar$s2LN) %>% 
      dplyr::select(-prop) %>% 
      median_hdci(s2LN, .simple_names = F) %>% 
      dplyr::select(-.width, -.point, -.interval) %>% 
      rename(s2LN.median = s2LN) %>%
      mutate(
        # Confidence interval size
        s2LN_ci_size = s2LN.upper - s2LN.lower,
      )
  # RMSE and ARB
	tt2_s1LN <- bind_rows(ar$s1LN)  %>% 
      tidyPERF(s1LN)
	tt2_s2LN <- bind_rows(ar$s2LN)  %>% 
      tidyPERF(s2LN)
    
  # true_prop
    agg_tp <- bind_rows(pr$true_prop) %>%
      group_by(area) %>% 
      summarise_at(vars(-c(missing, ps_area)), median, na.rm = TRUE) %>% 
      left_join(.,(bind_rows(pr$true_prop) %>%
                     dplyr::select(area, missing, ps_area) %>% 
                     group_by(area) %>% 
                     summarise(missing = mean(as.numeric(missing)))), 
                by = "area")
    
  # Combine all measures
  sim_list$cd_pa[[QaS]] <- list(tt_s2LN, tt2_s2LN,
                                tt_s1LN, tt2_s1LN,
                                agg_tp) %>% 
                                  reduce(full_join, by = "area")
  
  # CD: global measures
  sim_list$cd_global[[QaS]] <- sim_list$cd_pa[[QaS]] %>% 
    dplyr::select(contains("RRMSE"), contains("ARB")) %>% 
    summarise_all(mean, na.rm = T) %>% 
    mutate(x = 1) %>% 
    pivot_longer(-x,
                 names_to = c("model", "metric"),
                 names_pattern = "(.*)_(.*)") %>% 
    pivot_wider(values_from = value, 
                names_from = metric) %>% 
    dplyr::select(-x)

# Other objects
pr_all$m[[QaS]] <- unlist(pr$m)
sim_list$mbar[[QaS]] <- median(unlist(pr$m), na.rm = T)

pr_all$nbar[[QaS]] <- unlist(pr$nbar)
sim_list$nbarbar[[QaS]] <- median(unlist(pr$nbar), na.rm = T)

pr_all$Rhat_okay[[QaS]] <- bind_rows(pr$Rhat_okay, .id = "rep_counter")
sim_list$Rhat_okay[[QaS]] <- bind_rows(pr$Rhat_okay) %>% 
  summarise_all(median, na.rm = TRUE)
 
pr_all$Rhat_max[[QaS]] <- bind_rows(pr$Rhat_max, .id = "rep_counter")
sim_list$Rhat_max[[QaS]] <- bind_rows(pr$Rhat_max) %>% 
  summarise_all(median, na.rm = TRUE)

pr_all$num_divs[[QaS]] <- bind_rows(pr$num_divs, .id = "rep_counter")
sim_list$num_divs_max[[QaS]] <- bind_rows(pr$num_divs) %>% 
  summarise_all(max, na.rm = TRUE)
  
pr_all$Rhat_okay_ip[[QaS]] <- bind_rows(pr$Rhat_okay_ip, .id = "rep_counter")
pr_all$loo[[QaS]] <- bind_rows(pr$loo, .id = "rep_counter")
pr_all$loo_prop_bad_k[[QaS]] <- bind_rows(pr$loo_prop_bad_k, .id = "rep_counter")

pr_all$model_duration[[QaS]] <- bind_rows(pr$model_duration)
sim_list$median_duration[[QaS]] <- bind_rows(pr$model_duration) %>% 
  summarise_all(median, na.rm = TRUE)
 
# nominal confidence interval
pr_all$nom_ci[[QaS]] <- bind_rows(pr$nom_ci)

# effective sample size
pr_all$btess[[QaS]] <- bind_rows(pr$btess)

# smoothing ratio
pr_all$SR[[QaS]] <- bind_rows(pr$SR)
pr_all$ALC[[QaS]] <- unlist(pr$ALC)
pr_all$WOLSB[[QaS]] <- unlist(pr$WOLSB)
pr_all$WCOR[[QaS]] <- unlist(pr$WCOR)
pr_all$samp_var_sm[[QaS]] <- bind_rows(pr$samp_var_sm, .id = "JaS")
pr_all$bias_red[[QaS]] <- bind_rows(pr$bias_red)
pr_all$ratio_topten[[QaS]] <- bind_rows(pr$ratio_topten)

sim_list$prop_conv[[QaS]] <- data.frame(
  s1LN = mean(s1LN_convergence, na.rm = T),
  s2LN = mean(s2LN_convergence, na.rm = T)
)
sim_list$total_reps[[QaS]] <- data.frame(
  s1LN = length(s1LN_convergence),
  s2LN = length(s2LN_convergence)
)

# Get frequentist RRMSE and ARB for each area
sim_list$fp_pa[[QaS]] <- getfp(pr$spm_pa) %>% 
						left_join(., 
						          dplyr::select(ss_list$true_prop, area, ps_area, prop, area_ss)
						          , by = "area")

# Get Frequentist Mean RRMSE and ARB
sim_list$fp_global[[QaS]] <- getfp(pr$spm_pa) %>%
  pivot_longer(-c(area, missing, n_reps), 
               names_to = c("model", "metric"),
               names_pattern = "(.*)_(.*)") %>% 
  pivot_wider(values_from = value, 
              names_from = metric) %>% 
  group_by(missing, model) %>% 
  summarise(nu_points = sum(n_reps),
            fMRRMSE = mean(fRRMSE, na.rm = TRUE),
            fMARB = mean(fARB, na.rm = TRUE),
			      Bias = mean(biasminus, na.rm = TRUE),
			      Variance = mean(varInFrac, na.rm = TRUE),
            fRRMSE_repvar = var(fRRMSE, na.rm = TRUE),
            fARB_repvar = var(fARB, na.rm = TRUE),
			      .groups = "drop") %>% 
  mutate(MSE = Bias^2 + Variance)

