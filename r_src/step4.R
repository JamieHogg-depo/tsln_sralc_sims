##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                            STEP 4: PER REP LIST                          ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Lowering of Bias by using TSLN model
pr$bias_red[[JaS]] <- data.frame(
			Direct = with(ss_list$true_prop, mean(abs(HT_Direct - prop), na.rm = T)),
			TSLN = (s1LN_fit %>% 
							spread_draws(mu_pd[ps_area]) %>% 
							median_hdci() %>% 
							left_join(.,ss_list$true_prop, by = "ps_area") %>%
							summarise(TSLN = mean( abs(mu_pd - prop), na.rm = T) ) ) ) %>%
		mutate(ratio = TSLN/Direct)

# Add the median SR
pr$SR[[JaS]] <- s1LN_fit %>% 
					spread_draws(SR) %>% 
					median_hdci()
					
# Add WOLSB
	wols_data <- s1LN_fit %>% 
		spread_draws(theta_pd[ps_area]) %>% 
		median_hdci() %>% 
		left_join(.,ss_list$true_prop, by = "ps_area") %>% 
		filter(HT_stable) %>% 
		mutate(weights = (1/HT_el_VAR)*(nrow(.)/sum((1/HT_el_VAR))))
		
	# fit weighted OLS
	wols_fit <- lm(theta_pd ~ HT_el_Direct, data = wols_data, weights = weights)
	pr$WOLSB[[JaS]] <- unname(coef(wols_fit)[2])
	pr$WCOR[[JaS]] <- cov.wt(dplyr::select(wols_data, theta_pd, HT_el_Direct), wt = wols_data$weights, cor = T)$cor[1,2]
					
pr$samp_var_sm[[JaS]] <- s1LN_fit %>% 
							spread_draws(gamma_pd[ps_area]) %>% 
							median_hdci() %>% 
							left_join(.,ss_list$true_prop, by = "ps_area") %>% 
							# add sampling variance to posterior variance of theta_pd
							mutate( gamma = gamma_pd + apply(s1LN_its$theta_pd, 2, var) ) %>%
							dplyr::select(HT_el_VAR, gamma, HT_Direct, area_ss, area) %>% 
							mutate(ratio = gamma / HT_el_VAR) %>%
						relocate(area, area_ss, HT_Direct, HT_el_VAR, gamma, ratio)
							
# Reduction in sampling variance for direct sampling variances in the top 10%
pr$ratio_topten[[JaS]] <- pr$samp_var_sm[[JaS]] %>% 
							slice_max(HT_el_VAR, prop = 0.10) %>% 
							summarise(#median_ratio2 = median(ratio2),
									  median_ratio = median(ratio))

# Add the simulations parameters
pr$SPL[[JaS]] <- ss_list$SPL

# Add number of sampled areas
pr$m[[JaS]] <- ss_list$SPL$m

# Add median sample size
pr$nbar[[JaS]] <- median(ss_list$true_prop$area_ss, na.rm = T)

# Add divergence checks
pr$num_divs[[JaS]] <- data.frame(s1LN = s1LN_div,
                                 s2LN = s2LN_div)

# Add the Rhat_checks
pr$Rhat_okay[[JaS]] <- data.frame(s1LN = mean(as.numeric(ifelse(s1LN_c$Rhat<1.01, 1, 0)), na.rm = TRUE),
                                  s2LN = mean(as.numeric(ifelse(s2LN_c$Rhat<1.01, 1, 0)), na.rm = TRUE))

# Add the Rhat_checks
pr$Rhat_okay_ip[[JaS]] <- data.frame(s1LN = mean(as.numeric(ifelse(s1LN_c_ip$Rhat<1.01, 1, 0)), na.rm = TRUE),
                                     s2LN = mean(as.numeric(ifelse(s2LN_c_ip$Rhat<1.01, 1, 0)), na.rm = TRUE))

# Add the Rhat_max
pr$Rhat_max[[JaS]] <- data.frame(s2LN = max(s2LN_c$Rhat, na.rm = TRUE),
                                 s1LN = max(s1LN_c$Rhat, na.rm = TRUE))

# Add duration
pr$model_duration[[JaS]] <- data.frame(s1LN = s1LN_d, s2LN = s2LN_d)
								
# Add loocv estimate
pr$loo[[JaS]] <- data.frame(s1LN = s1LN_loo$estimates[1,1])
							
# check validity of loocv
pr$loo_prop_bad_k[[JaS]] <- 
data.frame(s1LN = round(mean(ifelse(s1LN_loo$diagnostics$pareto_k > 0.7, 1, 0))*100, 1) ) 

# Add effective sample size
pr$btess[[JaS]] <- cbind(
      s1LN_btess$prop_ESS %>% setNames(paste0("s1LN_", names(.))),
      s2LN_btess$prop_ESS %>% setNames(paste0("s2LN_", names(.)))
)

# Add _c objects
pr$s2LN_c[[JaS]] <- s2LN_c
pr$s1LN_c[[JaS]] <- s1LN_c

# Add MEASURES object
MEASURES <- list(
  ss_list$true_prop, 
  s1LN_measures,
  s2LN_measures) %>% 
  reduce(full_join, by = "ps_area") %>% 
  arrange(area) %>% 
  mutate(rep_counter = JaS)
pr$spm_pa[[JaS]] <- MEASURES

# Nominal COnfidence interval
pr$nom_ci[[JaS]] <- MEASURES %>% 
  dplyr::select(ends_with(c("median", "upper", "lower")), area, prop) %>% 
  pivot_longer(-c(area, prop)) %>% 
  separate(name, c("model", "metric")) %>% 
  pivot_wider(names_from = metric,
              values_from = value) %>% 
  mutate(in_ci = ifelse(prop > lower & prop < upper, 1, 0)) %>% 
  group_by(model) %>% 
  summarise(ci_coverage = mean(in_ci, na.rm = T)) %>% 
  pivot_wider(names_from = model, values_from = ci_coverage)

# Global measures
global <- MEASURES %>% 
  dplyr::select(missing, contains("RRMSE"), contains("ARB")) %>% 
  group_by(missing) %>% 
  summarise_all(mean, na.rm = T) %>% 
  pivot_longer(-missing,
               names_to = c("model", "metric"),
               names_pattern = "(.*)_(.*)") %>% 
  pivot_wider(values_from = value, 
              names_from = metric)

# Best Model performance for each area
arb <- MEASURES %>% 
  dplyr::select(area, missing, contains("ARB")) %>% 
  rename_with(~gsub("_ARB", "", .x)) %>%
  pivot_longer(-c(area, missing)) %>% 
  group_by(area, missing) %>% 
  arrange(area, value) %>% 
  mutate(V = ifelse(row_number()==1, 1, 0)) %>% 
  filter(V == 1) %>% 
  ungroup() %>%
  group_by(name, missing) %>%
  summarise(n = n(), .groups = "drop") %>% 
  group_by(missing) %>% 
  mutate(miss_n = sum(n)) %>% 
  ungroup() %>% 
  mutate(smallestARB = n/miss_n) %>% 
  dplyr::select(-n, -miss_n) %>% 
  rename(model = name)

rrmse <- MEASURES %>% 
  dplyr::select(area, missing, contains("RRMSE")) %>% 
  rename_with(~gsub("_RRMSE", "", .x)) %>%
  pivot_longer(-c(area, missing)) %>% 
  group_by(area, missing) %>% 
  arrange(area, value) %>% 
  mutate(V = ifelse(row_number()==1, 1, 0)) %>% 
  filter(V == 1) %>% 
  ungroup() %>%
  group_by(name, missing) %>%
  summarise(n = n(), .groups = "drop") %>% 
  group_by(missing) %>% 
  mutate(miss_n = sum(n)) %>% 
  ungroup() %>% 
  mutate(smallestRRMSE = n/miss_n) %>% 
  dplyr::select(-n, -miss_n) %>% 
  rename(model = name)

# create spm_global
pr$spm_global[[JaS]] <- list(arb, rrmse, global) %>% 
  reduce(full_join, by = c("missing", "model")) %>% 
  arrange(missing, model) %>% 
  mutate(rep_counter = JaS)
  
