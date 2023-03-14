## Create list

# define grand list object - overall performance
# each element is a performance metric
# element 1 of each element is the first scenario performance metrics 
sim_list <- list(scenario = list(),
                 # FP: per area
                 fp_pa = list(),# FP: frequentist point estimates
                 # FP: global
                 fp_global = list(), 
                 # SPM: global measures
                 spm_global = list(), # SPM: summarize_point-estimates_median
                 # SPM: per area
                 spm_pa = list(),
                 # CD: global measures
                 cd_global = list(), # CD: combine draws
                 # CD: per area 
                 cd_pa = list(),
                 # median model coefficients
                 mpl = list(),
                 # median number of sampled areas
                 mbar = list(),
                 # median of median area sampled sizes
                 nbarbar = list(),
                 # Convergence checks
                 Rhat_okay = list(),
                 Rhat_max = list(),
                 num_divs_max = list(),
                 prop_conv = list(),
                 total_reps = list(),
                 # how many do the models take
                 median_duration = list())

# initial update
sim_list$scenario[[QaS]] <- data.frame(sigma_res = sigma_res)

# Define list to store results from each rep
# pr: per rep
# 14 elements where each element will have `no_reps` number of elements
pr <-  list(m = list(),               # nu. of sampled areas
            nbar = list(),            # median sampled size of areas
            SPL = list(),             # simulation parameters
            s1LN_c = list(),        	# ALL parameter summary matrix from s1LN model
            s2LN_c = list(),       	# ALL parameter summary matrix from s2LN model
            spm_pa = list(),          # Summarize Point-estimate Median (spm) per area
            spm_global = list(),      # global measures for rep
            nom_ci = list(),          # nominal confidence interval coverage
            btess = list(),           # proportion of parameters with good bulk and tail ESS
            Rhat_okay = list(),       # proportion of parameters which Rhat<1.01
            Rhat_max = list(),        # max Rhat across ALL parameters
            Rhat_okay_ip = list(),    # proportion of IMPORTANT parameters which Rhat<1.01
            num_divs = list(),        # the number of divergence iterations in Stan
            true_prop = list(),       # current true_prop dataset
            model_duration = list(),  # how long for the model to fit  
            loo = list(),             # LOOCV for models
            loo_prop_bad_k = list(),  # LOOCV diagnostic
            SR = list(),				# Smoothing Ratio (SR)
            WOLSB = list(),			# Regression coefficient from WOLS of pd vs d
            WCOR = list(),			# Weighted pearson correlation between pd and d
            samp_var_sm = list(),		# dataset comparing the reduction in sampling variance under TSLN
            bias_red = list(),		# comparing the absolute bias between the HT_Direct and P_o from TSLN
            ratio_topten = list() )	# Median ratio of reduction in sampling variance for the largest HT_el_VAR values
pr_all <- pr

# Define list to collect posterior draws across all reps
# ar: all reps
ar <- list(s1LN = list(),
           s2LN = list())

# Define list to collect posterior medians of variance terms 
# and fixed coefficients across all reps
# mpl: model parameter list
mpl <- list(s1LN = list(),
            s2LN = list())

# Define convergence vectors
# each vector will be varying lengths based on how many times we believe they
# have converged during simulation
s1LN_convergence = as.numeric()
s2LN_convergence = as.numeric()