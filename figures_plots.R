# SR
SR <- bind_rows(lapply(pr_all$SR, 
                       FUN = function(x){cbind(x, rep_counter = 1:50)}), 
                .id = "QaS") %>% 
  group_by(QaS) %>% 
  summarise(SR = median(SR))

# ACL
foo <- function(x){cbind(ACL = x, rep_counter = 1:50)}
ACL <- bind_rows(lapply(pr_all$WOLSB, 
                       FUN = function(x){data.frame(ACL = x, rep_counter = 1:50)}), 
                .id = "QaS") %>% 
  group_by(QaS) %>% 
  summarise(ACL = median(ACL))

# Coverage
Coverage <- bind_rows(lapply(pr_all$nom_ci, FUN = function(x){cbind(x, rep_counter = 1:50)[,-1]}), .id = "QaS") %>% 
  group_by(QaS) %>% 
  summarise(coverage = mean(s2LN))

# Plot
new_df <- bind_rows(pr_all$spm_global, .id = "QaS") %>% 
  filter(model == "s2LN") %>% 
  group_by(QaS, missing) %>% 
  summarise(MRRMSE = median(RRMSE),
            MRRMSE_l = quantile(RRMSE, probs = 0.25),
            MRRMSE_u = quantile(RRMSE, probs = 0.75),
            MARB = median(ARB), .groups = "drop") %>% 
  left_join(.,SR, by = "QaS") %>% 
  left_join(.,ACL, by = "QaS") %>% 
  left_join(.,Coverage, by = "QaS")

# MRRMSE
new_df %>% 
  ggplot(aes(y = MRRMSE, ymin = MRRMSE_l, ymax = MRRMSE_u, 
             x = ACL))+
  geom_errorbar()+
  geom_point()+
  facet_grid(.~missing)

# Coverage
new_df %>% 
  ggplot(aes(y = coverage, 
             x = SR))+
  geom_point()




