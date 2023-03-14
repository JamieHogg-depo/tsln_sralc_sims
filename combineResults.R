##
# Combine results
##

.libPaths("r_lib")
library(tidyverse)

dir = "tsln_sralc_sims/outputs"

## sim_list ## -----------------------------------------------------------------

which_list = "sim_list"
  
r_loc <- paste0(dir, "/", cur_date, "/r")
file_names <- list.files(r_loc)
cur_files <- file_names[which(str_detect(file_names, which_list))]

# load lists into R environment
for(i in 1:length(cur_files)){
  assign(cur_files[i],
         readRDS(paste0(r_loc, "/", cur_files[i])))
}

assign("main_list", eval(parse(text=cur_files[1])))
nu_elements <- length(main_list)
nu_scens <- length(cur_files)

for(QaS in 1:nu_scens){
  for(i in 1:nu_elements){
    what_to_eval <- paste0(cur_files[QaS], "[[", i, "]][[", QaS, "]]")
    main_list[[i]][[QaS]] <- eval(parse(text=what_to_eval))
  }
}

saveRDS(
  main_list, 
  file=paste0(dir, "/", cur_date, "/", which_list, ".rds")
)

rm(main_list, cur_files, which_list)

## pr_all ## -------------------------------------------------------------------

which_list = "pr_all"

r_loc <- paste0(dir, "/", cur_date, "/r")
file_names <- list.files(r_loc)
cur_files <- file_names[which(str_detect(file_names, which_list))]

# load lists into R environment
for(i in 1:length(cur_files)){
  assign(cur_files[i],
         readRDS(paste0(r_loc, "/", cur_files[i])))
}

assign("main_list", eval(parse(text=cur_files[1])))
nu_elements <- length(main_list)
elements_id <- which(!names(main_list) %in% c("SPL", "true_prop"))
nu_scens <- length(cur_files)

for(QaS in 1:nu_scens){
  for(i in elements_id){
    what_to_eval <- paste0(cur_files[QaS], "[[", i, "]][[", QaS, "]]")
    main_list[[i]][[QaS]] <- eval(parse(text=what_to_eval))
  }
}

saveRDS(
  main_list, 
  file=paste0(dir, "/", cur_date, "/", which_list, ".rds")
)