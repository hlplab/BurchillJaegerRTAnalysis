# 'old' version
library(dplyr)
library(lme4)
library(zplyr)

path <- "~/Box Sync/Power simulations for RTs/data/"
unix_path <- "~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/data/"
path_on_server <- "/u/zburchil/workspace/launchpad/"

build_setal_fetal_files <- TRUE

if (build_setal_fetal_files==TRUE) {
  
  fetal_file <- "Fine-Exp2-after RT-based exclusions.RDS"
  # fetal_file <- "Fine-Exp2-before any exclusions.RDS"
  setal_file <- "Stack-Exp2-OSF-060818-after bad subject & items exclusions + RT-based exclusions.RDS"
  gamm_predictions_file <- "../PowerSimulations/gamm_predictions_all_data_right.RDS"
  lmer_predictions_file <- "../PowerSimulations/RT_lmer_predictions_setalstyle.RDS"
  
  gamm_predictions <- readRDS(paste0(path, gamm_predictions_file)) %>%
    filter(Structure=="Filler") %>%
    select(Subject, Trial.Order, Experiment_Group, Word.Order, Word, predictedLogRT) %>%
    mutate(Subject = stringr::str_match(as.character(Subject), "[0-9]+")[,1])
  
  lmer_predictions <- readRDS(paste0(path, lmer_predictions_file)) %>%
    filter(Structure=="Filler") %>%
    select(Subject, Trial.Order, Experiment_Group, Word.Order, Word, predictedRegRT)
  
  fetal <- readRDS(paste0(path, fetal_file)) %>%
    filter(Experiment == "Experiment 2") %>%
    filter(Structure == "Filler") %>%
    filter(Trial.Order > 0) %>%
    filter(!(Subject == "53" & Trial.Order == "35")) %>%
    filter(Subject != "23") %>% 
    zplyr::left_join(
      gamm_predictions %>% filter(Experiment_Group=="FETAL"),
      by = c("Subject", "Trial.Order", "Word.Order")
    ) %>%
    zplyr::left_join(
      lmer_predictions %>% filter(Experiment_Group=="FETAL"),
      by = c("Subject", "Trial.Order", "Word.Order")
    ) %>%
    select(Experimenters, Experiment, Group, Subject, Item, Block, Original.Condition, Structure,
           Trial.Order, Item.Order, Word, Word.Order, Word.Length, RT, predictedLogRT, predictedRegRT)
  
  setal <- readRDS(paste0(path, setal_file)) %>%
    filter(Trial.Order > 0) %>%
    filter(Structure == "Filler") %>%
    zplyr::left_join(
      gamm_predictions %>% filter(Experiment_Group=="SETAL"),
      by = c("Subject", "Trial.Order", "Word.Order")
    ) %>%
    zplyr::left_join(
      lmer_predictions %>% filter(Experiment_Group=="SETAL"),
      by = c("Subject", "Trial.Order", "Word.Order")
    ) %>%
    select(Experimenters, Experiment, Group, Subject, Item, Block, Original.Condition, Structure,
           Trial.Order, Item.Order, Word, Word.Order, Word.Length, RT, predictedLogRT, predictedRegRT)
  
  setal %>% filter(is.na(RT) | is.na(predictedRegRT) | is.na(predictedLogRT))
  fetal %>% filter(is.na(RT) | is.na(predictedRegRT) | is.na(predictedLogRT))
  
  saveRDS(setal, paste0(path, "../PowerSimulations/setal_sim_df.RDS"))
  saveRDS(fetal, paste0(path, "../PowerSimulations/fetal_sim_df.RDS"))
  
  # Move it to the server
  system(
    paste0(
      "rsync -avz ", paste0(unix_path, "../PowerSimulations/setal_sim_df.RDS"), " zburchil@cycle1.cs.rochester.edu:", path_on_server))
  system(
    paste0(
      "rsync -avz ", paste0(unix_path, "../PowerSimulations/fetal_sim_df.RDS"), " zburchil@cycle1.cs.rochester.edu:", path_on_server))
  
}



# # 'new' version:
# library(dplyr)
# library(lme4)
# library(zplyr)
# library(future)
# 
# path <- "~/Box Sync/Power simulations for RTs/data/"
# unix_path <- "~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/data/"
# path_on_server <- "/u/zburchil/workspace/launchpad/"
# 
# output_setal_file <- "../PowerSimulations/setal_sim_df_sfetal_only_with_items.RDS"
# output_fetal_file <- "../PowerSimulations/fetal_sim_df_sfetal_only_with_items.RDS"
# on_server_setal_name <- "setal_sim_df_sfetal_only_with_items.RDS"
# on_server_setal_reduced_name <- "setal_sim_df_sfetal_only_with_items_new_blocks.RDS"
# 
# build_setal_fetal_files <- TRUE
# make_extra_setal_file_on_server <- TRUE
# 
# if (build_setal_fetal_files==TRUE) {
# 
#   fetal_file <- "Fine-Exp2-after RT-based exclusions.RDS"
#   setal_file <- "Stack-Exp2-OSF-060818-after bad subject & items exclusions + RT-based exclusions.RDS"
#   # gamm_predictions_file <- "../PowerSimulations/gamm_predictions_all_data_right.RDS" -- old
#   # from "~/workspace/gamm_model_sfetal_onlynode92.RDS"
#   gamm_predictions_file <- "../PowerSimulations/data/prediction_data/all_predictions_gamm_sfetal_only_with_items.RDS"
#   #  lmer_predictions_file <- "../PowerSimulations/RT_lmer_predictions_setalstyle.RDS" -- old
#   lmer_predictions_file <- "../PowerSimulations/data/prediction_data/actually_all_predictions_lmer_sfetal_only_subjectslopes_only.RDS" 
# 
#   gamm_predictions <- readRDS(paste0(path, gamm_predictions_file)) %>%
#     filter(Structure=="Filler") %>%
#     select(Subject, Trial.Order, Experiment_Group, Word.Order, Word, predictedLogRT) %>%
#     mutate(Subject = stringr::str_match(as.character(Subject), "[0-9]+")[,1])
# 
#   lmer_predictions <- readRDS(paste0(path, lmer_predictions_file)) %>%
#     filter(Structure == "Filler") %>%
#     mutate(Experiment_Group = Experiment) %>%
#     select(Subject, Trial.Order, Experiment_Group, Word.Order, Word, predictedRegRT)
# 
#   fetal <- readRDS(paste0(path, fetal_file)) %>%
#     filter(Experiment == "Experiment 2") %>%
#     filter(Structure == "Filler") %>%
#     filter(Trial.Order > 0) %>%
#     filter(!(Subject == "53" & Trial.Order == "35")) %>%
#     filter(!(Subject == "53" & Item == "53")) %>% # NAs here!
#     zplyr::left_join(
#       gamm_predictions %>% filter(Experiment_Group=="FETAL"),
#       by = c("Subject", "Trial.Order", "Word.Order")
#     ) %>%
#     zplyr::left_join(
#       lmer_predictions %>% filter(Experiment_Group=="FETAL"),
#       by = c("Subject", "Trial.Order", "Word.Order")
#     ) %>%
#     select(Experimenters, Experiment, Group, Subject, Item, Block, Original.Condition, Structure,
#            Trial.Order, Item.Order, Word, Word.Order, Word.Length, RT, predictedLogRT, predictedRegRT)
# 
#   setal <- readRDS(paste0(path, setal_file)) %>%
#     filter(Trial.Order > 0) %>%
#     filter(Structure == "Filler") %>%
#     zplyr::left_join(
#       gamm_predictions %>% filter(Experiment_Group=="SETAL"),
#       by = c("Subject", "Trial.Order", "Word.Order")
#     ) %>%
#     zplyr::left_join(
#       lmer_predictions %>% filter(Experiment_Group=="SETAL"),
#       by = c("Subject", "Trial.Order", "Word.Order")
#     ) %>%
#     select(Experimenters, Experiment, Group, Subject, Item, Block, Original.Condition, Structure,
#            Trial.Order, Item.Order, Word, Word.Order, Word.Length, RT, predictedLogRT, predictedRegRT)
# 
#   setal %>% filter(is.na(RT) | is.na(predictedRegRT) | is.na(predictedLogRT))
#   fetal %>% filter(is.na(RT) | is.na(predictedRegRT) | is.na(predictedLogRT))
# 
#   saveRDS(setal, paste0(path, output_setal_file))
#   saveRDS(fetal, paste0(path, output_fetal_file))
# 
#   # Move it to the server
#   system(
#     paste0(
#       "rsync -avz ", paste0(unix_path, output_setal_file), " zburchil@cycle1.cs.rochester.edu:", path_on_server))
#   system(
#     paste0(
#       "rsync -avz ", paste0(unix_path, output_fetal_file), " zburchil@cycle1.cs.rochester.edu:", path_on_server))
# 
#   if (make_extra_setal_file_on_server==TRUE) {
#     plan(remote, workers="zburchil@cycle1.cs.rochester.edu")
#     res <- future::future({
#       setal <- readRDS(paste0(path_on_server, on_server_setal_name))
#       setal %>% 
#         filter(Structure=="Filler") %>%
#         mutate(NewBlock = dense_rank(Trial.Order)) %>%
#         mutate(NewBlock = ifelse(NewBlock < 17, 1L, 
#                           ifelse(NewBlock < 37, 2L, 
#                           ifelse(NewBlock < 52, 3L, 99L)))) %>% 
#         filter(NewBlock < 4L) %>% # Don't want no blocks beyond that
#         filter(!(NewBlock == 1 & Group == "RC-first")) %>% # Exclude so that Block 1 only has filler-first peeps
#         mutate(Block = NewBlock) %>%
#         saveRDS(paste0(path_on_server, on_server_setal_reduced_name))
#       "done!"
#     })
#     if (cs:::wait_and_check(resolved(res), 20) == FALSE) {
#       stop("Timed out! Check code!!!!!!!!!!")
#     } else {
#       value(res)
#     }
#   }
# }
