library(mgcv)
library(dplyr)
library(future)
library(cs)
library(scam)
library(lme4)


load_sfetal <- function(
  path = "~/data/data/",
  fetal_file  =    "Fine-Exp2-after RT-based exclusions.RDS",
  setal_file =     "Stack-Exp2-OSF-060818-after bad subject & items exclusions + RT-based exclusions.RDS") {
  suppressWarnings(
    bind_rows(
      "SETAL" =     readRDS(paste0(path, setal_file)) %>% select(-List),
      "FETAL" =     readRDS(paste0(path, fetal_file)) %>%
        filter(Block!="Practice") %>% mutate(Block=as.numeric(as.character(Block))),
      .id = "Experiment_Group"
    )
  ) %>% filter(Trial.Order > 0)
}
load_setal <- function(
  path = "~/data/data/Stack-Exp2-OSF-060818-after bad subject & items exclusions + RT-based exclusions.RDS") {
  readRDS(path) %>% filter(Trial.Order > 0) %>%
    mutate(Experiment_Group = "SETAL")
}
load_fetal <- function(
  path = "~/data/data/Fine-Exp2-after RT-based exclusions.RDS") {
  readRDS(path) %>% filter(Trial.Order > 0, Block!="Practice") %>%
    mutate(Experiment_Group = "FETAL",
           Block=as.numeric(as.character(Block)))
}

process_factors <- function(df) {
  df %>%
    mutate(Experiment_Group =      factor(Experiment_Group),
           Item =                  factor(paste0(Item, "_", Experiment_Group)),
           Subject =               factor(paste0(Subject, "_", Experiment_Group)),
           ItemType = factor(ifelse(Structure == "Filler", "Filler", "Critical")))
}
process_data <- function(df) {
  df %>% 
    process_factors() %>%
    mutate(Trial.Order.c = scale(Trial.Order, scale=FALSE)[,1] %>% rlang::set_attrs(NULL),
           Word.Length.c = scale(Word.Length, scale=FALSE)[,1] %>% rlang::set_attrs(NULL),
           logRT = log10(RT))
}

get_predictions <- function(m, loaded_df, processed_df,
                            no_subj = T, no_item = T) {
  if (class(m)=="list")
    m <- m[["value"]]
  recenter_length = mean(processed_df$Word.Length)
  recenter_trial =  mean(processed_df$Trial.Order)
  if (no_subj == TRUE)
    dummy_subj = processed_df[processed_df$Group=="Filler-first",]$Subject[[1]]
  else
    dummy_subj = process_factors(loaded_df)$Subject
  if (no_item==TRUE)
    dummy_item = processed_df[processed_df$Group=="Filler-first",]$Item[[1]]
  else
    dummy_item = process_factors(loaded_df)$Item
  
  # Get the rest of the data ready for predictions
  all_data_for_pred <- loaded_df %>%
    mutate(Subject = dummy_subj,
           Item = dummy_item,
           logRT = log10(RT)) %>%
    mutate(Trial.Order.c = Trial.Order - recenter_trial,
           Word.Length.c = Word.Length - recenter_length)
  
  if ("gam" %in% class(m)) {
    pred_df <- predict(m, all_data_for_pred, type="terms")
    
    intercept <- attributes(pred_df)$constant[[1]]
    pred_df <- as.data.frame(pred_df) %>%
      mutate(Intercept = intercept)
    # Remove REs if necessary
    if (no_subj==TRUE)
      pred_df <- pred_df %>% select(-contains("Subject"))
    if (no_item==TRUE)
      pred_df <- pred_df %>% select(-contains("Item")) 
    pred_df <- pred_df %>%
      mutate(predictedLogRT = rowSums(.))
    
    big_df <- cbind(loaded_df, pred_df)
  } else {
    #assume lmerMod
    predicted_rawRT <- predict(m, all_data_for_pred) %>% 
      rlang::set_attrs(NULL)
    big_df <- loaded_df %>% mutate(predictedRegRT = predicted_rawRT)
  }
  big_df
}

############################################################################################################

plan(list(
  # multisession,
  tweak(remote, workers="zburchil@cycle1.cs.rochester.edu"),
  tweak(cluster, workers="node92")
))

setal_gamm %2% {
  loaded_data <- load_setal()
  
  processed_data <- loaded_data %>%
    filter(Structure == "Filler") %>%
    process_data()
  
  cl <- parallel::makeForkCluster(6)
  tictoc::tic()
  saveRDS("a","~/tmp/starting")
  m <- zplyr::collect_all(
    bam(formula = logRT ~ 1 + Word.Length.c +
          s(Trial.Order.c, k = 7) +
          s(Subject, bs = "re") +
          s(Subject, Trial.Order.c, bs = "re") +
          s(Subject, Word.Length.c, bs = "re"),
        # Since we want fillers-only across all blocks
        data = processed_data,
        cluster = cl),
    catchErrors = TRUE
  )
  
  saveRDS("a","~/tmp/ended")
  saveRDS(m[["value"]], paste0("~/workspace/gamm_setal_only_no_items_in_fit_linear_wl.RDS"))
  time_count <- tictoc::toc()
  
  pred_res <- zplyr::collect_all({
    big_df <- get_predictions(m[["value"]], loaded_data, 
                              processed_data, 
                              no_subj = FALSE)
    saveRDS(big_df, "~/workspace/predictions_gamm_setal_only_no_items_in_fit_linear_wl.RDS")
    "done"},
    catchErrors = TRUE)
  saveRDS(list(m[2:5], pred_res[2:5]), paste0("~/workspace/gamm_setal_only_no_items_in_fit_linear_wl_errors.RDS"))
}
done(setal_gamm)

#-------------------------------------------

plan(list(
  # multisession,
  tweak(remote, workers="zburchil@cycle2.cs.rochester.edu"),
  tweak(cluster, workers="node91")
))

fetal_gamm %2% {
  loaded_data <- load_fetal()
  
  processed_data <- loaded_data %>%
    filter(Structure == "Filler") %>%
    process_data()
  
  cl <- parallel::makeForkCluster(12)
  tictoc::tic()
  saveRDS("a","~/tmp/starting")
  m <- zplyr::collect_all(
    bam(formula = logRT ~ 1 + Word.Length.c +
          s(Trial.Order.c, k = 7) +
          s(Subject, bs = "re") +
          s(Subject, Trial.Order.c, bs = "re") +
          s(Subject, Word.Length.c, bs = "re"),
        # Since we want fillers-only across all blocks
        data = processed_data,
        cluster = cl),
    catchErrors = TRUE
  )
  
  saveRDS("a","~/tmp/ended")
  saveRDS(m[["value"]], paste0("~/workspace/gamm_fetal_only_no_items_in_fit_linear_wl.RDS"))
  time_count <- tictoc::toc()
  
  pred_res <- zplyr::collect_all({
    big_df <- get_predictions(m[["value"]], loaded_data, 
                              processed_data, 
                              no_subj = FALSE)
    saveRDS(big_df, "~/workspace/predictions_gamm_fetal_only_no_items_in_fit_linear_wl.RDS")
    "done"},
    catchErrors = TRUE)
  saveRDS(list(m[2:5], pred_res[2:5]), paste0("~/workspace/gamm_fetal_only_no_items_in_fit_linear_wl_errors.RDS"))
  
}
done(fetal_gamm)

# ---------------------------------------------------

plan(list(
  # multisession,
  tweak(remote, workers="zburchil@cycle1.cs.rochester.edu"),
  tweak(cluster, workers="node92")
))

fetal_lmer %2% {
  loaded_data <- load_fetal()
  
  processed_data <- loaded_data %>%
    filter(Structure == "Filler") %>%
    process_data()
  
  m <- zplyr::collect_all(
    lmer(RT ~ 1 + Word.Length.c + (1 + Word.Length.c | Subject),
         data = processed_data,
         control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
    catchErrors = TRUE
  )
  
  saveRDS(m[["value"]], paste0("~/workspace/lmer_fetal_only_fillers_only.RDS"))
  
  pred_res <- zplyr::collect_all({
    big_df <- get_predictions(m[["value"]], loaded_data, 
                              processed_data, 
                              no_subj = FALSE)
    saveRDS(big_df, "~/workspace/predictions_lmer_fetal_only_fillers_only.RDS")
    "done"},
    catchErrors = TRUE)
  saveRDS(list(m[2:5], pred_res[2:5]), paste0("~/workspace/lmer_fetal_only_fillers_only_errors.RDS"))
  "a"
  
}
done(fetal_lmer)

#-------------------------------------------

plan(list(
  # multisession,
  tweak(remote, workers="zburchil@cycle2.cs.rochester.edu"),
  tweak(cluster, workers="node91")
))

setal_lmer %2% {
  loaded_data <- load_setal()
  
  processed_data <- loaded_data %>%
    filter(Structure == "Filler") %>%
    process_data()
  
  m <- zplyr::collect_all(
    lmer(RT ~ 1 + Word.Length.c + (1 + Word.Length.c | Subject),
         data = processed_data,
         control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
    catchErrors = TRUE
  )
  
  saveRDS(m[["value"]], paste0("~/workspace/lmer_setal_only_fillers_only.RDS"))
  
  pred_res <- zplyr::collect_all({
    big_df <- get_predictions(m[["value"]], loaded_data, 
                              processed_data, 
                              no_subj = FALSE)
    saveRDS(big_df, "~/workspace/predictions_lmer_setal_only_fillers_only.RDS")
    "done"},
    catchErrors = TRUE)
  saveRDS(list(m[2:5], pred_res[2:5]), paste0("~/workspace/lmer_setal_only_fillers_only_errors.RDS"))
  
}
done(setal_lmer)












good_gamm %2% {
  loaded_data <- load_sfetal()
  
  processed_data <- loaded_data %>%
    filter(Structure == "Filler",
           Group == "Filler-first") %>%
    process_data()
  
  tictoc::tic()
  saveRDS("a","~/tmp/starting")
  m <- zplyr::collect_all(
    scam(formula = logRT ~ 1 +
           s(Trial.Order.c, k = 7, bs = "mpd") + s(Word.Length.c, k = 4) +
           s(Subject, bs = "re") +
           s(Subject, Trial.Order.c, bs = "re") +
           s(Subject, Word.Length.c, bs = "re") +
           s(Item, bs = "re") +
           s(Item, Word.Length.c, bs = "re"),
         # Since we want fillers-only across all blocks
         data = processed_data),
    catchErrors = TRUE
  )
  
  saveRDS("a","~/tmp/ended")
  saveRDS(m[["value"]], paste0("~/workspace/scam_sfetal_nonlinear.RDS"))
  saveRDS(m[2:5],  paste0("~/workspace/scam_sfetal_nonlinear_errors.RDS"))
  time_count <- tictoc::toc()
  
  big_df <- get_predictions(m, loaded_data, processed_data)
  saveRDS(big_df, "~/workspace/scam_sfetal_nonlinear_predictions.RDS")
  
}
done(good_gamm)
#---------------------------------------------------------------------------------------

plan(list(
  # multisession,
  tweak(remote, workers="zburchil@cycle1.cs.rochester.edu"),
  tweak(cluster, workers="node91")
))
# for the scam that gets the 'traditional' fit:
traditional_model %2% {
  loaded_data <- load_sfetal()
  
  processed_data <- loaded_data %>%
    filter(Structure == "Filler",
           Group == "Filler-first") %>%
    process_data()
  
  tictoc::tic()
  saveRDS("a","~/tmp/starting")
  m <- zplyr::collect_all(
    # not totally sure if we need a 'scam' here, but using it anyway
    scam(formula = RT ~ 1 + Word.Length.c +
           s(Subject, bs = "re") +
           s(Subject, Word.Length.c, bs = "re"),
         # Since we want fillers-only across all blocks
         data = processed_data),
    catchErrors = TRUE
  )
  
  saveRDS("a","~/tmp/ended")
  saveRDS(m[["value"]], paste0("~/workspace/scam_sfetal_linear.RDS"))
  saveRDS(m[2:5],  paste0("~/workspace/scam_sfetal_linear_errors.RDS"))
  time_count <- tictoc::toc()
  
  big_df <- get_predictions(m, loaded_data, processed_data)
  saveRDS(big_df, "~/workspace/scam_sfetal_linear_predictions.RDS")
}
done(traditional_model)



#---------------------------------------------------------------------------------------
plan(list(
  # multisession,
  tweak(remote, workers="zburchil@cycle1.cs.rochester.edu"),
  tweak(cluster, workers="node89")
))
group_tester_gam %2% {
  loaded_data <- load_sfetal()
  
  processed_data <- loaded_data %>%
    filter(Structure == "Filler") %>%
    process_data()
  
  tictoc::tic()
  saveRDS("a","~/tmp/starting")
  m <- zplyr::collect_all(
    scam(formula = logRT ~ 1 + s(Word.Length.c, k = 4) +
           s(Trial.Order.c, k = 7, bs = "mpd") +
           s(Trial.Order.c, k = 7, by = Group, bs = "mpd") +
           s(Subject, bs = "re") +
           s(Subject, Trial.Order.c, bs = "re") +
           s(Subject, Word.Length.c, bs = "re") +
           s(Item, bs = "re") +
           s(Item, Word.Length.c, bs = "re"),
         # Since we want fillers-only across all blocks
         data = processed_data),
    catchErrors = TRUE
  )
  
  saveRDS("a","~/tmp/ended")
  saveRDS(m[["value"]], paste0("~/workspace/scam_sfetal_nonlinear_group_test.RDS"))
  saveRDS(m[2:5],  paste0("~/workspace/scam_sfetal_nonlinear_group_test_errors.RDS"))
  time_count <- tictoc::toc()
  
  big_df <- get_predictions(m, loaded_data, processed_data)
  saveRDS(big_df, "~/workspace/scam_sfetal_nonlinear_group_test_predictions.RDS")
  
}
done(group_tester_gam)





#---------------------------------------------------------------------------------------

plan(list(
  # multisession,
  tweak(remote, workers="zburchil@cycle1.cs.rochester.edu"),
  tweak(cluster, workers="node92")
))
crit_vs_filler_gam %2% {
  loaded_data <- load_sfetal()
  
  processed_data <- loaded_data %>%
    process_data()
  
  tictoc::tic()
  saveRDS("a","~/tmp/starting")
  m <- zplyr::collect_all(
    scam(formula = logRT ~ 1 + 
           ItemType +
           s(Word.Length.c, k = 4) +
           s(Trial.Order.c, bs = "mpd", k = 7) +
           s(Trial.Order.c, by = ItemType, bs = "mpd", k = 7) +
           s(Subject, bs = "re") +
           s(Subject, Trial.Order.c, bs = "re") +
           s(Subject, Word.Length.c, bs = "re") +
           s(Item, bs = "re") +
           s(Item, Word.Length.c, bs = "re"),
         # Since we want fillers-only across all blocks
         data = processed_data),
    catchErrors = TRUE
  )
  
  saveRDS("a","~/tmp/ended")
  saveRDS(m[["value"]], paste0("~/workspace/scam_sfetal_nonlinear_critvfiller_test.RDS"))
  saveRDS(m[2:5],  paste0("~/workspace/scam_sfetal_nonlinear_critvfiller_test_errors.RDS"))
  time_count <- tictoc::toc()
  
  big_df <- get_predictions(m, loaded_data, processed_data)
  saveRDS(big_df, "~/workspace/scam_sfetal_nonlinear_critvfiller_test_predictions.RDS")
  
}
done(crit_vs_filler_gam)

