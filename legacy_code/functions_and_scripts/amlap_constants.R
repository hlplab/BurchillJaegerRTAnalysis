library(dplyr)
library(purrr)

# Amlap constants

# Connectivity constants ####################################
username <- "zburchil"
server <- "cycle1.cs.rochester.edu"
backup_server <- "cycle2.cs.rochester.edu"
tertiary_server <- "cycle3.cs.rochester.edu"
remote_login <-   paste0(username, "@", server)
backup_login <-   paste0(username, "@", backup_server)
tertiary_login <- paste0(username, "@", tertiary_server) 
first_node = "node88"
num_nodes <- 27
janky_nodes <- c(73, 33, 64, 63, 57, 54, 46)  # c(46, 44, 62, 55) #73, 74, 78, 77, 76)

# constants
k_per_file = 2000
total_k = 10000
k_list = split(1:total_k, ceiling(seq_along(1:total_k)/k_per_file))

# Residualizations have a lot more data and will crash R if the file sizes are too big
residual_k_per_file <- 20
residual_total_files <- total_k/residual_k_per_file



new_bins <- list("bin1" = c(3, 7),   # each bin has 5 trials
                 "bin2" = c(44, 50),  # (ie there are 5 filler
                 "bin3" = c(108, 115))#  trials between 80 and 87)

# setal_file <- "setal_sim_df.RDS"
setal_file <- "amlap_setal_sim_df.RDS"
fetal_file <- "fetal_sim_df.RDS"
nsc_file   <- "nsc_sim_df.RDS"
predicted_RT_of_mean_wl <- 327.8723

# The way this is calculated is included below
setal_basis_info <- structure(list(
  Bins = c("bin1", "bin2", "bin3"),
  meanRT = c(394.77842501106, 
             337.828536831483, 299.347902097902), 
  sdRT = c(204.377823899962, 
           165.142903447751, 135.67006135807), 
  meanLogRT = c(2.55114930536008, 
                2.49000322290013, 2.44366306023715), 
  sdLogRT = c(0.192444824828533, 
              0.175779208951909, 0.160181080794851)), 
  row.names = c(NA, -3L
  ), class = c("tbl_df", "tbl", "data.frame"))

# setal_basis_info <- readRDS(paste0(data_path, "setal_sim_df.RDS")) %>%
#   mutate(Bins = case_when(Trial.Order  >= new_bins$bin1[1] & Trial.Order <= new_bins$bin1[2] ~ "bin1",
#                           Trial.Order  >= new_bins$bin2[1] & Trial.Order <= new_bins$bin2[2] ~ "bin2",
#                           Trial.Order  >= new_bins$bin3[1] & Trial.Order <= new_bins$bin3[2] ~ "bin3",
#                           TRUE ~ NA_character_)) %>%
#   filter(!is.na(Bins)) %>%
#   group_by(Bins) %>%
#   summarise(meanRT = mean(RT), sdRT = sd(RT),
#             meanLogRT = mean(log10(RT)), sdLogRT = sd(log10(RT)))

simple_models <- list(
  linear_power = function(x) lm(RT_with_Effect ~ 1 + SimCond, data = x,
                                control = lme4::lmerControl(optimizer="bobyqa")),
  log_power    = function(x) lm(log10(RT_with_Effect) ~ 1 + SimCond, data = x,
                                control = lme4::lmerControl(optimizer="bobyqa")),
  linear_type1 = function(x) lm(RT ~ 1 + SimCond, data = x,
                                control = lme4::lmerControl(optimizer="bobyqa")),
  log_type1    = function(x) lm(log10(RT) ~ 1 + SimCond, data = x,
                                control = lme4::lmerControl(optimizer="bobyqa"))
)

residual_models <- list(
  linear_power = function(x) lm(RT_Resid_with_Effect ~ 1 + SimCond, data = x,
                                control = lme4::lmerControl(optimizer="bobyqa")),
  log_power    = function(x) lm(LogRT_Resid_with_Effect ~ 1 + SimCond, data = x,
                                control = lme4::lmerControl(optimizer="bobyqa")),
  linear_type1 = function(x) lm(RT_Resid ~ 1 + SimCond, data = x,
                                control = lme4::lmerControl(optimizer="bobyqa")),
  log_type1    = function(x) lm(LogRT_Resid ~ 1 + SimCond, data = x,
                                control = lme4::lmerControl(optimizer="bobyqa"))
)

# Supposes bin1 = -1, bin3 = 1
interaction_models <- list(
  linear_all =         function(x) lm(RT_with_ALL ~ 1 + SimCond*Block, data = x,
                                      control = lme4::lmerControl(optimizer="bobyqa")),
  linear_main =        function(x) lm(RT_with_Main ~ 1 + SimCond*Block, data = x,
                                      control = lme4::lmerControl(optimizer="bobyqa")),
  linear_interaction = function(x) lm(RT_with_Interaction ~ 1 + SimCond*Block, data = x,
                                      control = lme4::lmerControl(optimizer="bobyqa")),
  linear_none =        function(x) lm(RT ~ 1 + SimCond*Block, data = x,
                                      control = lme4::lmerControl(optimizer="bobyqa")),
  # -------
  log_all =         function(x) lm(log10(RT_with_ALL) ~ 1 + SimCond*Block, data = x,
                                   control = lme4::lmerControl(optimizer="bobyqa")),
  log_main =        function(x) lm(log10(RT_with_Main) ~ 1 + SimCond*Block, data = x,
                                   control = lme4::lmerControl(optimizer="bobyqa")),
  log_interaction = function(x) lm(log10(RT_with_Interaction) ~ 1 + SimCond*Block, data = x,
                                   control = lme4::lmerControl(optimizer="bobyqa")),
  log_none =        function(x) lm(log10(RT) ~ 1 + SimCond*Block, data = x,
                                   control = lme4::lmerControl(optimizer="bobyqa"))
)

# Supposes bin1 = -1, bin3 = 1
# these are the same as the interaction models, but without the interaction term
fake_interaction_models <- list(
  linear_fake_all =         function(x) lm(RT_with_ALL ~ 1 + SimCond+Block, data = x,
                                           control = lme4::lmerControl(optimizer="bobyqa")),
  linear_fake_main =        function(x) lm(RT_with_Main ~ 1 + SimCond+Block, data = x,
                                           control = lme4::lmerControl(optimizer="bobyqa")),
  linear_fake_interaction = function(x) lm(RT_with_Interaction ~ 1 + SimCond+Block, data = x,
                                           control = lme4::lmerControl(optimizer="bobyqa")),
  linear_fake_none =        function(x) lm(RT ~ 1 + SimCond+Block, data = x,
                                           control = lme4::lmerControl(optimizer="bobyqa")),
  # -------
  log_fake_all =         function(x) lm(log10(RT_with_ALL) ~ 1 + SimCond+Block, data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_fake_main =        function(x) lm(log10(RT_with_Main) ~ 1 + SimCond+Block, data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_fake_interaction = function(x) lm(log10(RT_with_Interaction) ~ 1 + SimCond+Block, data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_fake_none =        function(x) lm(log10(RT) ~ 1 + SimCond+Block, data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa"))
)

exgaussian_models <- list(
  exgauss_power = function(x) {
    b1 <- readRDS(paste0(path_on_server, "exgauss_brms.RDS"))
    new_m <- update(
      b1, newdata = x,
      formula. = RT_with_Effect ~ 1 + SimCond,
      iter = 4000, future = T,
      control = list(adapt_delta = 0.99999))},
  exgauss_type1 = function(x) {
    b1 <- readRDS(paste0(path_on_server, "exgauss_brms.RDS"))
    new_m <- update(
      b1, newdata = x,
      formula. = RT ~ 1 + SimCond,
      iter = 4000, future = T,
      control = list(adapt_delta = 0.99999))}
)

mixed_models <- list(
  linear_power = function(x) lme4::lmer(RT_with_Effect ~ 1 + SimCond + 
                                          (1 + SimCond | UniqueSubject) + 
                                          (1 + SimCond | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_power    = function(x) lme4::lmer(log10(RT_with_Effect) ~ 1 + SimCond + 
                                          (1 + SimCond | UniqueSubject) + 
                                          (1 + SimCond | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  linear_type1 = function(x) lme4::lmer(RT ~ 1 + SimCond + 
                                          (1 + SimCond | UniqueSubject) + 
                                          (1 + SimCond | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_type1    = function(x) lme4::lmer(log10(RT) ~ 1 + SimCond + 
                                          (1 + SimCond | UniqueSubject) + 
                                          (1 + SimCond | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa"))
)

no_subject_slopes <- list(
  linear_power = function(x) lme4::lmer(RT_with_Effect ~ 1 + SimCond + 
                                          (1  | UniqueSubject) + 
                                          (1 + SimCond | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_power    = function(x) lme4::lmer(log10(RT_with_Effect) ~ 1 + SimCond + 
                                          (1  | UniqueSubject) + 
                                          (1 + SimCond | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  linear_type1 = function(x) lme4::lmer(RT ~ 1 + SimCond + 
                                          (1 | UniqueSubject) + 
                                          (1 + SimCond | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_type1    = function(x) lme4::lmer(log10(RT) ~ 1 + SimCond + 
                                          (1 | UniqueSubject) + 
                                          (1 + SimCond | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa"))
)
no_slopes <- list(
  linear_power = function(x) lme4::lmer(RT_with_Effect ~ 1 + SimCond + 
                                          (1  | UniqueSubject) + 
                                          (1  | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_power    = function(x) lme4::lmer(log10(RT_with_Effect) ~ 1 + SimCond + 
                                          (1  | UniqueSubject) + 
                                          (1  | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  linear_type1 = function(x) lme4::lmer(RT ~ 1 + SimCond + 
                                          (1 | UniqueSubject) + 
                                          (1 | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_type1    = function(x) lme4::lmer(log10(RT) ~ 1 + SimCond + 
                                          (1 | UniqueSubject) + 
                                          (1 | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa"))
)

no_slopes_wo_infs <- list(
  linear_power = function(x) lme4::lmer(RT_with_Effect ~ 1 + SimCond + 
                                          (1  | UniqueSubject) + 
                                          (1  | UniqueItem), 
                                        data = filter(x, 
                                                      RT_with_Effect > 0,
                                                      RT > 0), 
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_power    = function(x) lme4::lmer(log10(RT_with_Effect) ~ 1 + SimCond + 
                                          (1  | UniqueSubject) + 
                                          (1  | UniqueItem), 
                                        data = filter(x, 
                                                      RT_with_Effect > 0,
                                                      RT > 0), 
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  linear_type1 = function(x) lme4::lmer(RT ~ 1 + SimCond + 
                                          (1 | UniqueSubject) + 
                                          (1 | UniqueItem), 
                                        data = filter(x, 
                                                      RT_with_Effect > 0,
                                                      RT > 0), 
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_type1    = function(x) lme4::lmer(log10(RT) ~ 1 + SimCond + 
                                          (1 | UniqueSubject) + 
                                          (1 | UniqueItem), 
                                        data = filter(x, 
                                                      RT_with_Effect > 0,
                                                      RT > 0), 
                                        control = lme4::lmerControl(optimizer="bobyqa"))
)


residual_no_slopes <- list(
  linear_power = function(x) lme4::lmer(RT_Resid_with_Effect ~ 1 + SimCond + 
                                          (1  | UniqueSubject) + 
                                          (1  | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_power    = function(x) lme4::lmer(LogRT_Resid_with_Effect ~ 1 + SimCond + 
                                          (1  | UniqueSubject) + 
                                          (1  | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  linear_type1 = function(x) lme4::lmer(RT_Resid ~ 1 + SimCond + 
                                          (1 | UniqueSubject) + 
                                          (1 | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa")),
  log_type1    = function(x) lme4::lmer(LogRT_Resid ~ 1 + SimCond + 
                                          (1 | UniqueSubject) + 
                                          (1 | UniqueItem), data = x,
                                        control = lme4::lmerControl(optimizer="bobyqa"))
)

# # Checks to see that the BATA residualization model converged
# safe_residual_no_slopes <- list(
#   linear_power = function(x) {
#     if (all(is.na(x$BATA_PredRawRT_with_Effect)) || all(is.na(x$RT_Resid_with_Effect))) {
#       stop("Linear power BATA residualization model did not converge")
#     }
#     lme4::lmer(RT_Resid_with_Effect ~ 1 + SimCond + 
#                  (1  | UniqueSubject) + 
#                  (1  | UniqueItem), data = x,
#                control = lme4::lmerControl(optimizer="bobyqa"))
#   },
#   log_power = function(x) {
#     if (all(is.na(x$BATA_PredLogRT_with_Effect)) || all(is.na(x$LogRT_Resid_with_Effect))) {
#       stop("Log power BATA residualization model did not converge")
#     }
#     lme4::lmer(LogRT_Resid_with_Effect ~ 1 + SimCond + 
#                  (1  | UniqueSubject) + 
#                  (1  | UniqueItem), data = x,
#                control = lme4::lmerControl(optimizer="bobyqa"))},
#   linear_type1 = function(x) {
#     if (all(is.na(x$BATA_PredRawRT)) || all(is.na(x$RT_Resid))) {
#       stop("Linear type 1 BATA residualization model did not converge")
#     }
#     lme4::lmer(RT_Resid ~ 1 + SimCond + 
#                  (1  | UniqueSubject) + 
#                  (1  | UniqueItem), data = x,
#                control = lme4::lmerControl(optimizer="bobyqa"))
#   },
#   log_type1 = function(x) {
#     if (all(is.na(x$BATA_PredLogRT)) || all(is.na(x$LogRT_Resid))) {
#       stop("Log type 1 BATA residualization model did not converge")
#     }
#     lme4::lmer(LogRT_Resid ~ 1 + SimCond + 
#                  (1  | UniqueSubject) + 
#                  (1  | UniqueItem), data = x,
#                control = lme4::lmerControl(optimizer="bobyqa"))}
# )




safe_residual_no_slopes_no_subject <- list(
  linear_power = function(x) {
    if (all(is.na(x$RT_Resid_with_Effect)) | all(is.nan(x$RT_Resid_with_Effect))) {
      stop("Linear power BATA residualization model did not converge")
    }
    lme4::lmer(RT_Resid_with_Effect ~ 1 + SimCond + (1  | UniqueItem), data = x,
               control = lme4::lmerControl(optimizer="bobyqa"))
  },
  log_power = function(x) {
    if (all(is.na(x$LogRT_Resid_with_Effect)) | all(is.nan(x$LogRT_Resid_with_Effect))) {
      stop("Log power BATA residualization model did not converge")
    }
    lme4::lmer(LogRT_Resid_with_Effect ~ 1 + SimCond + (1  | UniqueItem), data = x,
               control = lme4::lmerControl(optimizer="bobyqa"))
  },
  linear_type1 = function(x) {
    if (all(is.na(x$RT_Resid)) | all(is.nan(x$RT_Resid))) {
      stop("Linear type 1 BATA residualization model did not converge")
    }
    lme4::lmer(RT_Resid ~ 1 + SimCond + (1 | UniqueItem), data = x,
               control = lme4::lmerControl(optimizer="bobyqa"))
  },
  log_type1 = function(x) {
    if (all(is.na(x$LogRT_Resid)) | all(is.nan(x$LogRT_Resid))) {
      stop("Log type 1 BATA residualization model did not converge")
    }
    lme4::lmer(LogRT_Resid ~ 1 + SimCond + (1 | UniqueItem), data = x,
               control = lme4::lmerControl(optimizer="bobyqa"))
  }
)


i_mm_block_slopes <- list(
  linear_all =         function(x) lme4::lmer(RT_with_ALL ~ 1 + SimCond * Block + 
                                                (1  + Block | UniqueSubject) + 
                                                (1  | UniqueItem), data = x,
                                              control = lme4::lmerControl(optimizer="bobyqa")),
  linear_main =        function(x) lme4::lmer(RT_with_Main ~ 1 + SimCond * Block + 
                                                (1  + Block | UniqueSubject) + 
                                                (1  | UniqueItem), data = x,
                                              control = lme4::lmerControl(optimizer="bobyqa")),
  linear_interaction = function(x) lme4::lmer(RT_with_Interaction ~ 1 + SimCond * Block + 
                                                (1  + Block | UniqueSubject) + 
                                                (1  | UniqueItem), data = x,
                                              control = lme4::lmerControl(optimizer="bobyqa")),
  linear_none =        function(x) lme4::lmer(RT ~ 1 + SimCond * Block + 
                                                (1  + Block | UniqueSubject) + 
                                                (1  | UniqueItem), data = x,
                                              control = lme4::lmerControl(optimizer="bobyqa")),
  
  log_all =         function(x) lme4::lmer(log10(RT_with_ALL) ~ 1 + SimCond * Block + 
                                             (1  + Block | UniqueSubject) + 
                                             (1  | UniqueItem), data = x,
                                           control = lme4::lmerControl(optimizer="bobyqa")),
  log_main =        function(x) lme4::lmer(log10(RT_with_Main) ~ 1 + SimCond * Block + 
                                             (1  + Block | UniqueSubject) + 
                                             (1  | UniqueItem), data = x,
                                           control = lme4::lmerControl(optimizer="bobyqa")),
  log_interaction = function(x) lme4::lmer(log10(RT_with_Interaction) ~ 1 + SimCond * Block + 
                                             (1  + Block | UniqueSubject) + 
                                             (1  | UniqueItem), data = x,
                                           control = lme4::lmerControl(optimizer="bobyqa")),
  log_none =        function(x) lme4::lmer(log10(RT) ~ 1 + SimCond * Block + 
                                             (1  + Block | UniqueSubject) + 
                                             (1  | UniqueItem), data = x,
                                           control = lme4::lmerControl(optimizer="bobyqa"))
)