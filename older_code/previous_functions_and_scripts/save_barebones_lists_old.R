library(dplyr)
library(purrr)
library(future)
library(furrr)
library(cs)
library(lme4)

path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
path_on_server <- "/u/zburchil/workspace/launchpad/"

source(paste0(path, "barebones_making_functions.R"))

# setal <- readRDS(paste0(path, "../setal_sim_df.RDS"))
# fetal <- readRDS(paste0(path, "../fetal_sim_df.RDS"))

# Connectivity constants ####################################
username <- "zburchil"
server <- "cycle1.cs.rochester.edu"
remote_login <- paste0(username, "@", server)
first_node = "node88"


# Start Plan
plan(list(
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = first_node)
))

# effect sizes
effect_sizes_df <- readRDS(paste0(path, "../stack_and_fine_effect_sizes_2018-07-06.RDS"))
q1_setal_effect_size <- effect_sizes_df %>%
  filter(Prediction_Number == "1 (MV in block 3)",
         Experiment == "Harrington Stack et al. (2018)") %>%
         {.$EffectValue[1]}
q2_setal_effect_size <- effect_sizes_df %>%
  filter(Prediction_Number == "2 (RC b/w 1 and 2)",
         Experiment == "Harrington Stack et al. (2018)") %>%
         {.$EffectValue[1]}
q3_setal_effect_size <- effect_sizes_df %>%
  filter(Prediction_Number == "3 (RC in block 2)",
         Experiment == "Harrington Stack et al. (2018)") %>%
         {.$EffectValue[1]}
q1_fetal_effect_size <- effect_sizes_df %>%
  filter(Prediction_Number == "1 (MV in block 3)",
         Experiment == "Fine et al. (2013)") %>%
         {.$EffectValue[1]}
q2_fetal_effect_size <- effect_sizes_df %>%
  filter(Prediction_Number == "2 (RC b/w 1 and 2)",
         Experiment == "Fine et al. (2013)") %>%
         {.$EffectValue[1]}
q3_fetal_effect_size <- effect_sizes_df %>%
  filter(Prediction_Number == "3 (RC in block 2)",
         Experiment == "Fine et al. (2013)") %>%
         {.$EffectValue[1]}

########################################################

# constants
k1 <- 1:2000
k2 <- 2001:4000
k3 <- 4001:6000
k4 <- 6001:8000
k5 <- 8001:10000
k_list <- list(k1, k2, k3, k4, k5)
if (testing==TRUE)
  k_list <- list(k1)

# node 88 for the 10k run
plan(list(
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = "node88")
))

setal_q1 %2% {
  for (k_e in k_list) {
    save_list(server_path = path_on_server,
              real_df_file = "setal_sim_df.RDS",
              thing_id = "setal_fullsize",
              question = "q1",
              k = k_e,
              effect_size = q1_setal_effect_size,
              block_filter = "3",
              trials_per_subject = 10,
              subjects_per_group = 208,
              seed_function = identity)
  }
}
resolved(futureOf(setal_q1))


plan(list(
  tweak(remote,  workers = remote_login),
  tweak(remote, workers = "node89")
))
setal_q2 %2% {
  for (k_e in k_list) {
  save_list(server_path = path_on_server,
            real_df_file = "setal_sim_df.RDS",
            thing_id = "setal_fullsize",
            question = "q2",
            k = k_e,
            effect_size = q2_setal_effect_size,
            trials_per_block1 = 16,
            trials_per_block2 = 9,
            n_subjects = 205,
            seed_function = identity)
  }
}
resolved(futureOf(setal_q2))


plan(list(
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = "node86")
))
setal_q3 %2% {
  for (k_e in k_list) {
  save_list(server_path = path_on_server,
            real_df_file = "setal_sim_df.RDS",
            thing_id = "setal_fullsize",
            question = "q3",
            k = k_e,
            effect_size = q3_setal_effect_size,
            block_filter = "2",
            trials_per_subject = 9,
            subjects_per_group = 208,
            seed_function = identity)
  }
}
resolved(futureOf(setal_q3))

#################################################################
# Fetal
plan(list(
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = "node88")
))
fetal_q1 %2% {
  for (k_e in k_list) {
    save_list(server_path = path_on_server,
              real_df_file = "fetal_sim_df.RDS",
              thing_id = "fetal",
              question = "q1",
              k = k_e,
              effect_size = q1_fetal_effect_size,
              block_filter = "3",
              trials_per_subject = 5,
              subjects_per_group = 39,
              seed_function = identity)
  }
}
resolved(futureOf(fetal_q1))


plan(list(
  tweak(remote,  workers = remote_login),
  tweak(remote, workers = "node89")
))
fetal_q2 %2% {
  for (k_e in k_list) {
    save_list(server_path = path_on_server,
              real_df_file = "fetal_sim_df.RDS",
              thing_id = "fetal",
              question = "q2",
              k = k_e,
              effect_size = q2_fetal_effect_size,
              trials_per_block1 = 8,
              trials_per_block2 = 5,
              n_subjects = 39,
              seed_function = identity)
  }
}
resolved(futureOf(fetal_q2))


plan(list(
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = "node92")
))
fetal_q3 %2% {
  for (k_e in k_list) {
    save_list(server_path = path_on_server,
              real_df_file = "fetal_sim_df.RDS",
              thing_id = "fetal",
              question = "q3",
              k = k_e,
              effect_size = q3_fetal_effect_size,
              block_filter = "2",
              trials_per_subject = 5,
              subjects_per_group = 78,
              seed_function = identity)
  }
}
resolved(futureOf(fetal_q3))

###################################################
# short setal

plan(list(
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = "node88")
))
setal_q1 %2% {
  for (k_e in k_list) {
    save_list(server_path = path_on_server,
              real_df_file = "setal_sim_df.RDS",
              thing_id = "setal_small",
              question = "q1",
              k = k_e,
              effect_size = q1_setal_effect_size,
              block_filter = "3",
              trials_per_subject = 5,
              subjects_per_group = 39,
              seed_function = identity)
  }
}
resolved(futureOf(setal_q1))


plan(list(
  tweak(remote,  workers = remote_login),
  tweak(remote, workers = "node89")
))
setal_q2 %2% {
  for (k_e in k_list) {
    save_list(server_path = path_on_server,
              real_df_file = "setal_sim_df.RDS",
              thing_id = "setal_small",
              question = "q2",
              k = k_e,
              effect_size = q2_setal_effect_size,
              trials_per_block1 = 8,
              trials_per_block2 = 5,
              n_subjects = 39,
              seed_function = identity)
  }
}
resolved(futureOf(setal_q2))


plan(list(
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = "node86")
))
setal_q3 %2% {
  for (k_e in k_list) {
    save_list(server_path = path_on_server,
              real_df_file = "setal_sim_df.RDS",
              thing_id = "setal_small",
              question = "q3",
              k = k_e,
              effect_size = q3_setal_effect_size,
              block_filter = "2",
              trials_per_subject = 5,
              subjects_per_group = 78,
              seed_function = identity)
  }
}
resolved(futureOf(setal_q3))


# l2 <- purrr::map(l, ~get_model_stats(., setal, model_func_list2))


# setal_q1 %2% {
#   zplyr::collect_all({
#     setal <- readRDS(paste0(path_on_server, "setal_sim_df.RDS"))
#     make_q1_q3_barebones_list(df = setal, k=k,
#                               block_filter="3",
#                               effect_size = q1_setal_effect_size,
#                               trials_per_subject = 10,
#                               subjects_per_group = 208,
#                               seed_function = identity) %>%
#       saveRDS(paste0(path_on_server,
#                      "saved_files/bb_setal_fullsize_q1_1_to_100.RDS"))
#     "all good"},
#     catchErrors = TRUE)
# }
# resolved(futureOf(setal_q1))

122761652
2087366322
2766206
