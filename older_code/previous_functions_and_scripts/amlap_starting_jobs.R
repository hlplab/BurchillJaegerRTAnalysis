library(dplyr)
library(purrr)
library(future)
library(furrr)
library(cs)

path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
path_on_server <- "/u/zburchil/workspace/launchpad/"

source(paste0(path, "modeling_data_functions.R"))
source(paste0(path, "amlap_modeling_data_functions.R"))



# Connectivity constants ####################################
username <- "zburchil"
server <- "cycle1.cs.rochester.edu"
backup_server <- "cycle2.cs.rochester.edu"
remote_login <- paste0(username, "@", server)
backup_login <- paste0(username, "@", backup_server)
first_node = "node88"
num_nodes <- 27
janky_nodes <- c(46, 44, 62, 55, 73, 74, 78, 77, 76)

ok_nodes <- get_n_best_nodes(remote_login, num_nodes,
                              !(number %in% janky_nodes),
                             # number != 92,
                              number > 72) %>%
  mutate(node = stringr::str_match(node, "node[0-9]+")[,1]) %>%
  { .$node } 
ok_nodes <- ok_nodes %>%
  test_nodes(remote_login, timeout_sec = 5, verbose=TRUE)

# 87, 82, 78, 

# ok_nodes <- paste0("node", c(75, 78:81,83:86,88:91)) %>%
#   test_nodes(remote_login, timeout_sec = 5, verbose = TRUE)
# beepr::beep()


# Jankier way:
# ok_nodes <- paste0("node", 72:91) %>%
  # test_nodes(remote_login, timeout_sec = 5, verbose=TRUE)
# ok_nodes %beep%
  # test_nodes(paste0("node", 72:91), remote_login, timeout_sec = 5) %plan%
  # list(multisession, multisession)
# resolved(futureOf(ok_nodes))
# ok_nodes <- c("node75", "node79", "node81", "node84", "node86", "node88", "node89", "node91", "node92")
# ok_nodes <- c("node75", "node79", "node80", "node84", "node88", "node89", "node92")

# if (length(node_list) > 120) {
#   node_list<-node_list[1:120]
#   stop("Node list needs to be less than 128")
# }

#####################################################################
# constants
k_per_file = 2000
total_k = 10000
k_list = split(1:total_k, ceiling(seq_along(1:total_k)/k_per_file))

new_bins <- list("bin1" = c(3, 7),   # each bin has 5 trials
                 "bin2" = c(44, 50),  # (ie there are 5 filler
                 "bin3" = c(108, 115))#  trials between 80 and 87)

# setal_file <- "setal_sim_df.RDS"
setal_file <- "amlap_setal_sim_df.RDS"
fetal_file <- "fetal_sim_df.RDS"
# predicted_RT_of_mean_wl <- 327.8723

setal_basis_info <- readRDS(paste0(path, "../", "setal_sim_df.RDS")) %>%
  mutate(Bins = case_when(Trial.Order  >= new_bins$bin1[1] & Trial.Order <= new_bins$bin1[2] ~ "bin1",
                          Trial.Order  >= new_bins$bin2[1] & Trial.Order <= new_bins$bin2[2] ~ "bin2",
                          Trial.Order  >= new_bins$bin3[1] & Trial.Order <= new_bins$bin3[2] ~ "bin3",
                          TRUE ~ NA_character_)) %>%
  filter(!is.na(Bins)) %>%
  group_by(Bins) %>%
  summarise(meanRT = mean(RT), sdRT = sd(RT),
            meanLogRT = mean(log10(RT)), sdLogRT = sd(log10(RT)))

simple_models <- list(
  linear_power = function(x) lm(RT_with_Effect ~ 1 + SimCond, data = x),
  log_power    = function(x) lm(log10(RT_with_Effect) ~ 1 + SimCond, data = x),
  linear_type1 = function(x) lm(RT ~ 1 + SimCond, data = x),
  log_type1    = function(x) lm(log10(RT) ~ 1 + SimCond, data = x)
)

residual_models <- list(
  linear_power = function(x) lm(RT_Resid_with_Effect ~ 1 + SimCond, data = x),
  log_power    = function(x) lm(LogRT_Resid_with_Effect ~ 1 + SimCond, data = x),
  linear_type1 = function(x) lm(RT_Resid ~ 1 + SimCond, data = x),
  log_type1    = function(x) lm(LogRT_Resid ~ 1 + SimCond, data = x)
)

# Supposes bin1 = -1, bin3 = 1
interaction_models <- list(
  linear_all =         function(x) lm(RT_with_ALL ~ 1 + SimCond*Block, data = x),
  linear_main =        function(x) lm(RT_with_Main ~ 1 + SimCond*Block, data = x),
  linear_interaction = function(x) lm(RT_with_Interaction ~ 1 + SimCond*Block, data = x),
  linear_none =        function(x) lm(RT ~ 1 + SimCond*Block, data = x),
  # -------
  log_all =         function(x) lm(log10(RT_with_ALL) ~ 1 + SimCond*Block, data = x),
  log_main =        function(x) lm(log10(RT_with_Main) ~ 1 + SimCond*Block, data = x),
  log_interaction = function(x) lm(log10(RT_with_Interaction) ~ 1 + SimCond*Block, data = x),
  log_none =        function(x) lm(log10(RT) ~ 1 + SimCond*Block, data = x)
)

# Supposes bin1 = -1, bin3 = 1
# these are the same as the interaction models, but without the interaction term
fake_interaction_models <- list(
  linear_all =         function(x) lm(RT_with_ALL ~ 1 + SimCond+Block, data = x),
  linear_main =        function(x) lm(RT_with_Main ~ 1 + SimCond+Block, data = x),
  linear_interaction = function(x) lm(RT_with_Interaction ~ 1 + SimCond+Block, data = x),
  linear_none =        function(x) lm(RT ~ 1 + SimCond+Block, data = x),
  # -------
  log_all =         function(x) lm(log10(RT_with_ALL) ~ 1 + SimCond+Block, data = x),
  log_main =        function(x) lm(log10(RT_with_Main) ~ 1 + SimCond+Block, data = x),
  log_interaction = function(x) lm(log10(RT_with_Interaction) ~ 1 + SimCond+Block, data = x),
  log_none =        function(x) lm(log10(RT) ~ 1 + SimCond+Block, data = x)
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

################################### killling jobs
killer <- kill_r_on_nodes(ok_nodes, backup_login, # "zburchil@cycle3.cs.rochester.edu"
                          wait_until_resolved = TRUE,
                          beep_on_end = TRUE)
resolved(killer)





#################################################################################
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = ok_nodes),
  # sequential # for the brms models that probably do better with... wait no, i bet not
  tweak(multiprocess, workers = function() {max(round(future::availableCores()/2), 4)})
))

# A flatter way of doing running things
crossed_df <- tidyr::crossing(names(new_bins),
                              seq_along(k_list),
                              c(16, 32, 64, 128, 526)) %>%
  set_names(c("Bins","Knum", "Nsubj")) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(iprefix = "full_random")

bigger %beep% {
  future_map(
    c(1:nrow(crossed_df)),
    function(i) {
      effect_size <- 14
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==crossed_df[i,]$Bins, ]$meanRT,
        effect_size)
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[crossed_df[i,]$Knum]],
        iprefix = paste0(crossed_df[i,]$iprefix, "_",
                         crossed_df[i,]$Nsubj),
        isuffix = crossed_df[i,]$Bins,
        ofile_appended_name = "residualized_14", # "normal_n_log_effect_size_"
        list_of_model_funcs = residual_models, # simple_models,
        effect_adder_function = function(df) {
          # mutate(df, RT_with_Effect = add_in_log_space(RT,SimCond,log_effect))}, #$PredRT_WL
          # mutate(df, RT_with_Effect = RT + SimCond * effect_size)},
          # mutate(df,
          #        RT_Resid_with_Effect =         add_in_log_space(RT, SimCond, log_effect) - PredRT_WL,
          #        LogRT_Resid_with_Effect =log10(add_in_log_space(RT, SimCond, log_effect))- PredLogRT_WL)},
          mutate(df,
                 RT_Resid_with_Effect =         (RT + SimCond * effect_size) - PredRT_WL,
                 LogRT_Resid_with_Effect = log10(RT + SimCond * effect_size) - PredLogRT_WL)},
        data_tidyer = broom::tidy)
    })
}
resolved(futureOf(bigger))



#### For the interaction models: ##############-------------#####################-----------------###############
interaction_crossed_df <- crossed_df %>%
  filter(Bins=="bin1") %>% select(-Bins) %>%
  mutate(iprefix = "interaction_full_random")

interaction_stuff %beep% {
  future_map(
    c(1:nrow(interaction_crossed_df)),
    function(i) {
      effect_size <- 56
      interaction_effect <- 28
      # log_effect <- log_space_effect_size_calculator(
      #   setal_basis_info[setal_basis_info$Bins==interaction_crossed_df[i,]$Bins, ]$meanRT,
      #   effect_size)
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[interaction_crossed_df[i,]$Knum]],
        iprefix = paste0(interaction_crossed_df[i,]$iprefix, "_",
                         interaction_crossed_df[i,]$Nsubj),
        isuffix = "bin1",
        ofile_appended_name = "normal_n_log_effect_size_56_28",
        list_of_model_funcs = interaction_models,
        effect_adder_function = function(df) {
          # mutate(df, RT_with_Effect = add_in_log_space(RT,SimCond,log_effect))},
          # mutate(df, RT_with_Effect = RT + SimCond * effect_size)},
          mutate(df, 
                 Block = factor(NewGroup, levels=c("First","Second")) %>% `contrasts<-`(value = c(-1, 1)),
                 BlockDbl = case_when(Block=="First" ~ -1, Block=="Second" ~ 1, TRUE ~ NA_real_),
                 RT_with_ALL = RT + SimCond * effect_size + BlockDbl*SimCond*interaction_effect,
                 RT_with_Interaction = RT + BlockDbl*SimCond*interaction_effect,
                 RT_with_Main =  RT + SimCond * effect_size)},
        data_tidyer = broom::tidy)
    })
}
resolved(futureOf(interaction_stuff))





####### Running the parametric data generation: ############################################

# Ehhhh, ehhhhhh
temp_df <- mega_cross_df %>% group_by_at(vars(everything(), -Space, -filename, -iprefix)) %>% summarise()

gaussian_parametric %<-% {
  future_walk(
    c(1:nrow(temp_df)),
    function(i) {
      effect_size <- temp_df[i,]$EffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins == temp_df[i,]$Bins,]$meanRT,
        effect_size)
      generate_parametric_list(info_row = temp_df[i,],
                               iprefix = "gaussian_para",
                               effect_name = paste0("normal_n_log_effect_size_",
                                                    
                                                    effect_size, ""),
                               server_path = path_on_server,
                               effect_adder_func = function(df) {
                                 # mutate(df, RT_with_Effect = add_in_log_space(RT, SimCond, log_effect))}, # change
                                 mutate(df, RT_with_Effect = RT + effect_size * SimCond)}, # change
                               list_of_model_funcs = simple_models,
                               # function stuff
                               parametric_func = rnorm, # function(x)
                               # happens before adding the sim condition/effect
                               # The generated predictor is `Y`
                               df_modifier_function = function(df) { mutate(df, RT = Y) },
                               n = temp_df[i,]$Nsubj,
                               mean=filter(setal_basis_info, Bins == temp_df[i,]$Bins)$meanRT,
                               sd = filter(setal_basis_info, Bins == temp_df[i,]$Bins)$sdRT)
    })
}
resolved(futureOf(gaussian_parametric))



######### To get the results: ####################################

mega_cross_df <- tidyr::crossing(names(new_bins),
                                 seq_along(k_list),
                                 c(32, 64, 128),
                                 c("full_random", "gaussian_para", "lognormal_para"),
                                 c(7,14,28,56),
                                 c("", "_log_space")) %>%
  set_names(c("Bins","Knum", "Nsubj", "iprefix","EffectSize","Space")) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]])),
         ModelNames = "residualized") %>%
  filter(EffectSize %in% c(28, 56), iprefix=="full_random") %>%
  mutate(filename = paste0(paste("modeldata", iprefix, Nsubj, Bins, ModelNames, EffectSize, sep="_"), Space, "_", MinK, "_to_", MaxK, ".RDS"))


# Are all the files there?
file_checker_df %beep% {
  future_map_dfr(
    c(1:nrow(mega_cross_df)),
    function(i) {
      there <- FALSE
      if (file.exists(paste0(path_on_server,
                             "amlap_saved_files/",
                             mega_cross_df[i,]$filename)))
        there <- TRUE
      return(mega_cross_df[i,] %>% mutate(there=there))
    })
}

# Load it up
big_df %beep% {
  biggie <- future_map_dfr(
    c(1:nrow(mega_cross_df)),
    function(i) {
      info_row <- mega_cross_df[i,]
      results_path <- paste0(path_on_server, "amlap_saved_files/results/", info_row$filename)
      if (file.exists(results_path)) {
        df <- readRDS(results_path)
        if (min(df$k) == 1 & info_row$MinK != 1)
          df <- mutate(df, k = k + info_row$MinK - 1)
      } else {
        df <- collate_data(paste0(path_on_server, "amlap_saved_files/",
                                  info_row$filename), .multicore = TRUE) %>%
          mutate(Nsubj = info_row$Nsubj,
                 SampleName = info_row$iprefix,
                 Bin = info_row$Bins,
                 EffectSize = info_row$EffectSize,
                 ModelNames = info_row$ModelNames,
                 Space = info_row$Space)
        # Ugh, somehow forgot to do this
        if (min(df$k) == 1 & info_row$MinK != 1)
          df <- mutate(df, k = k + info_row$MinK - 1)
        # Saving the files in case it doesn't load all the way
        saveRDS(df, results_path)
      }
      df
    }
  )
  saveRDS(biggie,  paste0(path_on_server, "amlap_saved_files/results/experimental_biggie.RDS"))
  paste0("rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/biggie.RDS ~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/PowerSimulations/data/experimental/all_data_round_1.RDS")
}
resolved(futureOf(big_df))


# The BRMS stuff ###################################################################
plan(list(
  # multisession, # for the beep!
  tweak(remote,  workers = backup_login),
  tweak(cluster, workers = "node92"),
  tweak(multiprocess, workers = function() {max(round(future::availableCores()/2), 4)})
))

# A flatter way of doing running things
brms_df <- tidyr::crossing(names(new_bins),
                           seq_along(k_list),
                           c(128)) %>%
  set_names(c("Bins","Knum", "Nsubj")) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(iprefix = "full_random")

brms_stuff %<-% {
  future_walk(
    c(1:1),
    function(i) {
      effect_size <- 56
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==brms_df[i,]$Bins, ]$meanRT,
        effect_size)
      # log_effect <- log_space_effect_size_calculator(
      #   abs(predicted_RT_of_mean_wl - setal_basis_info[setal_basis_info$Bins==brms_df[i,]$Bins, ]$meanRT),
      #   effect_size)
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[brms_df[i,]$Knum]],
        iprefix = paste0(brms_df[i,]$iprefix, "_",
                         brms_df[i,]$Nsubj),
        isuffix = brms_df[i,]$Bins,
        ofile_appended_name = "exgaussian_56", # "normal_n_log_effect_size_"
        list_of_model_funcs = exgaussian_models, # simple_models,
        effect_adder_function = function(df) {
          # mutate(df, RT_with_Effect = add_in_log_space(RT,SimCond,log_effect))}, #$PredRT_WL
          mutate(df, RT_with_Effect = RT + SimCond * effect_size)},
        # mutate(df,
        #        RT_Resid_with_Effect =         add_in_log_space(RT, SimCond, log_effect) - PredRT_WL,
        #        LogRT_Resid_with_Effect =log10(add_in_log_space(RT, SimCond, log_effect))- PredLogRT_WL)},
        # mutate(df,
        #        RT_Resid_with_Effect =         (RT + SimCond * effect_size) - PredRT_WL,
        #        LogRT_Resid_with_Effect = log10(RT + SimCond * effect_size) - PredLogRT_WL)},
        data_tidyer = tidy_brms_model)
    })
}
resolved(futureOf(brms_stuff))








# Monitoring cluster activity ###################################################

cs::monitor_cluster_resources("zburchil", "zburchil@cycle3.cs.rochester.edu", ok_nodes,
                              save_path="/u/zburchil/bb_maker_resources.RDS",
                              sleeping_time = 30,
                              total_checks = 90)

plan(remote, workers = "zburchil@cycle3.cs.rochester.edu")

hmm %<-% readRDS("/u/zburchil/bb_maker_resources.RDS")
resolved(futureOf(hmm))

