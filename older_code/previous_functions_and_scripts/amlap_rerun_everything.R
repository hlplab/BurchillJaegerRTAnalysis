# So this doesn't rerun EVERYTHING, I realize. It only reruns the MODELS, not the "bb" functions

library(dplyr)
library(purrr)
library(future)
library(furrr)
library(cs)

path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
path_on_server <- "/u/zburchil/workspace/launchpad/"

source(paste0(path, "modeling_data_functions.R"))
source(paste0(path, "amlap_modeling_data_functions.R"))
source(paste0(path, "amlap_constants.R"))
source(paste0(path, "utils.R"))

done <- cs::done

ok_nodes <- get_n_best_nodes(remote_login, num_nodes,
                             !(number %in% janky_nodes),
                             # number != 92,
                             number > 72) %>%
  mutate(node = stringr::str_match(node, "node[0-9]+")[,1]) %>%
  { .$node } 
ok_nodes <- ok_nodes %>%
  test_nodes(remote_login, timeout_sec = 5, verbose = TRUE)
# c("node85", "node84", "node87", "node86", "node79", "node80", "node83", "node81", "node92", "node89", "node88", "node75")

# flatter_nodes <- get_n_best_nodes(remote_login, num_nodes,
#                                   !(number %in% janky_nodes)) %>%
#   mutate(node = stringr::str_match(node, "node[0-9]+")[,1]) %>%
#   {map2(.$node, .$ncpus,
#         function(node_name, cpus) {
#           rep(node_name, round(cpus * 0.5))
#         })
#   } %>% unlist() %>% {
#     good_nodes <- test_nodes(unique(.), remote_login, 
#                              timeout_sec = 5, verbose=TRUE)
#     all <- .[. %in% good_nodes]
#     if (length(all) > 120) {
#       all<-all[1:120]
#     }
#     all
#   }



################################### killling jobs
killer <- kill_r_on_nodes(ok_nodes, backup_login, # "zburchil@cycle3.cs.rochester.edu"
                          wait_until_resolved = TRUE,
                          beep_on_end = TRUE)


# ############################## PLAN ###################################################################
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = ok_nodes),
  tweak(multiprocess, workers = function() {max(round(future::availableCores()/2), 4)})
))

# for flatter things
# plan(list(
#   multisession, # for the beep!
#   tweak(remote,  workers = remote_login),
#   tweak(cluster, workers = flatter_nodes)
# ))

####################### GETTING THE JOB DATA FRAMES ######################################

# Running the normal models -------------------------------
crossed_df <- tidyr::crossing(Bins = names(new_bins),
                              Knum = seq_along(k_list),
                              Nsubj = c(16, 32, 64, 128, 256),
                              EffectSize = c(7, 14, 28, 56),
                              Space = c("_log_space", ""),
                              RType = c("residualized", "unresidualized"),
                              iprefix = c("full_random")
) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  rowwise() %>%
  mutate(ofilename = amlap_file_coordinator(
    server_path = path_on_server,
    real_df_file = setal_file,
    k = k_list[[Knum]],
    iprefix = paste0(iprefix, "_", Nsubj),
    isuffix = Bins,
    ofile_appended_name = paste(RType, "size", EffectSize, sep="_") %>% paste0(Space)
  )$ofilename) %>% ungroup()  #%>% filter(Nsubj %in% c(16, 256))
# Running the interaction models -------------------------------
interaction_crossed_df <- tidyr::crossing(
  Bins = c("bin1"), #since interaction is between block1 and block3
  Knum = seq_along(k_list),
  Nsubj = c(32, 128, 256),
  MainEffectSize = c(56),
  InteractionEffectSize = c(28),
  Block1Code = c(-1, 1),
  Interaction = c("real", "fake"),
  Space = c(""),
  RType = c("unresidualized"),
  iprefix = c("interaction")
) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(ofilename = paste("modeldata", Interaction, "interaction",
                           ifelse(Block1Code<0, "neg", "pos"),
                           Nsubj, Bins, "size",
                           MainEffectSize, InteractionEffectSize, 
                           MinK, "to", MaxK,
                           sep = "_") %>% paste0(".RDS"))
# --- Running the parametric models ----------------------------
# Ehhhh, ehhhhhh
parametric_df <- crossed_df %>% 
  filter(RType=="unresidualized") %>%
  tidyr::crossing(Distribution = c("gaussian_para", "lognormal_para")) %>%
  mutate(ofilename = paste(
    "modeldata",
    Distribution, Nsubj, Bins,
    paste0("size_", EffectSize, Space),
    MinK, "to", MaxK, sep="_") %>% paste0(".RDS"))
# ------ The box cox models ---------------------------------
bc_df <- crossed_df %>%
  filter(EffectSize %in% c(56),
         RType=="unresidualized",
         Nsubj %in% c(16, 64, 256))
# Making the mixed models ----------------------
mixed_models_df <- tidyr::crossing(
  Bins = names(new_bins),
  Knum = seq_along(k_list),
  Nsubj =  c(8, 16),
  Nitems = c(8, 16),
  EffectSize = c(56, 80), # got rid of the 7*2^x
  Space = c("_log_space", ""),
  RType = c("residualized", "unresidualized")) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(iprefix = "same_across_subj",
         ofile = paste0(path_on_server, "amlap_saved_files/modeldata_",
                        iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
                        RType, "_size_", EffectSize, Space, "_", 
                        MinK, "_to_", MaxK, ".RDS")) %>%
  # filter(RType == "unresidualized") %>%
  # This is for the other stuff
                        { d <- filter(., Nsubj == Nitems, Nsubj == 8)
                        bind_rows(
                          .,
                          mutate(d, Nsubj=32, Nitems=32),
                          mutate(d, Nsubj=32, Nitems=16),
                          mutate(d, Nsubj=16, Nitems=32),
                          mutate(d, Nsubj=64, Nitems=64),
                          mutate(d, Nsubj=32, Nitems=64))
                        } %>%
  mutate(ofile = paste0(path_on_server, "amlap_saved_files/modeldata_",
                        iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
                        RType, "_size_", EffectSize, Space, "_", 
                        MinK, "_to_", MaxK, ".RDS"))
# ----- "Giant" mixed models
mixed_giant_df <- mixed_models_df %>%
  filter(Bins == "bin1", Space == "", RType == "unresidualized") %>%
  select(-RType, -Bins) %>%
  tidyr::crossing(Bins = c("giant", "threeword", "noexclusions"),
                  RType = c("residualized", "unresidualized")) %>%
  filter((RType == "unresidualized" & Bins %in% c("noexclusions", "giant")) |
         (RType == "residualized" & Bins %in% c("giant", "threeword"))) %>%
  mutate(ofile = paste0(path_on_server, "amlap_saved_files/modeldata_",
                        iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
                        RType, "_size_", EffectSize, Space, "_", 
                        MinK, "_to_", MaxK, ".RDS")) %>%
  filter(Bins == "giant")
# ----- Parametric mixed models
mixed_parametric_df <- mixed_models_df %>% 
  filter(RType == "unresidualized") %>%
  tidyr::crossing(Distribution = c("gaussian_para", "lognormal_para")) %>%
  select(-RType) %>% 
  tidyr::crossing(RType = c("unresidualized")) %>%
  mutate(ofile = paste0(path_on_server, "amlap_saved_files/modeldata_",
                        Distribution, "_", iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
                        RType, "_size_", EffectSize, Space, "_", 
                        MinK, "_to_", MaxK, ".RDS"))

# ----- The mixed models with interactions -----------------------------------
interaction_mixed_models_df <-  tidyr::crossing(
  Bins = c("bin1"), #since interaction is between block1 and block3
  Knum = seq_along(k_list),
  Nsubj =  c(32, 64),
  Nitems = c(32, 64),
  MainEffectSize = c(56),
  InteractionEffectSize = c(28),
  Block1Code = c(-1, 1),
  Interaction = c("real"), #   "fake"),
  Space = c("_log_space", ""),
  RType = c("unresidualized"),
  iprefix = c("same_across_subj_interaction")
) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(ofile = paste(
    paste0(path_on_server, "amlap_saved_files/modeldata"), 
    Interaction, iprefix,
    ifelse(Block1Code<0, "neg", "pos"),
    Nsubj, Nitems, Bins, RType, "size",
    MainEffectSize, InteractionEffectSize,
    sep = "_") %>%
      paste0(Space, "_", MinK, "_to_", MaxK, ".RDS"))






############## THE ACTUAL JOBS #######################################

# Running the normal models -------------------------------
bigger %beep% {
  tictoc::tic()
  df<-future_map(
    c(1:nrow(crossed_df)),
    function(i) {
      # Setup
      r <- crossed_df[i,]
      effect_size <- r$EffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==r$Bins, ]$meanRT,
        effect_size)
      
      if (r$RType=="residualized") {
        model_funcs <- residual_models
        # Establish adding functions
        if (r$Space=="")
          adder_f <- function(df) {
            mutate(df,
                   RT_Resid_with_Effect =         (RT + SimCond * effect_size) - PredRT_WL,
                   LogRT_Resid_with_Effect = log10(RT + SimCond * effect_size) - PredLogRT_WL)
          }
        else if (r$Space=="_log_space")
          adder_f <- function(df) {
            mutate(df,
                   RT_Resid_with_Effect =         add_in_log_space(RT, SimCond, log_effect) - PredRT_WL,
                   LogRT_Resid_with_Effect =log10(add_in_log_space(RT, SimCond, log_effect))- PredLogRT_WL)
          }
      }
      else if (r$RType=="unresidualized") {
        model_funcs <- simple_models
        # Establish adding functions
        if (r$Space=="")
          adder_f <- function(df) {
            mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
        else if  (r$Space=="_log_space")
          adder_f <- function(df) {
            mutate(df, RT_with_Effect = add_in_log_space(RT,SimCond,log_effect)) }
      }
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[r$Knum]],
        iprefix = paste0(r$iprefix, "_", r$Nsubj),
        isuffix = r$Bins,
        ofile_appended_name = paste(r$RType, "size", effect_size, sep="_") %>% paste0(r$Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = adder_f,
        data_tidyer = broom::tidy)
    })
  tictoc::tic()
  df
}
done(bigger)

## Saving the normal models

# Are all the files there?
file_checker_df %beep% {
  future_map_dfr(
    c(1:nrow(crossed_df)),
    function(i) {
      there <- FALSE
      if (file.exists(paste0(#path_on_server,
        #"amlap_saved_files/",
        crossed_df[i,]$ofilename)))
        there <- TRUE
      return(crossed_df[i,] %>% mutate(there=there))
    })
}
# crossed_df <- file_checker_df %>% filter(there==FALSE)

# Load it up
big_df %beep% {
  biggie <- future_map_dfr(
    c(1:nrow(crossed_df)),
    function(i) {
      info_row <- crossed_df[i,]
      
      results_path <- paste0(path_on_server, "amlap_saved_files/results/", 
                             sub(".*/","", info_row$ofilename))
      if (file.exists(results_path)) {
        df <- readRDS(results_path)
        if (min(df$k) == 1 & info_row$MinK != 1)
          df <- mutate(df, k = k + info_row$MinK - 1)
      } else {
        df <- collate_data(paste0(path_on_server, "amlap_saved_files/",
                                  sub(".*/","", info_row$ofilename)), 
                           .multicore = TRUE) %>%
          mutate(Nsubj = info_row$Nsubj,
                 SampleName = info_row$iprefix,
                 Bin = info_row$Bins,
                 EffectSize = info_row$EffectSize,
                 RType = info_row$RType,
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
  saveRDS(biggie,  paste0(path_on_server, "amlap_saved_files/results/regular_models.RDS"))
  paste0("rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/regular_models.RDS ~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/PowerSimulations/data/powerpaper/regular_models.RDS")
}
resolved(futureOf(big_df))

# Interaction stuff ---------------------------------------------------------

remaining_interaction_files <- remaining_files(interaction_crossed_df, file.exists(paste0("workspace/launchpad/amlap_saved_files/", ofilename)))

# Hasn't been tested wiht `iterate_over_df()` yet
interaction_stuff %beep% {
  iterate_over_df(
    remaining_interaction_files,
    future_map,
    function(zaza) {
      # Setup
      effect_size <- MainEffectSize
      interaction_effect <- InteractionEffectSize
      b1c <- Block1Code
      b3c <- -1 * b1c
      direction_name <- ifelse(b1c==-1,"neg", ifelse(b1c==1,"pos",stop("ERROR in direction")))
      
      effect_adder <- function(df) {
        mutate(df,
               # whether Block1 is coded as +1 or -1
               Block = factor(BlockName, levels=c("Block1","Block3")) %>% `contrasts<-`(value = c(b1c, b3c)),
               BlockDbl = case_when(Block=="Block1" ~ b1c, Block=="Block3" ~ b3c, TRUE ~ NA_real_),
               RT_with_ALL = RT + SimCond * effect_size + BlockDbl*SimCond*interaction_effect,
               RT_with_Interaction = RT + BlockDbl*SimCond*interaction_effect,
               RT_with_Main =  RT + SimCond * effect_size)}
      
      if (Interaction=="real")
        model_funcs <- interaction_models
      else if (Interaction=="fake")
        model_funcs <- fake_interaction_models
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[Knum]],
        iprefix = paste(Interaction, iprefix, direction_name, Nsubj, sep="_"),
        isuffix = "bin1",
        ofile_appended_name = paste("size", effect_size, interaction_effect, sep="_") %>% paste0(Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = effect_adder,
        data_tidyer = broom::tidy,
        ifilepath = paste0(path_on_server, "amlap_saved_files/bb_interaction") %>%
          paste(Nsubj, "bin1", min(k_list[[Knum]]), "to", max(k_list[[Knum]]), sep = "_") %>%
          paste0(".RDS"))
    })
}
done(interaction_stuff)

# loading them into one big thing
# hasn't been tested with the iterate_over_df funciton either
interaction_load %beep% {
  biggie <- iterate_over_df(
    interaction_crossed_df,
    future_map,
    function(zaza) {
      collate_data(paste0(path_on_server, "amlap_saved_files/", ofilename), 
                   .multicore = TRUE,
                   terms_to_harvest = c("SimCond", "Block1", "SimCond:Block1")) %>%
        mutate(Nsubj, # I believe this is the same thing as `Nsubj = Nsubj`
               Space, Block1Code, RType, Interaction, MainEffectSize,
               InteractionEffectSize)
    }) %>% bind_rows()
  saveRDS(biggie, paste0(path_on_server, "amlap_saved_files/results/interactions.RDS"))
  biggie
}
resolved(futureOf(interaction_load))

######################## parametric #########################



gaussian_parametric %beep% {
  future_walk(
    c(1:nrow(parametric_df)),
    function(i) {
      r <- parametric_df[i,]
      effect_size <- r$EffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins == r$Bins,]$meanRT,
        effect_size)
      
      if (r$Space == "")
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + effect_size * SimCond) }
      else if (r$Space == "_log_space")
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = add_in_log_space(RT, SimCond, log_effect)) }
      
      # happens before adding the sim condition/effect
      # The generated predictor is `Y`
      if (r$Distribution=="lognormal_para") {
        modifier_func <- function(df) { mutate(df, RT = 10^(Y)) }
        mu =  filter(setal_basis_info, Bins == r$Bins)$meanLogRT
        sig = filter(setal_basis_info, Bins == r$Bins)$sdLogRT
      }
      else if (r$Distribution=="gaussian_para") {
        modifier_func <- function(df) { mutate(df, RT = Y) }
        mu =  filter(setal_basis_info, Bins == r$Bins)$meanRT
        sig = filter(setal_basis_info, Bins == r$Bins)$sdRT
      }
      
      generate_parametric_list(
        info_row = r,
        iprefix = r$Distribution,
        effect_name = paste0("size_", effect_size, r$Space),
        server_path = path_on_server,
        effect_adder_func = adder_f,
        list_of_model_funcs = simple_models,
        # function stuff
        parametric_func = rnorm,
        df_modifier_function = modifier_func,
        n = r$Nsubj,
        mean = mu,
        sd = sig)
    })
}
resolved(futureOf(gaussian_parametric))

# Load it up
big_df %beep% {
  biggie <- future_map_dfr(
    c(1:nrow(parametric_df)),
    function(i) {
      info_row <- parametric_df[i,]
      
      results_path <- paste0(path_on_server, "amlap_saved_files/results/",
                             info_row$ofilename)
      # stop(results_path)
      if (file.exists(results_path)) {
        df <- readRDS(results_path)
        if (min(df$k) == 1 & info_row$MinK != 1)
          df <- mutate(df, k = k + info_row$MinK - 1)
      } else {
        df <- collate_data(paste0(path_on_server, "amlap_saved_files/",
                                  info_row$ofilename), 
                           .multicore = TRUE) %>%
          mutate(Nsubj = info_row$Nsubj,
                 SampleName = info_row$iprefix,
                 Bin = info_row$Bins,
                 EffectSize = info_row$EffectSize,
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
  saveRDS(biggie,  paste0(path_on_server, "amlap_saved_files/results/parametric_data.RDS"))
  paste0("rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/parametric_data.RDS ~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/PowerSimulations/data/powerpaper/parametric_data.RDS")
}
resolved(futureOf(big_df))


################################# BOX COX


boxcox %beep% {
  tictoc::tic()
  future_map(
    c(1:nrow(bc_df)),
    function(i) {
      # Setup
      r <- bc_df[i,]
      effect_size <- r$EffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==r$Bins, ]$meanRT,
        effect_size)
      
      if (r$RType=="residualized") {
        model_funcs <- residual_models
        # Establish adding functions
        if (r$Space=="")
          adder_f <- function(df) {
            mutate(df,
                   RT_Resid_with_Effect =         (RT + SimCond * effect_size) - PredRT_WL,
                   LogRT_Resid_with_Effect = log10(RT + SimCond * effect_size) - PredLogRT_WL)
          }
        else if (r$Space=="_log_space")
          adder_f <- function(df) {
            mutate(df,
                   RT_Resid_with_Effect =         add_in_log_space(RT, SimCond, log_effect) - PredRT_WL,
                   LogRT_Resid_with_Effect =log10(add_in_log_space(RT, SimCond, log_effect))- PredLogRT_WL)
          }
      }
      else if (r$RType=="unresidualized") {
        model_funcs <- simple_models
        # Establish adding functions
        if (r$Space=="")
          adder_f <- function(df) {
            mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
        else if  (r$Space=="_log_space")
          adder_f <- function(df) {
            mutate(df, RT_with_Effect = add_in_log_space(RT,SimCond,log_effect)) }
      }
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[r$Knum]],
        iprefix = paste0(r$iprefix, "_boxcox_", r$Nsubj),
        isuffix = r$Bins,
        ofile_appended_name = paste(r$RType, "size", effect_size, sep="_") %>% paste0(r$Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = adder_f,
        data_tidyer = boxcox_tidyer,
        ifilepath = paste0(path_on_server, "amlap_saved_files/bb_full_random") %>%
          paste(r$Nsubj, r$Bins, r$MinK, "to", r$MaxK, sep = "_") %>%
          paste0(".RDS"))
    })
  tictoc::toc()
}
done(boxcox)

boxcox_results %beep% {
  biggie <- future_map_dfr(
    c(1:1),
    function(i) {
      info_row <- bc_df[i,]
      
      results_path <- paste0(path_on_server, "amlap_saved_files/results/",
                             info_row$ofilename)
      # stop(results_path)
      if (file.exists(results_path)) {
        df <- readRDS(results_path)
        if (min(df$k) == 1 & info_row$MinK != 1)
          df <- mutate(df, k = k + info_row$MinK - 1)
      } else {
        df <- collate_data(paste0(path_on_server, "amlap_saved_files/",
                                  gsub("full_random", "full_random_boxcox", info_row$ofilename)), 
                           .multicore = TRUE) %>%
          mutate(Nsubj = info_row$Nsubj,
                 SampleName = info_row$iprefix,
                 Bin = info_row$Bins,
                 EffectSize = info_row$EffectSize,
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
  saveRDS(biggie,  paste0(path_on_server, "amlap_saved_files/results/boxcox.RDS"))
  paste0("rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/parametric_data.RDS ~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/PowerSimulations/data/powerpaper/parametric_data.RDS")
}
resolved(futureOf(boxcox_results))


###################################################################
## MIXED MODELS #################################################
#####################################################################

remaining_files_df <- remaining_files(mixed_models_df, file.exists(ofile))

bigger %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    remaining_files_df, 
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==Bins, ]$meanRT,
        effect_size)
      
      if (RType=="unresidualized") {
        model_funcs <- no_slopes
        # Establish adding functions
        if (Space=="")
          adder_f <- function(df) {
            mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
        else if  (Space=="_log_space")
          adder_f <- function(df) {
            mutate(df, RT_with_Effect = add_in_log_space(RT, SimCond, log_effect)) }
        
      } else if (RType=="residualized") {
        model_funcs <- residual_no_slopes
        # Establish adding functions
        if (Space=="")
          adder_f <- function(df) {
            mutate(df,
                   RT_Resid_with_Effect =         (RT + SimCond * effect_size) - PredRT_WL,
                   LogRT_Resid_with_Effect = log10(RT + SimCond * effect_size) - PredLogRT_WL) }
        else if (Space=="_log_space")
          adder_f <- function(df) {
            mutate(df,
                   RT_Resid_with_Effect =         add_in_log_space(RT, SimCond, log_effect) - PredRT_WL,
                   LogRT_Resid_with_Effect =log10(add_in_log_space(RT, SimCond, log_effect))- PredLogRT_WL) }
      }
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[Knum]],
        iprefix = paste0(iprefix, "_",
                         Nsubj, "_",
                         Nitems),
        isuffix = Bins,
        ofile_appended_name = paste(RType, "size", effect_size, sep="_") %>% paste0(Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = adder_f,
        data_tidyer = mixed_cleaner)
    })
  tictoc::tic()
  l
}
done(bigger)

load_mixed_models %beep% {
  iterate_over_df(
    mixed_models_df,
    future_map_dfr,
    function(i) {
      collate_data(ofile, 
                   # terms_to_harvest = c("sd_(Intercept).UniqueItem", 
                   #                      "sd_(Intercept).UniqueSubject"),
                   .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space)
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/mixed_models_lmerTest.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
}
done(load_mixed_models)



###################################################################
## PARAMETRIC MIXED MODELS #########################################
#####################################################################

# mixed_parametric_df <- filter(mixed_parametric_df, Nsubj %in% c(8,16), Nitems %in% c(8, 16, 32), Distribution == "gaussian_para", Space == "", RType == "noexclusions")

para_remaining_files_df <- remaining_files(mixed_parametric_df, file.exists(ofile))

bigger %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    para_remaining_files_df,
    # mixed_parametric_df, 
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==Bins, ]$meanRT,
        effect_size)
      
      if (RType %in% c("unresidualized", "noexclusions")) {
        model_funcs <- no_slopes
        # Establish adding functions
        if (Distribution == "gaussian_para") {
          inverse_f = identity
          if (Space=="")
            adder_f <- function(df) {
              mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
          else if  (Space=="_log_space")
            adder_f <- function(df) {
              mutate(df, RT_with_Effect = add_in_log_space(RT, SimCond, log_effect)) }
        } else if (Distribution == "lognormal_para") {
          inverse_f = function(x) 10^x
          if (Space=="")
            adder_f <- function(df) {
              mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
          else if  (Space=="_log_space")
            adder_f <- function(df) {
              mutate(df, RT_with_Effect = add_in_log_space(RT, SimCond, log_effect)) }
        }
      }
      if (RType == "unresidualized")
        has_exclusions <- TRUE
      else if (RType == "noexclusions")
        has_exclusions <- FALSE
      
      run_parametrics_from_model(
        server_path = path_on_server,
        model_file = paste0("amlap_saved_files/parametric_models/", Distribution, "_", Bins, ".RDS"),
        k = k_list[[Knum]],
        iprefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        isuffix = Bins,
        ofile_appended_name = paste(RType, "size", effect_size, sep="_") %>% paste0(Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = adder_f,
        data_tidyer = mixed_cleaner, #mixed_cleaner_Kenward
        ifilepath = ofile,
        backtrans_f = inverse_f,
        exclusions = has_exclusions)
      
    })
  tictoc::tic()
  l
}
done(bigger) 

big_df %beep% {
  biggie <- iterate_over_df(
    mixed_parametric_df,
    future_map,
    function(zaza) {
      results_path <- gsub("amlap_saved_files","amlap_saved_files/results", ofile)
      if (file.exists(results_path)) {
        df <- readRDS(results_path)
      } else {
        df <- collate_data(ofile, .multicore = TRUE) %>%
          mutate(Nsubj, Nitems, EffectSize, Space, RType, MinK,
                 Bin = Bins, SampleName = iprefix, Distribution)
        saveRDS(df, results_path)
      }
      df
    }) %>% bind_rows()
  
  saveRDS(biggie,  paste0(path_on_server, "amlap_saved_files/results/parametric_mixed_models_lmerTest.RDS"))
  cat("rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/parametric_mixed_models_lmerTest.RDS ~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/PowerSimulations/data/powerpaper/")
}
done(big_df)


#######################################################################################
######################## "GIANT" MIXED MODELS ########################################
######################################################################################

giant_mm_remaining_files_df <- remaining_files(mixed_giant_df, file.exists(ofile))

giant_maker %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    giant_mm_remaining_files_df, 
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==Bins, ]$meanRT,
        effect_size)
      
      if (RType=="unresidualized") {
        model_funcs <- no_slopes
        # Establish adding functions
        if (Space=="")
          adder_f <- function(df) {
            mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
        else if  (Space=="_log_space")
          stop("Log-space not implemented!")
          # adder_f <- function(df) {
            # mutate(df, RT_with_Effect = add_in_log_space(RT, SimCond, log_effect)) }
      } else if (RType=="residualized") {
        model_funcs <- residual_no_slopes
        # Establish adding functions
        if (Space=="")
          adder_f <- function(df) {
            mutate(df,
                   RT_Resid_with_Effect =         (RT + SimCond * effect_size) - PredRT_WL,
                   LogRT_Resid_with_Effect = log10(RT + SimCond * effect_size) - PredLogRT_WL) }
        else if (Space=="_log_space")
          stop("Log-space not implemented!")
          # adder_f <- function(df) {
          #   mutate(df,
          #          RT_Resid_with_Effect =         add_in_log_space(RT, SimCond, log_effect) - PredRT_WL,
          #          LogRT_Resid_with_Effect =log10(add_in_log_space(RT, SimCond, log_effect))- PredLogRT_WL) }
      }
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[Knum]],
        iprefix = paste0(iprefix, "_",
                         Nsubj, "_",
                         Nitems),
        isuffix = Bins,
        ofile_appended_name = paste(RType, "size", effect_size, sep="_") %>% paste0(Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = adder_f,
        data_tidyer = mixed_cleaner)
    })
  message(tictoc::tic())
  l
}
done(giant_maker)



load_giant_models %beep% {
  iterate_over_df(
    mixed_giant_df,
    future_map_dfr,
    function(i) {
      collate_data(ofile,
                   # terms_to_harvest = c("sd_(Intercept).UniqueItem",
                   #                      "sd_(Intercept).UniqueSubject"),
                   .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space)
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/giant_mixed_models_lmerTest.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
}
done(load_giant_models)




#######################################################################################
############################ Interaction mixed models #################################
#######################################################################################
i_mm_remaining_files_df <- remaining_files(interaction_mixed_models_df, file.exists(ofile))

bigger %beep% {
  l <- iterate_over_df(
    i_mm_remaining_files_df, 
    future_cnd_map,
    function(zaa) {
      # Set up th effect size and log_effectsize
      effect_size <- MainEffectSize
      interaction_effect <- InteractionEffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==Bins, ]$meanRT,
        effect_size)
      log_interaction_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==Bins, ]$meanRT,
        interaction_effect)
      # Set up the interaction direction
      direction_name <- ifelse(Block1Code==-1,"neg", ifelse(Block1Code==1,"pos", stop("ERROR in direction")))
      
      if (RType=="unresidualized") {
        model_funcs <- no_slopes
        # Establish adding functions
        if (Space=="")
          effect_adder <- function(df) {
            mutate(df, # whether Block1 is coded as +1 or -1
                   Block = factor(BlockName, levels=c("Block1","Block3")) %>% `contrasts<-`(value = c(Block1Code, -Block1Code)),
                   BlockDbl = case_when(Block=="Block1" ~ Block1Code, Block=="Block3" ~ -Block1Code, TRUE ~ NA_real_),
                   RT_with_ALL = RT + SimCond * effect_size + BlockDbl*SimCond*interaction_effect,
                   RT_with_Interaction = RT + BlockDbl*SimCond*interaction_effect,
                   RT_with_Main =  RT + SimCond * effect_size)}
        else if  (Space=="_log_space")
          effect_adder <- function(df) {
            mutate(df, # whether Block1 is coded as +1 or -1
                   Block = factor(BlockName, levels=c("Block1","Block3")) %>% `contrasts<-`(value = c(Block1Code, -Block1Code)),
                   BlockDbl = case_when(Block=="Block1" ~ Block1Code, Block=="Block3" ~ -Block1Code, TRUE ~ NA_real_),
                   RT_with_ALL = 10^(log10(RT) + SimCond * log_effect + BlockDbl*SimCond*log_interaction_effect),
                   RT_with_Interaction = 10^(log10(RT) + BlockDbl*SimCond*log_interaction_effect),
                   RT_with_Main =  10^(log10(RT) + SimCond * log_effect))}
      }
      
      if (Interaction=="real")
        model_funcs <- i_mm_block_slopes
      else if (Interaction=="fake")
        stop("Not coded yet!")
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = k_list[[Knum]],
        # iprefix = paste(Interaction, iprefix, direction_name, Nsubj, sep="_")
        iprefix = paste(Interaction, iprefix, direction_name, Nsubj, Nitems, sep="_"),
        isuffix = "bin1",
        ofile_appended_name = paste(RType, "size", effect_size, 
                                    interaction_effect,  sep="_") %>% paste0(Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = effect_adder,
        data_tidyer = mixed_cleaner,
        ifilepath = paste0(path_on_server, "amlap_saved_files/bb_same_across_subj_interaction") %>%
          paste(Nsubj, Nitems, "bin1", min(k_list[[Knum]]), "to", max(k_list[[Knum]]), sep = "_") %>%
          paste0(".RDS"))
    })
  l
}
done(bigger)

load_i_mms %beep% {
  iterate_over_df(
    interaction_mixed_models_df,
    future_map_dfr,
    function(i) {
      collate_data(ofile, .multicore = TRUE,
                   terms_to_harvest = c("SimCond", "Block1", "SimCond:Block1")) %>%
        mutate(SampleName = iprefix,
               Nsubj, Nitems, Bins, MainEffectSize,
               InteractionEffectSize, Interaction,
               Block1Code, RType, Space)
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/interaction_mixed_models_lmerTest.RDS"))
  "done!"
}
done(load_i_mms)


#############################################################################

# d <- readRDS("~/Box Sync/Power simulations for RTs/PowerSimulations/data/mixed_models_lmerTest.RDS")
# head(d)
# d %>% filter(!grepl("converge",warnings)) %>%
#   filter(Nsubj == Nitems) %>%
#   group_by(type, Bin, Nsubj, Nitems) %>%
#   summarise(n=n(),
#             pos = mean(ifelse(statistic > 1.96, 1, 0)),
#             pos2 = mean(ifelse(statistic > 0 & p.value < 0.05, 1, 0)),
#             # neg = mean(ifelse(statistic < -1.96, 1, 0)),
#             all = mean(ifelse(abs(statistic) > 1.96, 1, 0)),
#             all2 =  mean(ifelse(p.value < 0.05, 1, 0))) %>%
#   tidyr::separate(type, c("Analysis", "Type"), sep="_") %>%
#   group_by(Analysis, Bin, Nsubj, Nitems) %>%
#   summarise(Power = first(pos2[Type=="power"]),
#             `Type I` = first(all2[Type=="type1"]),
#             avg_n = mean(n),
#             ndiff = abs(n[1]-n[2])) %>%
#   arrange(-Nsubj, -Nitems, Bin, Analysis) %>%
#   group_by(Bin, Nsubj, Nitems) %>%
#   mutate(avg_n=mean(avg_n)) %>%
#   tidyr::unite(Temp,Power,`Type I`)





# the following is rank deficient and drops the interaction term:
#(Nsubj== 32, Nitems ==32, Space=="", Block1Code == -1, k == 3313)
# a few others do this as well

hmmm %beep% {
  readRDS("workspace/launchpad/amlap_saved_files/results/interaction_mixed_models_lmerTest.RDS") %>%
    as_tibble()
}

summarized_df <- hmmm %>% 
  filter(!(messages != "" | grepl("converg", warnings) | errors != "")) %>%
  group_by(type, Block1Code, Interaction, term, RType, Nsubj, Nitems, Space,
           MainEffectSize, InteractionEffectSize) %>%
  summarise(all = mean(ifelse(p.value < 0.05, 1, 0)),
            pos = mean(ifelse(p.value < 0.05 & estimate > 0, 1, 0)),
            neg = mean(ifelse(p.value < 0.05 & estimate < 0, 1, 0)),
            mean = mean(estimate),
            n=n()) %>%
  ungroup() %>%
  tidyr::separate(type, into=c("Analysis","Type","Temp"),
                  fill = "right") %>%
  mutate(Type = ifelse(Type=="fake",Temp,Type)) %>% 
  select(-Temp) %>%
  {counter_c <<- 0; .} %>% # this is just a weird system to throw 
  # an error if there are cases that don't fit in the case_when
  mutate(term = case_when(
    term == "SimCond" ~ "MainEffect",
    term == "Block1" ~ "Block",
    term == "SimCond:Block1" ~ "Interaction",
    TRUE ~ { counter_c <<- counter_c+1; stopifnot(counter_c<2); NA_character_})) %>% 
    {counter_c <<- 0; .} %>%
  mutate(right = case_when(
    term %in% c("MainEffect", "Interaction") ~ pos,
    term == "Block" & Block1Code < 0 ~ neg,
    term == "Block" & Block1Code > 0 ~ pos,
    TRUE ~ { counter_c <<- counter_c+1; stopifnot(counter_c<2); NA_real_})) %>%
    {counter_c <<- 0; .} %>%
  mutate(Space = case_when(
    Space == "" ~ "RawRTs",
    Space == "_log_space" ~ "log(RTs)",
    TRUE ~ { counter_c <<- counter_c+1; stopifnot(counter_c<2); NA_character_}))

# Compared to itself, log analysis sucks at finding certain interactions
summarized_df %>% 
  filter(Interaction=="real",
         Type=="all",
         term=="Interaction",
         Nsubj==64, Nitems==64
  ) %>%
  select(Analysis, Space, Block1Code, mean:right) %>%
  arrange(Analysis, Block1Code)


PandT <- summarized_df %>% 
  mutate(PorT = case_when(
    term == "Interaction" ~ case_when(
      Type %in% c("all", "interaction") ~ "power",
      Type %in% c("main", "none") ~ "type1",
      TRUE ~ NA_character_),
    term == "MainEffect" ~ case_when(
      Type %in% c("all", "main") ~ "power",
      Type %in% c("interaction", "none") ~ "type1",
      TRUE ~ NA_character_),
    term == "Block" ~ "power", # block is always
    TRUE ~ NA_character_))  %>%
  mutate(Stat = case_when(
    PorT == "power" ~ right,
    PorT == "type1" ~ all,
    TRUE ~ NA_real_)) %>%
  select(term, PorT, Stat, Analysis:InteractionEffectSize)



z_colorer <- function(x, named_color_list) {
  beginning = "<span style=\"display: block; padding: 0 4px; border-radius: 4px; text-align: center; background-color: "
  end <- "</span>"
  paste0(beginning, named_color_list[x], "\">", x, end)
}

z_barer <- function(x, named_color_list="NA",
                    default_color="lightgreen",
                    fun = formattable::proportion) {
  beginning = "<span style=\"display: inline-block; direction: rtl; border-radius: 4px; padding-right: 2px; background-color: "
  colors <- ifelse(is.na(named_color_list[x]),
                   default_color, named_color_list[x])
  widths <- paste0("width: ", formattable::percent(fun(as.numeric(x))))
  end <- "</span>"
  paste0(beginning, 
         colors, "; ", 
         widths, "\">", x, end)
}



PandT %>%
  filter(Interaction=="real") %>%
  filter(Block1Code==-1) %>%
  filter(Nsubj==64, Nitems==64) %>%
  filter(Space=="RawRTs") %>%
  # filter(term !="Block") %>%
  select(Analysis, Type ,Stat, term, PorT) %>%
  tidyr::unite(New, Type, term, PorT) %>%
  tidyr::spread(New, Stat) %>% 
  mutate_if(is.numeric, ~format(., digits=3)) %>%
  # This is where the magic starts happening
  mutate_at(vars(everything(),-Analysis),
            ~z_barer(.)) %>%
  set_names(c("", 
              names(.)[2:length(names(.))] %>%
              { actual_names <<- gsub("_[[:alnum:]]*$","", 
                                      gsub("^[[:alnum:]]*_", "", .))
              gsub("^.*_", "", .) } %>%
                map_chr(~ifelse(.=="power","Power","Type I")) %>%
                z_colorer(c("Power"="#E0E0E0", "Type I"="WhiteSmoke"))
  )) %>% 
  knitr::kable("html", 
               caption="Main effect = 56 ms, interaction with Block1 vs. Block3 = 28 ms.",
               booktabs = T,
               escape = F,
               align = "l") %>%
  kableExtra::add_header_above(c(" ", actual_names)) %>%
  # row_spec(0, background = "red")  %>%
  kableExtra::add_header_above(c(" ", 
                                 "With Both Effects" = n_distinct(actual_names), 
                                 "Without Main Effect" = n_distinct(actual_names),
                                 "Without Interaction" = n_distinct(actual_names),
                                 "Without Either Effects" = n_distinct(actual_names)))  %>%
  kableExtra::kable_styling("striped") %>%
  kableExtra::row_spec(2,bold = T)





