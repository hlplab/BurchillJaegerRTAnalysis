# So this doesn't rerun EVERYTHING, I realize. It only reruns the MODELS, not the "bb" functions

library(dplyr) # originally used dplyr_0.7.4
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

################################### killling jobs
kill_r_on_nodes(ok_nodes, backup_login, # "zburchil@cycle3.cs.rochester.edu"
                wait_until_resolved = TRUE,
                beep_on_end = TRUE)


# ############################## PLAN ###################################################################
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = ok_nodes),
  tweak(multiprocess, workers = function() {max(round(future::availableCores()/2), 4)})
))

####################### GETTING THE JOB DATA FRAMES ######################################

# Making the mixed models ----------------------
mixed_models_df <- tidyr::crossing(
  Bins = names(new_bins),
  Knum = seq_along(k_list),
  Nsubj =  c(8, 16),
  Nitems = c(8, 16),
  EffectSize = c(56, 80), # got rid of the 7*2^x
  Space = c("_log_space", ""),
  RType = c("unresidualized")) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]])),
         iprefix = "same_across_subj") %>%
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

# The residualized mixed model dfs are not really "based" off the default mixed model df
#################################################
resid_mixed_models_df <- tidyr::crossing(FillerItemRatio = 15, 
                                         MinK = 1, MaxK = 1,
                                         iprefix = "same_across_subj",
                                         Nsubj = c(8, 16, 32, 64),
                                         Bins = "giant", 
                                         RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         BBitems = (FillerItemRatio + 1) * Nitems) %>%
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(BBitems) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  arrange(BBitems, MinK) %>%
  # Turning it into the model df
  tidyr::crossing(EffectSize = c(15, 35),
                  Space = c("", "_log_space")) %>%
  mutate(
    ifile = paste0(
      path_on_server, "amlap_saved_files/bb_",
      iprefix, "_",  Nsubj, "_", BBitems, "_", Bins, "_",
      MinK, "_to_", MaxK, ".RDS"),
    ofile = paste0(
      path_on_server, "amlap_saved_files/modeldata_",
      iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
      RType, "_size_", EffectSize, Space, "_",
      MinK, "_to_", MaxK, ".RDS")
  )

region_mixed_models_df <- tidyr::crossing(FillerItemRatio = 15, 
                                          MinK = 1, MaxK = 1,
                                          iprefix = "same_across_subj_regionized",
                                          Nsubj = c(8, 16, 32, 64),
                                          Bins = "giant", 
                                          RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         BBitems = FillerItemRatio * Nitems) %>%
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(BBitems) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  arrange(BBitems, MinK) %>%
  # Turning it into the model df
  tidyr::crossing(EffectSize = c(15, 35),
                  Space = c("", "_log_space")) %>%
  mutate(
    ifile = paste0(
      path_on_server, "amlap_saved_files/bb_",
      iprefix, "_",  Nsubj, "_", Nitems, "_", BBitems, "_fillers_", Bins, "_",
      MinK, "_to_", MaxK, ".RDS"),
    ofile = paste0(
      path_on_server, "amlap_saved_files/modeldata_",
      iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
      RType, "_size_", EffectSize, Space, "_",
      MinK, "_to_", MaxK, ".RDS")
  )
################################################

# ----- "Giant" mixed models
mixed_giant_df <- mixed_models_df %>%
  filter(Bins == "bin1", Space == "", RType == "unresidualized") %>%
  select(-Bins, -iprefix) %>%
  tidyr::crossing(Bins = c("giant"),
                  iprefix = c("same_across_subj", "same_across_subj_noexcl")) %>%
  mutate(ofile = paste0(path_on_server, "amlap_saved_files/modeldata_",
                        iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
                        RType, "_size_", EffectSize, Space, "_",
                        MinK, "_to_", MaxK, ".RDS"))
# ----- Parametric mixed models
mixed_parametric_df <- mixed_models_df %>%
  filter(RType == "unresidualized") %>%
  tidyr::crossing(Distribution = c("gaussian_para", "lognormal_para")) %>%
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
} # %plan% elite_squad_plan
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

#############################################################
### RESIDUALIZED MIXED MODELS ##############################
#############################################################


resid_remaining_files_df <- remaining_files(resid_mixed_models_df, file.exists(ofile))

resid_mm_l %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    resid_remaining_files_df,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      # This is different because there aren't Bins for this one!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      log_effect <- log_space_effect_size_calculator(
        mean(setal_basis_info$meanRT), effect_size)
      
      if (RType %in% c("residualized")) {
        # Doesn't double-dip subjects, since they've been residualized out
        model_funcs <- safe_residual_no_slopes_no_subject
        
        # Establish adding functions
        if (Space=="")
          adder_f <- function(df) {
            set_filler_data(df, n_remaining_items = Nitems) %>%
              mutate(logRT = log10(RT),
                     RT_with_Effect =         (RT + SimCond * effect_size),
                     LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
              amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
              amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
              amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
              amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
              remove_filler_data(n_remaining_items = Nitems) %>%
              mutate(
                RT_Resid =           RT - BATA_PredRawRT,
                LogRT_Resid = log10(RT) - BATA_PredLogRT,
                # RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT,
                # LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT
                RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT_with_Effect,
                LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT_with_Effect
              )
          }
        else if (Space=="_log_space")
          adder_f <- function(df) {
            set_filler_data(df, n_remaining_items = Nitems) %>%
              mutate(logRT = log10(RT),
                     RT_with_Effect =         add_in_log_space(RT, SimCond, log_effect),
                     LogRT_with_Effect =log10(add_in_log_space(RT, SimCond, log_effect))) %>%
              amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
              amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
              amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
              amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
              remove_filler_data(n_remaining_items = Nitems) %>%
              mutate(
                RT_Resid =           RT - BATA_PredRawRT,
                LogRT_Resid = log10(RT) - BATA_PredLogRT,
                # RT_Resid_with_Effect =         RT_with_Effect - BATA_PredRawRT,
                # LogRT_Resid_with_Effect =   LogRT_with_Effect - BATA_PredLogRT
                RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT_with_Effect,
                LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT_with_Effect
              )
          }
      } else {
        stop(Space, " as a value for `Space` is not implemented")
      }
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = MinK:MaxK,
        iprefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        isuffix = Bins,
        ofile_appended_name = paste(RType, "size", effect_size, sep="_") %>% paste0(Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = adder_f,
        data_tidyer = mixed_cleaner,
        ifilepath = ifile)
    })
  tictoc::tic()
  l
} #%plan% elite_squad_plan
done(resid_mm_l)

load_resid_mixed_models %beep% {
  iterate_over_df(
    resid_mixed_models_df,
    future_map_dfr,
    function(i) {
      collate_data(ofile, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               FillerItemRatio = FillerItemRatio,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space)
    }) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/residualized_mixed_models_lmerTest.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} 
done(load_resid_mixed_models)


#############################################################
### REGIONIZED MIXED MODELS ##############################
#############################################################

region_remaining_files_df <- remaining_files(region_mixed_models_df, file.exists(ofile))
  
region_mm_l %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    region_remaining_files_df,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      # This is different because there aren't Bins for this one!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      log_effect <- log_space_effect_size_calculator(
        mean(setal_basis_info$meanRT), effect_size)
      
      if (RType %in% c("residualized")) {
        # Doesn't double-dip subjects, since they've been residualized out
        model_funcs <- safe_residual_no_slopes_no_subject
        
        # Establish adding functions
        if (Space=="")
          adder_f <- function(df) {
            mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
              mutate(logRT = log10(RT),
                     RT_with_Effect =         (RT + SimCond * effect_size),
                     LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
              amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
              amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
              amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
              amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
              # Filters out fillers
              filter(ItemType == "CriticalRegion") %>%
              mutate(UniqueItem = factor(UniqueItem)) %>%
              # Calculates residuals BEFORE averaging over the regions
              mutate(
                RT_Resid =           RT - BATA_PredRawRT,
                LogRT_Resid = log10(RT) - BATA_PredLogRT,
                # RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT,
                # LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT
                RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT_with_Effect,
                LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT_with_Effect
              ) %>%
              # Really important to have this right!!!!!!!!!
              group_by(Item, UniqueItem, Subject, UniqueSubject, SimCond, ItemAmbValue, SubjectAmbValue) %>%
              summarise(RT_Resid = mean(RT_Resid, na.rm = TRUE),
                        LogRT_Resid = mean(LogRT_Resid, na.rm = TRUE),
                        RT_Resid_with_Effect = mean(RT_Resid_with_Effect, na.rm = TRUE),
                        LogRT_Resid_with_Effect = mean(LogRT_Resid_with_Effect, na.rm = TRUE),
                        AllRTs = paste(RT, collapse=", "),
                        RT = mean(RT, na.rm=TRUE)) %>%
              ungroup()
          }
        else if (Space=="_log_space")
          adder_f <- function(df) {
            mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
              mutate(logRT = log10(RT),
                     RT_with_Effect =         add_in_log_space(RT, SimCond, log_effect),
                     LogRT_with_Effect =log10(add_in_log_space(RT, SimCond, log_effect))) %>%
              amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
              amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
              amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
              amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
              # Filters out fillers
              filter(ItemType == "CriticalRegion") %>%
              mutate(UniqueItem = factor(UniqueItem)) %>%
              # Calculates residuals BEFORE averaging over the regions
              mutate(
                RT_Resid =           RT - BATA_PredRawRT,
                LogRT_Resid = log10(RT) - BATA_PredLogRT,
                # RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT,
                # LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT
                RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT_with_Effect,
                LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT_with_Effect
              ) %>%
              # Really important to have this right!!!!!!!!!
              group_by(Item, UniqueItem, Subject, UniqueSubject, Words.Start, SimCond, ItemAmbValue, SubjectAmbValue) %>%
              summarise(RT_Resid = mean(RT_Resid, na.rm = TRUE),
                        LogRT_Resid = mean(LogRT_Resid, na.rm = TRUE),
                        RT_Resid_with_Effect = mean(RT_Resid_with_Effect, na.rm = TRUE),
                        LogRT_Resid_with_Effect = mean(LogRT_Resid_with_Effect, na.rm = TRUE),
                        AllRTs = paste(RT, collapse=", "),
                        RT = mean(RT, na.rm=TRUE)) %>%
              ungroup()
          }
      } else {
        stop(Space, " as a value for `Space` is not implemented")
      }
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = setal_file,
        k = MinK:MaxK,
        iprefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        isuffix = Bins,
        ofile_appended_name = paste(RType, "size", effect_size, sep="_") %>% paste0(Space),
        list_of_model_funcs = model_funcs,
        effect_adder_function = adder_f,
        data_tidyer = mixed_cleaner,
        ifilepath = ifile)
    })
  tictoc::tic()
  l
} #%plan% b_squad_plan
done(region_mm_l)

load_region_mixed_models %beep% {
  iterate_over_df(
    region_mixed_models_df,
    future_map_dfr,
    function(i) {
      collate_data(ofile, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               FillerItemRatio = FillerItemRatio,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space)
    }) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/regionized_mixed_models_lmerTest.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} 
done(load_region_mixed_models)





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
                 Bin = Bins, SampleName = iprefix, Distribution, ModelFile = ofile)
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
      # log_effect <- log_space_effect_size_calculator(
      #   setal_basis_info[setal_basis_info$Bins==Bins, ]$meanRT,
      #   effect_size)
      
      model_funcs <- no_slopes
      if (Space=="")
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      else 
        stop("Log-space not implemented!")
      
      
      if (!grepl("noexcl", iprefix)) {
        temp_setal_file <- setal_file
      } else if (grepl("noexcl", iprefix)) {
        temp_setal_file <- "setal_noexcl_fillers.RDS"
      } else {
        stop("Uh oh")
      }
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = temp_setal_file,
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
} #%plan% elite_squad_plan
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
               Space = Space,
               ModelFile = ofile)
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


######################################################################
########## FETAL ######################################################
######################################################################
mixed_fetal_df <- mixed_giant_df %>%
  mutate(iprefix = paste0("fetal_", iprefix)) %>%
  filter(Bins == "giant") %>%
  mutate(ofile = paste0(path_on_server, "amlap_saved_files/modeldata_",
                        iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
                        RType, "_size_", EffectSize, Space, "_",
                        MinK, "_to_", MaxK, ".RDS"))

fetal_mm_remaining_files_df <- remaining_files(mixed_fetal_df, file.exists(ofile))

fetal_maker %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    fetal_mm_remaining_files_df,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins==Bins, ]$meanRT,
        effect_size)
      
      model_funcs <- no_slopes
      if (Space=="")
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      else 
        stop("Log-space not implemented!")
      
      if (Bins=="giant") {
        
      } else {
        stop("Not implemented")
      }
      
      amlap_model_coordinator(
        server_path = path_on_server,
        real_df_file = "fetal_sim_df.RDS",
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
done(fetal_maker)


load_fetal_models %beep% {
  iterate_over_df(
    mixed_fetal_df,
    future_map_dfr,
    function(i) {
      collate_data(ofile, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space)
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/fetal_mixed_models_lmerTest.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
}
done(load_fetal_models)









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





