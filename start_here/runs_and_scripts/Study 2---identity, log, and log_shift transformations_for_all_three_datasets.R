# Loading all the boilerplate. CHANGE PATH TO YOUR OWN DIRECTORY
main_path = "/Users/tiflo/Library/CloudStorage/Box-Box/_Papers - Box/Power simulations for RTs/PowerSimulations/BurchillJaegerRTAnalysis/start_here/"
# main_path = "/Users/zburchill/Box Sync/Power simulations for RTs/PowerSimulations/BurchillJaegerRTAnalysis/start_here/"

# Run this to load all the right libraries and get all the right constants for the scripts
source(paste0(main_path, "functions/boilerplate.R"))


# Make the future plans ----------------------------------
plan(list(
  tweak(multisession, workers = 20)
))


# Save all the models and load in constants -----------------------------------
source(paste0(code_path, "saving_parametric_models.R"))


# Used in all the log-shift effect-adding functions
min_shift_setter <- function(df) {
  min_val = min(df$RT_with_Effect, na.rm = TRUE) - 1 # since log(0) = -Inf
  min_val2 = min(df$RT, na.rm = TRUE) - 1 # since log(0) = -Inf
  df %>% mutate(PowerShift = min_val,
                TypeIShift = min_val2)
}


# Attention: because the NSC data is different enough from the HS18 and F13 data, this script separates them into different runs. We first cover the HS18 and F13, then NSC.
# It also breaks things up into pre-residualized data, then residualized, then 3-word regions


# -----------------------------------------------------

# Make the big df you're going to use
sfetal_giant_mm_df <- tidyr::crossing(
  Bins = "giant",
  Knum = seq_along(k_list),
  Nitems = c(8, 16, 32, 64),
  Nsubj = c(8, 16, 32, 64),
  # iprefix = c("same_across_subj", "same_across_subj_noexcl"),
  iprefix = c("same_across_subj_new", "same_across_subj_noexcl_new"),
  EffectSize = c(56, 80), # got rid of the 7*2^x
  Space = c(""),
  RType = c("unresidualized"),
  Dataset = c("setal","fetal")) %>%
  filter(Nsubj == Nitems) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(
    # Add in redundancy
    iprefix = paste0(Dataset,"_",iprefix),
    # Are there exclusions?
    Exclusions = ifelse(grepl("noexcl", iprefix), FALSE, TRUE),
    exclname = ifelse(Exclusions,"excl","noexcl"),
    # Barebone data files: stores indices of sampled data for each generated
    # data set (BATA)
    bb_file = paste0(bb_files_path, "bb_", iprefix, "_",
                     Nsubj, "_", Nitems, "_",
                     Bins, "_", MinK, "_to_", MaxK, ".RDS"),
    # Model summary file: stores all information from each mixed-effect analysis
    # (e.g., direction of effect, p-value, etc.)
    mm_file = paste0(mm_files_path, "modeldata_",
                     iprefix, "_", Nsubj, "_", Nitems, "_",
                     Bins, "_",
                     RType, "_size_", EffectSize, Space, "_",
                     MinK, "_to_", MaxK, ".RDS"),
    # Source data file: the indices in the barebones data file refer to this
    # data set
    real_df_file = map2_chr(Dataset, exclname, ~og_data_paths[[.x]][[.y]]))

sfetal_giant_bb_df <- sfetal_giant_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_sfetal_giant_bb_files <- remaining_files(sfetal_giant_bb_df, file.exists(bb_file))


# Save BB files
sfetal_giant_bb <- {
  iterate_over_df(
    remaining_sfetal_giant_bb_files,
    future_cnd_map,
    function(zaza) {
      k_e <- MinK:MaxK
      filterers <- quos(Group == "Filler-first")

      zplyr::collect_all({
        real_df <- readRDS(real_df_file)
        l <- make_simple_sample_list(
          df = real_df, k = k_e,
          filter_quosures = filterers,
          # The function and its named additional arguments
          single_sample_df_function = same_items_per_subj,
          items_per_subj = Nitems,
          n_subj = Nsubj)
        saveRDS(l, bb_file)
        "good!"
      },
      catchErrors = TRUE)
    })
}

remaining_sfetal_giant_mm_files <- remaining_files(sfetal_giant_mm_df, file.exists(mm_file))
# Saves MM files
sfetal_giant_mm <- {
  iterate_over_df(
    remaining_sfetal_giant_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      
      model_funcs <- no_slopes_wo_infs_plus_logshift
      
      adder_f <- function(df) {
        mutate(df, RT_with_Effect = RT + SimCond * effect_size) %>%
          min_shift_setter()
      }
      
      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
}

# Save MM DFs -----------
sfetal_giant_load <- {
  iterate_over_df(
    sfetal_giant_mm_df,
    map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               Dataset = Dataset,
               Exclusions = Exclusions,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               MinK = MinK)
    }
  ) %>%
    saveRDS(paste0(complete_path, "logshift_sfetal_v01.RDS"))
}
sfetal_giant_load


######################################################################################################
# Doing the residualized models ######################################################################
######################################################################################################

# This data frame is for the residualized BATAs. Due to the increased size of the BATAs (because they have ~15x more data due to the fillers) and due to the fact that many more models are run per BATA (we have residualize so many different ways), we break down the number of BATAs per file. Feel free to adjust the `residual_total_files` variable to whatever you see fit.
sfetal_resid_mm_df <- tidyr::crossing(
  FillerItemRatio = 15,
  MinK = 1, MaxK = 1,
  iprefix = "same_across_subj_resid_new",
  Nsubj = c(8, 16, 32, 64),
  Bins = "giant",
  Dataset = c("setal","fetal"),
  RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         BBitems = (FillerItemRatio + 1) * Nitems) %>%
  # Here we break down the number of BATAs per file to be significantly fewer
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(BBitems, Dataset) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  tidyr::crossing(EffectSize=c(15, 35),
                  Space = c("")) %>%
  arrange(BBitems, MinK) %>%
  mutate(bb_file = paste0(bb_files_path, "bb_", iprefix, "_",
                          Nsubj, "_", Nitems, "_",
                          Bins, "_", MinK, "_to_", MaxK, ".RDS"),
         # Model summary file: stores all information from each mixed-effect analysis
         # (e.g., direction of effect, p-value, etc.)
         mm_file = paste0(mm_files_path, "modeldata_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         # Source data file: the indices in the barebones data file refer to this
         # data set
         real_df_file = map_chr(Dataset, ~og_data_paths[[.x]][["excl"]])
  ) 

sfetal_resid_bb_df <- sfetal_resid_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_sfetal_resid_bb_files <- remaining_files(sfetal_resid_bb_df, file.exists(bb_file))

# Save BB files
sfetal_resid_bb <- {
  iterate_over_df(
    remaining_sfetal_resid_bb_files,
    future_cnd_map,
    function(zaza) {
      k_e <- MinK:MaxK
      filterers <- quos(Group == "Filler-first")

      zplyr::collect_all({
        real_df <- readRDS(real_df_file)
        l <- make_simple_sample_list(
          df = real_df, k = k_e,
          filter_quosures = filterers,
          # The function and its named additional arguments
          single_sample_df_function = same_items_per_subj_with_fillers, # Note that we use a different sampling function here
          items_per_subj = BBitems,
          n_subj = Nsubj,
          n_critical_items = Nitems)
        saveRDS(l, bb_file)
        "good!"
      },
      catchErrors = TRUE)
    })
}

# Make MM files
remaining_sfetal_resid_mm_files <- remaining_files(sfetal_resid_mm_df, file.exists(mm_file))
resid_mm_l <- {
  iterate_over_df(
    remaining_sfetal_resid_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      model_funcs <- safe_residual_no_slopes_no_subject_w_logshift

      
      
      # Establish adding functions. This is a LOOOOT bigger than before because a LOT of pre-processing needs to be done to residualize things. We add the effects, then residualize for each time of analysis.
      adder_f <- function(df) {
        mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
          mutate(logRT = log10(RT),
                 RT_with_Effect =         (RT + SimCond * effect_size),
                 LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
          min_shift_setter() %>%
          mutate(LgShftRT_with_Effect = log10(RT_with_Effect-PowerShift),
                 LgShftRT = log10(RT-TypeIShift)) %>%
          amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
          amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
          amlap_residualize_bata_one_y(      RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
          amlap_residualize_bata_one_y(   LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LgShftRT, BATA_PredLgShftRT) %>%
          amlap_residualize_bata_one_y(LgShftRT_with_Effect, BATA_PredLgShftRT_with_Effect) %>%
          # We then filter out the filler data after residualizing
          filter(ItemType=="CriticalRegion") %>%
          mutate(
            RT_Resid =             RT - BATA_PredRawRT,
            LogRT_Resid =   log10(RT) - BATA_PredLogRT,
            LgShftRT_Resid = LgShftRT - BATA_PredLgShftRT, 
            RT_Resid_with_Effect =             RT_with_Effect - BATA_PredRawRT_with_Effect,
            LogRT_Resid_with_Effect =       LogRT_with_Effect - BATA_PredLogRT_with_Effect,
            LgShftRT_Resid_with_Effect = LgShftRT_with_Effect - BATA_PredLgShftRT_with_Effect
            )
      }

      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
}

# # Load the model results for the residualized data frames
load_sfetal_resid_mixed_models <- {
  iterate_over_df(
    sfetal_resid_mm_df,
    map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               FillerItemRatio = FillerItemRatio,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               real_df_file = real_df_file)
    }) %>%
    saveRDS(paste0(complete_path, "logshift_sfetal_residualized_v01.RDS"))
  "done"
}

######################################################################################################
# Doing the regionized models ########################################################################
######################################################################################################

sfetal_region_mm_df <- tidyr::crossing(
  FillerItemRatio = 15,
  MinK = 1, MaxK = 1,
  iprefix = "same_across_subj_regionized_new",
  Nsubj = c(8, 16, 32, 64),
  Bins = "giant",
  Dataset = c("setal","fetal"),
  RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         FillerItems = FillerItemRatio * Nitems) %>%
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(FillerItems, Dataset) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  tidyr::crossing(EffectSize=c(15, 35),
                  Space = c("")) %>%
  arrange(FillerItems, MinK) %>%
  mutate(bb_file = paste0(bb_files_path, "bb_", iprefix, "_",
                          Nsubj, "_", Nitems, "_",
                          Bins, "_", MinK, "_to_", MaxK, ".RDS"),
         # Model summary file: stores all information from each mixed-effect analysis
         # (e.g., direction of effect, p-value, etc.)
         mm_file = paste0(mm_files_path, "modeldata_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         # Source data file: the indices in the barebones data file refer to this
         # data set
         real_df_file = map_chr(Dataset, ~og_data_paths[[.x]][["excl"]])
  ) 

sfetal_region_bb_df <- sfetal_region_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_sfetal_region_bb_files <- remaining_files(sfetal_region_bb_df, file.exists(bb_file))

# Save BB files
sfetal_region_bb <- {
  iterate_over_df(
    remaining_sfetal_region_bb_files,
    future_cnd_map,
    function(zaza) {
      k_e <- MinK:MaxK
      filterers <- quos(Group == "Filler-first")
      
      zplyr::collect_all({
        real_df <- readRDS(real_df_file)
        l <- make_simple_sample_list(
          df = real_df,
          k = k_e,
          filter_quosures = filterers,
          # The function and its named additional arguments
          single_sample_df_function = same_items_per_subj_for_residualized_regions_free_order, #same_items_per_subj_for_residualized_regions,
          critical_regions_per_subj = Nitems,
          filler_items_per_subj = FillerItems,
          n_subj = Nsubj,
          words_per_region = 3#,
          # exclude_items_with_fewer_than = 10
        )
        saveRDS(l, bb_file)
        "good!"
      },
      catchErrors = TRUE)
    })
}

# MM
remaining_sfetal_region_mm_files <- remaining_files(sfetal_region_mm_df, file.exists(mm_file))
# Run MM
region_mm_l <- {
  iterate_over_df(
    remaining_sfetal_region_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      model_funcs <- safe_residual_no_slopes_no_subject_w_logshift
      
      adder_f <- function(df) {
        mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
          mutate(logRT = log10(RT),
                 RT_with_Effect =         (RT + SimCond * effect_size),
                 LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
          min_shift_setter() %>%
          mutate(LgShftRT_with_Effect = log10(RT_with_Effect - PowerShift),
                 LgShftRT = log10(RT - TypeIShift)) %>%
          amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
          amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
          amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LgShftRT, BATA_PredLgShftRT) %>%
          amlap_residualize_bata_one_y(LgShftRT_with_Effect, BATA_PredLgShftRT_with_Effect) %>%
          # Filters out fillers
          filter(ItemType == "CriticalRegion") %>%
          mutate(UniqueItem = factor(UniqueItem)) %>%
          # Calculates residuals BEFORE averaging over the regions
          mutate(
            RT_Resid =             RT - BATA_PredRawRT,
            LogRT_Resid =   log10(RT) - BATA_PredLogRT,
            LgShftRT_Resid = LgShftRT - BATA_PredLgShftRT,
            RT_Resid_with_Effect =             RT_with_Effect - BATA_PredRawRT_with_Effect,
            LogRT_Resid_with_Effect =       LogRT_with_Effect - BATA_PredLogRT_with_Effect,
            LgShftRT_Resid_with_Effect = LgShftRT_with_Effect - BATA_PredLgShftRT_with_Effect 
          ) %>%
          group_by(Item, UniqueItem, Subject, UniqueSubject, SimCond, ItemAmbValue, SubjectAmbValue) %>%
          summarise(      RT_Resid = mean(      RT_Resid, na.rm = TRUE),
                          LogRT_Resid = mean(   LogRT_Resid, na.rm = TRUE),
                          LgShftRT_Resid = mean(LgShftRT_Resid, na.rm = TRUE),
                          RT_Resid_with_Effect = mean(     RT_Resid_with_Effect, na.rm = TRUE),
                          LogRT_Resid_with_Effect = mean(   LogRT_Resid_with_Effect, na.rm = TRUE),
                          LgShftRT_Resid_with_Effect = mean(LgShftRT_Resid_with_Effect, na.rm = TRUE),
                          AllRTs = paste(RT, collapse=", "),
                          RT = mean(RT, na.rm=TRUE),
                          .groups="drop_last") %>%
          ungroup()
      }
      
      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
}

# Loads files
sfetal_region_load <- {
  iterate_over_df(
    sfetal_region_mm_df,
    map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               real_df_file = real_df_file)
    }
  ) %>%
    saveRDS(paste0(complete_path, "logshift_sfetal_regionized_v01.RDS"))
  "done"
}





######################################################################################################
# NSC data ###########################################################################################
######################################################################################################

# Now onto the NSC data

# Preresidualized data =============================================================
nsc_giant_mm_df <- tidyr::crossing(
  Bins = "giant",
  Knum = seq_along(k_list),
  Nitems = c(8, 16, 32, 64),
  Nsubj = c(8, 16, 32, 64),
  Nstories = c(4),
  iprefix = c("same_across_subj_new", "same_across_subj_noexcl_new"),
  EffectSize = c(56, 80), 
  Space = c(""),
  RType = c("unresidualized"),
  Dataset = c("nsc")) %>%
  filter(Nsubj == Nitems) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(
    # Add in redundancy
    iprefix = paste0(Dataset,"_",iprefix),
    # Are there exclusions?
    Exclusions = ifelse(grepl("noexcl", iprefix), FALSE, TRUE),
    exclname = ifelse(Exclusions,"excl","noexcl"),
    # Barebone data files: stores indices of sampled data for each generated
    # data set (BATA)
    bb_file = paste0(bb_files_path, "bb_", iprefix, "_",
                     Nsubj, "_", Nitems, "_",
                     Bins, "_", MinK, "_to_", MaxK, ".RDS"),
    # Model summary file: stores all information from each mixed-effect analysis
    # (e.g., direction of effect, p-value, etc.)
    mm_file = paste0(mm_files_path, "modeldata_",
                     iprefix, "_", Nsubj, "_", Nitems, "_",
                     Bins, "_",
                     RType, "_size_", EffectSize, Space, "_",
                     MinK, "_to_", MaxK, ".RDS"),
    # Source data file: the indices in the barebones data file refer to this
    # data set
    real_df_file = map2_chr(Dataset, exclname, ~og_data_paths[[.x]][[.y]]))

nsc_giant_bb_df <- nsc_giant_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_nsc_giant_bb_files <- remaining_files(nsc_giant_bb_df, file.exists(bb_file))

# Save BB files
nsc_giant_bb <- {
  iterate_over_df(
    remaining_nsc_giant_bb_files,
    future_cnd_map,
    function(zaza) {
      words_per_region = 1
      k <- MinK:MaxK
      k_e <- MinK:MaxK
      real_df <- readRDS(real_df_file)
      filterers <- quos(Group == "Filler-first")
      
      # Setting this manually
      Nstories <- 4
      
      bb_df <- real_df %>%
        group_by(Item) %>%
        summarise(maxWords = max(Word.Order)) %>%
        mutate(Words.Start = purrr::map_chr(maxWords, ~paste0(1:., collapse=","))) %>%
        tidyr::separate_rows(Words.Start, sep=",") %>%
        mutate(Words.Start = as.numeric(Words.Start)) %>%
        group_by(Item) %>%
        filter(Words.Start != maxWords,
               Words.Start != 1,
               Words.Start <= maxWords - words_per_region) %>%
        ungroup() %>%
        add_possible_subjects(real_df, words_per_region) %>%
        zplyr::left_join(distinct(real_df, Story, Item), by="Item")
      
      l <- map(
        k,
        function(k_e) {
          sample_NSC_bata(bb_df             = bb_df,
                          n_stories         = Nstories,
                          n_subj            = Nsubj,
                          n_items_per_story = Nitems/Nstories,
                          words_per_region  = words_per_region)
        }) %>%
        purrr::set_names(nm = k)
      
      saveRDS(l, bb_file)
      "done!"
    })
}
nsc_giant_bb


remaining_nsc_giant_mm_files <- remaining_files(nsc_giant_mm_df, file.exists(mm_file))
# Saves MM files
nsc_giant_mm <- {
  iterate_over_df(
    remaining_nsc_giant_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize

      bb_list <- readRDS(bb_file)
      real_df <- readRDS(real_df_file)
      
      model_funcs <- no_slopes_wo_infs_plus_logshift
      
      effect_size <- EffectSize
      data_tidyer <- mixed_cleaner
      
      adder_f <- function(df) {
        mutate(df, RT_with_Effect = RT + SimCond * effect_size) %>%
          min_shift_setter()
      }
      
      purrr::map(
        bb_list,
        function(bb_df, i) {
          fl_df <- bb_df %>%
            tidyr::unnest(Word.Order) %>%
            left_join(real_df,
                      by = c("Story", "Item",
                             "Subject", "Word.Order")) %>%
            mutate(SimCond = SimCond/2) %>%
            adder_f()
          
          purrr::map(
            model_funcs,
            ~zplyr::collect_all(data_tidyer(.(fl_df)),
                                catchErrors=TRUE))
        }
      ) %>%
        saveRDS(mm_file)
      "done!"
    })
}
nsc_giant_mm

# Load data into one file -----------
nsc_giant_load <- {
  iterate_over_df(
    nsc_giant_mm_df,
    map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               Dataset = Dataset,
               Exclusions = Exclusions,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               MinK = MinK)
    }
  ) %>%
    saveRDS(paste0(complete_path, "logshift_nsc_v01.RDS"))
}
nsc_giant_load




######################################################################################################
# Doing the residualized models ######################################################################
######################################################################################################

nsc_resid_mm_df <- tidyr::crossing(
  FillerItemRatio = 15,
  MinK = 1, MaxK = 1,
  iprefix = "same_across_subj_resid_new",
  Nsubj = c(8, 16, 32, 64),
  Nstories = c(4),
  Bins = "giant",
  Dataset = c("nsc"),
  RType = "residualized") %>%
  mutate(TotalItems = Nsubj,
         BBitems = (FillerItemRatio + 1) * TotalItems, # the actual total # of items sampled
         BBitemsPerStory = BBitems/Nstories,
         Nitems = TotalItems/Nstories) %>% # important
  {
    bind_rows(rep(list(.), residual_total_files))
  } %>%
  group_by(BBitems) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  tidyr::crossing(EffectSize=c(15, 35),
                  Space = c("_log_space", "")) %>%
  arrange(BBitems, MinK) %>%
  mutate(bb_file = paste0(bb_files_path, "bb_", iprefix, "_",
                          Nsubj, "_", Nitems, "_",
                          Bins, "_", MinK, "_to_", MaxK, ".RDS"),
         # Model summary file: stores all information from each mixed-effect analysis
         # (e.g., direction of effect, p-value, etc.)
         mm_file = paste0(mm_files_path, "modeldata_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         # Source data file: the indices in the barebones data file refer to this
         # data set
         real_df_file = map_chr(Dataset, ~og_data_paths[[.x]][["excl"]])
  )

nsc_resid_bb_df <- nsc_resid_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_nsc_resid_bb_files <- remaining_files(nsc_resid_bb_df, file.exists(bb_file))

# Make the bb dfs for the residualized data frames
bb_nsc_resid_l <- {
  iterate_over_df(
    remaining_nsc_resid_bb_files,
    future_cnd_map,
    function(zaza) {
      words_per_region = 1
      k = MinK:MaxK
      real_df <- readRDS(real_df_file)
      n_conditions <- 2
      
      bb_df <- real_df %>%
        group_by(Item) %>%
        summarise(maxWords = max(Word.Order)) %>%
        mutate(Words.Start = purrr::map_chr(maxWords,
                                            ~paste0(1:., collapse=","))) %>%
        tidyr::separate_rows(Words.Start, sep=",") %>%
        mutate(Words.Start = as.numeric(Words.Start)) %>%
        group_by(Item) %>%
        filter(Words.Start != maxWords,
               Words.Start != 1,
               Words.Start <= maxWords - words_per_region) %>%
        ungroup() %>%
        add_possible_subjects(real_df, words_per_region) %>%
        zplyr::left_join(distinct(real_df, Story, Item), by="Item")
      
      l <- purrr::map(
        k,
        function(k_e) {
          sample_NSC_bata(bb_df             = bb_df,
                          n_stories         = Nstories,
                          n_subj            = Nsubj,
                          n_items_per_story = BBitemsPerStory,
                          words_per_region  = words_per_region) %>%
            # Sets the criticalregion flag HERE rather than later
            group_by(UniqueStory, SimCond) %>%
            mutate(ItemType = ifelse(UniqueItem %in% unique(UniqueItem)[1:Nitems/n_conditions],
                                     "CriticalRegion", "FillerItems")) %>%
            ungroup()
        }) %>%
        purrr::set_names(nm = k)
      
      saveRDS(l, bb_file)
      "done!"
      bb_file
    })
}
bb_nsc_resid_l


# Make the models for the residualized data frames
remaining_nsc_resid_mm_files <- remaining_files(nsc_resid_mm_df, file.exists(mm_file))

nsc_resid_mm_l <-{
  iterate_over_df(
    remaining_nsc_resid_mm_files,
    future_cnd_map,
    function(zaa) {
      
      l <- readRDS(bb_file)
      real_df <- readRDS(real_df_file)
      
      effect_size <- EffectSize
      data_tidyer <- mixed_cleaner
      model_funcs <- safe_residual_no_slopes_no_subject_w_logshift
      
      # Establish adding functions
      adder_f <- function(df) {
        mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
          mutate(logRT = log10(RT),
                 RT_with_Effect =         (RT + SimCond * effect_size),
                 LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
          min_shift_setter() %>%
          mutate(LgShftRT_with_Effect = log10(RT_with_Effect-PowerShift),
                 LgShftRT = log10(RT-TypeIShift)) %>%
          amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
          amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
          amlap_residualize_bata_one_y(      RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
          amlap_residualize_bata_one_y(   LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LgShftRT, BATA_PredLgShftRT) %>%
          amlap_residualize_bata_one_y(LgShftRT_with_Effect, BATA_PredLgShftRT_with_Effect) %>%
          filter(ItemType=="CriticalRegion") %>%
          mutate(
            RT_Resid =             RT - BATA_PredRawRT,
            LogRT_Resid =   log10(RT) - BATA_PredLogRT,
            LgShftRT_Resid = LgShftRT - BATA_PredLgShftRT, 
            RT_Resid_with_Effect =             RT_with_Effect - BATA_PredRawRT_with_Effect,
            LogRT_Resid_with_Effect =       LogRT_with_Effect - BATA_PredLogRT_with_Effect,
            LgShftRT_Resid_with_Effect = LgShftRT_with_Effect - BATA_PredLgShftRT_with_Effect
          )
      }
      
      purrr::imap(
        l,
        function(bb_df, i) {
          fl_df <- bb_df %>%
            tidyr::unnest(Word.Order) %>%
            left_join(real_df, by = c("Story", "Item", "Subject", "Word.Order")) %>%
            mutate(SimCond = SimCond/2) %>%
            adder_f()
          
          purrr::map(model_funcs, ~zplyr::collect_all(data_tidyer(.(fl_df)),
                                                      catchErrors=TRUE))
        }) %>%
        saveRDS(mm_file)
      "done!"
    })
}
nsc_resid_mm_l

# # Load the model results for the residualized data frames
load_nsc_resid_mixed_models <- {
  iterate_over_df(
    nsc_resid_mm_df,
    map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               FillerItemRatio = FillerItemRatio,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               real_df_file = real_df_file)
    }) %>%
    saveRDS(paste0(complete_path, "logshift_nsc_residualized_v01.RDS"))
  "done"
}


######################################################################################################
# Doing the regionized models ########################################################################
######################################################################################################

nsc_region_mm_df <- tidyr::crossing(
  FillerItemRatio = 15,
  MinK = 1, MaxK = 1,
  iprefix = "same_across_subj_regionized_new",
  Nsubj = c(8, 16, 32, 64),
  Nstories = c(4),
  Bins = "giant",
  Dataset = c("nsc"),
  Space="",
  RType = "residualized") %>%
  mutate(TotalItems = Nsubj,
         FillerItems = FillerItemRatio * TotalItems, # Because critical items are diff
         FillerItemsPerStory = FillerItems/Nstories,
         Nitems = TotalItems/Nstories) %>%
  {
    bind_rows(rep(list(.), residual_total_files))
  } %>%
  group_by(FillerItems, Nstories) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  arrange(FillerItems, MinK) %>%
  tidyr::crossing(EffectSize = c(15, 35)) %>%
  arrange(FillerItems, MinK) %>%
  mutate(bb_file = paste0(bb_files_path, "bb_", iprefix, "_",
                          Nsubj, "_", Nitems, "_",
                          Bins, "_", MinK, "_to_", MaxK, ".RDS"),
         # Model summary file: stores all information from each mixed-effect analysis
         # (e.g., direction of effect, p-value, etc.)
         mm_file = paste0(mm_files_path, "modeldata_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         # Source data file: the indices in the barebones data file refer to this
         # data set
         real_df_file = map_chr(Dataset, ~og_data_paths[[.x]][["excl"]])
  ) 

nsc_region_bb_df <- nsc_region_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_nsc_region_bb_files <- remaining_files(nsc_region_bb_df, file.exists(bb_file))

# Make the bb dfs for the regionized data frames
nsc_region_bb_lists <- {
  iterate_over_df(
    remaining_nsc_region_bb_files,
    future_cnd_map,
    function(zaza) {
      # First, n is 1, then it will be changed to n=3, don't worry
      n = 1
      words_per_critical_region = 3
      k = MinK:MaxK
      real_df <- readRDS(real_df_file)
      
      bb_df <- real_df %>%
        group_by(Item) %>%
        summarise(maxWords = max(Word.Order)) %>%
        mutate(Words.Start = purrr::map_chr(maxWords, ~paste0(1:., collapse=","))) %>%
        tidyr::separate_rows(Words.Start, sep=",") %>%
        mutate(Words.Start = as.numeric(Words.Start)) %>%
        group_by(Item) %>%
        filter(Words.Start != maxWords,
               Words.Start != 1,
               Words.Start <= maxWords - n) %>%
        ungroup() %>%
        add_possible_subjects(real_df, n) %>%
        zplyr::left_join(distinct(real_df, Story, Item), by="Item")
      
      l <- purrr::map(
        k,
        function(k_e) {
          sample_regionized_NSC_bata(bb_df,
                                     critical_regions_per_subj_per_story = Nitems,
                                     n_stories                 = Nstories,
                                     filler_items_per_story    = FillerItemsPerStory,
                                     n_subj                    = Nsubj,
                                     words_per_region          = words_per_critical_region)
        }) %>%
        purrr::set_names(nm = k)
      
      saveRDS(l, bb_file)
      "done!"
      bb_file
    })
}
nsc_region_bb_lists


remaining_nsc_region_mm_files <- remaining_files(nsc_region_mm_df, file.exists(mm_file))
# Make the models for the regionized data frames
nsc_region_mm_l <- {
  iterate_over_df(
    remaining_nsc_region_mm_files,
    future_cnd_map,
    function(zaa) {
      l <- readRDS(bb_file)
      real_df <- readRDS(real_df_file)
      
      effect_size <- EffectSize
      data_tidyer = mixed_cleaner
      model_funcs <- safe_residual_no_slopes_no_subject_w_logshift
      
      adder_f <- function(df) {
        mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
          mutate(logRT = log10(RT),
                 RT_with_Effect =         (RT + SimCond * effect_size),
                 LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
          min_shift_setter() %>%
          mutate(LgShftRT_with_Effect = log10(RT_with_Effect - PowerShift),
                 LgShftRT = log10(RT - TypeIShift)) %>%
          amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
          amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
          amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LgShftRT, BATA_PredLgShftRT) %>%
          amlap_residualize_bata_one_y(LgShftRT_with_Effect, BATA_PredLgShftRT_with_Effect) %>%
          # Filters out fillers
          filter(ItemType == "CriticalRegion") %>%
          mutate(UniqueItem = factor(UniqueItem)) %>%
          # Calculates residuals BEFORE averaging over the regions
          mutate(
            RT_Resid =             RT - BATA_PredRawRT,
            LogRT_Resid =   log10(RT) - BATA_PredLogRT,
            LgShftRT_Resid = LgShftRT - BATA_PredLgShftRT,
            RT_Resid_with_Effect =             RT_with_Effect - BATA_PredRawRT_with_Effect,
            LogRT_Resid_with_Effect =       LogRT_with_Effect - BATA_PredLogRT_with_Effect,
            LgShftRT_Resid_with_Effect = LgShftRT_with_Effect - BATA_PredLgShftRT_with_Effect
          ) %>%
          group_by(Item, UniqueStory, Story, UniqueItem, Subject, UniqueSubject, SimCond) %>%
          summarise(      RT_Resid = mean(      RT_Resid, na.rm = TRUE),
                          LogRT_Resid = mean(   LogRT_Resid, na.rm = TRUE),
                          LgShftRT_Resid = mean(LgShftRT_Resid, na.rm = TRUE),
                          RT_Resid_with_Effect = mean(     RT_Resid_with_Effect, na.rm = TRUE),
                          LogRT_Resid_with_Effect = mean(   LogRT_Resid_with_Effect, na.rm = TRUE),
                          LgShftRT_Resid_with_Effect = mean(LgShftRT_Resid_with_Effect, na.rm = TRUE),
                          AllRTs = paste(RT, collapse=", "),
                          RT = mean(RT, na.rm=TRUE),
                          .groups="drop_last") %>%
          ungroup()
      }
      
      purrr::imap(
        l,
        function(bb_df, i) {
          fl_df <- bb_df %>%
            tidyr::unnest(Word.Order) %>%
            left_join(real_df, by = c("Story", "Item", "Subject", "Word.Order")) %>%
            mutate(SimCond = SimCond/2) %>%
            adder_f()
          purrr::map(model_funcs, ~zplyr::collect_all(data_tidyer(.(fl_df)),
                                                      catchErrors=TRUE))  }) %>%
        saveRDS(mm_file)
      "done!"
    })
}
nsc_region_mm_l


# Loads files
nsc_region_load <- {
  iterate_over_df(
    nsc_region_mm_df,
    map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               real_df_file = real_df_file)
    }
  ) %>%
    saveRDS(paste0(complete_path, "logshift_nsc_regionized_v01.RDS"))
  "done"
}




