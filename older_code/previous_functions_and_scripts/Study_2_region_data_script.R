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
source(paste0(path, "amlap_bb_making_functions.R"))
source(paste0(path, "get_plans.R"))

done <- cs::done

################################### killling jobs
kill_r_on_nodes(decent_nodes$nodename, tertiary_login, # "zburchil@cycle3.cs.rochester.edu",
                wait_until_resolved = TRUE,
                beep_on_end = TRUE)

# SETAL ###########################################################
setal_region_mm_df <- tidyr::crossing(FillerItemRatio = 15,
                                      MinK = 1, MaxK = 1,
                                      iprefix = "same_across_subj_regionized_new",
                                      Nsubj = c(8, 16, 32, 64),
                                      Bins = "giant",
                                      RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         FillerItems = FillerItemRatio * Nitems) %>%
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(FillerItems) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  tidyr::crossing(EffectSize=c(15, 35),
                  Space = c("_log_space", "")) %>%
  arrange(FillerItems, MinK) %>%
  mutate(bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          FillerItems, "_fillers_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_",  Nsubj, "_", Nitems, "_",
                          FillerItems, "_fillers_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file =  paste0(path_on_server, setal_file)) %>%
  filter(Space == "")

setal_region_bb_df <- setal_region_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_setal_region_bb_files <- remaining_files(setal_region_bb_df, file.exists(bb_file))

# Save BB files
setal_region_bb %beep% {
  iterate_over_df(
    remaining_setal_region_bb_files,
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
} %plan% decent_plan
done(setal_region_bb)

# MM
remaining_setal_region_mm_files <- remaining_files(setal_region_mm_df, file.exists(mm_file))
# Run MM
region_mm_l %beep% {
  iterate_over_df(
    remaining_setal_region_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      # Doesn't double-dip subjects, since they've been residualized out
      model_funcs <- safe_residual_no_slopes_no_subject

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

      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% decent_plan
done(region_mm_l)

# Loads files
setal_region_load %beep% {
  iterate_over_df(
    setal_region_mm_df,
    future_map_dfr,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/regionized_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/regionized_mixed_models_lmerTest_new.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(setal_region_load)


# FETAL ##################################################################
fetal_region_mm_df <- tidyr::crossing(FillerItemRatio = 15,
                                      MinK = 1, MaxK = 1,
                                      iprefix = "fetal_same_across_subj_regionized_new",
                                      Nsubj = c(8, 16, 32, 64),
                                      Bins = "giant",
                                      RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         FillerItems = FillerItemRatio * Nitems) %>%
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(FillerItems) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  tidyr::crossing(EffectSize=c(15, 35),
                  Space = c("_log_space", "")) %>%
  arrange(FillerItems, MinK) %>%
  mutate(bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          FillerItems, "_fillers_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_",  Nsubj, "_", Nitems, "_",
                          FillerItems, "_fillers_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file =  paste0(path_on_server, fetal_file)) %>%
  filter(Space == "")

fetal_region_bb_df <- fetal_region_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_fetal_region_bb_files <- remaining_files(fetal_region_bb_df, file.exists(bb_file))

# Save BB files
fetal_region_bb %beep% {
  iterate_over_df(
    remaining_fetal_region_bb_files,
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
          single_sample_df_function = same_items_per_subj_for_residualized_regions_free_order, # same_items_per_subj_for_residualized_regions,
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
} %plan% decent_plan
done(fetal_region_bb)

# MM
remaining_fetal_region_mm_files <- remaining_files(fetal_region_mm_df, file.exists(mm_file))
# Run MM
region_mm_l %beep% {
  iterate_over_df(
    remaining_fetal_region_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      # Doesn't double-dip subjects, since they've been residualized out
      model_funcs <- safe_residual_no_slopes_no_subject

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

      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% decent_plan
done(region_mm_l)

# Loads files
fetal_region_load %beep% {
  iterate_over_df(
    fetal_region_mm_df,
    future_map_dfr,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/fetal_regionized_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/fetal_regionized_mixed_models_lmerTest_new.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(fetal_region_load)


# NSC #############################################################
# Define the mixed models df first, since the bb is a subset
nsc_region_mm_df <- tidyr::crossing(
  FillerItemRatio = 15,
  MinK = 1, MaxK = 1,
  iprefix = "nsc_same_across_item_regionized_new",
  Nsubj = c(8, 16, 32, 64),
  Nstories = c(4),
  Bins = "giant",
  RType = "residualized",
  Space = c("")) %>%
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
  mutate(bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", TotalItems, "_", FillerItems,
                          "_fillers_", Nstories, "_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_", Nsubj, "_", TotalItems, "_", FillerItems,
                          "_fillers_", Nstories, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file = paste0(path_on_server, nsc_file)) %>%
  filter(Nitems > 1) %>% # Cuz we need at least two conditions per story
  filter(Space == "")

nsc_region_bb_df <- nsc_region_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_nsc_region_bb_files <- remaining_files(nsc_region_bb_df, file.exists(bb_file))

# Make the bb dfs for the regionized data frames
nsc_region_bb_lists %beep% {
  iterate_over_df(
    remaining_nsc_region_bb_files,
    future_cnd_map,
    function(zaza) {
      # First, n is 1, then it will be changed to n=3
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

      l <- furrr::future_map(
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
} %plan% decent_plan
done(nsc_region_bb_lists)

remaining_nsc_region_mm_files <- remaining_files(nsc_region_mm_df, file.exists(mm_file))
# Make the models for the regionized data frames
nsc_region_mm_l %sayname% {
  iterate_over_df(
    remaining_nsc_region_mm_files,
    future_cnd_map,
    function(zaa) {
      l <- readRDS(bb_file)
      real_df <- readRDS(real_df_file)

      effect_size <- EffectSize
      data_tidyer = mixed_cleaner
      model_funcs <- safe_residual_no_slopes_no_subject

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
            RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT_with_Effect,
            LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT_with_Effect
          ) %>%
          group_by(Item, UniqueStory, Story, UniqueItem, Subject, UniqueSubject, SimCond) %>%
          summarise(RT_Resid = mean(RT_Resid, na.rm = TRUE),
                    LogRT_Resid = mean(LogRT_Resid, na.rm = TRUE),
                    RT_Resid_with_Effect = mean(RT_Resid_with_Effect, na.rm = TRUE),
                    LogRT_Resid_with_Effect = mean(LogRT_Resid_with_Effect, na.rm = TRUE),
                    AllRTs = paste(RT, collapse=", "),
                    RT = mean(RT, na.rm=TRUE)) %>%
          ungroup()
      }

      furrr::future_imap(
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
} %plan% b_squad_plan
done(nsc_region_mm_l)

# Load the model results for the regionized data frames
load_nsc_region_mixed_models %beep% {
  iterate_over_df(
    nsc_region_mm_df,
    future_map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = TotalItems,
               Nstories = Nstories,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/nsc_regionized_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/nsc_regionized_mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% decent_plan
done(load_nsc_region_mixed_models)




# SETAL Free region ###########################################################
setal_region_mm_df <- tidyr::crossing(FillerItemRatio = 15,
                                      MinK = 1, MaxK = 1,
                                      iprefix = "same_across_subj_regionized_free",
                                      Nsubj = c(8, 16, 32, 64),
                                      Bins = "giant",
                                      RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         FillerItems = FillerItemRatio * Nitems) %>%
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(FillerItems) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  tidyr::crossing(EffectSize=c(15, 35),
                  Space = c("_log_space", "")) %>%
  arrange(FillerItems, MinK) %>%
  mutate(bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          FillerItems, "_fillers_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_",  Nsubj, "_", Nitems, "_",
                          FillerItems, "_fillers_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file =  paste0(path_on_server, setal_file)) %>%
  filter(Space == "") %>% filter(Nsubj==64)

setal_region_bb_df <- setal_region_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_setal_region_bb_files <- remaining_files(setal_region_bb_df, file.exists(bb_file))

# Save BB files
setal_region_bb %beep% {
  iterate_over_df(
    remaining_setal_region_bb_files,
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
          single_sample_df_function = same_items_per_subj_for_residualized_regions_free_order,
          critical_regions_per_subj = Nitems, 
          filler_items_per_subj = FillerItems,
          n_subj = Nsubj,
          words_per_region = 3)
        saveRDS(l, bb_file)
        "good!"
      },
      catchErrors = TRUE)
    })
} %plan% c_squad_plan
done(setal_region_bb)

# MM
remaining_setal_region_mm_files <- remaining_files(setal_region_mm_df, file.exists(mm_file))
# Run MM
region_mm_l %beep% {
  iterate_over_df(
    remaining_setal_region_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      # Doesn't double-dip subjects, since they've been residualized out
      model_funcs <- safe_residual_no_slopes_no_subject
      
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
      
      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% c_squad_plan
done(region_mm_l)

# Loads files
setal_region_load %beep% {
  iterate_over_df(
    setal_region_mm_df,
    future_map_dfr,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/regionized_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/regionized_mixed_models_lmerTest_new.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(setal_region_load)




