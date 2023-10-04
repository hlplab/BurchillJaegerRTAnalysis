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
kill_r_on_nodes(decent_nodes$nodename, backup_login, # "zburchil@cycle3.cs.rochester.edu",
                wait_until_resolved = TRUE,
                beep_on_end = TRUE)

# SETAL ###########################################################
setal_resid_mm_df <- tidyr::crossing(
  FillerItemRatio = 15,
  MinK = 1, MaxK = 1,
  iprefix = "same_across_subj_resid_new",
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
  tidyr::crossing(EffectSize=c(15, 35),
                  Space = c("_log_space", "")) %>%
  arrange(BBitems, MinK) %>%
  mutate(bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          BBitems, "_fillers_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_",  Nsubj, "_", Nitems, "_",
                          BBitems, "_fillers_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file =  paste0(path_on_server, setal_file)) %>%
  filter(Space == "")

setal_resid_bb_df <- setal_resid_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_setal_resid_bb_files <- remaining_files(setal_resid_bb_df, file.exists(bb_file))

# Save BB files
setal_resid_bb %beep% {
  iterate_over_df(
    remaining_setal_resid_bb_files,
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
          single_sample_df_function = same_items_per_subj_with_fillers,
          items_per_subj = BBitems,
          n_subj = Nsubj,
          n_critical_items = Nitems)
        saveRDS(l, bb_file)
        "good!"
      },
      catchErrors = TRUE)
    })
} %plan% decent_plan
done(setal_resid_bb)

remaining_setal_resid_mm_files <- remaining_files(setal_resid_mm_df, file.exists(mm_file))
# Make MM files
resid_mm_l %beep% {
  iterate_over_df(
    remaining_setal_resid_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      model_funcs <- safe_residual_no_slopes_no_subject

      # Establish adding functions
      adder_f <- function(df) {
        mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
          mutate(logRT = log10(RT),
                 RT_with_Effect =         (RT + SimCond * effect_size),
                 LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
          amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
          amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
          amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
          filter(ItemType=="CriticalRegion") %>%
          mutate(
            RT_Resid =           RT - BATA_PredRawRT,
            LogRT_Resid = log10(RT) - BATA_PredLogRT,
            RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT_with_Effect,
            LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT_with_Effect)
      }

      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% elite_squad_plan
done(resid_mm_l)

# Load the model results for the residualized data frames
load_setal_resid_mixed_models %beep% {
  iterate_over_df(
    setal_resid_mm_df,
    future_map_dfr,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/residualized_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/residualized_mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
}
done(load_setal_resid_mixed_models)

# FETAL ##################################################################
fetal_resid_mm_df <- tidyr::crossing(
  FillerItemRatio = 15,
  MinK = 1, MaxK = 1,
  iprefix = "fetal_same_across_subj_resid_new",
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
  tidyr::crossing(EffectSize=c(15, 35),
                  Space = c("_log_space", "")) %>%
  arrange(BBitems, MinK) %>%
  mutate(bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_",
                          BBitems, "_fillers_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_",  Nsubj, "_", Nitems, "_",
                          BBitems, "_fillers_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file =  paste0(path_on_server, fetal_file)) %>%
  filter(Space == "")

fetal_resid_bb_df <- fetal_resid_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_fetal_resid_bb_files <- remaining_files(fetal_resid_bb_df, file.exists(bb_file))

# Save BB files
fetal_resid_bb %beep% {
  iterate_over_df(
    remaining_fetal_resid_bb_files,
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
          single_sample_df_function = same_items_per_subj_with_fillers,
          items_per_subj = BBitems,
          n_subj = Nsubj,
          n_critical_items = Nitems)
        saveRDS(l, bb_file)
        "good!"
      },
      catchErrors = TRUE)
    })
} %plan% elite_squad_plan
done(fetal_resid_bb)

remaining_fetal_resid_mm_files <- remaining_files(fetal_resid_mm_df, file.exists(mm_file))
# Make MM files
resid_mm_l %beep% {
  iterate_over_df(
    remaining_fetal_resid_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      model_funcs <- safe_residual_no_slopes_no_subject

      # Establish adding functions
      adder_f <- function(df) {
        mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
          mutate(logRT = log10(RT),
                 RT_with_Effect =         (RT + SimCond * effect_size),
                 LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
          amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
          amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
          amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
          filter(ItemType=="CriticalRegion") %>%
          mutate(
            RT_Resid =           RT - BATA_PredRawRT,
            LogRT_Resid = log10(RT) - BATA_PredLogRT,
            RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT_with_Effect,
            LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT_with_Effect)
      }

      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% elite_squad_plan
done(resid_mm_l)

# Load the model results for the residualized data frames
load_fetal_resid_mixed_models %beep% {
  iterate_over_df(
    fetal_resid_mm_df,
    future_map_dfr,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/fetal_residualized_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/fetal_residualized_mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
}
done(load_fetal_resid_mixed_models)

# NSC #############################################################
# Define the mixed models df first, since the bb is a subset
nsc_resid_mm_df <- tidyr::crossing(
  FillerItemRatio = 15,
  MinK = 1, MaxK = 1,
  iprefix = "nsc_same_across_item_resid_new",
  Nsubj = c(8, 16, 32, 64),
  Nstories = c(4),
  Bins = "giant",
  RType = "residualized",
  EffectSize = c(15, 35),
  Space = c("")) %>%
  mutate(TotalItems = Nsubj,
         BBitems = (FillerItemRatio + 1) * TotalItems, # the actual total # of items sampled
         BBitemsPerStory = BBitems/Nstories,
         Nitems = TotalItems/Nstories) %>% # important
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(BBitems, Nstories) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  arrange(BBitems, MinK) %>%
  mutate(bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", BBitems, "_", Nstories, "_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_", Nsubj, "_", BBitems, "_", Nstories, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file = paste0(path_on_server, nsc_file))  %>%
  filter(Nitems > 1) %>% # Cuz we need at least two conditions per story
  filter(Space == "")

nsc_resid_bb_df <- nsc_resid_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_nsc_resid_bb_files <- remaining_files(nsc_resid_bb_df, file.exists(bb_file))

# Make the bb dfs for the residualized data frames
bb_nsc_resid_l %beep% {
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

      l <- furrr::future_map(
        k,
        function(k_e) {
          sample_NSC_bata(bb_df             = bb_df,
                          n_stories         = Nstories,
                          n_subj            = Nsubj,
                          n_items_per_story = BBitemsPerStory,
                          words_per_region  = words_per_region) %>%
            # Sets the criticalregion flag HERE rather than later
            group_by(UniqueStories, SimCond) %>%
            mutate(ItemType = ifelse(UniqueItem %in% unique(UniqueItem)[1:Nitems/n_conditions],
                                     "CriticalRegion", "FillerItems")) %>%
            ungroup()
        }) %>%
        purrr::set_names(nm = k)

      saveRDS(l, bb_file)
      "done!"
      bb_file
    })
} %plan% decent_nodes
done(bb_nsc_resid_l)

# Make the models for the residualized data frames
remaining_nsc_resid_mm_files <- remaining_files(nsc_resid_mm_df, file.exists(mm_file))
nsc_resid_mm_l %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    remaining_nsc_resid_mm_files,
    future_cnd_map,
    function(zaa) {

      l <- readRDS(bb_file)
      real_df <- readRDS(real_df_file)

      effect_size <- EffectSize
      data_tidyer = mixed_cleaner
      model_funcs <- safe_residual_no_slopes_no_subject

      # Establish adding functions
      adder_f <- function(df) {
        mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0)) %>%
          # set_filler_data_NSC(df, n_remaining_items_per_story = Nitems,  n_conditions=2) %>%
          mutate(logRT = log10(RT),
                 RT_with_Effect =         (RT + SimCond * effect_size),
                 LogRT_with_Effect = log10(RT + SimCond * effect_size)) %>%
          amlap_residualize_bata_one_y(   RT, BATA_PredRawRT) %>%
          amlap_residualize_bata_one_y(logRT, BATA_PredLogRT) %>%
          amlap_residualize_bata_one_y(   RT_with_Effect, BATA_PredRawRT_with_Effect) %>%
          amlap_residualize_bata_one_y(LogRT_with_Effect, BATA_PredLogRT_with_Effect) %>%
          # filter(SimCond != 0) %>%
          filter(ItemType=="CriticalRegion") %>%
          mutate(UniqueItem = factor(UniqueItem)) %>%
          mutate(
            RT_Resid =           RT - BATA_PredRawRT,
            LogRT_Resid = log10(RT) - BATA_PredLogRT,
            RT_Resid_with_Effect =       RT_with_Effect - BATA_PredRawRT_with_Effect,
            LogRT_Resid_with_Effect = LogRT_with_Effect - BATA_PredLogRT_with_Effect)
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
                                                      catchErrors=TRUE))
        }) %>%
        saveRDS(mm_file)
      "done!"
    })
  tictoc::tic()
  l
} %plan% elite_squad_plan
done(nsc_resid_mm_l)

# Load the model results for the residualized data frames
load_nsc_resid_mixed_models %beep% {
  iterate_over_df(
    nsc_resid_mm_df,
    future_map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nstories = Nstories,
               Nitems = TotalItems,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/nsc_residualized_mixed_models_lmerTest.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/nsc_residualized_mixed_models_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
}
done(load_nsc_resid_mixed_models)


