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


# SETAL #############################################################
# Define the mixed models df first, since the bb is a subset
setal_giant_mm_df <- tidyr::crossing(
  Bins = "giant",
  Knum = seq_along(k_list),
  Nitems = c(8, 16, 32, 64),
  Nsubj = c(8, 16, 32, 64),
  # iprefix = c("same_across_subj", "same_across_subj_noexcl"),
  iprefix = c("same_across_subj_new", "same_across_subj_noexcl_new"),
  EffectSize = c(56, 80), # got rid of the 7*2^x
  Space = c("_log_space", ""),
  RType = c("unresidualized")) %>%
  filter(Nsubj == Nitems) %>%
  { d <- filter(., Nsubj == Nitems, Nsubj == 8)
  bind_rows(
    .,
    mutate(d, Nsubj=8, Nitems=16),
    mutate(d, Nsubj=16, Nitems=8),
    mutate(d, Nsubj=16, Nitems=32),
    mutate(d, Nsubj=32, Nitems=16), # For testing
    mutate(d, Nsubj=32, Nitems=64))
  } %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]])),
         bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file = ifelse(grepl("noexcl", iprefix),
                               paste0(path_on_server, "setal_noexcl_fillers.RDS" ),
                               paste0(path_on_server, setal_file))) %>%
  filter(Space == "")

setal_giant_bb_df <- setal_giant_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_setal_giant_bb_files <- remaining_files(setal_giant_bb_df, file.exists(bb_file))

# Save BB files
setal_giant_bb %beep% {
  iterate_over_df(
    remaining_setal_giant_bb_files,
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
} %plan% elite_squad_plan
done(setal_giant_bb)

# Make MM files
remaining_setal_giant_mm_files <- remaining_files(setal_giant_mm_df, file.exists(mm_file))
# Saves MM files
setal_giant_mm %beep% {
  iterate_over_df(
    remaining_setal_giant_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      
      model_funcs <- no_slopes_wo_infs
      
      if (Space=="")
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      else
        stop("Log-space not implemented!")
      
      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% decent_plan
done(setal_giant_mm)

# Loads files
setal_giant_load %beep% {
  iterate_over_df(
    setal_giant_mm_df,
    future_map_dfr,
    function(i) {
      collate_data(mm_file,
                   .multicore = TRUE) %>%
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/giant_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/giant_mixed_models_lmerTest_new.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(setal_giant_load)





# FETAL #############################################################
# Define the mixed models df first, since the bb is a subset
fetal_giant_mm_df <- tidyr::crossing(
  Bins = "giant",
  Knum = seq_along(k_list),
  Nitems = c(8, 16, 32, 64),
  Nsubj = c(8, 16, 32, 64),
  iprefix = c("fetal_same_across_subj_new", "fetal_same_across_subj_noexcl_new"),
  EffectSize = c(56, 80), # got rid of the 7*2^x
  Space = c("_log_space", ""),
  RType = c("unresidualized")) %>%
  filter(Nsubj == Nitems) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]])),
         bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_",  Nsubj, "_", Nitems, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file = ifelse(grepl("noexcl", iprefix),
                               paste0(path_on_server, "fetal_noexcl_fillers.RDS" ),
                               paste0(path_on_server, fetal_file))) %>%
  filter(Space == "") #%>%
  # filter(Nsubj==64) #%>%
  # filter(Nsubj==64, iprefix=="fetal_same_across_subj_noexcl_new_times2")

fetal_giant_bb_df <- fetal_giant_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_fetal_giant_bb_files <- remaining_files(fetal_giant_bb_df, file.exists(bb_file))

# Save BB files
fetal_giant_bb %beep% {
  iterate_over_df(
    remaining_fetal_giant_bb_files,
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
} %plan% elite_squad_plan
done(fetal_giant_bb)

# Make MM files
remaining_fetal_giant_mm_files <- remaining_files(fetal_giant_mm_df, file.exists(mm_file))
# Saves MM files
fetal_giant_mm %beep% {
  iterate_over_df(
    remaining_fetal_giant_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      
      model_funcs <- no_slopes_wo_infs
      
      if (Space=="")
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      else
        stop("Log-space not implemented!")
      
      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% elite_squad_plan
done(fetal_giant_mm)

# Loads files
fetal_giant_load %beep% {
  iterate_over_df(
    fetal_giant_mm_df,
    future_map_dfr,
    function(i) {
      collate_data(mm_file,
                   .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               real_df_file = real_df_file,
               MMnode = Sys.info()[["nodename"]])
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/fetal_giant_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/fetal_giant_mixed_models_lmerTest_runagain2.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(fetal_giant_load)














# NSC ################################################################
# Define the mixed models df first, since the bb is a subset
nsc_giant_mm_df <- tidyr::crossing(
  Bins = "giant",
  Knum = seq_along(k_list),
  Nsubj =  c(8, 16, 32, 64),
  Nstories = c(2, 4, 8),
  iprefix = c("nsc_same_across_item_new", "nsc_same_across_item_noexcl_new"),
  EffectSize = c(56, 80),
  Space = c("_log_space", ""),
  RType = "unresidualized") %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(TotalItems = Nsubj,
         Nitems = TotalItems/Nstories) %>%
  mutate(bb_file = paste0(path_on_server, "amlap_saved_files/bb_", iprefix, "_",
                        Nsubj, "_", TotalItems, "_", Nstories, "_",
                        Bins, "_", MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_", Nsubj, "_", TotalItems, 
                          "_", Nstories, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file = ifelse(grepl("noexcl", iprefix),
                               paste0(path_on_server, "nsc_sim_df_no_RT_exclusions.RDS"),
                               paste0(path_on_server, nsc_file))) %>%
  filter(Nitems > 1)  %>% # Cuz we need at least two conditions per story
  filter(Space == "")

nsc_giant_bb_df <- nsc_giant_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_nsc_giant_bb_files <- remaining_files(nsc_giant_bb_df, file.exists(bb_file))

# Make the bb dfs ---
nsc_bb %beep% {
  iterate_over_df(
    remaining_nsc_giant_bb_files,
    future_cnd_map,
    function(zaza) {
      words_per_region = 1
      k <- MinK:MaxK
      
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
                          n_items_per_story = Nitems,
                          words_per_region  = words_per_region)
        }) %>%
        purrr::set_names(nm = k)
      
      saveRDS(l, bb_file)
      "done!"
    })
} %plan% decent_plan
done(nsc_bb)

# Make the models 
remaining_nsc_giant_mm_files <- remaining_files(nsc_giant_mm_df, file.exists(mm_file))
nsc_mm %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    remaining_nsc_giant_mm_files,
    future_cnd_map,
    function(zaa) {
      
      bb_list <- readRDS(bb_file)
      real_df <- readRDS(real_df_file)
      
      model_funcs <- no_slopes_wo_infs
      effect_size <- EffectSize
      data_tidyer <- mixed_cleaner
      
      adder_f <- function(df) {
        mutate(df, RT_with_Effect = RT + SimCond * effect_size)  }
      
      furrr::future_imap(
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
  tictoc::tic()
  l
} %plan% elite_squad_plan
done(nsc_mm)

# Load the model results into one place
load_nsc_mm %beep% {
  iterate_over_df(
    nsc_giant_mm_df,
    future_map_dfr,
    function(i) {
      collate_data(mm_file,
                   .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = TotalItems,
               Nstories = Nstories,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/nsc_giant_mixed_models_lmerTest_new.RDS"))
  "dooone"
} %plan% elite_squad_plan
done(load_nsc_mm)



# #############################################################
##############################################################
#######################################


# test_newregion2.RDS used `exclude_items_with_fewer_than=13`
# test_newregion3.RDS used `exclude_items_with_fewer_than=8`
# test_newregion4.RDS let region placement be free
# test_newregion5.RDS let region placement be free and made sure no critical items are in the fillers
# test_newregion6.RDS used `exclude_items_with_fewer_than=10` but made sure no critical items are in the fillers 
# test_newregion7.RDS used `exclude_items_with_fewer_than=10` but made sure that fillers were sampled from the same regions
# test_newregion8.RDS used `exclude_items_with_fewer_than=10` but averaged BEFORE residualization

