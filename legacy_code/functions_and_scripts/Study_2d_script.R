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
setal_binned_mm_df <- tidyr::crossing(
  Bins = names(new_bins),
  Knum = seq_along(k_list),
  Nitems = c(8, 16, 32, 64),
  Nsubj = c(8, 16, 32, 64),
  iprefix = c("same_across_subj_new"),
  EffectSize = c(56, 80), # got rid of the 7*2^x
  Space = c("_c_log_space", ""),
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
                               paste0(path_on_server, "setal_noexcl_fillers.RDS" ),
                               paste0(path_on_server, setal_file)))

setal_binned_bb_df <- setal_binned_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_setal_binned_bb_files <- remaining_files(setal_binned_bb_df, file.exists(bb_file))

# Save BB files
setal_binned_bb %beep% {
  iterate_over_df(
    remaining_setal_binned_bb_files,
    future_cnd_map,
    function(zaza) {
      k_e <- MinK:MaxK
      bin <- new_bins[[Bins]]
      low <- bin[1]
      high<- bin[2]
      filterers <- quos(Trial.Order  >= !!enquo(low) &
                          Trial.Order <= !!enquo(high) &
                          Group == "Filler-first")

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
done(setal_binned_bb)

# Make MM files
remaining_setal_binned_mm_files <- remaining_files(setal_binned_mm_df, file.exists(mm_file))
# Saves MM files
setal_binned_mm %beep% {
  iterate_over_df(
    remaining_setal_binned_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      # Set constant to be equal for the middle bin
      log_effect <- log_space_effect_size_calculator(
        setal_basis_info[setal_basis_info$Bins=="bin2", ]$meanRT,
        effect_size)

      model_funcs <- no_slopes_wo_infs

      if (Space=="")
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      else if (Space=="_c_log_space")
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = add_in_log_space(RT, SimCond, log_effect)) }

      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% elite_squad_plan
done(setal_binned_mm)

# Loads files
setal_binned_load %beep% {
  iterate_over_df(
    setal_binned_mm_df,
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
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/binned_mixed_models_lmerTest_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/binned_mixed_models_lmerTest_new.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(setal_binned_load)



