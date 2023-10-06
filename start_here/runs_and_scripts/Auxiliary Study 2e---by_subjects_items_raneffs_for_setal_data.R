remotes::install_github("burchill/zplyr")
remotes::install_github("burchill/cs")

library(tidyverse)
library(future)
library(furrr)
library(cs)
library(beepr)
library(broom.mixed) # need to load this for weird S3 reasons

# Make the future plans
plan(list(
  tweak(multisession, workers = 1), # for the beep!
  tweak(multisession, workers = 1), # for the initialization
  tweak(multisession, workers = 6) # five for my comp
))

done <- cs::done

# Constants ------------------------------------

# Loading all the boilerplate. CHANGE PATH TO YOUR OWN DIRECTORY
main_path = "/Users/tiflo/Library/CloudStorage/Box-Box/_Papers - Box/Power simulations for RTs/PowerSimulations/BurchillJaegerRTAnalysis/start_here/"
# main_path = "/Users/zburchill/Box Sync/Power simulations for RTs/PowerSimulations/BurchillJaegerRTAnalysis/start_here/"

data_path     = paste0(main_path, "og_data/")
models_path   = paste0(main_path, "models/")
bb_files_path = paste0(main_path, "bb_files/")
mm_files_path = paste0(main_path, "mm_files/")
complete_path = paste0(main_path, "collated_files/")
code_path     = paste0(main_path, "functions/")

nsc_file_path = data_path
nsc_df_file = "nsc_sim_df.RDS"
model_files_path = models_path

setal_file_path = data_path
setal_df_file = "setal_sim_df.RDS"

# Loading code -------------------------------------------

source(paste0(code_path, "utils.R"))
source(paste0(code_path, "amlap_constants.R"))
source(paste0(code_path, "amlap_modeling_data_functions.R"))
source(paste0(code_path, "modeling_data_functions.R"))
source(paste0(code_path, "amlap_bb_making_functions.R"))

add_gauss <- function(df, joincol, sd, name,
                      logspace_mean = NA) {
  if (!is.na(logspace_mean))
    sd <- log_space_effect_size_calculator(logspace_mean, sd)

  qname = rlang::enquo(name)
  qcol = rlang::enquo(joincol)

  tmp <- df %>%
    distinct(!!qcol) %>%
    mutate(!!qname := rnorm(!!qcol, mean=0, sd=sd))

  suppressMessages(left_join(df,tmp))
}

add_bb_filename <- function(df, name=NA) {
  if (is.na(name)) {
    df %>%
      mutate(bb_file = paste0(bb_files_path, "bb_", iprefix, "_",
                              Nsubj, "_", Nitems, "_",
                              Bins, "_",
                              MinK, "_to_", MaxK, ".RDS"))
  } else {
    df %>%
      mutate(bb_file = paste0(bb_files_path, "bb_", name, "_",
                              Nsubj, "_", Nitems, "_",
                              Bins, "_",
                              MinK, "_to_", MaxK, ".RDS"))
  }
}
add_mm_filename <- function(df, name=NA) {
  if (is.na(name)) {
    df %>%
      mutate(mm_file = paste0(mm_files_path, "modeldata_",
                              iprefix, "_", Nsubj, "_", Nitems,
                              "_", Bins, "_",
                              SubjSlopeSD, "_", ItemSlopeSD, "_",
                              SubjIntSD, "_", ItemIntSD, "_",
                              RType, "_size_", EffectSize, Space, "_",
                              MinK, "_to_", MaxK, ".RDS"))
  } else {
    df %>%
      mutate(mm_file = paste0(mm_files_path, "modeldata_",
                              name, "_", Nsubj, "_", Nitems,
                              "_", Bins, "_",
                              SubjSlopeSD, "_", ItemSlopeSD, "_",
                              SubjIntSD, "_", ItemIntSD, "_",
                              RType, "_size_", EffectSize, Space, "_",
                              MinK, "_to_", MaxK, ".RDS"))
  }
}
  
# SETAL #############################################################

# Define the mixed models df first, since the bb is a subset
starter_mm_df <- tidyr::crossing(
  Bins = names(new_bins),
  Knum = seq_along(k_list),
  Nitems = c(64, 32),
  Nsubj = c(64, 32),
  iprefix = c("withinboth_w_raneffs_sepremovals"),
  EffectSize = c(7,10), 
  Space = c("_c_log_space", ""),
  SubjIntSD = c(105, 0),
  ItemIntSD = c(42, 0),
  SubjSlopeSD = c(15),
  ItemSlopeSD = c(10),
  RType = c("unresidualized"),
  real_df_file = paste0(setal_file_path, setal_df_file)) %>%
  filter(Nitems==Nsubj) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  add_bb_filename() %>%
  add_mm_filename() %>%
  # Only looking at a single bin and space
  filter(Bins=="bin2") %>%
  # If one of the intercepts is 0, both should be
  filter(!(ItemIntSD==0 & SubjIntSD!=0),
         !(ItemIntSD!=0 & SubjIntSD==0))

# Add in the other effects
full_ranefs_mm_df <- starter_mm_df %>%
  bind_rows(
    .,
    mutate(., iprefix="withinboth_wo_raneffs_sepremovals") %>%
      mutate_at(vars(SubjIntSD:ItemSlopeSD), ~0) %>%
      add_mm_filename()
  )

full_ranefs_bb_df <- full_ranefs_mm_df %>% distinct(bb_file, .keep_all=TRUE)
remaining_full_ranefs_bb_files <- remaining_files(full_ranefs_bb_df, file.exists(bb_file))

# Save BB files
full_ranefs_bb %beep% {
  iterate_over_df(
    remaining_full_ranefs_bb_files,
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
          single_sample_df_function = within_subj_and_item, # same_items_per_subj,
          items_per_subj = Nitems,
          n_subj = Nsubj)
        saveRDS(l, bb_file)
        "good!"
      },
      catchErrors = TRUE)
    })
}
full_ranefs_bb

# Make MM files
remaining_full_ranefs_mm_files <- remaining_files(full_ranefs_mm_df, file.exists(mm_file))
# Saves MM files
full_ranefs_mm %beep% {
  iterate_over_df(
    remaining_full_ranefs_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      bin2meanRT <- setal_basis_info[setal_basis_info$Bins=="bin2", ]$meanRT

      # Set constant to be equal for the middle bin
      log_effect <- log_space_effect_size_calculator(bin2meanRT, effect_size)
      
      model_funcs <- slopes_no_corr_wo_sep_infs # slopes_no_corr_wo_infs

      if (Space=="") {
        adder_f <- function(df) {
          df %>%
            add_gauss(UniqueSubject, SubjIntSD,   SubjIntEffect) %>%
            add_gauss(UniqueItem,    ItemIntSD,   ItemIntEffect) %>%
            add_gauss(UniqueSubject, SubjSlopeSD, SubjSlopeEffect) %>%
            add_gauss(UniqueItem,    ItemSlopeSD, ItemSlopeEffect) %>%
            # IMPORTANTLY THIS ADDS BY-SUBJ & ITEM INTERCEPTS FOR ALLLLLLLL SITUATIONS
            # EVEN FOR TYPE1 STUFFF
            mutate(RT = RT + ItemIntEffect + SubjIntEffect +
                     SimCond * SubjSlopeEffect +
                     SimCond * ItemSlopeEffect,
                   RT_with_Effect = RT + SimCond * effect_size )
        }
      } else {
        adder_f <- function(df) {
          df %>%
            add_gauss(UniqueSubject, SubjIntSD,   SubjIntEffect, bin2meanRT) %>%
            add_gauss(UniqueItem,    ItemIntSD,   ItemIntEffect, bin2meanRT) %>%
            add_gauss(UniqueSubject, SubjSlopeSD, SubjSlopeEffect, bin2meanRT) %>%
            add_gauss(UniqueItem,    ItemSlopeSD, ItemSlopeEffect, bin2meanRT) %>%
            # IMPORTANTLY THIS ADDS BY-SUBJ & ITEM INTERCEPTS FOR ALLLLLLLL SITUATIONS
            # EVEN FOR TYPE1 STUFFF
            mutate(RT = 10^(log10(RT) + ItemIntEffect + SubjIntEffect),
                   RT = add_in_log_space(RT, SimCond, 
                                         (SubjSlopeEffect + ItemSlopeEffect)),
                   RT_with_Effect = add_in_log_space(RT, SimCond, log_effect))
        }
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
full_ranefs_mm

# Loads files
full_ranefs_load %beep% {
  iterate_over_df(
    full_ranefs_mm_df,
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
               real_df_file = real_df_file,
               SubjIntSD,
               ItemIntSD,
               SubjSlopeSD,
               ItemSlopeSD)
    }
  ) %>%
    saveRDS(paste0(complete_path, "setal_w_subjitem_effects.RDS"))
}
full_ranefs_load

