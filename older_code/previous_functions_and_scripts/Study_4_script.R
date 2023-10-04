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
setal_int_mm_df <- tidyr::crossing(
  Bins = c("interaction"), #since interaction is between block1 and block3
  Knum = seq_along(k_list),
  Nsubj =  c(8, 16, 32, 64),
  Nitems = c(8, 16, 32, 64),
  MainEffectSize = c(56),
  InteractionEffectSize = c(28),
  Block1Code = c(-1, 1),
  Interaction = c("real"), #   "fake"),
  Space = c("", "_c_log_space"),
  RType = c("unresidualized"),
  iprefix = c("same_across_subj_new")) %>%
  filter(Nsubj == Nitems) %>%
  mutate(direction_name = ifelse(Block1Code==-1, "neg",
                          ifelse(Block1Code==1,  "pos", NA_character_))) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]])),
         bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_", RType, "_",
                          Interaction, "_", direction_name, "_",
                          "_main_", MainEffectSize, "_second_", InteractionEffectSize,
                          Space, "_", MinK, "_to_", MaxK, ".RDS"),
         real_df_file = ifelse(grepl("noexcl", iprefix),
                               paste0(path_on_server, "setal_noexcl_fillers.RDS" ),
                               paste0(path_on_server, setal_file)))# %>%
  # filter(Space == "")

setal_int_bb_df <- setal_int_mm_df %>%
  select(-Space, -MainEffectSize, -Block1Code, -direction_name,
         -InteractionEffectSize, -RType, -mm_file) %>% distinct()
remaining_setal_int_bb_files <- remaining_files(setal_int_bb_df, file.exists(bb_file))

# Save BB files
setal_int_bb %beep% {
  iterate_over_df(
    remaining_setal_int_bb_files,
    future_cnd_map,
    function(zaza) {
      k_e <- MinK:MaxK
      bin1_low <- new_bins$bin1[1]
      bin1_high <- new_bins$bin1[2]
      bin3_low <- new_bins$bin3[1]
      bin3_high <-new_bins$bin3[2]

      filterers <- quos(
        Block1 =  Trial.Order  >= !!enquo(bin1_low) &
          Trial.Order <= !!enquo(bin1_high) &
          Group == "Filler-first",
        Block3 =  Trial.Order  >= !!enquo(bin3_low) &
          Trial.Order <= !!enquo(bin3_high) &
          Group == "Filler-first")

      zplyr::collect_all({
        real_df <- readRDS(real_df_file)
        l <- make_simple_sample_list(
          df = real_df, k = k_e,
          filter_quosures = filterers,
          # The function and its named additional arguments
          single_sample_df_function = same_items_per_subj_interaction,
          n_trials = Nitems,
          n_subj = Nsubj)
        saveRDS(l, bb_file)
        "good!"
      },
      catchErrors = TRUE)
    })
} %plan% elite_squad_plan
done(setal_int_bb)


# Make MM files
remaining_setal_int_mm_files <- remaining_files(setal_int_mm_df, file.exists(mm_file))
# Saves MM files
setal_int_mm %beep% {
  iterate_over_df(
    remaining_setal_int_mm_files,
    future_cnd_map,
    function(zaa) {
      # Set up th effect size and log_effectsize
      effect_size <- MainEffectSize
      interaction_effect <- InteractionEffectSize
      # Currently uses mean of Bin 2. COULD make it so that it uses different effects in different blocks
      log_effect <- log_space_effect_size_calculator(
        mean(setal_basis_info[setal_basis_info$Bins=="bin2", ]$meanRT),
        effect_size)
      log_interaction_effect <- log_space_effect_size_calculator(
        mean(setal_basis_info[setal_basis_info$Bins=="bin2", ]$meanRT),
        interaction_effect)

      # What kind of interactions are being looked at
      if (Interaction=="real")
        model_funcs <- i_mm_block_slopes
      else if (Interaction=="fake")
        stop("Not coded yet!")
      else stop("Invalid interaction type")

      # Establish adding functions
      if (Space=="")
        adder_f <- function(df) {
          df %>%
            # whether Block1 is coded as +1 or -1
            mutate(Block = factor(BlockName, levels=c("Block1","Block3")) %>%
                     `contrasts<-`(value = c(Block1Code, -Block1Code)),
                   BlockDbl = case_when(Block=="Block1" ~ Block1Code,
                                        Block=="Block3" ~ -Block1Code,
                                        TRUE ~ NA_real_),
                   RT_with_ALL = RT + SimCond * effect_size + BlockDbl * SimCond * interaction_effect,
                   RT_with_Interaction = RT + BlockDbl * SimCond * interaction_effect,
                   RT_with_Main =  RT + SimCond * effect_size) }
      if  (Space=="_c_log_space")
        adder_f <- function(df) {
          mutate(df, # whether Block1 is coded as +1 or -1
                 Block = factor(BlockName, levels=c("Block1","Block3")) %>%
                   `contrasts<-`(value = c(Block1Code, -Block1Code)),
                 BlockDbl = case_when(Block=="Block1" ~ Block1Code,
                                      Block=="Block3" ~ -Block1Code,
                                      TRUE ~ NA_real_),
                 RT_with_ALL = 10^(log10(RT) + SimCond * log_effect + BlockDbl * SimCond * log_interaction_effect),
                 RT_with_Interaction = 10^(log10(RT) + BlockDbl * SimCond * log_interaction_effect),
                 RT_with_Main =  10^(log10(RT) + SimCond * log_effect)) }

      # New way
      zplyr::collect_all(
        amlap_make_models(
          bb_file,
          mm_file,
          real_df_file,
          adder_f, model_funcs, mixed_cleaner),
        catchErrors = TRUE)
    })
} %plan% elite_squad_plan
done(setal_int_mm)

# Loads files
setal_int_load %beep% {
  iterate_over_df(
    setal_int_mm_df,
    future_map_dfr,
    function(i) {
      collate_data(mm_file,
                   terms_to_harvest = c("SimCond", "Block1", "SimCond:Block1"),
                   .multicore = TRUE) %>%
        mutate(Nsubj, Nitems, Bins, MainEffectSize,
               InteractionEffectSize, Interaction,
               Block1Code, RType, Space,
               ifile = bb_file,
               ofile = mm_file,
               real_df_file = real_df_file)
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/interaction_mixed_models_lmerTest_new.RDS"))
  cat("rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/interaction_mixed_models_lmerTest_new.RDS ~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/PowerSimulations/data/powerpaper/")
} %plan% elite_squad_plan
done(setal_int_load)

#
#
# temp_thing %beep% {
#   readRDS(setal_int_mm_df$mm_file[[1]])[[1]]
# }
#
# gurp <- list(temp_thing) %>%
#   df_list_to_df(terms_to_harvest = c("SimCond", "Block1", "SimCond:Block1"))
#
# df_list_to_df <- function(df,
#                           terms_to_harvest = c("SimCond"),
#                           .multicore = FALSE) {
#   # if (.multicore == FALSE) map_func <- purrr::imap_dfr
#   # else
#   map_func <- furrr::future_imap_dfr
#
#   df %>%
#     map_func(
#       function(dfl, iter) {
#         purrr::imap_dfr(
#           dfl,
#           function(single, name_type) {
#             if (is.na(single$value)) {
#               new_df <- data.frame(estimate=NA)
#             } else {
#               new_df <- single$value %>%
#                 filter(term %in% terms_to_harvest)
#             }
#             new_df %>%
#               mutate(messages = paste0(single$messages, collapse="--------"),
#                      warnings = paste0(single$warnings, collapse="--------"),
#                      errors = paste0(single$errors, collapse="--------"),
#                      type = name_type)
#           }) %>%
#           mutate(k = as.numeric(iter))
#       })
# }
#
#
#
#
#
#
#
#
#
#
