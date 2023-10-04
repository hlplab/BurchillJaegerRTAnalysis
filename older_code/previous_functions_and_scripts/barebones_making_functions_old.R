library(dplyr)
library(lme4)

# How this works (CONCEPTUALLY):
# # The data frame (SETAL or FETAL)
# df %>%
#   # This controls how things are saved on the servers
#   save_list %>%
#   # Depending on the arguments, it calls on of these (similar) functions
#     make_q2_barebones_list / make_q1_q3_barebones_list %>%
#       # These functions first get the small ddf of items and words to sample from
#       make_barebones_sampler_df %>%
#       # And which is used in one of these functions to get a df with no subjects
#       get_q2_single_sample_df / get_q1_q3_single_sample_df %>%
#       # And it samples from the `make_barebones_sampler_df` with
#         sample_from_barebones
#         # And add subjects sampled with:
#         sample_subjects

# Makes the 'barebones' data frame of Items and Words to sample from
make_barebones_sampler_df <- function(df, n=3) {
  df %>%
    group_by(Item) %>%
    summarise(maxWords = max(Word.Order)) %>%
    mutate(Words.Start = purrr::map_chr(maxWords, ~paste0(1:., collapse=","))) %>%
    tidyr::separate_rows(Words.Start, sep=",") %>%
    mutate(Words.Start = as.numeric(Words.Start)) %>%
    group_by(Item) %>%
    filter(Words.Start != maxWords,
           Words.Start != 1,
           Words.Start <= maxWords - n) %>%
    ungroup()
}

# This samples from the barebones df
sample_from_barebones <- function(df, size,
                                  n=3, replace = TRUE, ...) {
  sample_n(df, size = size, replace = replace, ...) %>%
    mutate(Word.Order = purrr::map(Words.Start,
                                   ~c(.:(.+n-1))),
           GroupingVal = seq_along(Word.Order),
           UniqueItem = paste0(Item, "_", GroupingVal))
}

# Samples subjects uniformly
sample_subjects <- function(df, size, replace=TRUE, id_adder_func = ~paste0(., paste0("_", seq_along(.))), ...) {
  uniq_subs <- unique(df$Subject)
  id_adder_func <- purrr::as_mapper(id_adder_func)
  id_adder_func(sample(uniq_subs, size = size, replace = replace, ...))
}

get_q1_q3_single_sample_df <- function(bb_df, real_df, trials_per_subject, subjects_per_group, seed_val) {
  set.seed(seed_val)
  single_block <- sample_from_barebones(bb_df, trials_per_subject)
  blocks_for_group <- do.call("rbind", rep(list(single_block), subjects_per_group))

  both_blocks <- bind_rows(mutate(blocks_for_group, SimCond = 1),
                           mutate(blocks_for_group, SimCond = -1))

  subjects <- sample_subjects(real_df, subjects_per_group*2)
  both_blocks$UniqueSubject <- rep(subjects, each=trials_per_subject)
  both_blocks$Subject <- purrr::map_chr(strsplit(both_blocks$UniqueSubject, "_"),
                                        ~.[1])

  both_blocks$GroupingVal <- seq_along(both_blocks$GroupingVal)
  return(both_blocks)
}

make_q1_q3_barebones_list <- function(df, k, block_filter,
                                      effect_size, # a list of constant value
                                      trials_per_subject = 10,
                                      subjects_per_group = 208,
                                      seed_function = identity) {

  if (rlang::is_bare_numeric(effect_size, n=1))
    effect_size <- rep(effect_size, length(k))

  bb_df <- df %>%
    filter(Block==block_filter) %>%
    make_barebones_sampler_df()

  purrr::map(k,
             ~get_q1_q3_single_sample_df(bb_df, df,
                                         trials_per_subject = trials_per_subject,
                                         subjects_per_group = subjects_per_group,
                                         seed_function(.))) %>%
    purrr::set_names(nm = k) %>%
    purrr::map2(effect_size,
                ~`attr<-`(.x, "effect_size", .y))
}

# First block is coded as positive
get_q2_single_sample_df <- function(bb_block1_df, trials_per_block1,
                                    bb_block2_df, trials_per_block2,
                                    n_subjects,
                                    real_df, seed_val) {
  set.seed(seed_val)
  block1_single  <- sample_from_barebones(bb_block1_df, trials_per_block1)
  block2_single <- sample_from_barebones(bb_block2_df, trials_per_block2)

  subjects <- sample_subjects(real_df %>% filter(Group=="Filler-first"), n_subjects)

  first_blocks <-  do.call("rbind", rep(list(block1_single), n_subjects))
  second_blocks <- do.call("rbind", rep(list(block2_single), n_subjects))

  first_blocks$UniqueSubject <-  rep(subjects, each = trials_per_block1)
  second_blocks$UniqueSubject <- rep(subjects, each = trials_per_block2)

  first_blocks$Subject <-  purrr::map_chr(strsplit(first_blocks$UniqueSubject,  "_"), ~.[1])
  second_blocks$Subject <- purrr::map_chr(strsplit(second_blocks$UniqueSubject, "_"), ~.[1])

  both_blocks  <- bind_rows(mutate(first_blocks, SimCond = 1),
                            mutate(second_blocks, SimCond = -1))

  both_blocks$GroupingVal <- seq_along(both_blocks$GroupingVal)
  return(both_blocks)
}

make_q2_barebones_list <- function(df, k,
                                   effect_size, # a list of constant value
                                   trials_per_block1 = 16,
                                   trials_per_block2 = 9,
                                   n_subjects = 205,
                                   seed_function = identity) {
  if (rlang::is_bare_numeric(effect_size, n=1))
    effect_size <- rep(effect_size, length(k))

  bb1_df <- df %>%
    filter(Block == "1") %>%
    make_barebones_sampler_df()

  bb2_df <- df %>%
    filter(Block == "2") %>%
    make_barebones_sampler_df()

  purrr::map(k,
             ~get_q2_single_sample_df(bb_block1_df = bb1_df,
                                      trials_per_block1 = trials_per_block1,
                                      bb_block2_df = bb2_df,
                                      trials_per_block2 = trials_per_block2,
                                      real_df = df,
                                      n_subjects = n_subjects,
                                      seed_function(.))) %>%
    purrr::set_names(nm = k) %>%
    purrr::map2(effect_size,
                ~`attr<-`(.x, "effect_size", .y))
}

save_list <- function(server_path,
                      real_df_file,
                      thing_id, question, k, effect_size, ...) {
  zplyr::collect_all({
    real_df <- readRDS(paste0(server_path, real_df_file))
    {
      if (question=="q2") {
        make_q2_barebones_list(df = real_df, k = k,
                               effect_size = effect_size, ...)
      }
      else if (question %in% c("q1", "q3")) {
        make_q1_q3_barebones_list(df = real_df, k = k,
                                  effect_size = effect_size, ...)
      }
    } %>%
      saveRDS(paste0(server_path, "saved_files/bb_", thing_id, "_", question, "_", min(k), "_to_", max(k), ".RDS"))
    paste0(server_path, "saved_files/bb_", thing_id, "_", question, "_", min(k), "_to_", max(k), ".RDS")},
    catchErrors = TRUE)
}


# # Adds the unique variables (used for Item and Subject)
# # uhhhhhhh, actually I wrote all this and realized that Items aren't added the way subjects are, so I'm just not going to use it
# add_unique_variables <- function(df, v, bare_var_name, separator="_") {
#   short = rlang::quo_name(rlang::enquo(bare_var_name))
#   long = paste0("Unique", short)
#   print(long)
#   if (long %in% names(df) | short %in% names(df))
#     stop(paste0("`", long, "` or `", short, "` already named in the data frame!"))
#   df %>%
#     mutate(!! long := v,
#            !! short :=  purrr::map_chr(strsplit(!! sym(long), "_"), ~.[1]))
# }


