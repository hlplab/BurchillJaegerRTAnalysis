library(dplyr)
library(lme4)

# current print and pass

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

remove_nas <- function(x) {
  x[!rlang::is_na(x)]
}

# Makes the 'barebones' data frame of Items and Words to sample from
make_barebones_sampler_df <- function(df, n=3, ..., keep_poss_subj = FALSE) {
  bb_df <- df %>%
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
  
  if (keep_poss_subj==TRUE)
    bb_df <- add_possible_subjects(bb_df, df, n)
  
  bb_df
}

# This adds the possible subjects into the data frame
add_possible_subjects <- function(bb_df, real_df, n) {
  orig_names <- names(bb_df)
  purrr::reduce(
    0:(n-1), 
    ~bind_rows(.x, mutate(bb_df, All.Words = Words.Start + .y)) %>%
      distinct(),
    .init = tibble()
    ) %>%
    left_join(select(real_df, Subject, Item, Word.Order),
              by = c("Item", All.Words = "Word.Order")) %>%
    group_by_at(vars(orig_names)) %>%
    summarise(PossibleSubjects = list(remove_nas(unique(Subject)))) %>%
    ungroup()
}

# This samples from the barebones df
sample_from_barebones <- function(df, size,
                                  n=3, replace = TRUE, ...,
                                  starting_groupval = 1 # For processes that sample more than once
                                  ) {
  sample_n(df, size = size, replace = replace, ...) %>%
    mutate(Word.Order = purrr::map(Words.Start,
                                   ~c(.:(.+n-1))),
           GroupingVal = seq_along(Word.Order) + starting_groupval - 1,
           UniqueItem = paste0(Item, "_", GroupingVal))
}

# Samples subjects uniformly
sample_subjects <- function(df, size, replace=TRUE, id_adder_func = ~paste0(., paste0("_", seq_along(.))), ...) {
  uniq_subs <- unique(df$Subject)
  id_adder_func <- purrr::as_mapper(id_adder_func)
  id_adder_func(sample(uniq_subs, size = size, replace = replace, ...))
}

# For future code. I don't know why I needed to pass the whole df rather than just unique subjects...
# Samples subjects uniformly
sample_subjects_from_list <- function(uniq_subs, size, replace=TRUE, id_adder_func = ~paste0(., paste0("_", seq_along(.))), ...) {
  stopifnot(!is.null(size))
  id_adder_func <- purrr::as_mapper(id_adder_func)
  id_adder_func(sample(uniq_subs, size = size, replace = replace, ...))
}



get_q1_q3_single_sample_df <- function(bb_df, real_df, trials_per_subject, subjects_per_group, seed_val) {
  set.seed(seed_val)
  # The order and way things are declared here is to keep the randomization comparable to previous runs.
  # Sorrrrry!
  # Sample items for subjects, such that the first half of the ITEMS for each subject will have -1
  items <- sample_from_barebones(bb_df, trials_per_subject) %>%
    mutate(ItemAmbValue = ifelse(seq_along(GroupingVal) <= length(GroupingVal)/2, -1, 1))
  # This is in this weird format to keep the randomization the same as it was before...
  # This first samples ONCE to get AL subjects
  subjects <- tibble(UniqueSubject = sample_subjects(real_df, subjects_per_group*2)) %>%
    mutate(Subject = purrr::map_chr(strsplit(UniqueSubject,  "_"), ~.[1]))

  # Sample subjects for groups
  group1_single  <- subjects[1:subjects_per_group,] %>%
    # Make the first half of subjects each GROUP have a -1 AmbVal and the second half +1
    mutate(SubjectAmbValue = ifelse(seq_along(UniqueSubject) <= subjects_per_group/2, -1, 1))
  group2_single  <- subjects[(subjects_per_group+1):(subjects_per_group*2),] %>%
    mutate(SubjectAmbValue = ifelse(seq_along(UniqueSubject) <= subjects_per_group/2, -1, 1))

  # Repeat the subject dataframe so there are the right number of trials per subject
  group1_expanded <- group1_single[rep(seq(subjects_per_group), each = trials_per_subject),]
  group2_expanded <- group2_single[rep(seq(subjects_per_group), each = trials_per_subject),]

  first_group <-  do.call("rbind", rep(list(items), subjects_per_group)) %>%
    bind_cols(group1_expanded)
  second_group <- do.call("rbind", rep(list(items), subjects_per_group)) %>%
    bind_cols(group2_expanded)

  both_groups <- bind_rows(mutate(first_group,  SimCond = 1),
                           mutate(second_group, SimCond = -1))

  both_groups$GroupingVal <- seq_along(both_groups$GroupingVal)
  return(both_groups)
}

make_q1_q3_barebones_list <- function(df, k, block_filter,
                                      effect_size, # a three-value vector
                                      trials_per_subject = 20,
                                      subjects_per_group = 208,
                                      seed_function = identity,
                                      multicore = FALSE) {
  if (length(effect_size) != 3)
    stop("Effect size needs to be a vector of length 3")
  effect_size <- rep(list(effect_size), length(k))
  # Lets you use the multicore version
  # if (multicore == TRUE)
  map_function <- furrr::future_map
  # else
  #   map_function <- purrr::map

  # Technically, worriable in the sense that SETAL's Blocks are ints, and FETAL's Blocks are chars. Fortunately, `2 == "2"2` is true in R. But I probably should change this...
  bb_df <- df %>%
    filter(Block==block_filter) %>%
    make_barebones_sampler_df()

  map_function(k,
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
  # Sample from blocks
  block1_single  <- sample_from_barebones(bb_block1_df, trials_per_block1) %>%
    # Make the first half of each block have a -1 AmbVal and the second half +1
    mutate(ItemAmbValue = ifelse(seq_along(GroupingVal) <= trials_per_block1/2, -1, 1))
  block2_single <- sample_from_barebones(bb_block2_df, trials_per_block2) %>%
    mutate(ItemAmbValue = ifelse(seq_along(GroupingVal) <= trials_per_block2/2, -1, 1))

  # Sample subjects only from the filler first condition
  subject_list <- sample_subjects(real_df %>% filter(Group=="Filler-first"), n_subjects)
  # Make the first half of the subjects have a -1 AmbVal and the second half have +1
  subjects <- tibble(UniqueSubject = subject_list) %>%
    mutate(SubjectAmbValue = ifelse(seq_along(UniqueSubject) <= length(UniqueSubject)/2, -1, 1),
           Subject = purrr::map_chr(strsplit(UniqueSubject,  "_"), ~.[1]))

  block1_subjects <- subjects[rep(seq(n_subjects), each = trials_per_block1),]
  block2_subjects <- subjects[rep(seq(n_subjects), each = trials_per_block2),]

  first_blocks <-  do.call("rbind", rep(list(block1_single), n_subjects)) %>%
    bind_cols(block1_subjects)
  second_blocks <- do.call("rbind", rep(list(block2_single), n_subjects)) %>%
    bind_cols(block2_subjects)

  both_blocks  <- bind_rows(mutate(first_blocks, SimCond = 1),
                            mutate(second_blocks, SimCond = -1))

  both_blocks$GroupingVal <- seq_along(both_blocks$GroupingVal)
  return(both_blocks)
}

make_q2_barebones_list <- function(df, k,
                                   effect_size, # a three-value vector
                                   trials_per_block1 = 32,
                                   trials_per_block2 = 18,
                                   n_subjects = 205,
                                   seed_function = identity,
                                   multicore = FALSE) {
  if (length(effect_size) != 3)
    stop("Effect size needs to be a vector of length 3")
  effect_size <- rep(list(effect_size), length(k))
  # Lets you use the multicore version
  if (multicore == TRUE)
    map_function <- furrr::future_map
  else
    map_function <- purrr::map

  bb1_df <- df %>%
    filter(Block == "1") %>%
    make_barebones_sampler_df()

  bb2_df <- df %>%
    filter(Block == "2") %>%
    make_barebones_sampler_df()

  map_function(k,
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


