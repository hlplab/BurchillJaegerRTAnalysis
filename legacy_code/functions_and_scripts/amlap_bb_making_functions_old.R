library(dplyr)
library(purrr)
library(future)
library(furrr)
library(cs)

path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
path_on_server <- "/u/zburchil/workspace/launchpad/"

source(paste0(path, "barebones_making_functions.R"))

# Gets a list arguments that are bare filterers and makes a list of bb_dfs based on those
# The arguments must all be named or be a single unnamed one
get_bb_df_list <- function(df, ..., .already_quosure = FALSE) {
  akw <- zplyr::args_and_kwargs(..., .already_quosure = .already_quosure)
  if (length(akw$kwargs) > 0) {
    stopifnot(length(akw$args) == 0)
    stopifnot(n_distinct(names(akw$kwargs)) == length(akw$kwargs))
    named <- TRUE
    zargs = akw$kwargs
  } else {
    stopifnot(length(akw$args) == 1)
    named <- FALSE
    zargs = akw$args
  }

  bb_df_list <- purrr::map(
    zargs,
    function(argv) {
      df %>%
        filter(!!argv) %>%
        make_barebones_sampler_df() })
  if (named==TRUE)
    bb_df_list <- purrr::map2(bb_df_list, names(zargs), ~mutate(.x, Name=.y))
  bb_df_list
}

# This samples for a 'single' bb_df
# Importantly, this takes only the unique subjects
# In the style of Q1/Q3
structured_simple_single_sample <- function(bb_df, unique_subjects, trials_per_subject, subjects_per_group, seed_val) {
  set.seed(seed_val)
  subjects <- tibble(UniqueSubject = sample_subjects_from_list(unique_subjects)) %>%
    mutate(Subject = purrr::map_chr(strsplit(UniqueSubject,  "_"), ~.[1]))

  # Sample items for subjects, such that the first half of the ITEMS for each subject will have -1
  items <- sample_from_barebones(bb_df, trials_per_subject) %>%
    mutate(ItemAmbValue = ifelse(seq_along(GroupingVal) <= length(GroupingVal)/2, -1, 1))

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

# A more flexible framework
make_simple_sample_list <- function(df, k, single_sample_df_function,
                                    ...,
                                    .seed_function = identity,
                                    .multicore = FALSE) {
  # Lets you use the multicore version
  if (.multicore == TRUE) map_function <- furrr::future_map
  else  map_function <- purrr::map

  bb_df <- get_bb_df_list(df, ...) %>%
    bind_rows()
  print(bb_df)
  print(unique(bb_df$Item))
  # Might want to change
  unique_subjects <- unique(df$Subject)

  map_function(k, ~single_sample_df_function(bb_df, unique_subjects, .seed_function(.))) %>%
    purrr::set_names(nm = k)
}

amlap_save_list <- function(server_path,
                            real_df_file,
                            prefix, suffix,
                            k, ...) {
  zplyr::collect_all({
    filename <- paste0(server_path, "amlap_saved_files/bb_", prefix, "_", suffix, "_", min(k), "_to_", max(k), ".RDS")
    real_df <- readRDS(paste0(server_path, real_df_file))

    make_simple_sample_list(df = real_df, k = k, ...) %>%
      saveRDS(filename)
    filename},
    catchErrors = TRUE)
}





