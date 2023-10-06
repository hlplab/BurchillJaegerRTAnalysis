library(dplyr)
library(purrr)
library(future)
library(furrr)
library(cs)

source(paste0(code_path, "barebones_making_functions.R"))

# This takes a data frame, a list of quosures to filter by, and a function to apply to those filtered data frames
function_after_filter_by_quosures <- function(df, quosures, fn) {
  if (is.null(names(quosures)))
    q_names <- as.list(rep("", length(quosures)))
  else
    q_names <- names(quosures)
    
  bb_df_list <- purrr::map(
    quosures,
    function(q) {
      filter(df, !!q) %>%
        fn() 
    })
  
  purrr::map2(bb_df_list, q_names, ~mutate(.x, Name=.y))
}



# Gets a list arguments that are bare filterers and makes a list of bb_dfs based on those
# The arguments must all be named or be a single unnamed one
# The `Name` column will refer to the name of the quosure if the quosure list is named
get_bb_df_list <- function(df, quosures, n=1, add_poss_subj = TRUE) {
  function_after_filter_by_quosures(
    df, 
    quosures, 
    function(x) make_barebones_sampler_df(x, n=n, keep_poss_subj = add_poss_subj))
}

# gets the largest word order in a set of items
get_order_range <- function(bb_df) {
  bb_df %>% group_by(Item) %>%
    summarise(min=min(Words.Start),
              max=max(Words.Start)) %>%
    summarise(min=max(min),
              max=min(max)) %>%
              {c(.$min:.$max)}
}


# add_possible_subjects <- function(bb_df, real_df, filter_quosures, n) {
#   orig_names <- names(bb_df)
#   real_df <- filter(real_df,)
#   # Should probably use 'reduce' to keep memory lower
#   purrr::map_dfr(0:(n-1), ~mutate(bb_df, All.Words = Words.Start + .x)) %>%
#     distinct() %>%
#     left_join(select(real_df, Subject, Item, Word.Order), 
#               by = c("Item", All.Words = "Word.Order")) %>% 
#     group_by_at(vars(orig_names)) %>%
#     summarise(PossibleSubjects = list(remove_nas(unique(Subject)))) %>%
#     ungroup()
# }

# Should only be used when the `n` of the initial bb_df is 1
# Essentially pools the possible subjects between regions 
#   (i.e., if a subject has only answered a subset of the region, they're still included)
readjust_possible_subjects <- function(bb_df, new_n,
                                       join_by_s = c("Item", "Name")) {
  stopifnot("PossibleSubjects" %in% names(bb_df))
  
  loading_bb_df <- bb_df %>% 
    mutate(Region = purrr::map_chr(Words.Start, ~paste0(.x:(.x + new_n - 1), collapse=","))) %>%
    tidyr::separate_rows(Region, sep=",") %>%
    mutate(Region = as.numeric(Region)) %>%
    select(one_of(join_by_s), maxWords, Words.Start, Region)
  
  loading_bb_df %>% 
    left_join(select(bb_df, -maxWords), 
              by=c(join_by_s, "Region"="Words.Start")) %>%
    filter(Words.Start <= maxWords - new_n) %>%
    group_by_at(vars(-PossibleSubjects, -Region)) %>%
    summarise(PossibleSubjects = list(unique(unlist(PossibleSubjects)))) %>%
    ungroup()
}




# A more flexible framework
# `single_sample_df_function` needs to have `bb_df`, `unique_subjects`, and `seed_val` as arguments
make_simple_sample_list <- function(df, k,
                                    filter_quosures,
                                    single_sample_df_function,
                                    ...,
                                    .seed_function = identity,
                                    .multicore = FALSE, #deprecated
                                    n = 1
) {
  
  bb_df <- get_bb_df_list(df, filter_quosures, n=n, add_poss_subj = TRUE) %>%
    bind_rows()
  
  # Might want to change
  unique_subjects <- unique(df$Subject)
  
  furrr::future_map(k, function(k_e) {
    single_sample_df_function(bb_df = bb_df,
                              unique_subjects = unique_subjects,
                              seed_val = .seed_function(k_e),
                              ...) }) %>%
    purrr::set_names(nm = k)
}

# Loads everything on the server, and runs the function that will make the samples
amlap_save_list <- function(server_path,
                            real_df_file,
                            prefix, suffix,
                            k, ...,
                            .toSave = TRUE
                            ) {
  if ("n" %in% rlang::call_args_names(sys.call()))
    warning("Remember, `n` has special meaning and is meant to be used exclusively for the number of words in an item")
  
  zplyr::collect_all({
    filename <- paste0(server_path, "amlap_saved_files/bb_", prefix, "_", suffix, "_", min(k), "_to_", max(k), ".RDS")
    real_df <- readRDS(paste0(server_path, real_df_file))
    l <- make_simple_sample_list(df = real_df, k = k, ...)
    if (.toSave == TRUE) {
      saveRDS(l, filename)
      filename 
    } else if (rlang::is_string(.toSave)) {
      saveRDS(l, .toSave)
      .toSave
    } else {
      l
    }
  },
  catchErrors = TRUE)
}


# DIFFERENT BEHAVIOR THAN `amlap_save_list`
# Makes a model from a data frame
run_parametric_model_maker <- function(df_file_name, filterers, formula, ofile,
                                       data_f = identity) {
  df <- readRDS(df_file_name)
  df <- df %>% filter(!!!filterers) %>% droplevels() %>%
    data_f()
  m <- lme4::lmer(formula = formula, data = df)
  saveRDS(m, ofile)
}



# This just combines the first two halves of each type of SimCond from two lists
# It assumes they have whole numbers and are equally divided
amlap_combine_two_bbs <- function(bb_1, bb_2) {
  f <- function(df) sample(1:(nrow(df)/2), size=(nrow(df)/4), replace=F)
  
  bind_rows("First" = bb_1 %>%
              group_by(SimCond) %>%
              filter(seq_along(Item) %in% f(bb_1)) %>%
              ungroup(),
            "Second" = bb_2 %>%
              group_by(SimCond) %>%
              filter(seq_along(Item) %in% f(bb_2)) %>%
              ungroup(),
            .id="NewGroup")
}

# Sim effect is constant within pseudo-subjects
# Not guaranteed to have full size
# Technically, the sampling population for subjects is any subject who has any response to ANY
#   item in the item sampling population.
# DEFAULT IS `n` = 1!!!!!!!!!!!!!!!!!!!!!
same_items_per_subj <- function(bb_df, items_per_subj, n_subj, seed_val, ..., words_per_region = 1) {
  # set.seed(seed_val)
  
  # Used to be just `n = 1`
  items_df <-  sample_from_barebones(bb_df, size = items_per_subj, n = words_per_region)
  
  # This is just a way of getting all the subjects who have any responses and sampling from them
  subjects <- unique(purrr::reduce(
    bb_df$PossibleSubjects,
    # items_df$PossibleSubjects, # Used to be this
    ~c(.x, as.character(.y)), .init=NULL)) %>%
    remove_nas() %>%
    sample(size = n_subj, replace = TRUE)
  
  purrr::map_dfr(
    1:n_subj, 
    ~mutate(items_df, Subject = subjects[[.]],
            UniqueSubject = paste0(Subject, "_", .),
            SimCond = ifelse(. <= n_subj/2, -1, 1))) %>%
    mutate(ItemAmbValue = 0,
           SubjectAmbValue = 0) %>%
    select(-PossibleSubjects)
}


# All subjects get the same items; 2 groups of subjects, alternate simcond for the items; both groups have half + half -
within_subj_and_item <- function(bb_df, items_per_subj, n_subj, seed_val, ..., words_per_region = 1) {
  # set.seed(seed_val)
  
  # Used to be just `n = 1`
  items_df <-  sample_from_barebones(bb_df, size = items_per_subj, n = words_per_region) %>%
    mutate(ItemGrouper = ifelse(seq_along(UniqueItem) <= items_per_subj/2,-1,1))
  
  # This is just a way of getting all the subjects who have any responses and sampling from them
  subjects <- unique(purrr::reduce(
    bb_df$PossibleSubjects,
    # items_df$PossibleSubjects, # Used to be this
    ~c(.x, as.character(.y)), .init=NULL)) %>%
    remove_nas() %>%
    sample(size = n_subj, replace = TRUE)
  
  purrr::map_dfr(
    1:n_subj, 
    ~mutate(items_df, Subject = subjects[[.]],
            UniqueSubject = paste0(Subject, "_", .),
            SimCond = ifelse(. <= n_subj/2, -1, 1))) %>%
    mutate(ItemAmbValue = 0,
           SubjectAmbValue = 0) %>%
    mutate(SimCond = SimCond*ItemGrouper) %>%
    select(-PossibleSubjects,-ItemGrouper)
}



# A wrapper function for `same_items_per_subj` that sets critical regions
# Only to be used for F13 and HS18 data
same_items_per_subj_with_fillers <- function(bb_df, items_per_subj, n_subj, seed_val, ..., words_per_region = 1, n_critical_items) {
  bb_df <- same_items_per_subj(
    bb_df = bb_df, items_per_subj = items_per_subj, n_subj = n_subj, 
    seed_val = seed_val, ..., words_per_region = words_per_region) %>%
    mutate(ItemType = ifelse(UniqueItem %in% unique(UniqueItem)[1:n_critical_items],
                             "CriticalRegion", "FillerItems"))
  
}

# Randomly picks a n-word region shared by ALL items and then samples items
# Probably only works with `n` = 1
region_then_items_same_items_per_subj <- function(bb_df, regions_per_subj, n_subj, seed_val, words_per_region, ..., exclude_items_with_fewer_than = NULL) {
  # set.seed(seed_val)
  
  if (is.null(exclude_items_with_fewer_than)) 
    exclude_items_with_fewer_than <- -Inf
  
  # Just some sanity checking
  stopifnot("maxWords" %in% names(bb_df))
  one_max_word_per_item <- bb_df %>% group_by(Item) %>% 
    summarise(n=n_distinct(maxWords)) %>% {all(.$n == 1)}
  stopifnot(one_max_word_per_item)
  
  # If you want to exclude items without enough words, you can
  bb_df <- filter(bb_df, maxWords >= exclude_items_with_fewer_than)
  
  # Readjusting `PossibleSubjects` from `n==1` to `n==words_per_region`
  new_bb_df <- readjust_possible_subjects(bb_df, new_n = words_per_region)
  
  # Randomly selecting the region and filtering
  possible_region_starts <- get_order_range(new_bb_df)
  warning("Possible item regions: ", length(possible_region_starts))
  region_start <- sample(possible_region_starts, size = 1)
  new_bb_df <- filter(new_bb_df, Words.Start == region_start)
  
  same_items_per_subj(bb_df = new_bb_df, 
                      items_per_subj = regions_per_subj, 
                      n_subj = n_subj, 
                      seed_val = seed_val,
                      words_per_region = words_per_region)
}

# Samples items and THEN randomly picks a n-word region shared by of the sampled items
# UNTESTED AND UNUSED SO FAR. ALMOST CERTAINLY GARBAGE. 
items_then_region_same_items_per_subj <- function(bb_df, regions_per_subj, n_subj, seed_val, words_per_region, ...) {
  # set.seed(seed_val)
  
  # Readjusting `PossibleSubjects` from `n==1` to `n==words_per_region`
  bb_df <- readjust_possible_subjects(bb_df, new_n = words_per_region)
  
  items <- sample(bb_df$Item, size = regions_per_subj, replace = TRUE)
  start_range <- get_order_range(filter(bb_df, Item %in% items))
  region <- sample(start_range, size=1)
  
  items_df <- items %>% 
    purrr::map_dfr(~filter(bb_df, Item==.x, Words.Start == region)) %>%
    mutate(Word.Order = purrr::map(Words.Start, ~c(.:(.+words_per_region-1))),
           GroupingVal = seq_along(Word.Order),
           UniqueItem = paste0(Item, "_", GroupingVal))
  
  subjects <- unique(purrr::reduce(
    bb_df$PossibleSubjects,
    # items_df$PossibleSubjects, # Used to be this
    ~c(.x, as.character(.y)), .init=NULL)) %>%
    remove_nas() %>%
    sample(size = n_subj, replace = TRUE)
  
  purrr::map_dfr(
    1:n_subj, 
    ~mutate(items_df, Subject = subjects[[.]],
            UniqueSubject = paste0(Subject, "_", .),
            SimCond = ifelse(. <= n_subj/2, -1, 1))) %>%
    mutate(ItemAmbValue = 0,
           SubjectAmbValue = 0) %>%
    select(-PossibleSubjects)
}



same_items_per_subj_for_residualized_regions <- function(bb_df, critical_regions_per_subj, filler_items_per_subj, n_subj, seed_val, words_per_region, ..., exclude_items_with_fewer_than = NULL) {
  
  region_sim_df <- region_then_items_same_items_per_subj(
    bb_df = bb_df, 
    regions_per_subj = critical_regions_per_subj, 
    n_subj = n_subj, 
    seed_val = seed_val, 
    words_per_region = words_per_region, 
    exclude_items_with_fewer_than = exclude_items_with_fewer_than)
  
  # HACK
  chosen_subjects <- distinct(region_sim_df, Subject, UniqueSubject)$Subject
  
  filler_items_df <-  sample_from_barebones(
    df = bb_df, 
    size = filler_items_per_subj, 
    n = 1,
    starting_groupval = critical_regions_per_subj+1) 

  filler_sim_df <- purrr::map_dfr(
    1:n_subj, 
    ~mutate(filler_items_df, Subject = chosen_subjects[[.]],
            UniqueSubject = paste0(Subject, "_", .),
            SimCond = 0)) %>%
    mutate(ItemAmbValue = 0, SubjectAmbValue = 0) %>%
    select(-PossibleSubjects)
  
  stopifnot(setequal(region_sim_df$UniqueSubject, filler_sim_df$UniqueSubject))
  
  bind_rows("CriticalRegion" = region_sim_df,
            "FillerItems" = filler_sim_df, 
            .id = "ItemType")
}
# same_items_per_subj_for_residualized_regions(
#   bb_df = bb_df, critical_regions_per_subj = 4, filler_items_per_subj = 2, n_subj = 2, seed_val = 1, words_per_region = 3, exclude_items_with_fewer_than = 10)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

# like `same_items_per_subj_for_residualized_regions` but lets you exclude from from
#   front and end of possible critical regions
same_items_per_subj_for_residualized_regions_custom <- function(bb_df, critical_regions_per_subj, filler_items_per_subj, n_subj, seed_val, words_per_region, ..., exclude_items_with_fewer_than = NULL, exclude_from_front = 0, exclude_from_back = 0) {
  
  region_sim_df <- region_then_items_same_items_per_subj_custom(
    bb_df = bb_df, 
    regions_per_subj = critical_regions_per_subj, 
    n_subj = n_subj, 
    seed_val = seed_val, 
    words_per_region = words_per_region, 
    exclude_items_with_fewer_than = exclude_items_with_fewer_than,
    exclude_from_front = exclude_from_front,
    exclude_from_back = exclude_from_back)
  
  # HACK
  chosen_subjects <- distinct(region_sim_df, Subject, UniqueSubject)$Subject
  
  filler_items_df <-  sample_from_barebones(
    df = bb_df, 
    size = filler_items_per_subj, 
    n = 1,
    starting_groupval = critical_regions_per_subj+1) 
  
  filler_sim_df <- purrr::map_dfr(
    1:n_subj, 
    ~mutate(filler_items_df, Subject = chosen_subjects[[.]],
            UniqueSubject = paste0(Subject, "_", .),
            SimCond = 0)) %>%
    mutate(ItemAmbValue = 0, SubjectAmbValue = 0) %>%
    select(-PossibleSubjects)
  
  stopifnot(setequal(region_sim_df$UniqueSubject, filler_sim_df$UniqueSubject))
  
  bind_rows("CriticalRegion" = region_sim_df,
            "FillerItems" = filler_sim_df, 
            .id = "ItemType")
}


# Randomly picks a n-word region shared by ALL items and then samples items
# Probably only works with `n` = 1
region_then_items_same_items_per_subj_custom <- function(bb_df, regions_per_subj, n_subj, seed_val, words_per_region, ..., exclude_items_with_fewer_than = NULL, exclude_from_front = 0, exclude_from_back = 0) {
  
  if (is.null(exclude_items_with_fewer_than)) 
    exclude_items_with_fewer_than <- -Inf
  
  # Just some sanity checking
  stopifnot("maxWords" %in% names(bb_df))
  one_max_word_per_item <- bb_df %>% group_by(Item) %>% 
    summarise(n=n_distinct(maxWords)) %>% {all(.$n == 1)}
  stopifnot(one_max_word_per_item)
  
  # If you want to exclude items without enough words, you can
  bb_df <- filter(bb_df, maxWords >= exclude_items_with_fewer_than)
  
  # Readjusting `PossibleSubjects` from `n==1` to `n==words_per_region`
  new_bb_df <- readjust_possible_subjects(bb_df, new_n = words_per_region)
  
  # Randomly selecting the region and filtering
  prs <- get_order_range(new_bb_df)
  # Not really error checked
  stopifnot(1+exclude_from_front <= length(prs)-exclude_from_back)
  possible_region_starts <- prs[(1+exclude_from_front):(length(prs)-exclude_from_back)]
  
  warning("Possible item regions: ", length(possible_region_starts))
  region_start <- sample(possible_region_starts, size = 1)
  new_bb_df <- filter(new_bb_df, Words.Start == region_start)
  
  same_items_per_subj(bb_df = new_bb_df, 
                      items_per_subj = regions_per_subj, 
                      n_subj = n_subj, 
                      seed_val = seed_val,
                      words_per_region = words_per_region)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# like `same_items_per_subj_for_residualized_regions` but samples items then regions
same_items_per_subj_for_residualized_regions_items_first <- function(bb_df, critical_regions_per_subj, filler_items_per_subj, n_subj, seed_val, words_per_region, ..., exclude_items_with_fewer_than = NULL) {

  region_sim_df <- items_then_region_same_items_per_subj(
    bb_df = bb_df, 
    regions_per_subj = critical_regions_per_subj, 
    n_subj = n_subj, 
    seed_val = seed_val, 
    words_per_region = words_per_region)
  
  # HACK
  chosen_subjects <- distinct(region_sim_df, Subject, UniqueSubject)$Subject
  
  filler_items_df <-  sample_from_barebones(
    df = bb_df, 
    size = filler_items_per_subj, 
    n = 1,
    starting_groupval = critical_regions_per_subj+1) 
  
  filler_sim_df <- purrr::map_dfr(
    1:n_subj, 
    ~mutate(filler_items_df, Subject = chosen_subjects[[.]],
            UniqueSubject = paste0(Subject, "_", .),
            SimCond = 0)) %>%
    mutate(ItemAmbValue = 0, SubjectAmbValue = 0) %>%
    select(-PossibleSubjects)
  
  stopifnot(setequal(region_sim_df$UniqueSubject, filler_sim_df$UniqueSubject))
  
  bind_rows("CriticalRegion" = region_sim_df,
            "FillerItems" = filler_sim_df, 
            .id = "ItemType")
}

# like `same_items_per_subj_for_residualized_regions` ~~~but makes sure critical regions don't show up in the fillers~~~ makes sure fillers come from the same region
same_items_per_subj_for_residualized_regions_clean_fillers <- function(bb_df, critical_regions_per_subj, filler_items_per_subj, n_subj, seed_val, words_per_region, ..., exclude_items_with_fewer_than = NULL) {
  
  bb_df <- filter(bb_df, maxWords >= exclude_items_with_fewer_than)
  
  region_sim_df <- region_then_items_same_items_per_subj(
    bb_df = bb_df, 
    regions_per_subj = critical_regions_per_subj, 
    n_subj = n_subj, 
    seed_val = seed_val, 
    words_per_region = words_per_region, 
    exclude_items_with_fewer_than = exclude_items_with_fewer_than)
  
  critical_items <- region_sim_df %>%
    tidyr::unnest(Word.Order) %>% 
    distinct(Item, Word.Order) %>% 
    mutate(PseudoItem = paste(Item, Word.Order))
  
  # # Remove critical items
  # bb_df <- bb_df %>% filter(!(paste(Item, Words.Start) %in% critical_items$PseudoItem))
  
  
  # HACK
  chosen_subjects <- distinct(region_sim_df, Subject, UniqueSubject)$Subject
  
  filler_items_df <-  sample_from_barebones(
    df = bb_df, 
    size = filler_items_per_subj, 
    n = 1,
    starting_groupval = critical_regions_per_subj+1) 
  
  filler_sim_df <- purrr::map_dfr(
    1:n_subj, 
    ~mutate(filler_items_df, Subject = chosen_subjects[[.]],
            UniqueSubject = paste0(Subject, "_", .),
            SimCond = 0)) %>%
    mutate(ItemAmbValue = 0, SubjectAmbValue = 0) %>%
    select(-PossibleSubjects)
  
  stopifnot(setequal(region_sim_df$UniqueSubject, filler_sim_df$UniqueSubject))
  
  bind_rows("CriticalRegion" = region_sim_df,
            "FillerItems" = filler_sim_df, 
            .id = "ItemType")
}


same_items_per_subj_for_residualized_regions_free_order <- function(bb_df, critical_regions_per_subj, filler_items_per_subj, n_subj, seed_val, words_per_region, ...) {
  
  new_bb_df <- readjust_possible_subjects(bb_df, new_n = words_per_region)
  
  region_sim_df <- same_items_per_subj(bb_df = new_bb_df, 
                      items_per_subj = critical_regions_per_subj, 
                      n_subj = n_subj, 
                      seed_val = seed_val,
                      words_per_region = words_per_region)
  
  critical_items <- region_sim_df %>%
    tidyr::unnest(Word.Order) %>% 
    distinct(Item, Word.Order) %>% 
    mutate(PseudoItem = paste(Item, Word.Order))
  
  # HACK
  chosen_subjects <- distinct(region_sim_df, Subject, UniqueSubject)$Subject
  
  filler_items_df <-  sample_from_barebones(
    df = bb_df, 
    size = filler_items_per_subj, 
    n = 1,
    starting_groupval = critical_regions_per_subj+1) 
  
  filler_sim_df <- purrr::map_dfr(
    1:n_subj, 
    ~mutate(filler_items_df, Subject = chosen_subjects[[.]],
            UniqueSubject = paste0(Subject, "_", .),
            SimCond = 0)) %>%
    mutate(ItemAmbValue = 0, SubjectAmbValue = 0) %>%
    select(-PossibleSubjects)
  
  stopifnot(setequal(region_sim_df$UniqueSubject, filler_sim_df$UniqueSubject))
  
  bind_rows("CriticalRegion" = region_sim_df,
            "FillerItems" = filler_sim_df, 
            .id = "ItemType")
}





# Get specific stuff
# Divided structure is asummed to be: `c(<Name> = n, <Name2> = n2)``
same_items_per_subj_divided <- function(bb_df, divided_struct, n_subj, seed_val, ...) {
  # set.seed(seed_val)
  
  # All names must match
  stopifnot(rlang::is_named(divided_struct))
  stopifnot(sort(names(divided_struct)) == sort(unique(bb_df$Name)))
  
  purrr::imap_dfr(
    divided_struct,
    function(items_per_subj, q_name) {
      short_bb_df <- filter(bb_df, Name == q_name)
      
      items_df <-  sample_from_barebones(short_bb_df, size = items_per_subj, n = 1)
      subjects <- unique(purrr::reduce(
        short_bb_df$PossibleSubjects,
        function(x, y) c(x, as.character(y)), .init=NULL)) %>%
        remove_nas() %>%
        sample(size = n_subj, replace = TRUE)
      
      purrr::map_dfr(
        1:n_subj, 
        function(x) {
          mutate(items_df, Subject = subjects[[x]],
                 UniqueSubject = paste0(Subject, "_", x),
                 SimCond = ifelse(x <= n_subj/2, -1, 1)) 
        }) %>%
        mutate(ItemAmbValue = 0,
               SubjectAmbValue = 0) %>%
        select(-PossibleSubjects) %>%
        mutate(Name = q_name)
    }
  )
}


# Assumes n_trials will be divided equally between blocks
# Sim effect is constant within pseudo-subjects
# Not guaranteed to have full size (ie if subject only has one response in only one block, still might get chosen for all blocks)
# DOES NOT ALLOW `n` to vary!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
same_items_per_subj_interaction <- function(bb_df, n_trials, n_subj, seed_val, ...) {
  n_trials <- n_trials / n_distinct(bb_df$Name)
  # set.seed(seed_val)
  
  # Sample the sub-bb-dfs for items
  bb_dfs <- unique(bb_df$Name) %>% 
    purrr::map(~sample_from_barebones(filter(bb_df, Name == .x), size=n_trials, n=1))
  
  # Sample subjects from all subjects who have ANY trials in EITHER bb_df
  subjects <- bb_dfs %>% 
    bind_rows() %>%
    purrr::pluck("PossibleSubjects") %>%
    purrr::reduce(~c(.x, as.character(.y)), 
                  .init=NULL) %>%
    unique() %>%
    remove_nas() %>% # NEW
    sample(size = n_subj, replace = TRUE)
  
  purrr::map_dfr(
    bb_dfs,
    function(items_df) {
      purrr::map_dfr(
        1:n_subj,
        function(i) {
          mutate(items_df, Subject = subjects[[i]],
                 UniqueSubject = paste0(Subject, "_", i),
                 SimCond = ifelse(i <= n_subj/2, -1, 1))
        }) %>%
        mutate(ItemAmbValue = 0,
               SubjectAmbValue = 0,
               BlockName = Name) %>%
        select(-PossibleSubjects)
    })
}


# Just totally random, baby!
setsize_random_sample <- function(bb_df, unique_subjects, # unused
                                  n_trials, seed_val) {
  # set.seed(seed_val)
  
  sim_df <-  sample_from_barebones(bb_df, size = n_trials, n = 1) %>%
    mutate(Subject = purrr::map_chr(PossibleSubjects, ~sample(remove_nas(.), size = 1)),
           UniqueSubject = paste0(Subject, "_", seq_along(Subject)),
           ItemAmbValue = 0,
           SubjectAmbValue = 0) %>%
    mutate(SimCond = ifelse(seq_along(GroupingVal) <= length(GroupingVal)/2, -1, 1)) %>%
    select(-PossibleSubjects)
  
  return(sim_df)
}

# Totally random, except that it focuses on a single word order across ALL trials
# Needs to have `n` set previously
setsize_single_order_sample <- function(bb_df, unique_subjects, n_trials, seed_val) {
  # set.seed(seed_val)
  order_num <- sample(get_order_range(bb_df), size=1)
  bb_df <- bb_df %>% filter(Words.Start == order_num)
  sim_df <- setsize_random_sample(bb_df, unique_subjects, n_trials, seed_val)
  return(sim_df)
}

# Simple random interaction. Assumes n_trials will be divided equally between blocks
#   Also, if participants in the Filler-/RC-group have different overall
#     speeds on, say, Block 3, then this interaction might be biased, in the sense
#     that only Filler-first subjects have responses for Block 1
setsize_interaction_random_sample <- function(bb_df, unique_subjects, n_trials, seed_val) {
  n_trials <- n_trials / n_distinct(bb_df$Name)
  # set.seed(seed_val)
  
  sim_df <- purrr::map_dfr(
    unique(bb_df$Name),
    function(x) {
      sample_from_barebones(filter(bb_df, Name == x), size=n_trials, n=1) %>%
        mutate(Subject = purrr::map_chr(PossibleSubjects, ~sample(remove_nas(.), size = 1)),
               UniqueSubject = paste0(Subject, "_", seq_along(Subject)),
               ItemAmbValue = 0,
               SubjectAmbValue = 0) %>%
        mutate(SimCond = ifelse(seq_along(GroupingVal) <= length(GroupingVal)/2, -1, 1)) %>%
        mutate(Name = x, BlockName = x) %>%
        select(-PossibleSubjects)
    })
  
  return(sim_df)
}


# First samples stories,
# Then sample subjects who have responses to at least half of those stories
# Then sample an equal number of items for each story
sample_NSC_bata <- function(bb_df, n_stories, n_subj, 
                            n_items_per_story, words_per_region,
                            stories_to_use=NULL,
                            seed = NULL) {
  if (n_items_per_story < 2)
    stop("You need to have at least two items per story (because there are two conditions)")
  
  if (!is.null(seed))
    set.seed(seed)
  
  # If a subject has at least responses for two stories, they're cool
  if (n_stories == 1)
    subj_qual_criteria = quo( n_distinct(UniqueStory) >= 1 )
  else
    subj_qual_criteria = quo( n_distinct(UniqueStory) >= n_stories/2 & n_distinct(UniqueStory) >= 2 )
  
  if (is.null(stories_to_use)) {
    # First sample stories
    available_stories <- unique(bb_df$Story)
    stories_to_use <- sample(available_stories, size = n_stories, replace=TRUE)
  }
  
  # Then get all the subjects who have at least one response to one of the stories
  subjs_per_story <- bb_df %>%
    filter(Story %in% stories_to_use) %>%
    select(Story, PossibleSubjects)
  # This is annoying, because of the Stories vs. pseudo-Stories distinction
  subjs_per_story <- purrr::imap_dfr(
    stories_to_use,
    ~filter(subjs_per_story, Story==.x) %>%
      mutate(UniqueStory = paste0(Story, "_", .y)))
  
  available_subj <- subjs_per_story %>%
    rowwise() %>%
    mutate(Subject = paste0(PossibleSubjects, collapse=",")) %>%
    select(-PossibleSubjects) %>%
    ungroup() %>%
    group_by(UniqueStory) %>%
    tidyr::separate_rows(Subject, sep=",") %>%
    group_by(Subject) %>%
    # Then filter them again by the below criteria
    filter(!!subj_qual_criteria) %>%
    {unique(.$Subject)}
  subjects_to_use <- sample(available_subj, size = n_subj, replace=TRUE)
  
  # Just get rid of PossibleSubjects since it's a little dangerous
  bb_df <- select(bb_df, -PossibleSubjects)
  # Then, sample the same number of items per study
  items_df <-  purrr::imap_dfr(
    stories_to_use,
    ~sample_from_barebones(filter(bb_df, Story==.x),
                           size = n_items_per_story, n = words_per_region) %>%
      mutate(SimCond = ifelse(seq_along(UniqueItem) <= n_items_per_story/2, -1, 1),
             Story = .x,
             UniqueStory = paste0(Story, "_", .y)))
  
  purrr::map_dfr(
    1:n_subj,
    ~mutate(items_df, Subject = subjects_to_use[[.]],
            UniqueSubject = paste0(Subject, "_", .x)))
}

sample_regionized_NSC_bata <- function(bb_df, critical_regions_per_subj_per_story, n_stories,
                                       filler_items_per_story, n_subj, words_per_region = 3) {
  
  # ---------------------------- this is the "region_sim_df" function analogue
  # Just some sanity checking
  stopifnot("maxWords" %in% names(bb_df))
  one_max_word_per_item <- bb_df %>% group_by(Item) %>%
    summarise(n = n_distinct(maxWords)) %>% {all(.$n == 1)}
  stopifnot(one_max_word_per_item)
  
  new_bb_df <- bb_df %>%
    # Exclude sentences with less than the number of words per_region
    filter(maxWords >= words_per_region) %>%
    mutate(Name = "") %>% # To satisfy the next function
    # Readjusting `PossibleSubjects` from `n==1` to `n==words_per_region`
    readjust_possible_subjects(new_n = words_per_region)
  
  # Unlike the SETAl/FETAL data, we don't keep the place of the region the same
  region_sim_df <- sample_NSC_bata(new_bb_df,
                                   n_stories = n_stories,
                                   n_subj = n_subj,
                                   n_items_per_story = critical_regions_per_subj_per_story,
                                   words_per_region = words_per_region)
  
  # Kind of hacky, but agnostic to how the unique values are specified
  stories_used <-  distinct(region_sim_df, Story,   UniqueStory  )$Story
  subjects_used <- distinct(region_sim_df, Subject, UniqueSubject) # used for something slightly different
  # Could use the following instead:
  # unique(region_sim_df$UniqueStory) %>%
  #   purrr::map_chr(~strsplit(., "_")[[1]][[1]])
  
  # Also hacky, but eh
  wrong_subject_filler_df <- sample_NSC_bata(bb_df,
                                             n_stories = n_stories,
                                             n_subj = n_subj,
                                             n_items_per_story = filler_items_per_story,
                                             words_per_region = 1,
                                             stories_to_use = stories_used) %>%
    mutate(SimCond = 0,
           GroupingVal = GroupingVal + max(region_sim_df$GroupingVal),
           UniqueItem  = paste0(Item, "_", GroupingVal))
  
  wrong_subjects <- distinct(wrong_subject_filler_df, Subject, UniqueSubject)$UniqueSubject
  
  filler_sim_df <- wrong_subject_filler_df %>%
    rowwise() %>%
    mutate(index = which(UniqueSubject == wrong_subjects)) %>%
    ungroup() %>%
    mutate(Subject = factor(subjects_used[index,]$Subject),
           UniqueSubject = factor(subjects_used[index,]$UniqueSubject))
  
  stopifnot(setequal(region_sim_df$UniqueSubject, filler_sim_df$UniqueSubject))
  stopifnot(setequal(region_sim_df$UniqueStory,   filler_sim_df$UniqueStory))
  
  suppressWarnings(
    bind_rows("CriticalRegion" = region_sim_df,
              "FillerItems" = filler_sim_df,
              .id = "ItemType")
  )
}






# # ---------------------------------------------------------------------------------
# # These are the old, OUTDATED AMLaP functions
# #   The problem with them is that they sampled item-word pairs independently of subjects.
# #   So when a particular subject didn't have a row for that item-word pair (as in RC-first group), when we add the data back in, we get lots of missing data
# # ---------------------------------------------------------------------------------
# # Totally random, except that it focuses on a single word order across trials
# simple_single_order_sample <- function(bb_df, unique_subjects, n_trials, seed_val) {
#   # set.seed(seed_val)
#   order_num <- sample(get_order_range(bb_df), size=1)
#   bb_df <- bb_df %>% filter(Words.Start == order_num)
#   sim_df <- simple_random_sample(bb_df, unique_subjects, n_trials, seed_val)
#   return(sim_df)
# }
# # Just totally random, baby!
# simple_random_sample <- function(bb_df, unique_subjects, n_trials, seed_val) {
#   # set.seed(seed_val)
#   
#   sim_df <- sample_from_barebones(bb_df, size=n_trials, n=1) %>%
#     mutate(UniqueSubject = sample_subjects_from_list(unique_subjects,
#                                                      size=n_trials),
#            Subject = purrr::map_chr(strsplit(UniqueSubject,  "_"), ~.[1]),
#            ItemAmbValue = 0,
#            SubjectAmbValue = 0) %>%
#     mutate(SimCond = ifelse(seq_along(GroupingVal) <= length(GroupingVal)/2, -1, 1))
#   
#   return(sim_df)
# }
# # Simple random interaction. Assumes n_trials will be divided equally between blocks
# simple_interaction_random_sample <- function(bb_df, unique_subjects, n_trials, seed_val) {
#   n_trials <- n_trials / n_distinct(bb_df$Name)
#   # set.seed(seed_val)
#   
#   sim_df <- purrr::map_dfr(
#     unique(bb_df$Name),
#     function(x) {
#       sample_from_barebones(filter(bb_df, Name == x), size=n_trials, n=1) %>%
#         mutate(UniqueSubject = sample_subjects_from_list(unique_subjects,
#                                                          size = n_trials),
#                Subject = purrr::map_chr(strsplit(UniqueSubject,  "_"), ~.[1]),
#                ItemAmbValue = 0,
#                SubjectAmbValue = 0) %>%
#         mutate(SimCond = ifelse(seq_along(GroupingVal) <= length(GroupingVal)/2, -1, 1)) %>%
#         mutate(Name = x, BlockName = x)
#     })
#   
#   return(sim_df)
# }
# 
# # Unused, as of now:
# # This samples for a 'single' bb_df
# # Importantly, this takes only the unique subjects
# # In the style of Q1/Q3
# structured_simple_single_sample <- function(bb_df, unique_subjects, trials_per_subject, subjects_per_group, seed_val) {
#   # set.seed(seed_val)
#   subjects <- tibble(UniqueSubject = sample_subjects_from_list(unique_subjects, size=subjects_per_group*2)) %>%
#     mutate(Subject = purrr::map_chr(strsplit(UniqueSubject,  "_"), ~.[1]))
#   
#   # Sample items for subjects, such that the first half of the ITEMS for each subject will have -1
#   items <- sample_from_barebones(bb_df, trials_per_subject, n=1) %>%
#     mutate(ItemAmbValue = ifelse(seq_along(GroupingVal) <= length(GroupingVal)/2, -1, 1))
#   
#   # Sample subjects for groups
#   group1_single  <- subjects[1:subjects_per_group,] %>%
#     # Make the first half of subjects each GROUP have a -1 AmbVal and the second half +1
#     mutate(SubjectAmbValue = ifelse(seq_along(UniqueSubject) <= subjects_per_group/2, -1, 1))
#   group2_single  <- subjects[(subjects_per_group+1):(subjects_per_group*2),] %>%
#     mutate(SubjectAmbValue = ifelse(seq_along(UniqueSubject) <= subjects_per_group/2, -1, 1))
#   
#   # Repeat the subject dataframe so there are the right number of trials per subject
#   group1_expanded <- group1_single[rep(seq(subjects_per_group), each = trials_per_subject),]
#   group2_expanded <- group2_single[rep(seq(subjects_per_group), each = trials_per_subject),]
#   
#   first_group <-  do.call("rbind", rep(list(items), subjects_per_group)) %>%
#     bind_cols(group1_expanded)
#   second_group <- do.call("rbind", rep(list(items), subjects_per_group)) %>%
#     bind_cols(group2_expanded)
#   
#   both_groups <- bind_rows(mutate(first_group,  SimCond = 1),
#                            mutate(second_group, SimCond = -1))
#   
#   both_groups$GroupingVal <- seq_along(both_groups$GroupingVal)
#   return(both_groups)
# }
