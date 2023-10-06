library(dplyr)
library(lme4)

# Combines a sampled barebones df with actual data
add_real_data <- function(bb_df, real_df) {
  bb_df %>%
    tidyr::unnest(Word.Order) %>%
    left_join(real_df, by = c("Item", "Subject", "Word.Order")) #%>%
    # filter(!is.na(RT)) %>% # leave this decision for others
}

# To make `add_in_effects` more flexible, i move the sampling to this
# TAKE IN STANDARD DEVIATION, NOT VARIANCE
sample_uncertainty <- function(stdev, seed_val) {
  set.seed(as.numeric(seed_val))
  rnorm(n=1, mean=0, sd = stdev)
}

# adding in effects
# Assumes the first value in the effect_size vector corresponds to when SimCond > 0
add_in_effects <- function(df,
                           effect_size,
                           sampled_error) {
  mean_effect_size = mean(effect_size[1:2])
  
  df %>%
    # Prepare the... "multipliers"?
    mutate(Ambiguity = case_when(SubjectAmbValue*ItemAmbValue ==  1 ~ "Ambiguous",
                                 SubjectAmbValue*ItemAmbValue == -1 ~ "Unambiguous",
                                 TRUE ~ NA_character_) %>%
             factor(levels=c("Ambiguous", "Unambiguous")) %>%
             `contrasts<-`(value = c(0.5, -0.5)),
           # If you change this to +1/-1, change the NullEffect formula
           AmbiguityCode = case_when(Ambiguity=="Ambiguous"   ~  0.5,
                                     Ambiguity=="Unambiguous" ~ -0.5,
                                     TRUE ~ NA_real_),
           SimCond = SimCond/2) %>%
    # generate the effects
    mutate(Uncertainty = sampled_error * AmbiguityCode * SimCond,
           NullEffect = mean_effect_size * AmbiguityCode, # Ambiguity code takes care of the other halving
           RealEffect = case_when(
             SimCond > 0 ~ effect_size[1] * AmbiguityCode, # Ambiguity code takes care of the other halving
             SimCond < 0 ~ effect_size[2] * AmbiguityCode, # Ambiguity code takes care of the other halving
             TRUE ~ NA_real_)
           )
}

# Gets the data from the three-consecutive-word data and averages across it
# Adds the effect in "ourWay" space, which currently is log10
avg_out_real_data <- function(df,
                              our_way_bare = predictedLogRT,
                              their_way_bare = predictedRegRT,
                              exclude_row_if_any_NAs) {
  our_way_bare   <- rlang::enquo(our_way_bare)
  their_way_bare <- rlang::enquo(their_way_bare)

  df %>%
  { if (exclude_row_if_any_NAs==TRUE)
     filter(., !is.na(RT),
            !is.na(!! our_way_bare),
             !is.na(!! their_way_bare))
    else
      .
  } %>%
    # Don't just use grouping val!
    group_by(GroupingVal,
             Subject, UniqueSubject, # Subject vals
             Item, UniqueItem,       # Item vals
             Ambiguity, SimCond,     # Coding vals
             Block, Words.Start) %>% # Etc
    # Add effects
    mutate(ourValue = log10(RT) - !!our_way_bare,
           ourWayNoEffect = ourValue + NullEffect + Uncertainty,
           ourWayEffect =   ourValue + RealEffect + Uncertainty,
           theirWayEffect =     10^(ourWayEffect + !!our_way_bare) - !!their_way_bare,
           theirWayNoEffect = 10^(ourWayNoEffect + !!our_way_bare) - !!their_way_bare) %>%
    summarise_at(
      vars(ourWayNoEffect, ourWayEffect, theirWayEffect, theirWayNoEffect, RT),
      funs("list"=paste(., collapse=","), "mean" = mean(., na.rm=TRUE))
    ) %>%
    ungroup() %>%
    purrr::set_names(gsub("_mean", "", names(.)))
}


get_model_stats <- function(el, el_name, real_df, list_of_model_funcs,
                            .require_effect = TRUE,
                            .extra_info = TRUE,
                            inside_tidy_function = broom::tidy) {
  # hack
  if (.extra_info==TRUE) {
    cleaner_f <- function(m) {
      y_term <- terms(formula(m))[[2]]
      if (mode(y_term) == "name")
        y_name <- as.character(y_term)
      else if (mode(y_term) == "call")
        y_name <- as.character(format(y_term))
      else stop(paste0("ya fucked ", y_term))
      m %>%
      {cbind(inside_tidy_function(.),
             "stdevY"=sd(model.frame(.)[, y_name]),
             "meanY"=mean(model.frame(.)[, y_name]),
             "n_subj"=length(unique(model.frame(.)[,"UniqueSubject"])),
             "n_rows"=nrow(model.frame(.)))}}
  }
  else
    cleaner_f <- inside_tidy_function
    
  possible_effect_size <- attributes(el)[["effect_size"]]
  if (.require_effect==TRUE)
    stopifnot(!is.null(possible_effect_size))
  
  df <- add_real_data(el, real_df) %>%
    add_in_effects(.,
                   effect_size = possible_effect_size,
                   sampled_error = sample_uncertainty(possible_effect_size[3], seed_val = el_name)) %>%
    avg_out_real_data(.,
                      our_way_bare = predictedLogRT,
                      their_way_bare = predictedRegRT,
                      exclude_row_if_any_NAs = FALSE)
  future_map(list_of_model_funcs,
             ~zplyr::collect_all(cleaner_f(.(df)),
                                 catchErrors=TRUE))
}

make_models <- function(ifile, ofile, 
                        real_df_file, list_of_model_funcs,
                        inside_tidy_function = broom::tidy) {
  l <- readRDS(ifile)
  real_df <- readRDS(real_df_file)
  model_list <- future_imap(l, ~get_model_stats(.x, .y, real_df, list_of_model_funcs, inside_tidy_function=inside_tidy_function))
  #mark_progress("E")
  saveRDS(model_list, ofile)
  if (!file.exists(ofile)) {
    Sys.sleep(10)
    saveRDS(model_list, ofile)
    if (!file.exists(ofile)) {
      warning(head(model_list))
      warning(ofile)
      
      stop(paste0("length=", length(model_list)))
    }
  }
  #mark_progress("F")
  "done!"
}

df_list_to_df <- function(df,
                          terms_to_harvest = c("SimCond"),
                          .multicore = FALSE) {
  # if (.multicore == FALSE) map_func <- purrr::imap_dfr
  # else 
  map_func <- furrr::future_imap_dfr
  
  df %>%
    map_func(
      function(dfl, iter) {
        purrr::imap_dfr(
          dfl,
          function(single, name_type) {
            if (is.na(single$value)) {
              new_df <- data.frame(estimate=NA)
            } else {
              new_df <- single$value %>%
                filter(term %in% terms_to_harvest)
            }
            new_df %>%
              mutate(messages = paste0(single$messages, collapse="--------"),
                     warnings = paste0(single$warnings, collapse="--------"),
                     errors = paste0(single$errors, collapse="--------"),
                     type = name_type)
          }) %>%
          mutate(k = as.numeric(iter))
      })
}



collate_data <- function(ifile,
                         terms_to_harvest = c("SimCond"),
                         .multicore = FALSE) {
  readRDS(ifile) %>% 
    df_list_to_df(terms_to_harvest = terms_to_harvest, .multicore =.multicore)
}

# For simplified running of commands in the format I've been doing.
# Note the real_df_function
run_command <- function(ititle="bb", otitle="modeldata", qs, exps, ks, model_funcs,
                        real_df_function = ~ifelse(grepl("[Ff]etal", .),"fetal_sim_df.RDS","setal_sim_df.RDS")) {
  message("Make sure the sim_df files correctly correspond to the experiment names")
  real_df_function <- purrr::as_mapper(real_df_function)
  bleh %<-% {
    walk(
      qs, function(q_i) {
        walk(
          exps, function(expname) {
            walk(
              ks, function(k_e) {
                min_k <- min(k_e)
                max_k <- max(k_e)
                make_models(
                  ifile = paste(
                    paste0(saver_path, ititle), expname, paste0("q", q_i), min_k, "to", paste0(max_k, ".RDS"),
                    sep="_"),
                  ofile = paste(
                    paste0(saver_path, otitle), expname, paste0("q", q_i), min_k, "to", paste0(max_k, ".RDS"),
                    sep="_"),
                  real_df_file = paste0(path_on_server, real_df_function(expname)),
                  list_of_model_funcs = model_funcs
                )})})})
  }
  return(futureOf(bleh))
}
