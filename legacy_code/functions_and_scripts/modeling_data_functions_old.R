library(dplyr)
library(lme4)

# Combines a sampled barebones df with actual data
add_real_data <- function(bb_df, real_df) {
  effect_size = attributes(bb_df)[["effect_size"]]
  if (is.null(effect_size))
    stop("Barebones dataframe has no attribute 'effect_size'")

  bb_df %>%
    tidyr::unnest(Word.Order) %>%
    left_join(real_df, by = c("Item", "Subject", "Word.Order")) %>%
    filter(!is.na(RT)) %>%
    `attr<-`("effect_size", effect_size)
}

# Gets the data from the three-consecutive-word data and averages across it
avg_out_real_data <- function(df,
                              effect_size,
                              exclude_row_if_any_NAs,
                              our_way_bare = predictedLogRT,
                              their_way_bare = predictedRegRT) {
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
             Subject, UniqueSubject,
             Item, UniqueItem,
             SimCond, Block, Words.Start) %>%
    mutate(ourWayNoEffect = log10(RT) - !!our_way_bare,
           ourWayEffect =   ourWayNoEffect + (SimCond * effect_size/2),
           theirWayEffect = 10^(ourWayEffect + !!our_way_bare) - !!their_way_bare,
           theirWayNoEffect = RT - !!their_way_bare) %>%
    summarise_at(
      vars(ourWayNoEffect, ourWayEffect, theirWayEffect, theirWayNoEffect, RT),
      funs("list"=paste(., collapse=","), "mean" = mean(., na.rm=TRUE))
    ) %>%
    ungroup() %>%
    purrr::set_names(gsub("_mean", "", names(.)))
}


get_model_stats <- function(el, real_df, list_of_model_funcs) {
  df <- add_real_data(el, real_df) %>%
    avg_out_real_data(.,
                      effect_size = attributes(.)[["effect_size"]],
                      our_way_bare = predictedLogRT,
                      their_way_bare = predictedRegRT,
                      exclude_row_if_any_NAs = TRUE)
  purrr::map(list_of_model_funcs,
             ~zplyr::collect_all(broom::tidy(.(df)),
                                 catchErrors=TRUE))
}

make_models <- function(ifile, ofile, real_df_file, list_of_model_funcs) {
  l <- readRDS(ifile)
  real_df <- readRDS(real_df_file)
  model_list <- future_map(l, ~get_model_stats(., real_df, list_of_model_funcs))
  saveRDS(model_list, ofile)
  "done!"
}
collate_data <- function(ifile,
                         term_to_harvest="SimCond") {
  df <- readRDS(ifile) %>%
    purrr::imap_dfr(
      function(dfl, iter) {
        purrr::imap_dfr(dfl,
                        function(single, name_type) {
                          single$value %>%
                            filter(term == term_to_harvest) %>%
                            mutate(messages = paste0(single$messages, collapse="--------"),
                                   warnings = paste0(single$warnings, collapse="--------"),
                                   errors = paste0(single$errors, collapse="--------"),
                                   type = name_type)}) %>%
          mutate(k = as.numeric(iter))
      })
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
