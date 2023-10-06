convergence_handler <- function(exclude=TRUE, no_lmerTest_warnings=FALSE, return_quosures=FALSE,
                                df = NA) {
  filterers <- rlang::exprs(
    errors != "",
    grepl("unable to evaluate scaled", warnings, ignore.case = TRUE),
    grepl("rank deficient|unable to evaluate scaled", messages, ignore.case = TRUE))

  if (no_lmerTest_warnings==TRUE)
    filterers <- rlang::exprs(!!!filterers, grepl("converg", warnings, ignore.case = TRUE))
  else
    filterers <- rlang::exprs(!!!filterers,
                              (grepl("converg", warnings, ignore.case = TRUE) & !grepl("as_lmerModLT", warnings)),
                              (grepl("Model failed to converge", warnings) & grepl("as_lmerModLT", warnings)))

  if (exclude==TRUE)
    filterers <- filterers %>% map(function(q) rlang::expr(!(!!q)))
  else
    filterers <- filterers %>% reduce(function(q1, q2) rlang::expr(!!q1 | !!q2))

  if (return_quosures==TRUE)
    return(filterers)
  else
    return(filter(df, !!!filterers))
}

cond_df_to_df <- function(df_w_conds, terms_to_harvest="SimCond") {
  df_w_conds$value %>%
    filter(term %in% terms_to_harvest) %>%
    mutate(messages = paste0(df_w_conds$messages, collapse="--------"),
           warnings = paste0(df_w_conds$warnings, collapse="--------"),
           errors = paste0(df_w_conds$errors, collapse="--------"))
}

delete_this <- function(cond_carrier) {
  tibble() %>%
    mutate(messages = paste0(cond_carrier$messages, collapse="--------"),
           warnings = paste0(cond_carrier$warnings, collapse="--------"),
           errors = paste0(cond_carrier$errors, collapse="--------"))
}


model_converged <- function(df_w_conds, custom_converge_test=NA) {
  df <- tibble(
    messages = paste0(df_w_conds$messages, collapse="--------"),
    warnings = paste0(df_w_conds$warnings, collapse="--------"),
    errors   = paste0(df_w_conds$errors, collapse="--------"))

  if (!rlang::is_na(custom_converge_test))
    return(nrow(custom_converge_test(df)) > 0)
  else
    return(nrow(convergence_handler(df=df)) > 0)
}

# PRE-RESIDUALIZATION FALLBACK MODELING FUNCTIONS -------------------

# requires top level of `fallback_model_list` to be named
fall_back_modeler <- function(df, model_type, cleaner_list, fallback_model_list) {
  stopifnot(length(cleaner_list)==length(fallback_model_list))

  fall_back_id <- names(fallback_model_list[1])[[1]]
  model_funcs <- fallback_model_list[[1]]
  model_cleaner <- cleaner_list[[1]]
  
  message("Trying ", fall_back_id," models\n")

  model_results <- model_funcs %>%
    purrr::map(~zplyr::collect_all(model_cleaner(.[[model_type]](df)), catchErrors=TRUE))

  kept_results <- model_results %>% keep(model_converged)

  if (is_empty(kept_results)) {
    if (length(fallback_model_list)==1) stop("You've run out of fallback models")
    # Doing some cleaning
    rm(model_results, model_funcs)
    # Poppin'
    cleaner_list        <-        cleaner_list[2:length(cleaner_list)]
    fallback_model_list <- fallback_model_list[2:length(fallback_model_list)]
    # Recurse
    return(fall_back_modeler(df, model_type, cleaner_list, fallback_model_list))
  } else {
    return(list(results=kept_results, id=fall_back_id))
  }
}

# Starts everything out
fallback_checker <- function(df, model_type, cleaner_list, fallback_model_list, terms_to_harvest="SimCond") {
  # Name giver
  fallback_model_list <- rlang::names2(fallback_model_list) %>%
    ifelse(.=="", 1:length(fallback_model_list), .) %>%
    {rlang::set_names(fallback_model_list, .)}

  res <- fall_back_modeler(df, model_type, cleaner_list, fallback_model_list)
  model_level <- res$id
  # stop(length(res$results), "= length, typeof ", typeof(res$results[[1]]))
  # return(res$results)
  
  models_df <- res$results %>%
    map(~cond_df_to_df(., terms_to_harvest=terms_to_harvest)) %>% 
    bind_rows(.id="displayed_model")
  number_converged <- nrow(models_df)
  
  filter(models_df, p.value==min(p.value)) %>% 
    slice(1) %>%
    mutate(converged_id = model_level,
           number_converged = number_converged)
}

pre_resid_fallback_models <- list(
  one_intercept = list(
    no_item_int = list(
      linear_power = function(x) lme4::lmer(RT_with_Effect ~ 1 + SimCond + (1 | UniqueSubject),
                                            data = filter(x, RT_with_Effect > 0),
                                            control = lme4::lmerControl(optimizer="bobyqa")),
      log_power = function(x) lme4::lmer(log10(RT_with_Effect) ~ 1 + SimCond + (1 | UniqueSubject),
                                         data = filter(x, RT_with_Effect > 0),
                                         control = lme4::lmerControl(optimizer="bobyqa")),
      linear_type1 = function(x) lme4::lmer(RT ~ 1 + SimCond + (1 | UniqueSubject),
                                            data = filter(x, RT_with_Effect > 0),
                                            control = lme4::lmerControl(optimizer="bobyqa")),
      log_type1 = function(x) lme4::lmer(log10(RT) ~ 1 + SimCond + (1 | UniqueSubject),
                                         data = filter(x, RT_with_Effect > 0),
                                         control = lme4::lmerControl(optimizer="bobyqa"))
    ),
    no_subj_int = list(
      linear_power = function(x) lme4::lmer(RT_with_Effect ~ 1 + SimCond + (1 | UniqueItem),
                                            data = filter(x, RT_with_Effect > 0),
                                            control = lme4::lmerControl(optimizer="bobyqa")),
      log_power = function(x) lme4::lmer(log10(RT_with_Effect) ~ 1 + SimCond + (1 | UniqueItem),
                                         data = filter(x, RT_with_Effect > 0),
                                         control = lme4::lmerControl(optimizer="bobyqa")),
      linear_type1 = function(x) lme4::lmer(RT ~ 1 + SimCond + (1 | UniqueItem),
                                            data = filter(x, RT_with_Effect > 0),
                                            control = lme4::lmerControl(optimizer="bobyqa")),
      log_type1 = function(x) lme4::lmer(log10(RT) ~ 1 + SimCond + (1 | UniqueItem),
                                         data = filter(x, RT_with_Effect > 0),
                                         control = lme4::lmerControl(optimizer="bobyqa"))
    )),
  no_rand_effects = list(
    no_rand_effect_model = list(
      linear_power = function(x) lm(RT_with_Effect ~ 1 + SimCond, data = filter(x, RT_with_Effect > 0)),
      log_power = function(x) lm(log10(RT_with_Effect) ~ 1 + SimCond, data = filter(x, RT_with_Effect > 0)),
      linear_type1 = function(x) lm(RT ~ 1 + SimCond, data = filter(x, RT_with_Effect > 0)),
      log_type1 = function(x) lm(log10(RT) ~ 1 + SimCond, data = filter(x, RT_with_Effect > 0))
    ))
)
pre_resid_fallback_cleaners <- list(
  mixed_cleaner_fallback = mixed_cleaner,
  lm_cleaner_fallback = function(x) broom::tidy(x)
)



# POST-RESIDUALIZATION FALLBACK MODELING FUNCTIONS -------------------

# Residualization funcs should be a single set per second level in the list
# Also, the model funcs should NOT be a full set, since the predvar is transformed to the same
two_model_fall_back_modeler <- function(df, model_type, cleaner_list, fallback_model_list, fallback_resid_list) {
  stopifnot(length(cleaner_list)==length(fallback_model_list))
  exclude_resid_nonconv <- function(x) filter(x, !grepl("CUSTOM", errors))

  fall_back_id <- names(fallback_model_list[1])[[1]]
  resid_id     <- names(fallback_resid_list[1])[[1]]
  model_funcs <- fallback_model_list[[1]]
  resid_funcs <- fallback_resid_list[[1]]
  model_cleaner <- cleaner_list[[1]]

  model_results <- model_funcs %>%
    purrr::map(function(mf) {
      zplyr::collect_all(
        model_cleaner( mf( resid_funcs[[model_type]](df) )), 
        catchErrors=TRUE) })
  
  message("Trying ", fall_back_id," models\n")
  
  # Gets both types of convergence results
  resid_conv_results  <- model_results %>% keep(~model_converged(., custom_converge_test=exclude_resid_nonconv))
  convergence_results <- model_results %>% keep(~model_converged(.))
  
  
  # if (fall_back_id=="no_re") {
    # warning(model_results[[1]]$messages)
    # warning(model_results[[1]]$errors)
    # warning(model_results[[1]]$warning)
    # warning("^^^^^^")
  # }
  
  stopifnot(length(resid_conv_results)==0 | length(resid_conv_results) == length(model_funcs))

  # Since the resid model not converging implies general non-convergence, this means the resid conv check didn't work
  if (is_empty(resid_conv_results)) {
    if (length(fallback_resid_list)==1) stop("You've run out of fallback residualizing models")
    # Doing some cleaning
    rm(model_results, model_funcs)
    # Poppin'
    fallback_resid_list <- fallback_resid_list[2:length(fallback_resid_list)]
    # Recurse
    return(two_model_fall_back_modeler(df, model_type, cleaner_list, fallback_model_list, fallback_resid_list))
  } else if (is_empty(convergence_results)) {
    if (length(fallback_model_list)==1) stop("You've run out of fallback models")
    # Doing some cleaning
    rm(model_results, model_funcs)
    # Poppin'
    cleaner_list        <-        cleaner_list[2:length(cleaner_list)]
    fallback_model_list <- fallback_model_list[2:length(fallback_model_list)]
    # Recurse
    return(two_model_fall_back_modeler(df, model_type, cleaner_list, fallback_model_list, fallback_resid_list))
  } else {
    return(list(results=convergence_results, id=fall_back_id, resid_id = resid_id))
  }
}

# Starts everything out
two_model_fallback_checker <- function(df, model_type, cleaner_list, fallback_model_list, fallback_resid_list, terms_to_harvest="SimCond") {
  # Name giver
  fallback_model_list <- rlang::names2(fallback_model_list) %>%
    ifelse(.=="", 1:length(fallback_model_list), .) %>%
    {rlang::set_names(fallback_model_list, .)}
  # Name giver x2
  fallback_resid_list <- rlang::names2(fallback_resid_list) %>%
    ifelse(.=="", 1:length(fallback_resid_list), .) %>%
    {rlang::set_names(fallback_resid_list, .)}
  
  res <- two_model_fall_back_modeler(df, model_type, cleaner_list, fallback_model_list, fallback_resid_list)
  model_level <- res$id
  resid_level <- res$resid_id
  
  models_df <- res$results %>%
    map(~cond_df_to_df(., terms_to_harvest=terms_to_harvest)) %>% 
    bind_rows(.id="displayed_model")
  number_converged <- nrow(models_df)
  
  filter(models_df, p.value==min(p.value)) %>% 
    slice(1) %>%
    mutate(converged_id = model_level, resid_id = resid_level,
           number_converged = number_converged)
}

# The models
post_resid_fallback_models <- list(
  with_items = list(
    with_item_int = function(x) lme4::lmer(PredVar_Resid ~ 1 + SimCond + (1 | UniqueItem),
                                           data = filter(x, RT_with_Effect > 0),
                                           control = lme4::lmerControl(optimizer="bobyqa"))),
  no_re = list(
    no_re_model = function(x) lm(PredVar_Resid ~ 1 + SimCond,
                                 data = filter(x, RT_with_Effect > 0))))
# The cleaners
post_resid_fallback_cleaners <- list(
  mixed_cleaner_fallback = mixed_cleaner,
  lm_cleaner_fallback = function(x) broom::tidy(x)
)
# Assumes all resids will converge with by-subj REs
fallback_residualizers <- list(
  full_residualization = list(
      linear_power = residualizer_generator(
        RT_with_Effect,    1 + Word.Length + (1 + Word.Length | UniqueSubject)),
      log_power    = residualizer_generator(
        LogRT_with_Effect, 1 + Word.Length + (1 + Word.Length | UniqueSubject)),
      linear_type1 = residualizer_generator(
        RT,                1 + Word.Length + (1 + Word.Length | UniqueSubject)),
      log_type1    = residualizer_generator(
        LogRT,             1 + Word.Length + (1 + Word.Length | UniqueSubject))
  ),
  no_slopes_resid = list(
    linear_power = residualizer_generator(
      RT_with_Effect,    1 + Word.Length + (1 | UniqueSubject)),
    log_power    = residualizer_generator(
      LogRT_with_Effect, 1 + Word.Length + (1 | UniqueSubject)),
    linear_type1 = residualizer_generator(
      RT,                1 + Word.Length + (1 | UniqueSubject)),
    log_type1    = residualizer_generator(
      LogRT,             1 + Word.Length + (1 | UniqueSubject))
  )
)

# Used to generate residualization functions. 
residualizer_generator <- function(col_name, # ie RT, LogRT, RT_with_Effect, LogRT_with_Effect
                                   predictors#, # effect struct of resid model
                                   # has_regions # TRUE for regionized data
) {
  col_name <- rlang::enquo(col_name)
  predictors <- rlang::enquo(predictors)
  # has_regions <- identity(has_regions)
  function(df) {
    
    # This is a super hack to trick `future`
    f <- lmerTestBroom
    
    has_regions <- df %>% distinct(UniqueItem, Word.Order) %>% 
      group_by(UniqueItem) %>% filter(n()>1) %>% {nrow(.) > 0}
    
    if (has_regions==FALSE) {
      # HACK
      number_of_critical_items <- ceiling(n_distinct(df$UniqueItem)/16)
      df <- set_filler_data(df, n_remaining_items = number_of_critical_items)
    } else
      df <- mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0))
    
    df <- df %>%
      mutate(LogRT = log10(RT),
             LogRT_with_Effect = log10(RT_with_Effect)) %>%
      amlap_residualize_bata_one_y_generalized(!!col_name, BATA_PredVar, !!predictors)
    
    if (has_regions==FALSE) 
      df <- remove_filler_data(df, n_remaining_items = number_of_critical_items)
    else
      df <- filter(df, ItemType == "CriticalRegion")
    
    df <- df %>%
      mutate(UniqueItem = factor(UniqueItem)) %>%
      # Calculates residuals BEFORE averaging over the regions
      mutate(PredVar_Resid = !!col_name - BATA_PredVar) %>%
      # Really important to have this right!!!!!!!!!
      group_by(Item, UniqueItem, Subject, UniqueSubject, SimCond, ItemAmbValue, SubjectAmbValue) %>%
      summarise(PredVar_Resid =  mean(PredVar_Resid,  na.rm = TRUE),
                RT_with_Effect = mean(RT_with_Effect, na.rm = TRUE)) %>%
      ungroup()
    
    if (all(is.na(df$PredVar_Resid)) | all(is.nan(df$PredVar_Resid)))
      stop("CUSTOM: Residualization model did not converge")
    
    df
  }
}
  
  

# Temp testing functions ----------------------------------------


all_things <- readRDS("~/Desktop/delete.RDS")

# Fall
setal_giant_pre_resid_fallbacks %beep% {
  res <- iterate_over_df(
    all_things,
    future_cnd_map,
    function(zaa) {
      
      real_df <- readRDS(real_df_file)
      bad_batas <- BadBATAs[[1]] %>%
      {split(., 1:nrow(.))} %>% unname()
      
      bb_list <- readRDS(ifile)[map_chr(bad_batas, ~as.character(.$k))] 
      
      furrr::future_map2( #future_cnd_map2(
        bb_list, bad_batas,
        function(bb_df, bata_info) {
          stopifnot(nrow(bata_info)==1)
          effect_size <- as.numeric(as.character(bata_info$EffectSize))
          adder_f <- function(df) mutate(df, RT_with_Effect = RT + SimCond * effect_size) 
          # Flesh out data frame
          fl_df <- bb_df %>%
            tidyr::unnest(Word.Order) %>%
            { if ("Story" %in% names(real_df)) #ie is a NSC data frame
              left_join(., real_df, by = c("Story", "Item", "Subject", "Word.Order"))
              else 
                left_join(., real_df, by = c("Item", "Subject", "Word.Order"))} %>%
            mutate(SimCond = SimCond/2) %>%
            adder_f()
          
          fallback_checker(fl_df, bata_info$type, 
                           cleaner_list = pre_resid_fallback_cleaners,
                           fallback_model_list = pre_resid_fallback_models,
                           terms_to_harvest = bata_info$term) %>%
            bind_cols(select(bata_info, -term)) %>%
            mutate(ifile, SampleName, Experiment, real_df_file)
        }) %>%
        bind_rows()
    })
  res %>% map_dfr(~.$value) %>% saveRDS(paste0(path_on_server, "amlap_saved_files/results/sfetal_giant_pre_resid_fallbacks.RDS"))
  res %>% map(~.[2:4])
} %plan% elite_squad_plan
done(setal_giant_pre_resid_fallbacks)




# RESIDUALIZED #####################################################################

post_resid_all_things <- readRDS("~/Desktop/delete.RDS") %>% mutate(Experiment="HS18")
single_thing <- post_resid_all_things
single_thing$BadBATAs[[1]] <- single_thing$BadBATAs[[1]][1,]

# Fall
setal_giant_post_resid_fallbacks %beep% {
  res <- iterate_over_df(
    post_resid_all_things,
    future_cnd_map,
    function(zaa) {
      
      real_df <- readRDS(real_df_file)
      bad_batas <- BadBATAs[[1]] %>%
      {split(., 1:nrow(.))} %>% unname()
      
      bb_list <- readRDS(ifile)[map_chr(bad_batas, ~as.character(.$k))] 
      
      furrr::future_map2( #future_cnd_map2(
        bb_list, bad_batas,
        function(bb_df, bata_info) {
          stopifnot(nrow(bata_info)==1)
          effect_size <- as.numeric(as.character(bata_info$EffectSize))
          adder_f <- function(df) mutate(df, RT_with_Effect = RT + SimCond * effect_size) 
          # Flesh out data frame
          fl_df <- bb_df %>%
            tidyr::unnest(Word.Order) %>%
            { if ("Story" %in% names(real_df)) #ie is a NSC data frame
              left_join(., real_df, by = c("Story", "Item", "Subject", "Word.Order"))
              else 
                left_join(., real_df, by = c("Item", "Subject", "Word.Order"))} %>%
            mutate(SimCond = SimCond/2) %>%
            adder_f()
          
          residualizer_generator <- function(col_name, # ie RT, LogRT, RT_with_Effect, LogRT_with_Effect
                                             predictors#, # effect struct of resid model
                                             # has_regions # TRUE for regionized data
          ) {
            col_name <- rlang::enquo(col_name)
            predictors <- rlang::enquo(predictors)
            # has_regions <- identity(has_regions)
            function(df) {
              
              has_regions <- df %>% distinct(UniqueItem, Word.Order) %>% 
                group_by(UniqueItem) %>% filter(n()>1) %>% {nrow(.) > 0}
              
              if (has_regions==FALSE) {
                # HACK
                number_of_critical_items <- ceiling(n_distinct(df$UniqueItem)/16)
                df <- set_filler_data(df, n_remaining_items = number_of_critical_items)
              } else
                df <- mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0))
              
              df <- df %>%
                mutate(LogRT = log10(RT),
                       LogRT_with_effect = log10(RT_with_Effect)) %>%
                amlap_residualize_bata_one_y_generalized(!!col_name, BATA_PredVar, !!predictors)
              
              if (has_regions==FALSE) 
                df <- remove_filler_data(df, n_remaining_items = number_of_critical_items)
              else
                df <- filter(df, ItemType == "CriticalRegion")
              f <- lmerTestBroom
              df <- df %>%
                mutate(UniqueItem = factor(UniqueItem)) %>%
                # Calculates residuals BEFORE averaging over the regions
                mutate(PredVar_Resid = !!col_name - BATA_PredVar) %>%
                # Really important to have this right!!!!!!!!!
                group_by(Item, UniqueItem, Subject, UniqueSubject, SimCond, ItemAmbValue, SubjectAmbValue) %>%
                summarise(PredVar_Resid =  mean(PredVar_Resid,  na.rm = TRUE),
                          RT_with_Effect = mean(RT_with_Effect, na.rm = TRUE)) %>%
                ungroup()
              
              if (all(is.na(df$PredVar_Resid)) | all(is.nan(df$PredVar_Resid)))
                stop("CUSTOM: Residualization model did not converge")
              
              df
            }
          }
          
          two_model_fallback_checker(fl_df, bata_info$type, 
                                     cleaner_list = post_resid_fallback_cleaners,
                                     fallback_model_list = post_resid_fallback_models,
                                     fallback_resid_list = fallback_residualizers) %>%
            bind_cols(select(bata_info, -term)) %>%
            mutate(ifile, SampleName, Experiment, real_df_file)
        }) %>%
        bind_rows()
    })
  res %>% map_dfr(~.$value) %>% saveRDS(paste0(path_on_server, "amlap_saved_files/results/setal_giant_resid_fallbacks.RDS"))
  res %>% map(~.[2:4])
} %plan% elite_squad_plan
done(setal_giant_post_resid_fallbacks)

# REGIONIZED ########################################################################3

post_resid_all_things <- readRDS("~/Desktop/delete.RDS") %>% mutate(Experiment="HS18")
single_thing <- post_resid_all_things

setal_giant_post_resid_fallbacks %beep% {
  res <- iterate_over_df(
    post_resid_all_things,
    future_cnd_map,
    function(zaa) {
      
      real_df <- readRDS(real_df_file)
      bad_batas <- BadBATAs[[1]] %>%
      {split(., 1:nrow(.))} %>% unname()
      
      bb_list <- readRDS(ifile)[map_chr(bad_batas, ~as.character(.$k))] 
      
      furrr::future_map2( #future_cnd_map2(
        bb_list, bad_batas,
        function(bb_df, bata_info) {
          stopifnot(nrow(bata_info)==1)
          effect_size <- as.numeric(as.character(bata_info$EffectSize))
          adder_f <- function(df) mutate(df, RT_with_Effect = RT + SimCond * effect_size) 
          # Flesh out data frame
          fl_df <- bb_df %>%
            tidyr::unnest(Word.Order) %>%
            { if ("Story" %in% names(real_df)) #ie is a NSC data frame
              left_join(., real_df, by = c("Story", "Item", "Subject", "Word.Order"))
              else 
                left_join(., real_df, by = c("Item", "Subject", "Word.Order"))} %>%
            mutate(SimCond = SimCond/2) %>%
            adder_f()
          
          residualizer_generator <- function(col_name, # ie RT, LogRT, RT_with_Effect, LogRT_with_Effect
                                             predictors#, # effect struct of resid model
                                             # has_regions # TRUE for regionized data
          ) {
            col_name <- rlang::enquo(col_name)
            predictors <- rlang::enquo(predictors)
            # has_regions <- identity(has_regions)
            function(df) {
              
              has_regions <- df %>% distinct(UniqueItem, Word.Order) %>% 
                group_by(UniqueItem) %>% filter(n()>1) %>% {nrow(.) > 0}
              
              if (has_regions==FALSE) {
                # HACK
                number_of_critical_items <- ceiling(n_distinct(df$UniqueItem)/16)
                df <- set_filler_data(df, n_remaining_items = number_of_critical_items)
              } else
                df <- mutate(df, SimCond = ifelse(ItemType=="CriticalRegion", SimCond, 0))
              
              df <- df %>%
                mutate(LogRT = log10(RT),
                       LogRT_with_effect = log10(RT_with_Effect)) %>%
                amlap_residualize_bata_one_y_generalized(!!col_name, BATA_PredVar, !!predictors)
              
              if (has_regions==FALSE) 
                df <- remove_filler_data(df, n_remaining_items = number_of_critical_items)
              else
                df <- filter(df, ItemType == "CriticalRegion")
              f <- lmerTestBroom
              df <- df %>%
                mutate(UniqueItem = factor(UniqueItem)) %>%
                # Calculates residuals BEFORE averaging over the regions
                mutate(PredVar_Resid = !!col_name - BATA_PredVar) %>%
                # Really important to have this right!!!!!!!!!
                group_by(Item, UniqueItem, Subject, UniqueSubject, SimCond, ItemAmbValue, SubjectAmbValue) %>%
                summarise(PredVar_Resid =  mean(PredVar_Resid,  na.rm = TRUE),
                          RT_with_Effect = mean(RT_with_Effect, na.rm = TRUE)) %>%
                ungroup()
              
              if (all(is.na(df$PredVar_Resid)) | all(is.nan(df$PredVar_Resid)))
                stop("CUSTOM: Residualization model did not converge")
              
              df
            }
          }
          
          two_model_fallback_checker(fl_df, bata_info$type, 
                                     cleaner_list = post_resid_fallback_cleaners,
                                     fallback_model_list = post_resid_fallback_models,
                                     fallback_resid_list = fallback_residualizers) %>%
            bind_cols(select(bata_info, -term)) %>%
            mutate(ifile, SampleName, Experiment, real_df_file)
        }) %>%
        bind_rows()
    })
  res %>% map_dfr(~.$value) %>% saveRDS(paste0(path_on_server, "amlap_saved_files/results/setal_giant_regionized_fallbacks.RDS"))
  res %>% map(~.[2:4])
} %plan% elite_squad_plan
done(setal_giant_post_resid_fallbacks)

