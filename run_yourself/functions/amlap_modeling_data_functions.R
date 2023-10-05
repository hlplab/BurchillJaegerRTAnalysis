library(dplyr)
library(lme4)

# This function replaces the `catchr::catch_expr` calls I used previously
catch_warning <- function(expr, warning) {
  my_warnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(
    expr, warning=function(z) {warning(z);invokeRestart("muffleWarning")})
  val
}


# RT + effect = 10^(log(RT) + bleh)
# log(RT + effect) = log(RT) + bleh
# bleh = log10(meanRT + effect) - log10(meanRT)
log_space_effect_size_calculator <- function(meanRT, effect_size) {
  log10(meanRT + effect_size) - log10(meanRT)
}

add_in_log_space <- function(rtcol, simcol, effect) {
  10^(log10(rtcol)+effect*simcol)
}

# adds the data
amlap_add_real_data <- function(bb_df, real_df) {
  bb_df %>%
    tidyr::unnest(Word.Order) %>%
    left_join(real_df, by = c("Item", "Subject", "Word.Order"))
  # # Used to be:
  # bb_df %>%
  #   select(-Word.Order) %>%
  #   mutate(Word.Order = Words.Start) %>%
  #   left_join(real_df, by = c("Item", "Subject", "Word.Order"))
}

# Gets the brms HPDI
tidy_brms_model <- function(br_m) {
  br_m %>%
    brms::as.mcmc() %>%
    HDInterval::hdi(credMass = 0.95) %>%
    t() %>% `attributes<-`(NULL)
}
# Functions for the brms models so it doesn't have to recompile each time
load_and_use_precompiled_brms <- function(brms_file, df, formula, ...) {
  use_precompiled_brms(readRDS(brms_file), df, formula, ...)
}
use_precompiled_brms <- function(brms_obj, df, formula, ...) {
  new_m <- update(brms_obj, formula. = formula, newdata = df, ...)
  summary(new_m) # just as a precaution, this should give warning messages that will be collected (probably)
  new_m
}

#adds effects and data
amlap_fleshed_out <- function(el, el_name, real_df,
                              effect_adder_function) {
  amlap_add_real_data(el, real_df) %>%
    add_in_effects(.,
                   effect_size = 0,
                   sampled_error = 0) %>%
    effect_adder_function()
}

# gets a data frame and turns it into model output
amlap_get_stats <- function(el, el_name, real_df,
                            effect_adder_function,
                            list_of_model_funcs,
                            data_tidyer=broom::tidy) {
  df <- amlap_fleshed_out(el, el_name, real_df,
                          effect_adder_function)
  purrr::map(list_of_model_funcs,
             ~zplyr::collect_all(data_tidyer(.(df)),
                                 catchErrors=TRUE))
}

# Gets a file to a list of data frames and saves the model output
amlap_make_models <- function(ifile, ofile, real_df_file,
                              effect_adder_function, list_of_model_funcs,
                              data_tidyer=broom::tidy,
                              max_global_size_MiB=NULL) {
  if (!is.null(max_global_size_MiB))
    options("future.globals.maxSize" = max_global_size_MiB *1024^2 )
  l <- readRDS(ifile)
  real_df <- readRDS(real_df_file)
  model_list <- imap(l,
                            ~amlap_get_stats(el=.x, el_name=.y, real_df=real_df,
                                             effect_adder_function = effect_adder_function,
                                             list_of_model_funcs = list_of_model_funcs,
                                             data_tidyer=data_tidyer))
  saveRDS(model_list, ofile)
  "done!"
}

# retuns a list of the filenames, given certain variables
amlap_file_coordinator <- function(server_path,
                                   real_df_file,
                                   k,
                                   iprefix,
                                   isuffix,
                                   ofile_appended_name) {
  ifilename <- paste0(server_path, "amlap_saved_files/bb_", iprefix, "_", isuffix, "_", min(k), "_to_", max(k), ".RDS")
  ofilename <- paste0(server_path, "amlap_saved_files/modeldata_", iprefix, "_", isuffix, "_", ofile_appended_name, "_", min(k), "_to_", max(k), ".RDS")
  real_df_file <- paste0(server_path, real_df_file)
  list("ifilename" = ifilename, "ofilename" = ofilename, "real_df_file"=real_df_file)
}

# saves the fleshed out samples
amlap_bb_file_saver <- function(server_path,
                                real_df_file,
                                k,
                                iprefix,
                                isuffix,
                                ofile_appended_name,
                                effect_adder_function) {
  info <- amlap_file_coordinator(server_path, real_df_file, k, iprefix, isuffix, ofile_appended_name)
  l <- readRDS(info[["ifilename"]])
  real_df <- readRDS(info[["real_df_file"]])
  
  dfl <- future_imap(l,
                     ~amlap_fleshed_out(el=.x, el_name=.y, real_df=real_df,
                                        effect_adder_function = effect_adder_function))
  saveRDS(dfl, info[["ofilename"]])
}


# Basically sets some presents for `amlap_make_models`
amlap_model_coordinator <- function(server_path,
                                    real_df_file,
                                    k,
                                    iprefix,
                                    isuffix,
                                    ofile_appended_name,
                                    list_of_model_funcs,
                                    effect_adder_function,
                                    data_tidyer=broom::tidy,
                                    ifilepath = NULL) {
  info <- amlap_file_coordinator(server_path, real_df_file, k, iprefix, isuffix, ofile_appended_name)
  
  if (!is.null(ifilepath))  info$ifilename <- ifilepath
  
  zplyr::collect_all(
    amlap_make_models(info[["ifilename"]], info[["ofilename"]], info[["real_df_file"]],
                      effect_adder_function, list_of_model_funcs, data_tidyer),
    catchErrors = TRUE)
}

# turns the model lists into a dataframe
collect_amlap_results <- function(info_row, effect_name, server_path,
                                  .save_path = NULL,
                                  .load_temps = FALSE) {
  iname <- paste(info_row$iprefix,
                 info_row$Nsubj,
                 info_row$Bins,
                 effect_name,
                 info_row$MinK, "to",
                 info_row$MaxK, sep="_")
  # message(iname)
  if (.load_temps == TRUE)
    return(readRDS(paste0(.save_path, "temp_", iname, ".RDS")))
  
  df <- paste0(server_path,
               "amlap_saved_files/modeldata_",
               iname, ".RDS") %>%
    collate_data(.multicore = TRUE) %>%
    mutate(Nsubj = info_row$Nsubj,
           SampleName = info_row$iprefix,
           EffectName = effect_name,
           Bin = info_row$Bins)
  if (!is.null(.save_path))
    saveRDS(df, paste0(.save_path, "temp_", iname, ".RDS"))
  return(df)
}

# used to generate the parametric data
make_parametric_df <- function(parametric_func, df_modifier_function, seed_val, ...) {
  set.seed(seed_val)
  df <- data.frame(Y = parametric_func(...)) %>%
    df_modifier_function()
  df %>% mutate(SimCond = ifelse(seq_along(RT) <= nrow(df)/2, -0.5, 0.5))
}

# Generates parametric data and then applies the functions (no need to save that data)
generate_parametric_list <- function(info_row, iprefix, effect_name, server_path,
                                     effect_adder_func, list_of_model_funcs,
                                     parametric_func, df_modifier_function,
                                     ..., data_tidyer = broom::tidy) {
  oname <- paste(iprefix,
                 info_row$Nsubj,
                 info_row$Bins,
                 effect_name,
                 info_row$MinK, "to",
                 info_row$MaxK, sep="_") %>%
    zplyr::print_and_pass(message)
  l <- future_map(
    c(info_row$MinK:info_row$MaxK),
    function(el) {
      df <- make_parametric_df(parametric_func,
                               df_modifier_function,
                               seed_val = as.numeric(el),
                               ...) %>%
        effect_adder_func()
      purrr::map(list_of_model_funcs,
                 ~zplyr::collect_all(data_tidyer(.(df)),
                                     catchErrors=TRUE)) })
  if (!is.null(server_path))
    saveRDS(l, paste0(server_path, "amlap_saved_files/modeldata_", oname, ".RDS"))
  l
}


# # -------------The new way I calculated word length corrections
# hmmm %2% {
#   zplyr::collect_all({
#     library(tidyverse)
#     library(lme4)
# 
#     df_Data <- readRDS("~/workspace/launchpad/setal_sim_df.RDS")
#     # df_Data <- readRDS("~/workspace/launchpad/fetal_sim_df.RDS")
#     sub_df <- df_Data %>%
#       filter(Group == "Filler-first") %>%
#       mutate(Subject = factor(Subject))
# 
#     mRT <- lmer(RT ~ 1 + Word.Length + (1 + Word.Length | Subject),
#                 data = sub_df, control = lmerControl(optimizer="bobyqa"))
#     summary(mRT)
# 
#     mLogRT <- lmer(log10(RT) ~ 1 + Word.Length + (1 + Word.Length | Subject),
#                    data = sub_df, control = lmerControl(optimizer="bobyqa"))
#     summary(mLogRT)
# 
#     sub_df$RT_Resid <- residuals(mRT)
#     sub_df$LogRT_Resid <- residuals(mLogRT)
#     sub_df$PredRT_WL <- predict(mRT)
#     sub_df$PredLogRT_WL <- predict(mLogRT)
# 
#     new_df <- left_join(df_Data, sub_df)
#     new_df$Subject <- df_Data$Subject
# 
#     saveRDS(new_df, "~/workspace/launchpad/amlap_setal_sim_df.RDS")
#     # saveRDS(new_df, "~/workspace/launchpad/amlap_fetal_sim_df.RDS")
#     "done!"
#   }, catchErrors = TRUE)
# 
# }


# # -------------The previous way I calculated total residualized word lengths
# hmmm %2% {
#   zplyr::collect_all({
#     library(lme4)
#     df_Data <- readRDS("~/workspace/launchpad/setal_sim_df.RDS")
#     mRT <- lmer(RT ~ 1 + Word.Length + (1 + Word.Length | Subject),
#                 data = df_Data, control = lmerControl(optimizer="bobyqa"))
#     summary(mRT)
# 
#     mLogRT <- lmer(log10(RT) ~ 1 + Word.Length + (1 + Word.Length | Subject),
#                    data = df_Data, control = lmerControl(optimizer="bobyqa"))
#     summary(mLogRT)
# 
#     df_Data$RT_Resid <- residuals(mRT)
#     df_Data$LogRT_Resid <- residuals(mLogRT)
#     df_Data$PredRT_WL <- predict(mRT)
#     df_Data$PredLogRT_WL <- predict(mLogRT)
# 
#     saveRDS(df_Data, "~/workspace/launchpad/amlap_setal_sim_df.RDS")
#   }, catchErrors = TRUE)
# 
# }
# 
# mean(df_Data$Word.Length)
# predict(mRT, newdata=data.frame(Word.Length = c(4.8)), re.form=NA)
# # 327.8723
# # mean for bin 3 is 299
# # 327.8723-299
# # 28.8723
# 
# log_space_effect_size_calculator(28.8723, 56)
# 


alpha_correcter <- function(type1_p_vals,
                            desired_type_1 = 0.05,
                            tolerance = 0.003) {
  get_num <- function(ps, alph) length(ps[ps < alph])
  n <- length(type1_p_vals)
  desired_num <- n * desired_type_1
  s <- sort(type1_p_vals)
  potential_alpha <- s[round(desired_num)+1]
  num <- get_num(type1_p_vals, potential_alpha)
  # Edge cases
  if (num < desired_num) {
    alt_num <-  get_num(type1_p_vals, s[round(desired_num)+2])
    if (abs(alt_num - desired_num) < abs(num - desired_num))
      potential_alpha <- s[round(desired_num) + 2]
  } else if (num > desired_num) {
    alt_num <-  get_num(type1_p_vals, s[round(desired_num)])
    if (abs(alt_num - desired_num) < abs(num - desired_num))
      potential_alpha <- s[round(desired_num)]
  }
  new_type1 <- get_num(type1_p_vals, potential_alpha)/n
  if (!between(new_type1, desired_type_1-tolerance, desired_type_1 + tolerance))
    warning(paste0("New type1 (", new_type1, ") not in tolerance of ", desired_type_1))
  potential_alpha
}



# alpha_correcter <- function(type1_p_vals,
#                             desired_type_1 = 0.05,
#                             tolerance = 0.003) {
#   n = length(type1_p_vals)
#   desired_sig = n * desired_type_1
#   temp_l <- sort(type1_p_vals)
#   alpha = 0.5
#   prev_alpha=NULL
#   f <- function(x, y) mean(ifelse(x < y, 1, 0))
#   g <- function(x, y) abs(f(x, y) - desired_type_1) <= tolerance
#   if (!g(type1_p_vals, alpha))
# 
# 
#   while (!g(type1_p_vals, alpha)) {
#     if (f(type1_p_vals, alpha) > desired_type_1)
# 
#   }
# 
# }


# WORD LENGTH CORRECTION  --log model for log, etc.
# don't need to do it for EVERYTHING, eg not the parametric stuff
#     only for outliers removed  -- FIRST remove outliers, THEN residualize
# maybe just do one effect size (you have to make new ones)

# LogRT ~ Word.Length  ---> Resid ~ SimCond


#do residuals
# use ACTUAL normal distribution
# bring up idea that it OFTEN doesn't matter, but sometimes it does
#  to the extent that the question been explored, it's almost exclusively been done with paraetric data--it seems cryptic, but the data has been generated under the sassmmption of normality, and then analyzed with those assumptions
# 'standard procedure' 'often unquestioned'
# get 'stats' from google scholar, 'reading time' 'ANOVA'
# "power is the new type1"
# for outliers, best practice is to make a model that incorporates them, but if you can't model the latent, different processes, like going to the toilet, remove them
# ' the correct solution we have to analyze RTs is unknwon, and there are lots of complex, blah blha, what we care about here is finding out how bad things actually are.'
# make sure that 2.3, etc., the effect size doesn't change ordering, etc, and then you can just use one

# for the satellite chart, difference in normal odds, the small ones might need log "SCALE"

# for the parametric plots, put MODELS on the same facets (ie color)
# "parametrically genreated data can lead you to the wrong conclusion"  (one plot)  (only use linearly generated data) really want to reference other paper
# put outliers on the backburner
# we want the residual stuff
# also, the interaction stuff
#   but also, there are situations where it isn't try

# really need to aknowledge bruno nicenboom, waggenmaker, etc., peter gordon is the


# applies a boxcox transformation to the data, reruns the model
boxcox_tidyer <- function(m, lamda_vals=NULL,
                          secondary_tidyer = broom::tidy) {
  l <- box_cox_transformation(m, lamda_vals)
  secondary_tidyer(l$m) %>%
    mutate(lambda = l$lambda)
}

# Returns a list(m=<model>, lambda=<lambda>)
box_cox_transformation <- function(m, lamda_vals=NULL) {
  # Taken from https://stackoverflow.com/questions/33999512/how-to-use-the-box-cox-power-transformation-in-r
  powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
    boxcoxTrans <- function(x, lam1, lam2 = NULL) {
      # if we set lambda2 to zero, it becomes the one parameter transformation
      lam2 <- ifelse(is.null(lam2), 0, lam2)
      
      if (lam1 == 0L) {
        log(x + lam2)
      } else {
        (((x + lam2)^lam1) - 1) / lam1
      }
    }
    switch(method, 
           boxcox = boxcoxTrans(y, lambda1, lambda2), 
           tukey = y^lambda1
    )
  }
  
  # fuck me, this boxcox shit is awful
  formula1 <- formula(m)
  y_term <- terms(formula1)[[2]]
  if (mode(y_term) == "name")
    y_name <- as.character(y_term)
  else if (mode(y_term) == "call")
    y_name <- format(y_term)
  else stop(paste0("the y name is fucked", y_term))
  
  # FUCK NSE FUCK IT ALL
  new_df <- model.frame(m) %>% 
    mutate(NewY = !!rlang::sym(y_name))
  newy_formula <- update(formula1, NewY ~ .)
  
  if (is.null(lamda_vals))
    bc <- MASS::boxcox(newy_formula, plotit = FALSE, data=new_df)
  else
    bc <- MASS::boxcox(newy_formula, lambda=lamda_vals, plotit = FALSE, data=new_df)
  
  lambda <- bc$x[which.max(bc$y)]
  new_df <- new_df %>% 
    mutate(BoxCoxY = powerTransform(NewY, lambda))
  
  m2 <- lm(update(formula1, BoxCoxY ~ .), data = new_df)
  
  list(m = m2, lambda = lambda)
}
# 
# df <- readRDS("~/Desktop/delete2.RDS")
# m1 <- simple_models$linear_power(df)
# m2 <- simple_models$log_power(df)
# boxcox_tidyer(m1)
# boxcox_tidyer(m2)

# 
# boxcox_tidyer(m)
# 
# sa <- data.frame(x = runif(100,1,100),
#            x2 = scale(runif(100,1,100)^2),
#            fac = factor(c(rep(c("a","b","c"),33), NA))) %>%
#   mutate(y= 2*x+3*x2+rnorm(100))
# m2<-lm(cbind(log(y),x) ~ 1 + x2 + fac, data = sa)
# model.frame(m2) %>% as.data.frame() %>% head() %>% tbl_df()

# Requires lmerTest >= 3.0
lmerTestBroom <- function(m, ddf = "Satterthwaite") {
  fixed_effects <- summary(lmerTest::as_lmerModLmerTest(m), ddf = ddf)$coefficients %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "term") %>%
    purrr::set_names(c("term","estimate","std.error","df","statistic","p.value")) %>%
    select(-df) %>%
    mutate(group="fixed")
  random_effects <- broom.mixed::tidy(m) %>%
    filter(group != "fixed")
  x <- suppressWarnings(bind_rows(fixed_effects, random_effects))
}

mixed_cleaner <- function(m) {
  y_term <- terms(formula(m))[[2]]
  if (mode(y_term) == "name")
    y_name <- as.character(y_term)
  else if (mode(y_term) == "call")
    y_name <- as.character(format(y_term))
  else stop(paste0("ya fucked ", y_term))
  
  m %>%
  {cbind(lmerTestBroom(.),
         "stdevY"=sd(model.frame(.)[, y_name]),
         "meanY"=mean(model.frame(.)[, y_name]),
         "skew" = moments::skewness(model.frame(.)[, y_name]),
         "kurtosis" = moments::kurtosis(model.frame(.)[, y_name]),
         # "n_subj"=length(unique(model.frame(.)[,"UniqueSubject"])),
         "n_rows"=nrow(model.frame(.)),
         "node" = Sys.info()[["nodename"]])}
}

by_item_skewfinder <- function(df, y_name, top_n) {
  y_col <- rlang::sym(y_name)
  df <- df %>%
    as.data.frame() %>%
    group_by(UniqueItem) %>%
    summarise(skew = moments::skewness(!!y_col, na.rm=TRUE),
              kurtosis = moments::kurtosis(!!y_col, na.rm=TRUE)) 
  df %>%
    summarise(mean_skew = mean(skew, na.rm=TRUE),
              mean_kurt = mean(kurtosis, na.rm=TRUE),
              top_skew = paste(round(sort(df$skew, decreasing=TRUE)[1:top_n], digits=1),
                               collapse=","),
              top_kurt = paste(round(sort(df$kurtosis, decreasing=TRUE)[1:top_n], digits=1),
                               collapse=",")
    )
  
}

mixed_cleaner_byItem <- function(m) {
  y_term <- terms(formula(m))[[2]]
  if (mode(y_term) == "name")
    y_name <- as.character(y_term)
  else if (mode(y_term) == "call")
    y_name <- as.character(format(y_term))
  else stop(paste0("ya fucked ", y_term))
  
  df <- m %>%
  {cbind(lmerTestBroom(.),
         "stdevY"=sd(model.frame(.)[, y_name]),
         "meanY"=mean(model.frame(.)[, y_name]),
         "skew" = moments::skewness(model.frame(.)[, y_name]),
         "kurtosis" = moments::kurtosis(model.frame(.)[, y_name]),
         "n_item"=length(unique(model.frame(.)[,"Item"])),
         "n_rows"=nrow(model.frame(.)),
         "node" = Sys.info()[["nodename"]],
         by_item_skewfinder(model.frame(.), y_name, 10))}
  bad_rows  <- which(!(df$term %in% c("SimCond", "(Intercept)")))
  good_rows <- which(df$term %in% c("SimCond", "(Intercept)"))
  top <- slice(df,good_rows)
  bot <- slice(df, bad_rows)
  bot <- bot %>% mutate_at(vars(std.error:top_kurt),
                           ~ifelse(is.numeric(.), NA_real_, NA_character_))
  bind_rows(top,bot)
}



mixed_cleaner_Kenward <- function(m) {
  y_term <- terms(formula(m))[[2]]
  if (mode(y_term) == "name")
    y_name <- as.character(y_term)
  else if (mode(y_term) == "call")
    y_name <- as.character(format(y_term))
  else stop(paste0("ya fucked ", y_term))
  
  m %>%
  {cbind(lmerTestBroom(., ddf = "Kenward-Roger"),
         "stdevY"=sd(model.frame(.)[, y_name]),
         "meanY"=mean(model.frame(.)[, y_name]),
         "n_subj"=length(unique(model.frame(.)[,"UniqueSubject"])),
         "n_rows"=nrow(model.frame(.)),
         "node" = Sys.info()[["nodename"]])}
}


# Simulates a data fram parametrically
# HACK
# Is hard coded to only sample for Subject and Item
simulate_para_df_from_model <- function(df, m, backtrans_f) {
  num_subj <- n_distinct(df$UniqueSubject)
  num_item <- n_distinct(df$UniqueItem)
  num_rows <- nrow(df)
  
  coefs <- suppressWarnings(broom.mixed::tidy(m))
  subject_stdev <- coefs[coefs$term=="sd__(Intercept)" & coefs$group=="Subject",]$estimate
  item_stdev <- coefs[coefs$term=="sd__(Intercept)"  & coefs$group=="Item",]$estimate
  residual_stdev <- coefs[coefs$term=="sd__Observation" & coefs$group=="Residual",]$estimate
  
  intercept_val <- coefs[coefs$term=="(Intercept)",]$estimate
  v_subj_effects <- rnorm(num_subj, mean=0, sd = subject_stdev)
  v_item_effects <- rnorm(num_item, mean=0, sd = item_stdev)
  v_resid_effects<- rnorm(num_rows, mean=0, sd = residual_stdev)
  
  
  subj_df <- df %>% 
    group_by(UniqueSubject) %>% 
    summarise() %>% 
    bind_cols(data.frame(SubjectEffect = v_subj_effects))
  item_df <- df %>% 
    group_by(UniqueItem) %>% 
    summarise() %>% 
    bind_cols(data.frame(ItemEffect = v_item_effects))
  
  df %>% 
    left_join(subj_df, by = "UniqueSubject") %>% 
    left_join(item_df, by = "UniqueItem") %>% 
    mutate(InterceptEffect = v_resid_effects + intercept_val) %>%
    mutate(RT = InterceptEffect + SubjectEffect + ItemEffect) %>%
    mutate(RT = backtrans_f(RT))
}
# simulate_para_df_from_model <- function(df, m, backtrans_f) {
#   df <- df %>% mutate(Subject = factor(paste0("XXX", UniqueSubject)),
#                          Item = factor(paste0("XXX", UniqueItem)))
#   df$RT <- simulate(m, newdata = df, allow.new.levels = TRUE)$sim_1 %>%
#     backtrans_f()
#   df
# }
simulate_para_df_from_model_w_wl <- function(df, m, backtrans_f) {
  num_subj <- n_distinct(df$UniqueSubject)
  num_item <- n_distinct(df$UniqueItem)
  num_rows <- nrow(df)
  
  coefs <- suppressWarnings(broom::tidy(m))
  item_stdev <- coefs[coefs$term=="sd_(Intercept).Item",]$estimate
  residual_stdev <- coefs[coefs$term=="sd_Observation.Residual",]$estimate
  main_wl_effect <- coefs[coefs$term=="Word.Length",]$estimate
  subject_varcor <- VarCorr(m)$Subject
  word_lengths <- sample(m@frame$Word.Length, size = num_item,replace = TRUE)
  
  re_subj_effects <- MASS::mvrnorm(num_subj, mu=c(0,0), Sigma=subject_varcor,  empirical=TRUE) %>%
    as.data.frame() %>%
    rename(SubjIntercept = `(Intercept)`, SubjWLSlope = Word.Length)
  
  intercept_val <- coefs[coefs$term=="(Intercept)",]$estimate
  v_item_effects <- rnorm(num_item, mean=0, sd = item_stdev)
  v_resid_effects<- rnorm(num_rows, mean=0, sd = residual_stdev)
  
  subj_df <- df %>% 
    group_by(UniqueSubject) %>% 
    summarise() %>% 
    bind_cols(re_subj_effects)
  item_df <- df %>% 
    group_by(UniqueItem) %>% 
    summarise() %>% 
    bind_cols(data.frame(ItemEffect = v_item_effects)) %>%
    mutate(ActualWordLengths = word_lengths)
  
  df %>% 
    left_join(subj_df, by = "UniqueSubject") %>% 
    left_join(item_df, by = "UniqueItem") %>% 
    mutate(InterceptEffect = v_resid_effects + intercept_val) %>%
    mutate(RT = InterceptEffect + SubjIntercept + ItemEffect) %>%
    mutate(RT = RT + (SubjWLSlope + main_wl_effect) * word_lengths) %>%
    mutate(RT = backtrans_f(RT))
}





# Does NOT set seed!
run_parametrics_from_model <- function(server_path,
                                       model_file,
                                       k,
                                       iprefix,
                                       isuffix,
                                       ofile_appended_name,
                                       list_of_model_funcs,
                                       effect_adder_function,
                                       data_tidyer=broom::tidy,
                                       ifilepath = NULL,
                                       backtrans_f,
                                       exclusions = TRUE) {
  info <- amlap_file_coordinator(server_path, "", k, iprefix, isuffix, ofile_appended_name)
  
  # CURRENTLY LOADS BIG DATA ***BEFORE*** SENDING OUT TO SUBPROCESSES
  m <- readRDS(paste0(server_path, model_file))
  bb_dfs <- readRDS(info$ifilename)
  
  l <- future_map(
    bb_dfs,
    function(bb_df) {
      df <- simulate_para_df_from_model(bb_df, m, backtrans_f) %>%
      {if (exclusions == TRUE)
        filter(., RT > 100, RT < 2000)
        else . } %>%
        effect_adder_function()
      purrr::map(list_of_model_funcs,
                 ~zplyr::collect_all(data_tidyer(.(df)),
                                     catchErrors=TRUE)) })
  if (!is.null(ifilepath))
    saveRDS(l, ifilepath)
  l
}




annoying_test <- function(server_path,
                          model_file,
                          k,
                          iprefix,
                          isuffix,
                          ofile_appended_name,
                          list_of_model_funcs,
                          effect_adder_function,
                          data_tidyer=broom::tidy,
                          ifilepath = NULL,
                          backtrans_f,
                          para_func,
                          exclusions = TRUE) {
  info <- amlap_file_coordinator(server_path, "", k, iprefix, isuffix, ofile_appended_name)
  
  # CURRENTLY LOADS BIG DATA ***BEFORE*** SENDING OUT TO SUBPROCESSES
  m <- readRDS(paste0(server_path, model_file))
  bb_dfs <- readRDS(info$ifilename)
  
  l <- future_map(
    bb_dfs,
    function(bb_df) {
      df <- para_func(bb_df, m, backtrans_f) %>%  
      {if (exclusions == TRUE)
        filter(., RT > 100, RT < 2000)
        else . } %>%
        effect_adder_function()
      purrr::map(list_of_model_funcs,
                 ~zplyr::collect_all(data_tidyer(.(df)),
                                     catchErrors=TRUE)) })
  if (!is.null(ifilepath))
    saveRDS(l, ifilepath)
  l
}

# bb_l <- readRDS("/u/zburchil/workspace/launchpad/amlap_saved_files/bb_same_across_subj_8_16_bin1_1_to_2000.RDS")
# fuller_bb <- amlap_add_real_data(bb_df, real_df) %>% add_in_effects(effect_size=0,sampled_error=0)
# effect_adder_function(fuller_bb) %>% select(BATA_PredLogRT:LogRT_Resid_with_Effect)
# 
# 
# istop <- 0
# tester <- imap(bb_l, function(x,y) {
#   amlap_add_real_data(x, real_df) %>% 
#     add_in_effects(effect_size=0, sampled_error=0) %>% 
#     effect_adder_function()
#   message(y)
#   istop <<- y
#   })


##### Ugh, here's how the residualization procedure works.
# First, it marks the "Filler" data with a SimCond of 0.
# Then, you need to add the effects
# THEN, you need to remove the filler data...
# AND THEN YOU NEED TO MAKE THE RESIDUALS

set_filler_data <- function(bb_df, n_remaining_items) {
  unique_items <- unique(bb_df$UniqueItem)
  unique_items <- unique_items[1:n_remaining_items]
  
  bb_df %>%
    mutate(SimCond = ifelse(UniqueItem %in% unique_items, SimCond, 0))
}

# For the NSC data, where we want to keep a certain number of items per story
#  We also can't just keep the first n, because simcond is by item rather than subject
set_filler_data_NSC <- function(bb_df, n_remaining_items_per_story, n_conditions=2) {
  good_items <- bb_df %>%
    group_by(UniqueStory, SimCond) %>%
    filter(UniqueItem %in% unique(UniqueItem)[ 1:(n_remaining_items_per_story/n_conditions) ]) %>%
    { unique(.$UniqueItem) }
  
  bb_df %>%
    mutate(SimCond = ifelse(UniqueItem %in% good_items, SimCond, 0))
}




remove_filler_data <- function(bb_df, n_remaining_items) {
  unique_items <- unique(bb_df$UniqueItem)
  unique_items <- unique_items[1:n_remaining_items]
  
  bb_df %>%
    filter(UniqueItem %in% unique_items) %>%
    mutate(UniqueItem = factor(UniqueItem))
}

amlap_residualize_bata_one_y <- function(bb_df, outcome_var, new_col_name) {
  outcome_sym <- rlang::ensym(outcome_var)
  new_col_name <- rlang::enquo(new_col_name)
  
  predictors <- quote(1 + Word.Length + (1 + Word.Length | UniqueSubject))
  
  effect_formula <- rlang::new_formula(outcome_sym, predictors)
  
  convergence_bool <- FALSE
  m <- catch_warning(
    lme4::lmer(effect_formula, data = bb_df,
               control = lme4::lmerControl(optimizer="bobyqa")),
    warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {convergence_bool <<- TRUE})
  
  if (convergence_bool)
    bb_df[[rlang::as_name(new_col_name)]] <- NA
  else
    bb_df[[rlang::as_name(new_col_name)]] <- predict(m, newdata = bb_df, allow.new.levels = TRUE) %>% 
    rlang::set_names(NULL)
  
  # `allow.new.levels` is set to true in case an entire level (ie Subject) is comprised of NAs, which effectively removes them from the model
  # Using predict on the same data now throws an error otherwise, so to keep these predictions based on NAs from the final dataset we remove them below
  bb_df %>%
    mutate(!!new_col_name := ifelse(is.na(!!outcome_sym), NA, !!new_col_name))
}
amlap_residualize_bata_one_y_generalized <- function(bb_df, outcome_var, new_col_name, 
                                                     predictors) {
  outcome_sym <- rlang::ensym(outcome_var)
  new_col_name <- rlang::enquo(new_col_name)
  predictors <- rlang::quo_get_expr(rlang::enquo(predictors))
  effect_formula <- rlang::new_formula(outcome_sym, predictors)
  
  convergence_bool <- FALSE
  m <- catch_warning(
    lme4::lmer(effect_formula, data = bb_df,
               control = lme4::lmerControl(optimizer="bobyqa")),
    warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {convergence_bool <<- TRUE})
  
  if (convergence_bool)
    bb_df[[rlang::as_name(new_col_name)]] <- NA
  else
    bb_df[[rlang::as_name(new_col_name)]] <- predict(m, newdata = bb_df, allow.new.levels = TRUE) %>% 
    rlang::set_names(NULL)
  
  # `allow.new.levels` is set to true in case an entire level (ie Subject) is comprised of NAs, which effectively removes them from the model
  # Using predict on the same data now throws an error otherwise, so to keep these predictions based on NAs from the final dataset we remove them below
  bb_df %>%
    mutate(!!new_col_name := ifelse(is.na(!!outcome_sym), NA, !!new_col_name))
}



amlap_residualize_bata <- function(bb_df, rawEffect, logEffect, n_remaining_items) {
  rawEffect <- rlang::ensym(rawEffect)
  logEffect <- rlang::ensym(logEffect)
  predictors <- quote(1 + Word.Length + (1 + Word.Length | UniqueSubject))
  
  raw_effect_formula <- rlang::new_formula(rawEffect, predictors)
  log_effect_formula <- rlang::new_formula(logEffect, predictors)
  
  
  # These are the models without the effects:---------------------------------------------------
  logbool <- FALSE
  m_log_resid <- catch_warning(
    lme4::lmer(log10(RT) ~ 1 + Word.Length + (1 + Word.Length | UniqueSubject), data = bb_df,
               control = lme4::lmerControl(optimizer="bobyqa")),
    warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {logbool <<- TRUE})
  
  rawbool <- FALSE
  m_raw_resid <- catch_warning(
    lme4::lmer(RT ~ 1 + Word.Length + (1 + Word.Length | UniqueSubject), data = bb_df,
               control = lme4::lmerControl(optimizer="bobyqa")),
    warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {rawbool <<- TRUE})
  
  if (logbool)
    bb_df$BATA_PredLogRT <- NA
  else
    bb_df$BATA_PredLogRT <- predict(m_log_resid, newdata = bb_df) %>% 
    rlang::set_names(NULL)
  
  if (rawbool)
    bb_df$BATA_PredRawRT <- NA
  else
    bb_df$BATA_PredRawRT <- predict(m_raw_resid, newdata = bb_df) %>% 
    rlang::set_names(NULL)
  
  # This is the same thing but for the data WITH the effects included: ------------------------------
  logbool <- FALSE
  m_log_resid2 <- catch_warning(
    lme4::lmer(log_effect_formula, data = bb_df,
               control = lme4::lmerControl(optimizer="bobyqa")),
    warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {logbool <<- TRUE})
  
  rawbool <- FALSE
  m_raw_resid2 <- catch_warning(
    lme4::lmer(raw_effect_formula, data = bb_df,
               control = lme4::lmerControl(optimizer="bobyqa")),
    warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {rawbool <<- TRUE})
  
  if (logbool)
    bb_df$BATA_PredLogRT_with_Effect <- NA
  else
    bb_df$BATA_PredLogRT_with_Effect <- predict(m_log_resid2, newdata=bb_df) %>% 
    rlang::set_names(NULL)
  
  if (rawbool)
    bb_df$BATA_PredRawRT_with_Effect <- NA
  else
    bb_df$BATA_PredRawRT_with_Effect <- predict(m_raw_resid2, newdata=bb_df) %>% 
    rlang::set_names(NULL)
  
  unique_items <- unique(bb_df$UniqueItem)
  unique_items <- unique_items[1:n_remaining_items]
  
  bb_df %>%
    filter(UniqueItem %in% unique_items) %>%
    mutate(UniqueItem = factor(UniqueItem))
}


amlap_residualize_bata_no_effects <- function(bb_df, end_after, divisor=2, filter_via_zero=FALSE) {
  logbool <- FALSE
  m_log_resid <- catch_warning(
    lme4::lmer(log10(RT) ~ 1 + Word.Length + (1 + Word.Length | UniqueSubject), data = bb_df,
               control = lme4::lmerControl(optimizer="bobyqa")),
    warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {logbool <<- TRUE})
  
  rawbool <- FALSE
  m_raw_resid <- catch_warning(
    lme4::lmer(RT ~ 1 + Word.Length + (1 + Word.Length | UniqueSubject), data = bb_df,
               control = lme4::lmerControl(optimizer="bobyqa")),
    warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {rawbool <<- TRUE})
  
  if (logbool)
    bb_df$BATA_PredLogRT <- NA
  else
    bb_df$BATA_PredLogRT <- predict(m_log_resid, newdata=bb_df) %>% 
    rlang::set_names(NULL)
  
  if (rawbool)
    bb_df$BATA_PredRawRT <- NA
  else
    bb_df$BATA_PredRawRT <- predict(m_raw_resid, newdata=bb_df) %>% 
    rlang::set_names(NULL)
  
  if (filter_via_zero) {
    bb_df %>%
      filter(SimCond != 0) %>%
      mutate(UniqueItem = factor(UniqueItem))
  } else {
    unique_items <- unique(bb_df$UniqueItem)
    unique_items <- unique_items[1:(length(unique_items)/divisor)]
    
    bb_df %>%
      filter(UniqueItem %in% unique_items) %>%
      mutate(UniqueItem = factor(UniqueItem))
  }
  
}


# The old way
# If residualization models don't converge, it replaces them with NAs
# amlap_residualize_bata <- function(bb_df, rawEffect, logEffect, divisor=2, filter_via_zero=FALSE) {
#   rawEffect <- rlang::ensym(rawEffect)
#   logEffect <- rlang::ensym(logEffect)
#   predictors <- quote(1 + Word.Length + (1 + Word.Length | UniqueSubject))
#   
#   raw_effect_formula <- rlang::new_formula(rawEffect, predictors)
#   log_effect_formula <- rlang::new_formula(logEffect, predictors)
#   
#   
#   # These are the models without the effects:---------------------------------------------------
#   logbool <- FALSE
#   m_log_resid <- catchr::catch_expr(
#     lme4::lmer(log10(RT) ~ 1 + Word.Length + (1 + Word.Length | UniqueSubject), data = bb_df,
#                control = lme4::lmerControl(optimizer="bobyqa")),
#     warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {logbool <<- TRUE})
#   
#   rawbool <- FALSE
#   m_raw_resid <- catchr::catch_expr(
#     lme4::lmer(RT ~ 1 + Word.Length + (1 + Word.Length | UniqueSubject), data = bb_df,
#                control = lme4::lmerControl(optimizer="bobyqa")),
#     warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {rawbool <<- TRUE})
#   
#   if (logbool)
#     bb_df$BATA_PredLogRT <- NA
#   else
#     bb_df$BATA_PredLogRT <- predict(m_log_resid, newdata=bb_df) %>% 
#     rlang::set_names(NULL)
#   
#   if (rawbool)
#     bb_df$BATA_PredRawRT <- NA
#   else
#     bb_df$BATA_PredRawRT <- predict(m_raw_resid, newdata=bb_df) %>% 
#     rlang::set_names(NULL)
#   
#   # This is the same thing but for the data WITH the effects included: ------------------------------
#   logbool <- FALSE
#   m_log_resid2 <- catchr::catch_expr(
#     lme4::lmer(log_effect_formula, data = bb_df,
#                control = lme4::lmerControl(optimizer="bobyqa")),
#     warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {logbool <<- TRUE})
#   
#   rawbool <- FALSE
#   m_raw_resid2 <- catchr::catch_expr(
#     lme4::lmer(raw_effect_formula, data = bb_df,
#                control = lme4::lmerControl(optimizer="bobyqa")),
#     warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {rawbool <<- TRUE})
#   
#   if (logbool)
#     bb_df$BATA_PredLogRT_with_Effect <- NA
#   else
#     bb_df$BATA_PredLogRT_with_Effect <- predict(m_log_resid2, newdata=bb_df) %>% 
#     rlang::set_names(NULL)
#   
#   if (rawbool)
#     bb_df$BATA_PredRawRT_with_Effect <- NA
#   else
#     bb_df$BATA_PredRawRT_with_Effect <- predict(m_raw_resid2, newdata=bb_df) %>% 
#     rlang::set_names(NULL)
#   
#   if (filter_via_zero) {
#     bb_df %>%
#       filter(SimCond != 0) %>%
#       mutate(UniqueItem = factor(UniqueItem))
#   } else {
#     unique_items <- unique(bb_df$UniqueItem)
#     unique_items <- unique_items[1:(length(unique_items)/divisor)]
#     
#     bb_df %>%
#       filter(UniqueItem %in% unique_items) %>%
#       mutate(UniqueItem = factor(UniqueItem))
#   }
# }


# Used to turn bb_dfs into a real dfs and return them
return_fleshed_out <- function(server_path,
                               real_df_file,
                               real_k,
                               min_k, max_k,
                               iprefix,
                               isuffix,
                               ofile_appended_name) {
  ifilename <- paste0(server_path, "amlap_saved_files/bb_", 
                      iprefix, "_", isuffix, "_", min_k, "_to_", max_k, ".RDS")
  real_df_file <- paste0(server_path, real_df_file)
  
  l <- readRDS(ifilename)[real_k]
  real_df <- readRDS(real_df_file)
  
  furrr::future_map(1:length(real_k), ~amlap_add_real_data(l[[.x]], real_df)  %>%
                      add_in_effects(effect_size = 0,
                                     sampled_error = 0))
}

# Takes a row from a data frame and returns list of fleshedout data frames
row_to_fleshed_out <- function(row, k,
                               server_path = path_on_server, 
                               real_df_file) {
  stopifnot(nrow(row)==1)
  with(
    row,
    return_fleshed_out(server_path = path_on_server, real_df_file = real_df_file,
                       real_k=k,
                       min_k = MinK, max_k = MaxK, 
                       iprefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
                       isuffix = Bins,
                       ofile_appended_name = paste(RType, "size", effect_size, sep="_") %>% paste0(Space))
  )
}



# amlap_residualize_bata <- function(bb_df, divisor=2, filter_via_zero=FALSE) {
#   logbool <- FALSE
#   m_log_resid <- catchr::catch_expr(
#     lme4::lmer(log10(RT) ~ 1 + Word.Length + (1 + Word.Length | UniqueSubject), data = bb_df),
#     warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {logbool <<- TRUE})
#   
#   rawbool <- FALSE
#   m_raw_resid <- catchr::catch_expr(
#     lme4::lmer(RT ~ 1 + Word.Length + (1 + Word.Length | UniqueSubject), data = bb_df),
#     warning = function(x) if (grepl("converg", x$message, ignore.case = TRUE)) {rawbool <<- TRUE})
#   
#   if (logbool)
#     bb_df$BATA_PredLogRT <- NA
#   else
#     bb_df$BATA_PredLogRT <- predict(m_log_resid, newdata=bb_df) %>% 
#     rlang::set_names(NULL)
#   
#   if (rawbool)
#     bb_df$BATA_PredRawRT <- NA
#   else
#     bb_df$BATA_PredRawRT <- predict(m_raw_resid, newdata=bb_df) %>% 
#     rlang::set_names(NULL)
#  
#   if (filter_via_zero) {
#     bb_df %>%
#       filter(SimCond != 0) %>%
#       mutate(UniqueItem = factor(UniqueItem))
#   } else {
#   unique_items <- unique(bb_df$UniqueItem)
#   unique_items <- unique_items[1:(length(unique_items)/divisor)]
#   
#   bb_df %>%
#     filter(UniqueItem %in% unique_items) %>%
#     mutate(UniqueItem = factor(UniqueItem))
#   }
# }










# 
# manually_make_sum_subj_df <- function(df, m, backtrans_f) {
#   num_subj <- n_distinct(df$UniqueSubject)
#   num_item <- n_distinct(df$UniqueItem)
#   subjects <- paste0("SOSOS", 1:num_subj)
#   items <- paste0("SOSOS", 1:num_item)
#   df <- tidyr::crossing(Subject = subjects,
#                         Item = items) %>%
#     arrange(Subject) %>%
#     mutate_all(factor) %>%
#     mutate(Ya = seq_along(Subject)) %>%
#     mutate(SimCond = ifelse(seq_along(Subject) <= num_subj*num_item/2, -1, 1))
#   
#   df$RT <- simulate(m, newdata = df, allow.new.levels = TRUE)$sim_1 %>%
#     backtrans_f()
#   df %>% mutate(UniqueItem = Item, UniqueSubject = Subject)
# }
# 
# alter_distros <- function(df, m, val) {
#   # Old way of making df --------------------------
#   num_subj <- n_distinct(df$UniqueSubject)
#   num_item <- n_distinct(df$UniqueItem)
#   subjects <- paste0("SOSOS", 1:num_subj)
#   items <- paste0("SOSOS", 1:num_item)
#   
#   df <- tidyr::crossing(Subject = subjects,
#                         Item = items) %>%
#     arrange(Subject) %>%
#     mutate_all(factor) %>%
#     mutate(Ya = seq_along(Subject)) %>%
#     mutate(SimCond = ifelse(seq_along(Subject) <= num_subj*num_item/2, -1, 1))
#   
#   # New way of generating effects --------------------------------------------
#   
#   coefs <- suppressWarnings(broom::tidy(m))
#   subject_stdev <- coefs[coefs$term=="sd_(Intercept).Subject",]$estimate
#   item_stdev <- coefs[coefs$term=="sd_(Intercept).Item",]$estimate
#   if (val=="subj")
#     item_stdev <- subject_stdev
#   else if (val=="item")
#     subject_stdev <- item_stdev
#   else if (val == "both") {
#     subject_stdev <- (item_stdev+subject_stdev)/2.0
#     item_stdev <- subject_stdev
#   }
#   subject_effects <- rnorm(num_subj, mean=0, sd=subject_stdev)
#   item_effects <-    rnorm(num_item, mean=0, sd=item_stdev)
#   
#   intercepts <- simulate(m, newdata = df, re.form = ~0, allow.new.levels = TRUE)$sim_1
#   
#   effects <- expand.grid(Item_effects = item_effects, 
#                          Subject_effects = subject_effects) %>%
#     mutate(Intercepts = intercepts,
#            RT = Item_effects + Intercepts + Subject_effects) %>%
#     select(RT)
#   
#   stopifnot(nrow(effects) == nrow(df))
#   df <- bind_cols(df, effects) %>%
#     mutate(UniqueItem = Item, UniqueSubject = Subject)
# }
# 
# 
# why_god_why <- function(df, m, val) {
#   num_subj <- n_distinct(df$UniqueSubject)
#   num_item <- n_distinct(df$UniqueItem)
#   
#   df <- df %>% 
#     mutate(Subject = factor(paste0("ITEM", UniqueSubject)),
#            Item = factor(paste0("SUBJ", UniqueItem)))
#   
#   # -----------------------------------------------------------------------
#   coefs <- suppressWarnings(broom::tidy(m))
#   subject_stdev <- coefs[coefs$term=="sd_(Intercept).Subject",]$estimate
#   item_stdev <- coefs[coefs$term=="sd_(Intercept).Item",]$estimate
#   # ....
#   subject_effects <- rnorm(num_subj, mean=0, sd=subject_stdev)
#   item_effects <-    rnorm(num_item, mean=0, sd=item_stdev)
#   
#   intercepts <- simulate(m, newdata = df, re.form = ~0, allow.new.levels = TRUE)$sim_1
#   # -----------------------------------------------------
#   
#   subj_df <- df %>% 
#     group_by(Subject) %>% 
#     summarise() %>% 
#     bind_cols(data.frame(subjeff = subject_effects))
#   item_df <- df %>% 
#     group_by(Item) %>% 
#     summarise() %>% 
#     bind_cols(data.frame(itemeff = item_effects))
#   
#   df %>% 
#     left_join(subj_df) %>% 
#     left_join(item_df) %>% 
#     mutate(Intercepts = intercepts) %>%
#     mutate(RT = Intercepts + subject_effects + item_effects) %>%
#     mutate(UniqueItem = Item, UniqueSubject = Subject)
# }
# 
# 
# # Makes the data frame the same way as alter_distros,
# #   but attaches the effects via left_join
# why_god_why2 <- function(df, m, val) {
#   num_subj <- n_distinct(df$UniqueSubject)
#   num_item <- n_distinct(df$UniqueItem)
#   
#   # df <- df %>% 
#   #   mutate(Subject = factor(paste0("ITEM", UniqueSubject)),
#   #          Item = factor(paste0("SUBJ", UniqueItem)))
#   
#   subjects <- paste0("SOSOS", 1:num_subj)
#   items <- paste0("SOSOS", 1:num_item)
#   
#   df <- tidyr::crossing(Subject = subjects,
#                         Item = items) %>%
#     arrange(Subject) %>%
#     mutate_all(factor) %>%
#     mutate(Ya = seq_along(Subject)) %>%
#     mutate(SimCond = ifelse(seq_along(Subject) <= num_subj*num_item/2, -1, 1))
#   
#   # -----------------------------------------------------------------------
#   coefs <- suppressWarnings(broom::tidy(m))
#   subject_stdev <- coefs[coefs$term=="sd_(Intercept).Subject",]$estimate
#   item_stdev <- coefs[coefs$term=="sd_(Intercept).Item",]$estimate
#   # ....
#   subject_effects <- rnorm(num_subj, mean=0, sd=subject_stdev)
#   item_effects <-    rnorm(num_item, mean=0, sd=item_stdev)
#   
#   intercepts <- simulate(m, newdata = df, re.form = ~0, allow.new.levels = TRUE)$sim_1
#   # -----------------------------------------------------
#   
#   subj_df <- df %>% 
#     group_by(Subject) %>% 
#     summarise() %>% 
#     bind_cols(data.frame(subjeff = subject_effects))
#   item_df <- df %>% 
#     group_by(Item) %>% 
#     summarise() %>% 
#     bind_cols(data.frame(itemeff = item_effects))
#   
#   df %>% 
#     left_join(subj_df) %>% 
#     left_join(item_df) %>% 
#     mutate(Intercepts = intercepts) %>%
#     mutate(RT = Intercepts + subject_effects + item_effects) %>%
#     mutate(UniqueItem = Item, UniqueSubject = Subject)
# }
# 
# 
# # Makes the data frame exactly the
# why_god_why3 <- function(df, m, val) {
#   # Old way of making df --------------------------
#   num_subj <- n_distinct(df$UniqueSubject)
#   num_item <- n_distinct(df$UniqueItem)
#   subjects <- paste0("SOSOS", 1:num_subj)
#   items <- paste0("SOSOS", 1:num_item)
#   
#   df <- tidyr::crossing(Subject = subjects,
#                         Item = items) %>%
#     arrange(Subject) %>%
#     mutate_all(factor) %>%
#     mutate(Ya = seq_along(Subject)) %>%
#     mutate(SimCond = ifelse(seq_along(Subject) <= num_subj*num_item/2, -1, 1))
#   
#   # New way of generating effects --------------------------------------------
#   df$RT <- simulate(m, newdata = df, re.form = ~0, allow.new.levels = TRUE)$sim_1
#   df %>%
#     mutate(UniqueItem = Item, UniqueSubject = Subject)
# }
# 
# new_hope <- function(df, m, val) {
#   # Old way of making df --------------------------
#   num_subj <- n_distinct(df$UniqueSubject)
#   num_item <- n_distinct(df$UniqueItem)
#   subjects <- paste0("SOSOS", 1:num_subj)
#   items <- paste0("SOSOS", 1:num_item)
#   
#   df <- tidyr::crossing(Subject = subjects,
#                         Item = items) %>%
#     arrange(Subject) %>%
#     mutate_all(factor) %>%
#     mutate(Ya = seq_along(Subject)) %>%
#     mutate(SimCond = ifelse(seq_along(Subject) <= num_subj*num_item/2, -1, 1))
#   
#   # New way of generating effects --------------------------------------------
#   
#   coefs <- suppressWarnings(broom::tidy(m))
#   subject_stdev <- coefs[coefs$term=="sd_(Intercept).Subject",]$estimate
#   item_stdev <- coefs[coefs$term=="sd_(Intercept).Item",]$estimate
#   if (val=="subj") {
#     item_stdev <- subject_stdev
#   } else if (val=="item") {
#     subject_stdev <- item_stdev
#   } else if (val == "both") {
#     subject_stdev <- (item_stdev+subject_stdev)/2.0
#     item_stdev <- subject_stdev
#   }
#   subject_effects <- rnorm(num_subj, mean=0, sd=subject_stdev)
#   item_effects <-    rnorm(num_item, mean=0, sd=item_stdev)
#   
#   intercepts <- simulate(m, newdata = df, re.form = NULL, allow.new.levels = TRUE)$sim_1
#   
#   effects <- expand.grid(Item_effects = item_effects, 
#                          Subject_effects = subject_effects) %>%
#     mutate(Intercepts = intercepts,
#            RT = Item_effects + Intercepts + Subject_effects) %>%
#     select(RT)
#   
#   stopifnot(nrow(effects) == nrow(df))
#   df <- bind_cols(df, effects) %>%
#     mutate(UniqueItem = Item, UniqueSubject = Subject)
# }
