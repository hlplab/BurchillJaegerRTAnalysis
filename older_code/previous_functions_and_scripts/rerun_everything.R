library(dplyr)
library(purrr)
library(future)
library(furrr)
library(cs)
library(lme4)

path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
path_on_server <- "/u/zburchil/workspace/launchpad/"
saver_path <- "/u/zburchil/workspace/launchpad/saved_files/"

source(paste0(path, "barebones_making_functions.R"))
source(paste0(path, "modeling_data_functions.R"))
source(paste0(path, "amlap_modeling_data_functions.R"))
source(paste0(path, "utils.R"))

# Connectivity constants ####################################
username <- "zburchil"
server <- "cycle1.cs.rochester.edu"
backup_server <- "cycle2.cs.rochester.edu"
remote_login <- paste0(username, "@", server)
backup_login <- paste0(username, "@", backup_server)
first_node = "node88"
num_nodes <- 27
janky_nodes <- c(46, 44, 62, 55, 73, 74, 78, 77, 76)

ok_nodes <- get_n_best_nodes(remote_login, num_nodes,
                             !(number %in% janky_nodes),
                             # number != 92,
                             number > 72) %>%
  mutate(node = stringr::str_match(node, "node[0-9]+")[,1]) %>%
  { .$node }
ok_nodes <- ok_nodes %>%
  test_nodes(remote_login, timeout_sec = 5, verbose=TRUE) %>%
  sort(decreasing=T)


#####################################################################
# constants
k_per_file = 2000
total_k = 10000
k_list = split(1:total_k, ceiling(seq_along(1:total_k)/k_per_file))
if (testing==TRUE) {
  k_per_file = 10
  total_k = 10
  k_list = split(1:total_k, ceiling(seq_along(1:total_k)/k_per_file))
}

setal_file <- "setal_sim_df_sfetal_only_with_items.RDS"
fetal_file <- "fetal_sim_df_sfetal_only_with_items.RDS"
setal_scaled_file <- "setal_sim_df_sfetal_only_with_items_new_blocks.RDS"


# Models ---------

SETAL_q1_model_funcs <- list(
  ourwaypower = function(x) lme4::lmer(
    ourWayEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity | UniqueSubject) + 
      (1 | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  ourwaytype1 = function(x) lme4::lmer(
    ourWayNoEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity | UniqueSubject) + 
      (1 | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaypower = function(x) lme4::lmer(
    theirWayEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity | UniqueSubject) + 
      (1 | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaytype1 = function(x) lme4::lmer(
    theirWayNoEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity | UniqueSubject) + 
      (1 | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)))
)

SETAL_q2_model_funcs <- list(
  ourwaypower = function(x) lme4::lmer(
    ourWayEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity * SimCond | UniqueSubject) + 
      (1 + Ambiguity | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  ourwaytype1 = function(x) lme4::lmer(
    ourWayNoEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity * SimCond | UniqueSubject) + 
      (1 + Ambiguity | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaypower = function(x) lme4::lmer(
    theirWayEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity * SimCond | UniqueSubject) + 
      (1 + Ambiguity | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaytype1 = function(x) lme4::lmer(
    theirWayNoEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity * SimCond | UniqueSubject) + 
      (1 + Ambiguity | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)))
)

SETAL_q3_model_funcs <- list(
  ourwaypower = function(x) lme4::lmer(
    ourWayEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity | UniqueSubject) + 
      (1 + Ambiguity | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  ourwaytype1 = function(x) lme4::lmer(
    ourWayNoEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity | UniqueSubject) + 
      (1 + Ambiguity | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaypower = function(x) lme4::lmer(
    theirWayEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity | UniqueSubject) + 
      (1 + Ambiguity | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaytype1 = function(x) lme4::lmer(
    theirWayNoEffect ~ 1 + Ambiguity * SimCond + 
      (1 + Ambiguity | UniqueSubject) + 
      (1 + Ambiguity | UniqueItem), 
    data = x,  
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)))
)

full_RE_mixed_model_funcs <- list(
  ourwaypower = function(x) lme4::lmer(
    ourWayEffect ~ 1 + Ambiguity * SimCond +
      (1 + Ambiguity * SimCond | UniqueSubject) +
      (1 + Ambiguity * SimCond | UniqueItem),
    data = x,
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  ourwaytype1 = function(x) lme4::lmer(
    ourWayNoEffect ~ 1 + Ambiguity * SimCond +
      (1 + Ambiguity * SimCond | UniqueSubject) +
      (1 + Ambiguity * SimCond | UniqueItem),
    data = x,
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaypower = function(x) lme4::lmer(
    theirWayEffect ~ 1 + Ambiguity * SimCond +
      (1 + Ambiguity * SimCond | UniqueSubject) +
      (1 + Ambiguity * SimCond | UniqueItem),
    data = x,
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaytype1 = function(x) lme4::lmer(
    theirWayNoEffect ~ 1 + Ambiguity * SimCond +
      (1 + Ambiguity * SimCond | UniqueSubject) +
      (1 + Ambiguity * SimCond | UniqueItem),
    data = x,
    control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)))
)

# ---------

# Get the effect sizes
# effect_size_filename <- paste0(path, "../stack_and_fine_GP_effects2018-08-06.RDS")
effect_size_filename <- paste0(path, "../stack_and_fine_GP_effects_GAMM_SFETAL_only_2018-08-15.RDS")

effect_sizes_df <- readRDS(effect_size_filename) %>%
{counter_c <<- 0; .} %>% # this is just a weird system to throw 
  # an error if there are cases that don't fit in the case_when
  mutate(Experiment = case_when(
    Experiment == "Harrington Stack et al. (2018)" ~ "SETAL",
    Experiment == "Fine et al. (2013)" ~ "FETAL",
    TRUE ~ { counter_c <<- counter_c+1; stopifnot(counter_c<2); NA_character_})) %>%
    {counter_c <<- 0; .} %>%
  mutate(Question = case_when(
    grepl("MV GP Effect for", Prediction_Number) ~ "q1",
    grepl("RC GP Effect Block", Prediction_Number) ~ "q2",
    grepl("Group in Block 2", Prediction_Number) ~ "q3",
    TRUE ~ { counter_c <<- counter_c+1; stopifnot(counter_c<2); NA_character_})) %>%
  group_by(Experiment, Question) %>%
  tidyr::nest(.key = "EffectSize") %>%
  ungroup() %>%
  mutate(EffectSize = purrr::map(EffectSize,
                                 ~c(.$EffectValue[1], .$EffectValue[2], .$StdDeviation[1])))




# Getting the setal data frame
sfetal_cross_df <- tidyr::crossing(
  Knum = seq_along(k_list),
  Question = c("q1","q2","q3"),
  DataSet = c("setal_fullsize", # the regular setal database
              "fetal", 
              "setal_small", 
              "fetal_big"
              )
) %>%
  mutate(MinK = map_dbl(Knum, ~min(k_list[[.]])),
         MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         Experiment = ifelse(grepl("fetal", DataSet), "FETAL", "SETAL")) %>%
  rowwise() %>%
  mutate(ifile = paste0(saver_path, "bb_", DataSet, "_", Question,
                        "_", MinK, "_to_", MaxK, ".RDS"),
         ofile = gsub("bb_", "modeldata_", ifile),
         datafiletouse = case_when(
           DataSet == "fetal" ~ fetal_file,
           DataSet == "fetal_big" ~ setal_scaled_file,
           DataSet == "setal_fullsize" ~ setal_file,
           DataSet == "setal_small" ~ setal_file,
           TRUE ~ "NOTAFILE")) %>%
  ungroup() %>%
  zplyr::left_join(effect_sizes_df) #%>%
  # For the new stuff
  # filter((DataSet == "fetal" & Question == "q3") |
  #          (DataSet == "setal_fullsize" & Question == "q1") |  
  #          DataSet == "setal_small" | DataSet=="fetal_big") 




##############################################################


################################### killling jobs
killer <- kill_r_on_nodes(ok_nodes, backup_login, # "zburchil@cycle3.cs.rochester.edu",
                          wait_until_resolved = TRUE,
                          beep_on_end = TRUE)

#################################################################################
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = remote_login),
  # sequential,
  tweak(cluster, workers = ok_nodes),
  # sequential # for the brms models that probably do better with... wait no, i bet not
  tweak(multiprocess, workers = function() {max(round(future::availableCores()/2), 4)})
))




# Save the barebones data frames ##########################################################
##########################################################
##########################################################
##########################################################
good_bbs %beep% {
  iterate_over_df(
    sfetal_cross_df,
    future_cnd_map,
    function(i) {
      k_e <- k_list[[Knum]]
      
      # full size setal parameters: ----------------------------------
      if (DataSet == "setal_fullsize") {
        if (Question == "q1")
          partial_f <- purrr::partial(save_list,
                                      block_filter = "3",
                                      trials_per_subject = 16,
                                      subjects_per_group = 208)
        else if (Question == "q2")
          partial_f <- purrr::partial(save_list,
                                      trials_per_block1 = 32,
                                      trials_per_block2 = 18,
                                      n_subjects = 205)
        else if (Question == "q3")
          partial_f <- purrr::partial(save_list,
                                      block_filter = "2",
                                      trials_per_subject = 18,
                                      subjects_per_group = 208)
        else
          stop("Impossible question")
        # Scaled-up FETAL (uses SETAL data): lots of participants, less trials
      } else if (DataSet == "fetal_big") {
        if (Question == "q1")
          partial_f <- purrr::partial(save_list,
                                      block_filter = "3",
                                      trials_per_subject = 10,
                                      subjects_per_group = 208)
        else if (Question == "q2")
          partial_f <- purrr::partial(save_list,
                                      trials_per_block1 = 16,
                                      trials_per_block2 = 10,
                                      n_subjects = 205)
        else if (Question == "q3")
          partial_f <- purrr::partial(save_list,
                                      block_filter = "2",
                                      trials_per_subject = 10,
                                      subjects_per_group = 208)
        else
          stop("Impossible question")
        # -- Regular-sized FETAL ----------------------------------------
      } else if (DataSet %in% c("fetal")) {
        if (Question == "q1")
          partial_f <- purrr::partial(save_list,
                                      block_filter = "3",
                                      trials_per_subject = 10,
                                      subjects_per_group = 39)
        else if (Question == "q2")
          partial_f <- purrr::partial(save_list,
                                      trials_per_block1 = 16,
                                      trials_per_block2 = 10,
                                      n_subjects = 39)
        else if (Question == "q3")
          partial_f <- purrr::partial(save_list,
                                      block_filter = "2",
                                      trials_per_subject = 10,
                                      subjects_per_group = 39)
        else
          stop("Impossible question")
        # -- SETAL, but with FETAL-sized participant numbers --------
      } else if (DataSet %in% c("setal_small")) {
        if (Question == "q1")
          partial_f <- purrr::partial(save_list,
                                      block_filter = "3",
                                      trials_per_subject = 16,
                                      subjects_per_group = 39)
        else if (Question == "q2")
          partial_f <- purrr::partial(save_list,
                                      trials_per_block1 = 32,
                                      trials_per_block2 = 18,
                                      n_subjects = 39)
        else if (Question == "q3")
          partial_f <- purrr::partial(save_list,
                                      block_filter = "2",
                                      trials_per_subject = 18,
                                      subjects_per_group = 39)
        else
          stop("Impossible question")
      } else
        stop("Impossible dataset")
      
      partial_f(server_path = path_on_server,
                real_df_file = datafiletouse,
                thing_id = DataSet,
                question = Question,
                k = k_e,
                effect_size = EffectSize[[1]],
                seed_function = identity,
                multicore = TRUE)
    })
}
done(good_bbs)

# Run the models ##########################################################
##########################################################
##########################################################
##########################################################

# SETAL_q1_model_funcs

model_sfetal_df <- tidyr::crossing(
  sfetal_cross_df,
  # Models = c("full_re", "no_itemslopes", "no_subjectslopes")
  Models = c("SETAL_re")
) %>% rowwise() %>%
  mutate(ofile = gsub("modeldata_", paste0("modeldata_", Models, "_"), ofile)) %>%
  ungroup()

files_to_run_df <- remaining_files(model_sfetal_df, file.exists(ofile))

setal_models3 %beep% {
  iterate_over_df(
    filter(files_to_run_df, Models == "SETAL_re"), # files_to_run_df,
    future_cnd_map,
    function(i) {

      if (Models == "SETAL_re")
        model_funcs <- switch(Question,
                              "q1" = SETAL_q1_model_funcs,
                              "q2" = SETAL_q2_model_funcs,
                              "q3" = SETAL_q3_model_funcs,
                              stop("Incorrect question"))
      else if (Models == "full_re")
        model_funcs <- full_RE_mixed_model_funcs
      else stop("Impossible model")

      make_models(ifile = ifile,
                  ofile = ofile,
                  real_df_file = paste0(path_on_server, datafiletouse),
                  list_of_model_funcs = model_funcs, #model_funcs
                  inside_tidy_function = lmerTestBroom
      )
    }) %>%
    {saveRDS(., paste0("~/tmp/tempres_",
                       strftime(Sys.time(), format = "%m-%d_%Hh%Mm"),
                       ".RDS"))
      . }
}
done(setal_models3)

# setal_models4 %beep% {
#   future_cnd_map(
#     1:nrow(files_to_run_df),
#     function(i) {
#       mark_progress("A", i)
# 
#       r <- files_to_run_df[i,]
#       
#       if (r$Models == "full_re")
#         model_funcs <- mixed_model_func_list
#       else if (r$Models == "no_subjectslopes")
#         model_funcs <- no_subject_slopes_mixed_model_func_list
#       else if (r$Models == "no_itemslopes")
#         model_funcs <- no_item_slopes_mixed_model_func_list
#       else stop("Impossible model")
#       mark_progress("B", i)
#       
#       make_models(ifile = r$ifile,
#                   ofile = r$ofile,
#                   real_df_file = paste0(path_on_server, r$datafiletouse),
#                   list_of_model_funcs = model_funcs #model_funcs
#       )
#     }, .options=future_options(scheduling = 1)) %>%
#     { mark_progress("Z")
#       saveRDS(., paste0("~/tmp/tempres_",  
#                        strftime(Sys.time(), format = "%m-%d_%Hh%Mm"), 
#                        ".RDS"))
#       . }
# }
# done(setal_models4)


################# collecting stuff

single_file %beep% {
  iterate_over_df(
    filter(model_sfetal_df, Models == "SETAL_re"),
    future_map_dfr,
    function(za) {
      readRDS(ofile) %>%
        future_imap_dfr(
          function(dfl, iter) {
            future_imap_dfr(
              dfl,
              function(single, name_type) {
                single$value %>%
                  tibble::as_tibble() %>%
                  mutate(type = name_type) %>%
                  tidyr::nest(-type) %>%
                  mutate(messages = paste0(single$messages, collapse="--------"),
                         warnings = paste0(single$warnings, collapse="--------"),
                         errors = paste0(single$errors, collapse="--------"))
              }) %>%
              mutate(k = as.numeric(iter))
          }) %>%
        mutate(Models = Models,
               Question = Question,
               DataSet = DataSet)
    }
  ) %>%
    saveRDS(paste0(saver_path, "now_with_big_fetalv2.RDS"))
  "done!"
}
done(single_file)
