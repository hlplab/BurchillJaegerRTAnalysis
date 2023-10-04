library(dplyr)
library(purrr)
library(future)
library(furrr)
library(cs)
library(lme4)


path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
path_on_server <- "/u/zburchil/workspace/launchpad/"
saver_path <- "/u/zburchil/workspace/launchpad/saved_files/"

source(paste0(path, "modeling_data_functions.R"))

# setal <- readRDS(paste0(path, "../PowerSimulations/setal_sim_df.RDS"))
# fetal <- readRDS(paste0(path, "../PowerSimulations/setal_sim_df.RDS"))

# Connectivity constants ####################################
username <- "zburchil"
server <- "cycle1.cs.rochester.edu"
backup_server <- "cycle2.cs.rochester.edu"
remote_login <- paste0(username, "@", server)
backup_login <- paste0(username, "@", backup_server)
first_node = "node88"
num_nodes <- 27
janky_nodes <- c(46, 44, 62, 55, 73, 74, 78, 77, 76)

node_list <- get_n_best_nodes(remote_login, num_nodes,
                              !(number %in% janky_nodes),
                              number > 72) %>%
  mutate(node = stringr::str_match(node, "node[0-9]+")[,1]) %>%
  {map2(.$node, .$ncpus,
        function(node_name, cpus) {
          rep(node_name, round(cpus * 3/4))
        })
  } %>% unlist() %>%
  test_nodes(remote_login, timeout_sec = 2)
if (length(node_list) > 120) {
  node_list<-node_list[1:120]
  stop("Node list needs to be less than 128")
}


#############################################

# effect sizes:
# effect_sizes_df <- readRDS(paste0(path, "../stack_and_fine_effect_sizes_2018-07-06.RDS"))
# q1_setal_effect_size <- effect_sizes_df %>%
#   filter(Prediction_Number == "1 (MV in block 3)",
#          Experiment == "Harrington Stack et al. (2018)") %>%
#          {.$EffectValue[1]}
# q2_setal_effect_size <- effect_sizes_df %>%
#   filter(Prediction_Number == "2 (RC b/w 1 and 2)",
#          Experiment == "Harrington Stack et al. (2018)") %>%
#          {.$EffectValue[1]}



#####################################################################
# constants
k1 <- 1:2000
k2 <- 2001:4000
k3 <- 4001:6000
k4 <- 6001:8000
k5 <- 8001:10000
k_list <- list(k1, k2, k3, k4, k5)
if (testing==TRUE)
  k_list <- list(k1)

# formula_predictors = " ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 | UniqueItem)"
mixed_model_func_list <- list(
  ourwaypower = function(x)       lme4::lmer(ourWayEffect ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  ourwaytype1 = function(x)     lme4::lmer(ourWayNoEffect ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaypower = function(x)   lme4::lmer(theirWayEffect ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaytype1 = function(x) lme4::lmer(theirWayNoEffect ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)))
)

# Testing scaling
scaled_mixed_model_func_list <- list(
  ourwaypower = function(x)       lme4::lmer(scale(ourWayEffect)[,1] ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  ourwaytype1 = function(x)     lme4::lmer(scale(ourWayNoEffect)[,1] ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaypower = function(x)   lme4::lmer(scale(theirWayEffect)[,1] ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaytype1 = function(x) lme4::lmer(scale(theirWayNoEffect)[,1] ~ 1 + SimCond + (1 + SimCond | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)))
)

# Testing removing slopes
no_subject_slopes_mixed_model_func_list <- list(
  ourwaypower = function(x)     lme4::lmer(ourWayEffect ~ 1 + SimCond + (1 | UniqueSubject) + (1 + SimCond  | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  ourwaytype1 = function(x)     lme4::lmer(ourWayNoEffect ~ 1 + SimCond + (1 | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaypower = function(x)   lme4::lmer(theirWayEffect ~ 1 + SimCond + (1 | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))),
  theirwaytype1 = function(x) lme4::lmer(theirWayNoEffect ~ 1 + SimCond + (1 | UniqueSubject) + (1 + SimCond | UniqueItem), data = x,  control = lme4::lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)))
)





############################
# to kill all the jobs:
plan(list(
  tweak(remote, workers = backup_login),
  tweak(cluster, workers = unique(node_list))
))
killer %<-% {
  xxx<-future_walk(c(1:length(unique(node_list))),
              function(x) { names <- system("pidof R", intern=TRUE); system(paste0("kill ",names))})
  xxx
}
resolved(futureOf(killer))

###################################################

plan(list(
  tweak(remote, workers = remote_login),
  tweak(cluster, workers = node_list)
))

# no_sslopes <- run_command(ititle="bb", otitle="modeldata_no_subject_slopes",
#                          qs = c(1,2,3),
#                          exps=c("setal_fullsize"),
#                          ks = k_list,
#                          model_funcs = no_subject_slopes_mixed_model_func_list)
# resolved(no_sslopes)

test_resolver <- run_command(ititle="bb", otitle="modeldata_scaled_DV",
            qs = c(1,2,3),
            exps=c("setal_fullsize"),
            ks = k_list,
            model_funcs = scaled_mixed_model_func_list)
resolved(test_resolver)

# no_slopes <- run_command(ititle="bb", otitle="modeldata_no_slopes",
#                              qs = c(1,2,3),
#                              exps=c("setal_fullsize"),
#                              ks = k_list,
#                              model_funcs = no_slopes_mixed_model_func_list)
# resolved(no_slopes)





cs::monitor_cluster_resources("zburchil", backup_login, node_list,
                              save_path="/u/zburchil/resources_full_hour.RDS",
                              sleeping_time = 30, 120)

plan(remote, workers = "zburchil@cycle3.cs.rochester.edu")

hmm %<-% readRDS("/u/zburchil/resources_full_hour.RDS")
resolved(futureOf(hmm))


#
# testing_command %<-% {
#   walk(c(1,2,3),
#        function(q_i) {
#          walk(k_list,
#               function(k_e) {
#                 min_k <- min(k_e)
#                 max_k <- max(k_e)
#                 make_models(
#                   ifile = paste0(saver_path, "bb_",        "setal_fullsize_q", q_i, "_", min_k, "_to_", max_k, ".RDS"),
#                   ofile = paste0(saver_path, "modeldata_", "setal_fullsize_q", q_i, "_", min_k, "_to_", max_k, ".RDS"),
#                   real_df_file = paste0(path_on_server, "setal_sim_df.RDS"),
#                   list_of_model_funcs = mixed_model_func_list
#                 )
#
#               })
#        })
# }
# resolved(futureOf(testing_command))



command %<-% {

  walk(c(1,2,3),
       function(q_i) {
         walk(k_list,
              function(k_e) {
                min_k <- min(k_e)
                max_k <- max(k_e)
                make_models(
                  ifile = paste0(saver_path, "bb_",        "fetal_q", q_i, "_", min_k, "_to_", max_k, ".RDS"),
                  ofile = paste0(saver_path, "modeldata_", "fetal_q", q_i, "_", min_k, "_to_", max_k, ".RDS"),
                  real_df_file = paste0(path_on_server, "fetal_sim_df.RDS"),
                  list_of_model_funcs = mixed_model_func_list
                )

              })
       })
  walk(c(1,2,3),
       function(q_i) {
         walk(k_list,
              function(k_e) {
                min_k <- min(k_e)
                max_k <- max(k_e)
                make_models(
                  ifile = paste0(saver_path, "bb_",        "setal_fullsize_q", q_i, "_", min_k, "_to_", max_k, ".RDS"),
                  ofile = paste0(saver_path, "modeldata_", "setal_fullsize_q", q_i, "_", min_k, "_to_", max_k, ".RDS"),
                  real_df_file = paste0(path_on_server, "setal_sim_df.RDS"),
                  list_of_model_funcs = mixed_model_func_list
                )

              })
       })
  walk(c(1,2,3),
       function(q_i) {
         walk(k_list,
              function(k_e) {
                min_k <- min(k_e)
                max_k <- max(k_e)
                make_models(
                  ifile = paste0(saver_path, "bb_",        "setal_small_q", q_i, "_", min_k, "_to_", max_k, ".RDS"),
                  ofile = paste0(saver_path, "modeldata_", "setal_small_q", q_i, "_", min_k, "_to_", max_k, ".RDS"),
                  real_df_file = paste0(path_on_server, "setal_sim_df.RDS"),
                  list_of_model_funcs = mixed_model_func_list
                )

              })
       })

}
resolved(futureOf(command))




check %<-% {
  tic()
  make_models(
    ifile = "/u/zburchil/workspace/launchpad/saved_files/bb_setal_fullsize_q1_1_to_2000.RDS",
    ofile = "/u/zburchil/workspace/launchpad/saved_files/modeldata_setal_fullsize_q1_1_to_2000.RDS",
    real_df_file = paste0(path_on_server, "setal_sim_df.RDS"),
    list_of_model_funcs = mixed_model_func_list
  )
  toc()
}
resolved(futureOf(check))



plan(list(
  tweak(remote, workers = "zburchil@cycle2.cs.rochester.edu"),
  tweak(cluster, workers = node_list)
))
setal_q1_get_data %<-% {
  exp_group <- "setal_fullsize_q"
  q_i <- 2
  l<-future_map_dfr(k_list,
          function(k_e) {
            min_k <- min(k_e)
            max_k <- max(k_e)
            collate_data(paste0(saver_path, "modeldata_", exp_group, q_i, "_", min_k, "_to_", max_k, ".RDS"))})
}
resolved(futureOf(setal_q1_get_data))




check2 %<-% {
  tic()
  l <- readRDS(paste0(path_on_server,
                      "saved_files/bb_setal_fullsize_q2_1_to_100.RDS"))
  setal <- readRDS(paste0(path_on_server, "setal_sim_df.RDS"))
  model_list <- future_map(l[1:2], ~get_model_stats(., setal, mixed_model_func_list))
  saveRDS(model_list,
          paste0(path_on_server,
                 "saved_files/modeldata_setal_fullsize_q2_1_to_100.RDS"))
  toc()
}
resolved(futureOf(check2))



# l2 <- purrr::map(l, ~get_model_stats(., setal, model_func_list2))
#
#
# ###################################################################
#
# l <- make_q1_q3_barebones_list(df=setal, k=c(1:2), block_filter="3",
#                                effect_size = q1_setal_effect_size, # a list or constant value
#                                trials_per_subject = 10,
#                                subjects_per_group = 208,
#                                seed_function = identity)
# #
# l2 <- l %>%
#   purrr::map(
#     function(el) {
#       l1 <- add_real_data(el, setal)
#       l1$Dummy <- 0
#       l1_df <- avg_out_real_data(l1,
#                                  effect_size = attributes(l1)[["effect_size"]],
#                                  our_way_bare = predictedLogRT,
#                                  their_way_bare = predictedRegRT,
#                                  exclude_row_if_any_NAs = TRUE)
#     }
#   )
#
#
#
# l <- make_q1_q3_barebones_list(df=setal, k=c(1:3), block_filter="3",
#                                effect_size = 100, # a list of constant value
#                                trials_per_subject = 10,
#                                subjects_per_group = 208,
#                                seed_function = identity)
# # effect_size = l[[1]]
# l1 <- add_real_data(l[[1]], setal)
# l1$Dummy <- 0
# l1_df <- avg_out_real_data(l1,
#                            effect_size = attributes(l1)[["effect_size"]],
#                            our_way_bare = predictedLogRT,
#                            their_way_bare = predictedRegRT,
#                            exclude_row_if_any_NAs = TRUE)



#############################################################################




############
# bb_df <- make_barebones_sampler_df(setal)
# single_block <- sample_from_barebones(bb_df, 4)
# single_block$Subject <- sample_subjects(setal, nrow(single_block))
# single_block %>%
#   tidyr::unnest(Word.Order) %>%
#   left_join(setal)





# setal_full_file = "Stack et al 2018 data shared via OSF (after 1st revision)/Stack-Exp2-final-before ANY exclusions.RDS"
#
# deez <- readRDS(paste0(path,setal_full_file))
# deez %>%
#   filter(Structure != "Filler") %>%
#   filter(Detailed.Region == "ambiguous NP/PP") %>%
#   group_by(Block,Group,Ambiguity,Structure,Cond,Subject) %>%
#   summarise(n=n_distinct(Item)) %>%
#   summarise(n=max(n))
# setal %>%
#   filter(Structure != "Filler") %>%
#   group_by(Block,Group,Original.Condition,Ambiguity,Structure,Subject) %>%
#   summarise(n=n_distinct(Item)) %>%
#   summarise(n=max(n))


# " Q1:  block 3: mv stimuli slower in the rc-first GROUP (between subjects)
#   Q2: block 1 vs block 2 (rc-first group ONLY): RCs in the two blocks (within subjects?)
#   Q3: Block 2: rcs slower in rc-first group





# same_blocks=TRUE
# data_per_block = 10
#
# bb_df <- setal %>%
#   filter(Block == "3") %>%
#   make_barebones_sampler_df()
# single_block <- sample_from_barebones(bb_df, data_per_block)
# single_block$Subject <- sample_subjects(setal, nrow(single_block))
# single_block %>%
#   tidyr::unnest(Word.Order) %>%
#   left_join(setal)
#
#
#
#
# purrr::imap(1:3,~.) %>%
#   purrr::map(~`attr<-`(., "effect", c(1,2,3)))
#

saver_path <- "/u/zburchil/workspace/launchpad/saved_files/"
exp_group <- "setal_fullsize_q"

q_i <- 1
q1_df <- collate_data(paste0(saver_path, "modeldata_", exp_group, q_i, "_1_to_2000", ".RDS"))

q_i <- 2
q2_df <- collate_data(paste0(saver_path, "modeldata_", exp_group, q_i, "_1_to_2000", ".RDS"))

q_i <- 3
q3_df <- collate_data(paste0(saver_path, "modeldata_", exp_group, q_i, "_1_to_2000", ".RDS"))

q1_df %>%
  filter(!grepl("converg|Hessian", warnings) & errors == "") %>%
  group_by(type) %>%
  summarise(n=n(),
            effect = mean(ifelse(abs(statistic) > 1.96, 1, 0)),
            right = mean(ifelse(statistic > 1.96, 1, 0)))
q2_df %>%
  filter(!grepl("converg|Hessian", warnings) & errors == "") %>%
  group_by(type) %>%
  summarise(n=n(),
            effect = mean(ifelse(abs(statistic) > 1.96, 1, 0)),
            right = mean(ifelse(statistic > 1.96, 1, 0)))
q3_df %>%
  filter(!grepl("converg|Hessian", warnings) & errors == "") %>%
  group_by(type) %>%
  summarise(n=n(),
            effect = mean(ifelse(abs(statistic) > 1.96, 1, 0)),
            right = mean(ifelse(statistic < -1.96, 1, 0)))




q2_all_converge_ks <- q2_df %>%
  filter(!grepl("converg|Hessian", warnings) & errors == "") %>%
  group_by(k) %>%
  summarise(n=n()) %>% arrange(-n) %>% head(10) %>% {.$k}





