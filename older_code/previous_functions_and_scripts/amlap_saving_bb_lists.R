library(dplyr)
library(purrr)
library(future)
library(furrr)
library(cs)

path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
path_on_server <- "/u/zburchil/workspace/launchpad/"
source(paste0(path, "amlap_bb_making_functions.R"))
source(paste0(path, "amlap_constants.R"))
source(paste0(path, "utils.R"))

done <- cs::done

ok_nodes <- get_n_best_nodes(remote_login, num_nodes,
                             !(number %in% janky_nodes),
                             number >= 79) %>%
  mutate(node = stringr::str_match(node, "node[0-9]+")[,1]) %>%
  { .$node } 
ok_nodes <- ok_nodes %>%
  test_nodes(remote_login, timeout_sec = 5, verbose = TRUE)

################################### killling jobs
killer <- kill_r_on_nodes(ok_nodes, backup_login, # "zburchil@cycle3.cs.rochester.edu"
                          wait_until_resolved = TRUE,
                          beep_on_end = TRUE)


#################################################################################
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = remote_login),
  tweak(cluster, workers = ok_nodes),
  tweak(multiprocess, workers = function() { max(round(future::availableCores()/2), 4) })
))



crossed_df <- tidyr::crossing(
  Bins = names(new_bins),
  Knum = seq_along(k_list),
  Nsubj= c(16, 32, 64, 128, 256)) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(iprefix = "full_random") ############-----------------------

made_lists %beep% {
  iterate_over_df(
    crossed_df,
    future_map,
    function(zaza) {
      bin <- new_bins[[Bins]]
      low <- bin[1]
      high<- bin[2]
      k_e <- k_list[[Knum]]
      filterers <- quos(Trial.Order  >= !!enquo(low) & 
                          Trial.Order <= !!enquo(high) & 
                          Group == "Filler-first")
      
      amlap_save_list(server_path = path_on_server,
                      real_df_file = setal_file,
                      prefix = paste0(iprefix, "_", Nsubj),
                      suffix = Bins,
                      k = k_e,
                      filter_quosures = filterers,
                      .multicore = TRUE,
                      # The function and its named additional arguments
                      single_sample_df_function = setsize_random_sample,
                      # single_sample_df_function = simple_single_order_sample,
                      n_trials = Nsubj)
    })
}
done(made_lists)



# ######### Interactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE NEEEEEEEEEEWWWWWWWWW WAY

# A flatter way of doing it
real_int_df <- tidyr::crossing(
  Knum = seq_along(k_list),
  Nsubj = c(16, 32, 64, 128, 256)) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(iprefix = "interaction") ############-----------------------

made_lists %beep% {
  future_map(
    c(1:nrow(real_int_df)),
    function(i) {
      k_e <- k_list[[real_int_df[i,]$Knum]]
      bin1_low <- new_bins$bin1[1]
      bin1_high <- new_bins$bin1[2] 
      bin3_low <- new_bins$bin3[1]
      bin3_high <-new_bins$bin3[2] 
      filterers <- quos(
        Block1 =  Trial.Order  >= !!enquo(bin1_low) & 
          Trial.Order <= !!enquo(bin1_high) & 
          Group == "Filler-first",
        Block3 =  Trial.Order  >= !!enquo(bin3_low) & 
          Trial.Order <= !!enquo(bin3_high) & 
          Group == "Filler-first")
      amlap_save_list(server_path = path_on_server,
                      real_df_file = setal_file,
                      prefix = paste0(real_int_df[i,]$iprefix, "_",
                                      real_int_df[i,]$Nsubj),
                      suffix = "bin1",
                      k = k_e,
                      filter_quosures = filterers,
                      .multicore = TRUE,
                      # The function and its named additional arguments
                      single_sample_df_function = setsize_interaction_random_sample,
                      n_trials = real_int_df[i,]$Nsubj)
    })
}
done(made_lists)



##################### Mixed models! #################################

mixed_models_df <- tidyr::crossing(
  Bins = names(new_bins),
  Knum = seq_along(k_list),
  Nsubj =  c(8, 16),
  Nitems = c(8, 16)) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(iprefix = "same_across_subj") %>%
  # To get the other configurations we want to run:
  { d <- filter(., Nsubj == Nitems, Nsubj == 8)
  bind_rows(
    .,
    mutate(d, Nsubj=32, Nitems=32),
    mutate(d, Nsubj=32, Nitems=16),
    mutate(d, Nsubj=16, Nitems=32),
    mutate(d, Nsubj=32, Nitems=64), # For testing
    mutate(d, Nsubj=64, Nitems=64))
  } %>%
  mutate(ofile = paste0(path_on_server, "amlap_saved_files/bb_",
                        paste0(iprefix, "_", Nsubj, "_", Nitems),
                        "_", Bins, "_", MinK, "_to_", MaxK, ".RDS"))

remaining_mixed_models_df <- remaining_files(mixed_models_df, file.exists(ofile))

mm_listsx %beep% {
  iterate_over_df(
    remaining_mixed_models_df,
    future_map,
    function(zaza) {
      bin <- new_bins[[Bins]]
      low <- bin[1]
      high<- bin[2]
      k_e <- k_list[[Knum]]
      filterers <- quos(Trial.Order  >= !!enquo(low) & 
                          Trial.Order <= !!enquo(high) & 
                          Group == "Filler-first")
      amlap_save_list(
        server_path = path_on_server,
        real_df_file = setal_file,
        prefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        suffix = Bins,
        k = k_e,
        filter_quosures = filterers,
        # The function and its named additional arguments
        single_sample_df_function = same_items_per_subj,
        items_per_subj = Nitems, 
        n_subj = Nsubj)
    })
} #%plan% b_squad_plan
done(mm_listsx)

###### Make the models that will be used for the parametric mixed model data #############
para_model_df <- tidyr::crossing(Bins = names(new_bins),
                                 Distribution = c("gaussian_para", "lognormal_para")) %>%
  mutate(forms = case_when(
    Distribution == "gaussian_para" ~ "RT ~ 1 + (1|Subject) + (1|Item)",
    Distribution == "lognormal_para" ~ "log10(RT) ~ 1 + (1|Subject) + (1|Item)",
    TRUE ~ NA_character_))

para_models %beep% {
  iterate_over_df(
    para_model_df,
    future_map,
    function(zaza) {
      bin <- new_bins[[Bins]]
      low <- bin[1]
      high<- bin[2]
      filterers <- quos(Trial.Order  >= !!enquo(low) & 
                          Trial.Order <= !!enquo(high) & 
                          Group == "Filler-first")
      remove_first_last_and_change_items <- function(df) {
        df %>% group_by(Item, Trial.Order) %>%
          filter(Word.Order != 1,
                 Word.Order != max(Word.Order)) %>%
          ungroup() %>%
          mutate(Item = as.factor(paste0(Word.Order, "-", Trial.Order)))
      }
      
      run_parametric_model_maker(
        df_file_name = paste0(path_on_server, setal_file),
        filterers = filterers,
        formula = as.formula(forms),
        ofile = paste0(path_on_server, "amlap_saved_files/parametric_models/", Distribution, "_", Bins, ".RDS"),
        data_f = remove_first_last_and_change_items)
    })
}
done(para_models)




########## INTERACTION MIXED MODELS!!!!!!!!!! BOTHHHHHHHHH!!!!!!!!!!!!!!!!!!!! #####################

# A flatter way of doing it
interaction_mm_df <- tidyr::crossing(
  Knum = seq_along(k_list),
  Nsubj =  c(32, 64),
  Nitems = c(32, 64)) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(iprefix = "same_across_subj_interaction")

made_lists %beep% {
  iterate_over_df(
    interaction_mm_df,
    future_map,
    function(zaza) {
      k_e <- k_list[[Knum]]
      bin1_low <- new_bins$bin1[1]
      bin1_high <- new_bins$bin1[2] 
      bin3_low <- new_bins$bin3[1]
      bin3_high <-new_bins$bin3[2] 
      filterers <- quos(
        Block1 =  Trial.Order  >= !!enquo(bin1_low) & 
          Trial.Order <= !!enquo(bin1_high) & 
          Group == "Filler-first",
        Block3 =  Trial.Order  >= !!enquo(bin3_low) & 
          Trial.Order <= !!enquo(bin3_high) & 
          Group == "Filler-first")
      
      amlap_save_list(
        server_path = path_on_server,
        real_df_file = setal_file,
        prefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        suffix = "bin1",
        k = k_e,
        filter_quosures = filterers,
        # The function and its named additional arguments
        single_sample_df_function = same_items_per_subj_interaction,
        n_trials = Nitems, 
        n_subj = Nsubj)
    })
}
done(made_lists)


################ Make the FULL-experiment data ###########################################

giant_mm_df <- tidyr::crossing(Knum = seq_along(k_list),
                               Nsubj =  c(8, 16),
                               Nitems = c(8, 16)) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(iprefix = "same_across_subj") %>%
  # To get the other configurations we want to run:
  { d <- filter(., Nsubj == Nitems, Nsubj == 8)
  bind_rows(
    .,
    mutate(d, Nsubj=32, Nitems=32),
    mutate(d, Nsubj=32, Nitems=16),
    mutate(d, Nsubj=16, Nitems=32),
    mutate(d, Nsubj=32, Nitems=64), # For testing
    mutate(d, Nsubj=64, Nitems=64))
  } %>%
  mutate(ofile = paste0(path_on_server, "amlap_saved_files/bb_",
                        paste0(iprefix, "_", Nsubj, "_", Nitems),
                        "_", "giant", "_", MinK, "_to_", MaxK, ".RDS"))

remaining_giant_models_df <- remaining_files(giant_mm_df, file.exists(ofile))

giant_mm_lists %beep% {
  iterate_over_df(
    giant_mm_df,
    future_map,
    function(zaza) {
      k_e <- k_list[[Knum]]
      filterers <- quos(Group == "Filler-first")
      amlap_save_list(
        server_path = path_on_server,
        real_df_file = setal_file,
        prefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        suffix = "giant",
        k = k_e,
        filter_quosures = filterers,
        # The function and its named additional arguments
        single_sample_df_function = same_items_per_subj,
        items_per_subj = Nitems, 
        n_subj = Nsubj)
    })
} # %plan% elite_squad_plan
done(giant_mm_lists)


giant_noexcl_mm_df <- giant_mm_df %>%
  mutate(iprefix = "same_across_subj_noexcl") %>%
  mutate(ofile = paste0(path_on_server, "amlap_saved_files/bb_",
                        paste0(iprefix, "_", Nsubj, "_", Nitems),
                        "_", "giant", "_", MinK, "_to_", MaxK, ".RDS"))

giant_noexcl_lists %beep% {
  iterate_over_df(
    giant_noexcl_mm_df,
    future_map,
    function(zaza) {
      k_e <- k_list[[Knum]]
      filterers <- quos(Group == "Filler-first")
      amlap_save_list(
        server_path = path_on_server,
        real_df_file = "setal_noexcl_fillers.RDS",
        prefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        suffix = "giant",
        k = k_e,
        filter_quosures = filterers,
        # The function and its named additional arguments
        single_sample_df_function = same_items_per_subj,
        items_per_subj = Nitems, 
        n_subj = Nsubj)
    })
} %plan% elite_squad_plan
done(giant_noexcl_lists)


######### Residualization stuff ######################################
# TO-DO: make sure residualization for bins doesn't just use the data FROM the bins
# TO-DO: make the residualization df be a function of the default mm df

bb_resid_mm_df <- tidyr::crossing(FillerItemRatio = 15, 
                                  MinK = 1, MaxK = 1,
                                  iprefix = "same_across_subj",
                                  Nsubj = c(8, 16, 32, 64),
                                  Bins = "giant", 
                                  RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         BBitems = (FillerItemRatio + 1) * Nitems) %>%
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(BBitems) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  arrange(BBitems, MinK) %>%
  mutate(ifile = paste0(
    path_on_server, "amlap_saved_files/bb_",
    iprefix, "_",  Nsubj, "_", BBitems, "_", Bins, "_",
    MinK, "_to_", MaxK, ".RDS"))

remaining_resid_bb_models_df <- remaining_files(bb_resid_mm_df, file.exists(ifile))

resid_mm_lists %beep% {
  iterate_over_df(
    remaining_resid_bb_models_df,
    future_map,
    function(zaza) {
      filterers <- quos(Group == "Filler-first")
      amlap_save_list(
        server_path = path_on_server,
        real_df_file = setal_file,
        prefix = paste0(iprefix, "_", Nsubj, "_", BBitems),
        suffix = Bins,
        k = MinK:MaxK,
        filter_quosures = filterers,
        # The function and its named additional arguments
        single_sample_df_function = same_items_per_subj,
        items_per_subj = BBitems, 
        n_subj = Nsubj)
    })
} #%plan% elite_squad_plan
done(resid_mm_lists)


############ Residualized 3-word region barebones ##################################
bb_region_mm_df <- tidyr::crossing(FillerItemRatio = 15, 
                                   MinK = 1, MaxK = 1,
                                   iprefix = "same_across_subj_regionized",
                                   Nsubj = c(8, 16, 32, 64),
                                   Bins = "giant", 
                                   RType = "residualized") %>%
  mutate(Nitems = Nsubj,
         FillerItems = FillerItemRatio * Nitems) %>%
         {
           bind_rows(rep(list(.), residual_total_files))
         } %>%
  group_by(FillerItems) %>%
  mutate(Seq = seq_along(MinK),
         MaxK = Seq * residual_k_per_file,
         MinK = MaxK - residual_k_per_file + 1) %>%
  ungroup() %>%
  arrange(FillerItems, MinK) %>%
  mutate(ifile = paste0(
    path_on_server, "amlap_saved_files/bb_",
    iprefix, "_",  Nsubj, "_", Nitems, "_", FillerItems, "_fillers_", Bins, "_",
    MinK, "_to_", MaxK, ".RDS"))

remaining_region_bb_models_df <- remaining_files(bb_region_mm_df, file.exists(ifile))

region_mm_lists %beep% {
  iterate_over_df(
    remaining_region_bb_models_df,
    future_map,
    function(zaza) {
      filterers <- quos(Group == "Filler-first")
      
      amlap_save_list(
        server_path = path_on_server,
        real_df_file = setal_file,
        prefix = NA,
        suffix = NA,
        k = MinK:MaxK,
        filter_quosures = filterers,
        .toSave = ifile,
        # The function and its named additional arguments
        single_sample_df_function = same_items_per_subj_for_residualized_regions,
        critical_regions_per_subj = Nitems, 
        filler_items_per_subj = FillerItems,
        n_subj = Nsubj,
        words_per_region = 3, 
        exclude_items_with_fewer_than = 10)
    })
}
done(region_mm_lists)


###################### FETAL DATA ######################################

fetal_giant <- giant_mm_df %>%
  mutate(iprefix = paste0("fetal_", iprefix))

fetal_mm_lists %beep% {
  iterate_over_df(
    fetal_giant,
    future_map,
    function(zaza) {
      k_e <- k_list[[Knum]]
      filterers <- quos(Group == "Filler-first")
      amlap_save_list(
        server_path = path_on_server,
        real_df_file = "fetal_sim_df.RDS",
        prefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        suffix = "giant",
        k = k_e,
        filter_quosures = filterers,
        # The function and its named additional arguments
        single_sample_df_function = same_items_per_subj,
        items_per_subj = Nitems, 
        n_subj = Nsubj)
    })
}
done(fetal_mm_lists)








###################### FRANKETAL DATA ######################################

franketal_giant <- giant_mm_df %>%
  mutate(iprefix = paste0("franketal_", iprefix)) %>%
  select(-ofile)

franketal_mm_lists %beep% {
  iterate_over_df(
    franketal_giant,
    future_map,
    function(zaza) {
      k_e <- k_list[[Knum]]
      filterers <- quos(Trial.Order > -Inf) # a dummy
      amlap_save_list(
        server_path = path_on_server,
        real_df_file = "franketal_sim_df.RDS",
        prefix = paste0(iprefix, "_", Nsubj, "_", Nitems),
        suffix = "giant",
        k = k_e,
        filter_quosures = filterers,
        # The function and its named additional arguments
        single_sample_df_function = same_items_per_subj,
        items_per_subj = Nitems, 
        n_subj = Nsubj)
    })
} %plan% elite_squad_plan
done(franketal_mm_lists)






######################################

cs::monitor_cluster_resources("zburchil", "zburchil@cycle3.cs.rochester.edu", ok_nodes,
                              save_path="/u/zburchil/wait_and_see2.RDS",
                              sleeping_time = 1,
                              total_checks = 1)

plan(remote, workers = "zburchil@cycle3.cs.rochester.edu")

hmm %<-% readRDS("/u/zburchil/wait_and_see2.RDS")
resolved(futureOf(hmm))

# hmm <- hmm %>%
# mutate(Temp=CMD, CMD = PID, PID=Temp)

hmm %>% 
  group_by(Nodename, SampleTime) %>% 
  filter(CMD =="R") %>% 
  summarise(RSS = sum(as.numeric(RSS)),
            CPU = sum(as.numeric(`%CPU`))) %>% 
  summarise(RSS = last(RSS), CPU = last(CPU)) %>%
  ggplot(aes(x=Nodename, y=RSS, fill=Nodename)) + 
  geom_bar(stat="identity")

hmm %>% 
  group_by(Nodename, SampleTime) %>% 
  filter(CMD =="R") %>% 
  summarise(RSS = sum(as.numeric(RSS)),
            CPU = sum(as.numeric(`%CPU`))) %>% 
  # summarise(RSS = last(RSS), CPU = last(CPU)) %>%
  ggplot(aes(x=SampleTime, y=RSS, color=Nodename)) + 
  geom_line()




# Saving the basic lists
# hmm %<-% {
#   future_iwalk(
#     new_bins,
#     function(bin, binname) {
#       low <- bin[1]
#       high<- bin[2]
#       filterers <- quos(Trial.Order  >= !!enquo(low) & Trial.Order <= !!enquo(high))
#       for (k_e in k_list) {
#       amlap_save_list(server_path = path_on_server,
#                       real_df_file = setal_file,
#                       prefix = "simple_start",
#                       suffix = binname,
#                       k = k_e,
#                       filter_quosures = filterers,
#                       .multicore = TRUE,
#                       # The function and its named additional arguments
#                       single_sample_df_function = structured_simple_single_sample,
#                       trials_per_subject = setal_basis_info$words_per_subject_per_trial * 5,
#                       subjects_per_group =  setal_basis_info$n_subjects / 2)
#                    }
#     })
# }
# resolved(futureOf(hmm))


#
# binname <- crossed_df[1,]$Bins
# bin <- new_bins[[binname]]
# low <- bin[1]
# high<- bin[2]
# k_e <- c(1)
# filterers <- quos(Trial.Order  >= !!enquo(low) & Trial.Order <= !!enquo(high))
# hmmm
# amlap_save_list(server_path = path_on_server,
#                 real_df_file = setal_file,
#                 prefix = "simple_start",
#                 suffix = binname,
#                 k = k_e,
#                 filter_quosures = filterers,
#                 .multicore = TRUE,
#                 # The function and its named additional arguments
#                 single_sample_df_function = structured_simple_single_sample,
#                 trials_per_subject = setal_basis_info$words_per_subject_per_trial * 5,
#                 subjects_per_group =  setal_basis_info$n_subjects / 2)

