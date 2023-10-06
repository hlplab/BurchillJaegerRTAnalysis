library(dplyr) # originally used dplyr_0.7.4
library(purrr)
library(future)
library(furrr)
library(cs)

path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
path_on_server <- "/u/zburchil/workspace/launchpad/"

source(paste0(path, "modeling_data_functions.R"))
source(paste0(path, "amlap_modeling_data_functions.R"))
source(paste0(path, "amlap_constants.R"))
source(paste0(path, "utils.R"))
source(paste0(path, "amlap_bb_making_functions.R"))
source(paste0(path, "get_plans.R"))

done <- cs::done

################################### killling jobs
kill_r_on_nodes(decent_nodes$nodename, backup_login, # "zburchil@cycle3.cs.rochester.edu"
                wait_until_resolved = TRUE,
                beep_on_end = TRUE)

# Make model files
# library(lme4)
# save_model_files %beep% {
#   df <- readRDS(paste0(path_on_server, setal_file)) %>%
#     filter(Group == "Filler-first")
#
# m_raw <- lme4::lmer(RT ~ 1 + (1 | Subject) + (1 | Item), data=df,
#                     control = lme4::lmerControl(optimizer="bobyqa"))
# saveRDS(m_raw, "/u/zburchil/workspace/launchpad/amlap_saved_files/parametric_models/gaussian_para_giant.RDS")
#
# m_log <- lme4::lmer(log10(RT) ~ 1 + (1 | Subject) + (1 | Item), data=df,
#                     control = lme4::lmerControl(optimizer="bobyqa"))
# saveRDS(m_log, "/u/zburchil/workspace/launchpad/amlap_saved_files/parametric_models/lognormal_para_giant.RDS")

# m_wl_raw <- lme4::lmer(RT ~ 1 + Word.Length + (1 + Word.Length | Subject) + (1 | Item), data=df,
#                        control = lme4::lmerControl(optimizer="bobyqa"))
# saveRDS(m_wl_raw, "/u/zburchil/workspace/launchpad/amlap_saved_files/parametric_models/gaussian_para_w_wl_giant.RDS")
#
# m_wl_log <- lme4::lmer(log10(RT) ~ 1 + Word.Length + (1 + Word.Length | Subject) + (1 | Item), data=df,
#                        control = lme4::lmerControl(optimizer="bobyqa"))
# saveRDS(m_wl_log, "/u/zburchil/workspace/launchpad/amlap_saved_files/parametric_models/lognormal_para_w_wl_giant.RDS")
#
#   "done"
# } %plan% node92_plan
# done(save_model_files)

para_giant_mm_df <- tidyr::crossing(
  Bins = "giant",
  Knum = seq_along(k_list),
  Nitems = c(8, 16, 32, 64),
  Nsubj = c(8, 16, 32, 64),
  iprefix = c("same_across_subj_new"),
  Distribution = c("gaussian_para", "lognormal_para"),
  EffectSize = c(56, 80),
  Space = c("_log_space", ""),
  RType = c("unresidualized")) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]])),
         model_file = paste0(path_on_server, "amlap_saved_files/parametric_models/",
                             Distribution, "_", Bins, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_", Distribution, "_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS")) %>%
  filter(Nsubj == Nitems) %>%
  filter(Space == "") #%>% filter(Distribution != "gaussian_para")

# Check bb_dfs
stopifnot(nrow(remaining_files(para_giant_mm_df, file.exists(bb_file)))   ==0)
stopifnot(nrow(remaining_files(para_giant_mm_df, file.exists(model_file)))==0)

# Make MM files
para_remaining_files_df <- remaining_files(para_giant_mm_df, file.exists(mm_file))
para_giant_mm_files %beep% {
  tictoc::tic()
  l <- iterate_over_df(
    para_remaining_files_df,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      data_tidyer <- mixed_cleaner
      model_funcs <- no_slopes
      
      m <- readRDS(model_file)
      bb_dfs <- readRDS(bb_file)
      
      # Establish adding functions
      if (Distribution == "gaussian_para") {
        inverse_f = identity
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      } else if (Distribution == "lognormal_para") {
        inverse_f = function(x) 10^x
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      }
      
      future_map(
        bb_dfs,
        function(bb_df) {
          df <- simulate_para_df_from_model(bb_df, m, inverse_f) %>% adder_f()
          purrr::map(model_funcs,
                     ~zplyr::collect_all(data_tidyer(.(df)),
                                         catchErrors=TRUE))
        }) %>%
        saveRDS(mm_file)
    })
  tictoc::tic()
  l
} %plan% elite_squad_plan
done(para_giant_mm_files)

# Loads files
para_giant_load %beep% {
  iterate_over_df(
    para_giant_mm_df,
    future_map_dfr,
    function(i) {
      collate_data(mm_file,
                   .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               Distribution = Distribution,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               MinK = MinK)
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/parametric_giant_lmerTest.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/parametric_giant_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(para_giant_load)

#############################################################################
########## Get distributional stats #########################################
#############################################################################

sample_stats <- list(
  stdev = function(x) sd(x$RT, na.rm = TRUE),
  mu = function(x) mean(x$RT, na.rm = TRUE),
  skewness = function(x) moments::skewness(x$RT, na.rm = TRUE),
  kurtosis = function(x) moments::kurtosis(x$RT, na.rm = TRUE),
  n_rows = function(x) nrow(filter(x, !is.na(RT)))
)


para_distro_df <- para_giant_mm_df %>%
  mutate(mm_file = paste0(path_on_server, "amlap_saved_files/statdistro_", Distribution, "_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"))

# Get the stats of the parametric stuff
para_distro_remaining_files <- remaining_files(para_distro_df, file.exists(mm_file))
para_distro_files %beep% {
  iterate_over_df(
    para_distro_remaining_files,
    future_cnd_map,
    function(zaa) {
      model_funcs <- sample_stats
      m <- readRDS(model_file)
      bb_dfs <- readRDS(bb_file)
      
      if (grepl("gaussian_para", Distribution)) 
        inverse_f = identity
      else if (grepl("lognormal_para", Distribution)) 
        inverse_f = function(x) 10^x
      else stop("Unknown distribution!")
      
      future_map(
        bb_dfs,
        function(bb_df) {
          df <- simulate_para_df_from_model(bb_df, m, inverse_f)
          purrr::map(model_funcs, ~.(df)) %>% as.data.frame()
        }) %>%
        saveRDS(mm_file)
      "done"
    })
} %plan% elite_squad_plan
done(para_distro_files)


# Make the 'natural' distribution df
natch_distro_df <- para_distro_df %>% 
  filter(Distribution == first(Distribution)) %>%
  mutate(Distribution = "Natural") %>%
  mutate(mm_file = paste0(path_on_server, "amlap_saved_files/statdistro_", Distribution, "_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         real_df_file = paste0(path_on_server, setal_file))

# Get the stats of the parametric stuff
natch_distro_remaining_files <- remaining_files(natch_distro_df, file.exists(mm_file))
natch_distro_files %beep% {
  iterate_over_df(
    natch_distro_remaining_files,
    future_cnd_map,
    function(zaa) {
      model_funcs <- sample_stats
      real_df <- readRDS(real_df_file)
      bb_dfs <- readRDS(bb_file)
      
      if (!grepl("Natural", Distribution)) 
        stop("Unknown distribution!")
      
      future_map(
        bb_dfs,
        function(bb_df) {
          df <- amlap_add_real_data(bb_df, real_df)
          purrr::map(model_funcs, ~.(df)) %>% as.data.frame()
        }) %>%
        saveRDS(mm_file)
      "done"
    })
} %plan% elite_squad_plan
done(natch_distro_files)

# Loads files
distro_load %beep% {
  iterate_over_df(
    bind_rows(para_distro_df, natch_distro_df),
    future_map_dfr,
    function(i) {
      readRDS(mm_file) %>% 
        future_imap_dfr(~mutate(.x, k=.y)) %>%
        mutate(Nsubj, Nitems, Distribution, RType, Space,
               Bin = Bins, SampleName = iprefix, ifile = bb_file,
               ofile = mm_file, MinK = MinK)
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/stats_moments_new.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/stats_moments_new.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(distro_load)





#############################################################################
########## Parametric models with word length ###############################
##############################################################################

para_wl_mm_df <- tidyr::crossing(Bins = "giant",
                                 Knum = seq_along(k_list),
                                 Nitems = c(8, 16, 32, 64),
                                 Nsubj = c(8, 16, 32, 64),
                                 iprefix = c("same_across_subj_new"),
                                 Distribution = c("gaussian_para_w_wl", "lognormal_para_w_wl"),
                                 EffectSize = c(56, 80),
                                 Space = c("_log_space", ""),
                                 RType = c("unresidualized")) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]])),
         model_file = paste0(path_on_server, "amlap_saved_files/parametric_models/",
                             Distribution, "_", Bins, ".RDS"),
         mm_file = paste0(path_on_server, "amlap_saved_files/modeldata_", Distribution, "_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         bb_file = paste0(path_on_server, "amlap_saved_files/bb_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS")) %>%
  filter(Nsubj == Nitems) %>%
  filter(Space == "") #%>% filter(Distribution == "gaussian_para_w_wl")

# Check bb_dfs
stopifnot(nrow(remaining_files(para_wl_mm_df, file.exists(bb_file)))   ==0)
stopifnot(nrow(remaining_files(para_wl_mm_df, file.exists(model_file)))==0)

# Make MM files
para_remaining_files_df <- remaining_files(para_wl_mm_df, file.exists(mm_file))
para_wl_mm_files %beep% {
  iterate_over_df(
    para_remaining_files_df,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      data_tidyer <- mixed_cleaner
      model_funcs <- no_slopes
      
      m <- readRDS(model_file)
      bb_dfs <- readRDS(bb_file)
      
      # Establish adding functions
      if (grepl("gaussian_para", Distribution)) {
        inverse_f = identity
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      } else if (grepl("lognormal_para", Distribution)) {
        inverse_f = function(x) 10^x
        adder_f <- function(df) {
          mutate(df, RT_with_Effect = RT + SimCond * effect_size) }
      }
      
      future_map(
        bb_dfs,
        function(bb_df) {
          df <- simulate_para_df_from_model_w_wl(bb_df, m, inverse_f) %>% adder_f()
          purrr::map(model_funcs,
                     ~zplyr::collect_all(data_tidyer(.(df)),
                                         catchErrors=TRUE))
        }) %>%
        saveRDS(mm_file)
      "done"
    })
} %plan% elite_squad_plan
done(para_wl_mm_files)

# Loads files
para_wl_load %beep% {
  iterate_over_df(
    para_wl_mm_df,
    future_map_dfr,
    function(i) {
      collate_data(mm_file,
                   .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               Distribution = Distribution,
               SampleName = iprefix,
               Bin = Bins,
               EffectSize = EffectSize,
               RType = RType,
               Space = Space,
               ifile = bb_file,
               ofile = mm_file,
               MinK = MinK)
    }
  ) %>%
    saveRDS(paste0(path_on_server, "amlap_saved_files/results/parametric_wl_lmerTest.RDS"))
  "rsync -avzh zburchil@cycle1.cs.rochester.edu:/u/zburchil/workspace/launchpad/amlap_saved_files/results/parametric_wl_lmerTest.RDS ~/Box\ Sync/Power\ simulations\ for\ RTs/PowerSimulations/data/powerpaper/"
} %plan% elite_squad_plan
done(para_wl_load)




