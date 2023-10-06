# Loading all the boilerplate. CHANGE PATH TO YOUR OWN DIRECTORY
main_path = "/Users/tiflo/Library/CloudStorage/Box-Box/_Papers - Box/Power simulations for RTs/PowerSimulations/BurchillJaegerRTAnalysis/start_here/"
# main_path = "/Users/zburchill/Box Sync/Power simulations for RTs/PowerSimulations/BurchillJaegerRTAnalysis/start_here/"

# Run this to load all the right libraries and get all the right constants for the scripts
source(paste0(main_path, "functions/boilerplate.R"))


# Make the future plans ----------------------------------
plan(list(
  tweak(multisession, workers = 5)
))

# Save all the models -----------------------------------
source(paste0(code_path, "saving_parametric_models.R"))



# NATURAL DF ##########################################################################
all_natch_mm_df <- tidyr::crossing(
  Bins     = "giant",
  Knum     = seq_along(k_list), # The batch
  Nsubj    = c(64),
  Nitems   = c(64),
  Nstories = c(4),
  iprefix  = c("same_across_item_new", "noexcl_same_across_item_new"),
  EffectSize = c(56),
  Space    = c(""),
  Dataset  = c("nsc", "setal", "fetal"),
  RType = "unresidualized"
) %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]]))) %>%
  mutate(
    # Add in redundancy
    iprefix = paste0(Dataset,"_",iprefix),
    # Are there exclusions?
    Exclusions = ifelse(grepl("noexcl", iprefix), FALSE, TRUE),
    exclname = ifelse(Exclusions,"excl","noexcl"),
    # Barebone data files: stores indices of sampled data for each generated
    # data set (BATA)
    bb_file = paste0(bb_files_path, "bb_", iprefix, "_",
                     Nsubj, "_", Nitems, "_",
                     Bins, "_", MinK, "_to_", MaxK, ".RDS"),
    # Model summary file: stores all information from each mixed-effect analysis
    # (e.g., direction of effect, p-value, etc.)
    mm_file = paste0(mm_files_path, "modeldata_",
                     iprefix, "_", Nsubj, "_", Nitems, "_",
                     Bins, "_",
                     RType, "_size_", EffectSize, Space, "_",
                     MinK, "_to_", MaxK, ".RDS"),
    # Source data file: the indices in the barebones data file refer to this
    # data set
    real_df_file = map2_chr(Dataset, exclname, ~og_data_paths[[.x]][[.y]]))


# Natural BB dfs -----------------------------------------------------------------
all_natch_bb_df <- all_natch_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_all_natch_bb_files <- remaining_files(all_natch_bb_df, file.exists(bb_file))

all_bb <- {
  iterate_over_df(  # Go through all the runs in the table
    remaining_all_natch_bb_files,
    future_cnd_map,
    function(zaza) {  # 'zaza' is just a dummy kwarg for NSE
      words_per_region = 1
      k <- MinK:MaxK
      k_e <- MinK:MaxK
      real_df <- readRDS(real_df_file)
      filterers <- quos(Group == "Filler-first")

      if (Dataset == "nsc") {
        # Setting this manually
        Nstories <- 4

        bb_df <- real_df %>%
          group_by(Item) %>%
          summarise(maxWords = max(Word.Order)) %>%
          mutate(Words.Start = purrr::map_chr(maxWords, ~paste0(1:., collapse=","))) %>%
          tidyr::separate_rows(Words.Start, sep=",") %>%
          mutate(Words.Start = as.numeric(Words.Start)) %>%
          group_by(Item) %>%
          filter(Words.Start != maxWords,
                 Words.Start != 1,
                 Words.Start <= maxWords - words_per_region) %>%
          ungroup() %>%
          add_possible_subjects(real_df, words_per_region) %>%
          zplyr::left_join(distinct(real_df, Story, Item), by="Item")

        l <- map(
          k,
          function(k_e) {
            sample_NSC_bata(bb_df             = bb_df,
                            n_stories         = Nstories,
                            n_subj            = Nsubj,
                            n_items_per_story = Nstories/Nitems,
                            words_per_region  = words_per_region)
          }) %>%
          purrr::set_names(nm = k)

        saveRDS(l, bb_file)
        "done!"
      } else {
        zplyr::collect_all({
          l <- make_simple_sample_list(
            df = real_df, k = k_e,
            filter_quosures = filterers,
            # The function and its named additional arguments
            single_sample_df_function = same_items_per_subj,
            items_per_subj = Nitems,
            n_subj = Nsubj)
          saveRDS(l, bb_file)
          "good!"
        },
        catchErrors = TRUE)
      }
    })
}
all_bb


# PARAMETRIC DFS ###############################################################
para_all_mm_df <- all_natch_mm_df %>%
  tidyr::crossing(Distribution = distributions)  %>%
  mutate(MaxK = map_dbl(Knum, ~max(k_list[[.]])),
         MinK = map_dbl(Knum, ~min(k_list[[.]])),
         format_file = bb_file,
         model_file = pmap_chr(., function(Distribution, Dataset, Exclusions, ...) {
           save_if_necessary(Distribution, Dataset, Exclusions, TRUE)
         }),
         bb_file = paste0(bb_files_path, "bb_", Distribution, "_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          MinK, "_to_", MaxK, ".RDS"),
         mm_file = paste0(mm_files_path, "modeldata_", Distribution, "_",
                          iprefix, "_", Nsubj, "_", Nitems, "_", Bins, "_",
                          RType, "_size_", EffectSize, Space, "_",
                          MinK, "_to_", MaxK, ".RDS"))

# Parametric BB dfs -----------------------------------------------------------------
para_all_bb_df <- para_all_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_para_all_bb_files <- remaining_files(para_all_bb_df, file.exists(bb_file))

para_all_bb <- {
  iterate_over_df(
    remaining_para_all_bb_files,
    future_cnd_map,
    function(dontusethis) {
      inverse_f <- get_inverse_function(Distribution, Dataset, Exclusions)
      inverse_f

      zplyr::collect_all(
        {
          m  <- readRDS(model_file)
          format_dfs <- readRDS(format_file)
          l <- format_dfs %>%
            map(function(b_df)
              simulate_para_df_from_model(df = b_df,
                                          m = m,
                                          backtrans_f = inverse_f))
          saveRDS(l, bb_file)
          "good!"
        },
        catchErrors = TRUE)

    })
}
para_all_bb


# Parametric MM dfs -----------------------------------------------------------------
remaining_para_all_mm_files <- remaining_files(para_all_mm_df, file.exists(mm_file))
# Saves MM files
para_all_mm <- {
  iterate_over_df(
    remaining_para_all_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      data_tidyer <- mixed_cleaner
      model_funcs <- no_slopes

      bb_dfs <- readRDS(bb_file)
      adder_f <- function(df)
        mutate(df, RT_with_Effect = RT + SimCond * effect_size)

      map(
        bb_dfs,
        function(bb_df) {
          df <- bb_df %>% adder_f()
          purrr::map(model_funcs,
                     ~zplyr::collect_all(data_tidyer(.(df)),
                                         catchErrors=TRUE))
        }) %>%
        saveRDS(mm_file)
      "GOOD!"
    })
}
para_all_mm


# Save Parametric MM DFs -----------------------------------------------------
para_giant_load <- {
  iterate_over_df(
    para_all_mm_df,
    map_dfr,
    function(i) {
      collate_data(mm_file, .multicore = TRUE) %>%
        mutate(Nsubj = Nsubj,
               Nitems = Nitems,
               Distribution = Distribution,
               Dataset = Dataset,
               Exclusions = Exclusions,
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
    saveRDS(paste0(complete_path, "parametric_all_v02.RDS"))
}
para_giant_load





#############################################################################
########## Get distributional stats #########################################
#############################################################################

plan(list(
  tweak(multisession, workers = 12)
))

sample_stats_w_lambda <- list(
  stdev = function(x) sd(x$RT, na.rm = TRUE),
  mu = function(x) mean(x$RT, na.rm = TRUE),
  skewness = function(x) moments::skewness(x$RT, na.rm = TRUE),
  kurtosis = function(x) moments::kurtosis(x$RT, na.rm = TRUE),
  lambda = function(x) { x2 <- dplyr::filter(x,!is.na(x$RT),x$RT>0); b <- MASS::boxcox(x2$RT~1, plotit = F, lambda=seq(-2, 2, 1/50)); rm(x2); b$x[which.max(b$y)] },
  n_rows = function(x) nrow(dplyr::filter(x, !is.na(x$RT), x$RT>0))
)


# Get the distributional info of Parametric MMs -----------------------------
para_distro_df <- para_all_mm_df %>%
  mutate(mm_file = gsub("modeldata_","statdistro_", mm_file))
para_distro_remaining_files <- remaining_files(para_distro_df, file.exists(mm_file))

para_distro_files <- {
  iterate_over_df(
    para_distro_remaining_files,
    future_cnd_map,
    function(zaa) {
      model_funcs <- sample_stats_w_lambda
      bb_dfs <- readRDS(bb_file)

      map(
        bb_dfs,
        function(bb_df) {  
          if (Exclusions==TRUE) # Repeat post-generation exclusion
            bb_df <- filter(bb_df, RT>100, RT<2000)
          purrr::map(model_funcs, ~.(bb_df)) %>% 
            as.data.frame()
        }) %>%
        saveRDS(mm_file)
      gc()
      "done"
    })
}
para_distro_files


# Get the distro info on the natural MM dfs --------------------------
natch_distro_df <- all_natch_mm_df %>%
  mutate(Distribution = "Natural",
         mm_file = gsub("modeldata_","statdistro_", mm_file))

natch_distro_remaining_files <- remaining_files(natch_distro_df, file.exists(mm_file))
natch_distro_files <- {
  iterate_over_df(
    natch_distro_remaining_files,
    future_cnd_map,
    function(zaa) {
      model_funcs <- sample_stats_w_lambda
      real_df <- readRDS(real_df_file)
      bb_dfs <- readRDS(bb_file)

      if (!grepl("Natural", Distribution))
        stop("Unknown distribution!")

      map(
        bb_dfs,
        function(bb_df) {
          df <- amlap_add_real_data(bb_df, real_df)
          purrr::map(model_funcs, ~.(df)) %>% 
            as.data.frame()
        }) %>%
        saveRDS(mm_file)
      "done"
    })
}
natch_distro_files


# Loads files
distro_load <- {
  iterate_over_df(
    bind_rows(para_distro_df, natch_distro_df),
    map_dfr,
    function(i) {
      readRDS(mm_file) %>%
        future_imap_dfr(~mutate(.x, k=.y)) %>%
        mutate(Nsubj, Nitems, Distribution, RType, Space,
               Dataset, Exclusions,
               Bin = Bins, SampleName = iprefix, ifile = bb_file,
               ofile = mm_file, MinK = MinK)
    }
  ) %>%
    saveRDS(paste0(complete_path, "statdistro_all_v02.RDS"))

}
distro_load

