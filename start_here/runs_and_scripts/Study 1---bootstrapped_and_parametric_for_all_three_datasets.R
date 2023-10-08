# Loading all the boilerplate. CHANGE PATH TO YOUR OWN DIRECTORY
main_path = "/Users/YYY/BurchillJaegerRTAnalysis/start_here/"

# Run this to load all the right libraries and get all the right constants for the scripts.
# This includes two libraries developed for this purpose. For details, see boilerplate.R
source(paste0(main_path, "functions/boilerplate.R"))
# Please note that, as we have done above, we provide comments in the lines preceding the 
# the code the comments apply to.


# Make the future plans to parallelize work ----------------------------------
# The optimal number of workers depends on the memory demands of the script, and
# the available memory and cores on your computer. The worker numbers present in 
# our code were chosen for a 2013 MacPro with 24 cores and 64GB RAM.
plan(list(
  tweak(multisession, workers = 5)
))

# Save all model specifications -----------------------------------
# (e.g., LMEs over raw RTs, log-RTs, log-shift RTs, reciprical RTs, etc.)
source(paste0(code_path, "saving_parametric_models.R"))


# NATURAL DF ##########################################################################
# The way these scripts are structured is that each set of BATAs is generated from a row 
# of parameters in data frames, which the columns of the data frame corresponding to 
# parameters used in creating the BATAs. This lets us easily control, edit, and monitor 
# large numbers of runs, using R's df editing code (eg dplyr)

# Below we are making the data frame that will control all the `mm` files (see top-level README)
all_natch_mm_df <- tidyr::crossing(
  # Not important here, but important to the legacy code
  Bins     = "giant",   
  Knum     = seq_along(k_list), # The batch
  Nsubj    = c(64),  
  Nitems   = c(64),
  # For the NSC data
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

# Because there are multiple mm files per BB file, we start with the mm_df and then work 
# backwards to create the corresponding bb_df so we don't need to make bb files multiple 
# times redundantly
all_natch_bb_df <- all_natch_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()

# `remaining_files()` is used often in these scripts: it detects whether certain rows of 
# the data frames already have output files saved, and if so, removes them from the data 
# frame so they aren't rerun. This is useful if a run crashes or you need to restart 
# something without redoing all the work
remaining_all_natch_bb_files <- remaining_files(all_natch_bb_df, file.exists(bb_file))

# Here we have the call that will create all the BB dfs for the previous data frame
all_bb <- {
  # Go through all the runs in the table
  iterate_over_df(  
    remaining_all_natch_bb_files,
    # `future_map()`` is from `furrr` and runs these jobs in paralell; `future_cnd_map()` 
    # is from `cs` and does the same thing but catches raised conditions so we can see errors, 
    # messages, etc. that were raised. when work first began on this project, such code was 
    # necessary, but `future` and `furrr` have come a long way since then and now likely have 
    # ways of doing something similar
    future_cnd_map,  
    # 'zaza' is just a dummy kwarg for NSE. You can ignore it.
    function(zaza) {  
      words_per_region = 1
      k <- MinK:MaxK
      k_e <- MinK:MaxK
      real_df <- readRDS(real_df_file)
      # Filtering that will be used for HS18 and F13 BATAs
      filterers <- quos(Group == "Filler-first") 

      # The HS18 and F13 share similar structure (e.g., they are both the outcome of factorial
      # experiments that cross repeated measures of participants and items). We thus developed 
      # a number of functions that contain the code that would otherwise be shared between the 
      # simulations for these two source data. The NSC data differs from the F13 and HS18 data 
      # (e.g., it is a corpus withe different repeated measures structure, see paper for 
      # details). We hence capture this special case, and code the treatement of NSC data here,
      # rather than in separate functions.
      if (Dataset == "nsc") {
        # Setting this manually, but we could have used the parameter in the df
        Nstories <- 4
        
        
        # This is not really a 'bb_df', although the variable is called that. This is essentially a 
        # cleaned up, trimmed down version of all possible sampling indices and possible subjects for 
        # those sampling indices. It's also flexible enough to accomodate word regions that are more 
        # than a single word.
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
        
        # This creates the listed of sampled BATA indices (the bb_dfs)
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
        
        #`zplyr::collect_all()` essentially just collects all the errors, messages, and 
        # warnings raised in executing the code it contains, making sure a single error doesn't 
        # derail the entire set of runs. This means that when all_bb is executed below, it will
        # return a long slew of messages and warnings, collected from the different workers. 
        # These outputs can be used for debugging.
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

# The code above was for sampling the natural BATAs. Most of the code below here is 
# now for constructing parametrically generated BATAs. I believe a quirk of the code 
# for parametrically generated BATAs requires there to be natural bb_df BATAs already 
# existing to take as an example of the structure of the BATAs in order for it to work. 
# This was a legacy hack from code where much more elaborate parametrically generated 
# BATAs might be envisioned being made. 
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

# Again, we do the same thing as we have seen above. 
para_all_bb_df <- para_all_mm_df %>%
  select(-Space, -EffectSize, -RType, -mm_file) %>% distinct()
remaining_para_all_bb_files <- remaining_files(para_all_bb_df, file.exists(bb_file))

# Note for the parametrically generated BATAs, since there ARE no indices to link to, we 
# include the actual RTs in the bb_dfs, unlike elsewhere.
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

# Here we start generating the `mm` files that contain model results.
remaining_para_all_mm_files <- remaining_files(para_all_mm_df, file.exists(mm_file))
# Saves MM files
para_all_mm <- {
  iterate_over_df(
    remaining_para_all_mm_files,
    future_cnd_map,
    function(zaa) {
      effect_size <- EffectSize
      # The function that extracts the important metrics from the model and data
      data_tidyer <- mixed_cleaner
      # The regression functions that make the models; note that we break out Type 1 and 
      # power analyses as separate analyses here. This means that there are actually some 
      # redundant Type1 analyses across different effect sizes. But, this is a very flexible 
      # way of doing it.
      model_funcs <- no_slopes 
      bb_dfs <- readRDS(bb_file)
      # This function adds the artificial effects to the BATAs
      adder_f <- function(df)
        mutate(df, RT_with_Effect = RT + SimCond * effect_size)
      
      # Here we add the effects to the BATAs, and apply each of the model funcs to the data
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

# Here we create the 'collated_file`, which is the results file that combines all 
# results of the mm_files into a single file and adds in the parameter metadata.
# We actually didn't end up using this for Study 1: this code is a bit of a relic 
# for what we did for Study 3.
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

# This section gets the distributional statistics we presented in Study 1.
# Here, we increase the amount of parallel workers because this part of the 
# process requires less RAM, allowing more cores to be used to speed up the
# computation.
plan(list(
  tweak(multisession, workers = 12)
))

# This function extracts the distributional information we want from the BATAs
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

# Since the natural BATAs need to be "fleshed out" before we can get the 
# distro stats for them, we do so below.
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


# Loads files into a single results file
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

