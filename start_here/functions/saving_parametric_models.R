## ---------------------------------------------------------------------------------------------
# This script saves the specifications for the LMEs that will be used to generate parametric data 
## ---------------------------------------------------------------------------------------------

# This lists stores for each of the three data (nsc, setal, fetal) some of the parameters that
# are used in the transformations of RTs prior to fitting an LMM to the (transformed) source data.
# These LMMs are then used to parametrically generate data sets in our simulations (for details,
# see text). The quantities in the code below (e.g., 15centile_shifted_para) were calculated 
# manually from the source data.
#
# This parameter is required for the log-shift transformation.
mins_and_percentiles <- list(
  "nsc" = list(
    # For each source data, we fit parametric models to the pre- and post-exclusion data (noexcl 
    # and excl). We originally considered two versions of the log-shift model: one that used the 
    # minimum RT - 1 (but at least 0) as the lower bound, and one that used the lower 15th 
    # percentile of RTs as the lower bound. We only reported the former in the paper since the 
    # latter was primarily intended as a sanity check (and indeed never performed better). We 
    # include the latter here, so that potential users of the script can see how different lower
    # bounds can be specified for the log-shift transformation.
    "noexcl" = list(
      "min_shifted_para" = 0,
      "15centile_shifted_para" = 35
    ),
    "excl"  = list(
      "min_shifted_para" = 101,
      "15centile_shifted_para" = 164
    )
  ),
  "setal" = list(
    "noexcl" = list(
      "min_shifted_para" = 1,
      "15centile_shifted_para" = 147
    ),
    "excl"  = list(
      "min_shifted_para" = 101,
      "15centile_shifted_para" = 151
    )
  ),
  "fetal" = list(
    "noexcl" = list(
      "min_shifted_para" = 45,
      "15centile_shifted_para" = 169
    ),
    "excl"  = list(
      "min_shifted_para" = 101,
      "15centile_shifted_para" = 169
    )
  )
)

# This specifies where the source data files (pre- and post-exclusion can be found)
og_data_paths <- list(
  "nsc" = list(
    "excl"   = paste0(nsc_file_path, nsc_df_file),
    "noexcl" = paste0(nsc_file_path, nsc_noexcl_df_file)
  ),
  "setal" = list(
    "excl"   = paste0(setal_file_path, setal_df_file),
    "noexcl" = paste0(setal_file_path, setal_noexcl_df_file)
  ),
  "fetal" = list(
    "excl"   = paste0(fetal_file_path, fetal_df_file),
    "noexcl" = paste0(fetal_file_path, fetal_noexcl_df_file)
  )
)

# This specifies the parametric distributions that we are considering in Study 1 to 
# generate data.
distributions <- c(
  "gaussian_para",
  "lognormal_para",
  "invsqrt_para",
  "min_shifted_para",
  "15centile_shifted_para",
  "reciprocal_para"
)

# Inverse functions for each distribution are required since our studies---following the 
# most common parametric analysis and simulation approaches in the field---do NOT actually
# use GLMMs with the above outcome distributions but rather use LMMs over RTs that have 
# been transformed by the link function of the GLMM. In order to transform the generated 
# transformed RTs back into raw RTs (e.g., for subsequent analysis by various analysis 
# approaches), we specify the inverse link functions for each parametric model.
get_inverse_function <- function(distro, dataset, w_excl) {
  excl_name <- ifelse(w_excl, "excl", "noexcl")
  possible_shifts <-mins_and_percentiles[[dataset]][[excl_name]]
  fn <- switch(
    distro,
    "gaussian_para"          = identity,
    "lognormal_para"         = function(x) 10^x,
    "invsqrt_para"           = function(x) (1/x)^2,
    "min_shifted_para"       = function(x) 10^x + possible_shifts[[distro]],
    "15centile_shifted_para" = function(x) 10^x + possible_shifts[[distro]],
    "reciprocal_para"        = function(x) 1/x
  )
  fn # should touch it to cause it to be evaluated
  return(fn)
}

# Fit each model to the source data. These models are later used to generate new RT data 
# based on the parametric assumptions they represent.
save_if_necessary <- function(distro, dataset, w_excl,
                              return_filename = FALSE) {
  excl_name <- ifelse(w_excl, "excl", "noexcl")
  filename <- paste0(model_files_path, distro, "_", dataset, "_", excl_name, "_giant.RDS")
  if (return_filename==TRUE) return(filename)
  # End if already exists
  if (file.exists(filename)) return(NULL)

  # Standard effects
  effects <- rlang::expr(1 + (1 | Subject) + (1 | Item))
  shift <- 0
  # If it has a shift!
  if (distro %in% c("min_shifted_para", "15centile_shifted_para"))
    shift <- mins_and_percentiles[[dataset]][[excl_name]][[distro]]

  df <- readRDS(og_data_paths[[dataset]][[excl_name]]) %>%
    filter(RT > shift)

  rhs <- switch(
    distro,
    "gaussian_para"          = rlang::expr(RT),
    "lognormal_para"         = rlang::expr(log10(RT)),
    "invsqrt_para"           = rlang::expr(1/sqrt(RT)),
    "min_shifted_para"       = rlang::expr(log10(RT - shift)),
    "15centile_shifted_para" = rlang::expr(log10(RT - shift)),
    "reciprocal_para"        = rlang::expr(1/RT)
    )

  formula <- rlang::new_formula(rhs, effects)
  model <- lme4::lmer(formula, data=df,
                      control = lme4::lmerControl(optimizer="bobyqa"))
  saveRDS(model, filename)

}

# RUNNING EVERYTHING ------------------------
for (distro in distributions) {
  for (w_excl in c(TRUE, FALSE)) {
    for (dataset in c("nsc","setal","fetal")) {
      save_if_necessary(distro, dataset, w_excl)
    }
  }
}
