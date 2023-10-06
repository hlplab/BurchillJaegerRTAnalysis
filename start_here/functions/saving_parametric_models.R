# Saving the models that will be used to generate parametric data ----------------

# Calculated manually
mins_and_percentiles <- list(
  "nsc" = list(
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

distributions <- c(
  "gaussian_para",
  "lognormal_para",
  "invsqrt_para",
  "min_shifted_para",
  "15centile_shifted_para",
  "reciprocal_para"
)
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






