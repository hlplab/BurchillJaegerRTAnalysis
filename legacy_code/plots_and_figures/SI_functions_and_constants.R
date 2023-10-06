library(tidyverse)
library(rlang)
library(cowplot)
library(ggrepel)
library(grid)
library(gridExtra)

# Constants --------------------------------------------------------------

source(paste0(cur_dir, "../functions_and_scripts/amlap_constants.R"))

# --------------------------------------------------------------


# Plotting constants --------------------------------------------------------------
# Plotting text
space_header="Effects\nlinear in:" # With linebreak
space_text = "Effects linear in:"
odds_ratio_text = "Power advantage over\nuntransformed analysis (odds ratio)"
corrected_odds_ratio_text = "Corrected power advantage over\nuntransformed analysis (odds ratio)"
log_odds_diff_text = "Power advantage over\nuntransformed analysis (log-odds)"
corrected_log_odds_diff_text = "Corrected power advantage over\nuntransformed analysis (log-odds)"
corrected_power_text = "Power (corrected for Type I error)"
power_text = "Power (uncorrected)"

item_text = "# of items"
stories_text = "Subsampled stories"
subject_text = "Number of subjects and items"

jaeger_top_text = "Mean RTs (by trial)" #"Mean RTs (binned over trials)"
jaeger_bottom_text = "Mean RTs (by trial)"

convergence_y_text = "Proportion of convergence failures\n(out of 10,000)"
convergence_y_text_si = "Proportion of convergence failures"
singfit_y_text = "Proportion singular fits"

intdeptext1 = "super-additive"
intdeptext2 = "more scale-dependent"
intgoodtext1 = "sub-additive"
intgoodtext2 = "less scale-dependent"

intdeptext_full  = paste0(intdeptext1, "\n", "(",intdeptext2,")")
intgoodtext_full = paste0(intgoodtext1, "\n", "(",intgoodtext2,")")

shared_power_type_yrange_exp1 <- c(-0.45, 5)

# Theme stuff and dimensions:

pic_width = 6.50
pic_height = 3.10
def_dodge_width = 0.5

knitr::opts_chunk$set(fig.width = pic_width, fig.height=pic_height)

theme_set(theme_bw())
theme_update(strip.background = element_blank())
theme_update(panel.grid.minor = element_blank())
theme_update(text=element_text(size=11))
theme_update(legend.position = "bottom")

smaller_legends_theme <- theme(legend.spacing.y = unit(0, "inches"),
                               legend.key.size  = unit(0.2, "inches"))
smallest_legends_theme<- theme(legend.spacing.y = unit(0, "inches"),
                               legend.key.size  = unit(0.2, "inches"),
                               legend.text  = element_text(size=rel(0.8)),
                               legend.title = element_text(size=rel(0.8)))

# Palettes and values:

analysis_colors <- c("linear" = "#28B463",  "log" = "#A569BD",
                     "RT ~ ..." = "#28B463", "log(RT) ~ ..." = "#A569BD",
                     "lgshift" = "#FF8C00", "log(RT-x) ~ ..." = "#FF8C00") 
analysis_labels <- c("linear" = "RT ~ ...", "log" = "log(RT) ~ ...",
                     "lgshift"= "log(RT-x) ~ ...")
residual_labels <- c("linear" = "residual RT ~ ...", "log" = "residual log(RT) ~ ...",
                     "lgshift"="residual log(RT-x) ~ ...")
type_alphas <- c("power" = 1, "type1" = 0.50)
type_labels <- c("power" = "Power", "type1" = "Type I")


sampling_colors <- c(
  "Natural"               = "#E74C3C", 
  "Parametric normal"     = "#F39C12", "gaussian_para"  = "#F39C12",
  "Parametric log-normal" = "#5DADE2", "lognormal_para" = "#5DADE2"
)
sampling_labels_parens <- c(
  "Natural"               = "Natural (bootstrap)", 
  "Parametric normal"     = "Normal (parametric)",
  "gaussian_para"         = "Normal (parametric)",
  "Parametric log-normal" = "Log-normal (parametric)", 
  "lognormal_para"        = "Log-normal (parametric)"
)
sampling_labels_parens_newline <- c(
  "Natural"               = "Natural\n(bootstrap)", 
  "Parametric normal"     = "Normal\n(parametric)",
  "gaussian_para"         = "Normal\n(parametric)",
  "Parametric log-normal" = "Log-normal\n(parametric)", 
  "lognormal_para"        = "Log-normal\n(parametric)"
)

effect_size_colors <- c(
  "null"  = "#00BA38",
  "medium" = "#F8766D",
  "large" = "#619CFF"
)
effect_size_color_scale <- scale_color_manual(
  name = "Effect size",
  values = effect_size_colors
)

# For the interaction plots
effect_type_colors <- c(
  "none" = "#F8766D",
  "main" = "#7CAE00",
  "interaction" = "#00BFC4",
  "both" = "#C77CFF"
)
effect_type_color_scale <- scale_color_manual(
  name = "Effects:",
  values = effect_type_colors
)

mean_RT_colors <- c(
  high = "#f6b48f",
  medium = "#e13342",
  low = "#701f57"
)
mean_RT_color_scale <- scale_color_manual(
  name = "Mean RT",
  values = mean_RT_colors
)


# For non-power plots, we only plot one effect size in the paper
temp_paper_plot_filterer <- function(df) {
  if ("Interaction" %in% df$term | "Block" %in% df$term) {
    df %>%
      filter(term=="Block", Type %in% c("all","both"))
  } else {
    df %>%
      filter(EffectSizeName == "medium", Type=="power")
  }
}

# --------------------------------------------------------------


# Caption constants --------------------------------------------------------------
setal_items_text = "Note that there are simulations here that are not in the main paper (i.e., BATAs where the number of items per subject and the number of subjects are not equal)."
binom_ci_text = "The CIs are 95% binomial score-test-based confidence intervals."
small_ci_text = "The small CIs in some simulations are partially hidden by the mean points."
nsc_story_text = "Note that there are simulations with numbers of stories sampled shown here that are not shown in the main paper: BATAs sampling from two stories and eight."
boot_ci_text = "The CIs are 95% bootstrapped confidence intervals."
inflated_text= "The inflated Type I errors for the three-word region condition result from dictating that the critical region be in the same sentence position across sampled items."
effect_size_text = "Effect sizes have been collapsed into two categories: 'medium' (15 and 56 ms) and 'large' (35 and 80 ms)."

# For SI_nonexperimental:
jaeger_cap_intro = "Correlation between means and standard deviations (SD) of outlier-excluded word-by-word reading times (RTs, cf. Figure 1b), adapted from Figure 14 in [REFS]."
# partial_panel_a = "Panel A: Randomly drawn participants from F13. 
jaeger_cap_a = "Points reflect mean sentence RTs. Note that axis limits vary across panels." # "RTs were binned by trial order into 20 bins with equally many trials. Note that axis limits vary across panels."
jaeger_cap_b = "Panel B: Correlation across sentence items. Each data point represents one sentence. Shape and color show the sentence structure. RT means are a highly significant linear predictor of SDs (p<0.001)."

# For the omnibus tables
omnibus_table_intro <- "The results of the nested model omnibus tests for"
omnibus_table_mid <- "comparing a full interaction model"
omnibus_table_end <- "to models with progressively fewer higher-order analysis approach terms (as i increases)."


# --------------------------------------------------------------
# helper fn
unipluck <- function(df, col) {
  q <- enexpr(col)
  if (class(q) != "character")
    q <- rlang::as_name(q)
  unique(pluck(df, q))
}

# Type I functions ----------------------------------------------------------
chisq_checker <- function(col1, col2) {
  m <- as.table(rbind(col1, col2))
  chisq.test(m)$p.value
}

zetween <- function(x, range) {
  x >= range[[1]] & x <= range[[2]]
}

# Checks to see if type1s vary significantly and if they differ from 0.05
# dots are the grouping vars
type1_checker <- function(summarized_df, ...) {
  grouping_vars <- enquos(...)
  
  if (n_distinct(summarized_df$term) != 1)
    stop("You can't use interaction data here!")
  
  summarized_df %>%
    filter(Type == "type1") %>%
    mutate(NumSignificant = round(all*n)) %>%
    # select(-Type, -pos:-mean, -right) %>%
    group_by(!!!grouping_vars) %>%
    mutate(good = purrr::map2_lgl(
      NumSignificant, n, 
      ~zetween(0.05, binom.test(.x, .y, p = 0.05, alternative="t")$conf.int))) %>%
    zplyr::zummarise(
      Type1sDiffer.p.v = chisq_checker(NumSignificant, n),
      LowerType1 = paste(Analysis[all==min(all)], collapse = " AND "),
      # oops: `all` function and `all` variable named the same...
      both_near_05 = all(good),
      differs_from_05 = list(set_names(!good, Analysis)),
      type1_vals = list(set_names(Stat, Analysis)),
      logit_diff = diff(qlogis(all[order(all)]))
    )
}


# dots are to group by, eg
#   Analysis, Nsubj, SampleName, ModelNames, `mean RT`, EffectSize, Space
get_power_type1 <- function(df, ...) {
  if (n_distinct(df$term) != 1)
    stop("You can't use interaction data here!")
  
  groupers <- rlang::enquos(...)
  
  df %>%
    # Groups without type first
    group_by(!!!groupers) %>% 
    # This gets the new alpha
    mutate(new_alpha = alpha_correcter(p.value[Type=="type1"])) %>%
    # Adds the Type back in
    group_by(Type, !!!groupers) %>% 
    zplyr::zummarise(
      right_direction = mean(ifelse(p.value < 0.05 & estimate > 0, 1, 0)),
      corrected_right_direction = mean(ifelse(p.value < new_alpha & estimate > 0, 1, 0)),
      all = mean(ifelse(p.value < 0.05, 1, 0)),
      corrected_all = mean(ifelse(p.value < new_alpha, 1, 0))) %>%
    mutate(Stat = ifelse(Type=="power", right_direction, all),
           CorrectedStat = ifelse(Type=="power", corrected_right_direction, corrected_all)) %>%
    select(-right_direction,-all,-corrected_right_direction,-corrected_all)
}

alpha_correcter <- function(type1_p_vals,
                            desired_type_1 = 0.05,
                            tolerance = 0.003) {
  if (!(length(type1_p_vals) > 0)) {
    warning("No type1 p values supplied! Returning 0.05")
    return(0.05)
  }
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
# --------------------------------------------------------------


# Functions --------------------------------------------------------------
# format_p_val <- function(p.value, 
#                          thresholds = c(0.001, 0.01, 0.05)) {
#   purrr::reduce(
#     sort(thresholds),
#     function(res, next_thresh) {
#       if (p.value < next_thresh) {
#         print("XXXX")
#         rlang::done(next_thresh)
#       } else {
#         res
#       }
#     }, .init = p.value) %>% paste0()
# }


collapse_effect_sizes <- function(v, include_vals=FALSE) {
  v <- as.character(v)
  names <- c("medium", "large")
  if (include_vals==TRUE)
    names <- paste(names, c("(15/56)", "(35/80)"))
  v <- case_when(v %in% c("15", "56") ~ names[[1]],
                 v %in% c("35", "80") ~ names[[2]],
                 TRUE ~ "ERROR IN `collapse_effect_sizes`!") %>%
    factor(levels=names)
  v
}
# # General, but too slow
# collapse_effect_sizes <- function(v, ...) {
#   qs <- enquos(...)
#   if (rlang::is_empty(qs))
#     vals <- list(small = c("15","56"), large=c("35","80"))
#   else
#     vals <- map(qs, rlang::eval_tidy)
#   v <- as.character(v)
#   map_chr(v, function(x) 
#     names(vals)[(map_lgl(vals, ~x %in% .x))]) %>%
#     factor(levels=names(vals))
# }
exclude_nonconverged <- function(df, lmerTest_warnings=FALSE, opposite=FALSE) { 
  # This is the regular exclusions
  reg_excl <- df %>% 
    filter(!grepl("converg|unable to evaluate scaled", warnings), 
           errors == "") %>%
    filter(!(grepl("rank deficient|unable to evaluate scaled", messages))) 
  
  if (lmerTest_warnings == FALSE) {
    lmerTest_ones <- df %>%
      filter(!grepl("rank deficient|unable to evaluate scaled", messages), 
             errors =="") %>%
      filter(grepl("as_lmerModLT", warnings),
             !grepl("Model failed to converge|unable to evaluate scaled", warnings))
    reg_excl<- bind_rows(reg_excl, lmerTest_ones)
  }
  
  if (opposite==TRUE) {
    anti_join(df, reg_excl)
  } else {
    reg_excl
  }
}


# Supposed to sum up the mixed models
summarise_mixed_dfs <- function(df, ..., no_messages = F, no_sing_fit = F) {
  # Deal with the grouping variables ----------------------------
  groupers <- rlang::enquos(...)
  group_wo_type <- quos(term, Analysis, p_sign, !!!groupers)
  full_grouping <- quos(term, Analysis, p_sign, !!!groupers, Type, PorT)
  
  # Stuff for when we have the interaction effects. Hacky, kinda
  is_study4 <- any("Block1" %in% df$term)
  # `group_wo_type` only used for power correction
  # `type1grouper` will be defined below
  if (is_study4)
    group_wo_type <- quos(term, type1grouper, Analysis, p_sign, !!!groupers)
  
  df <- df %>%
    # Gets Analysis and Type right ------------------------------
  tidyr::separate(type, into=c("Analysis","Type","Temp"),
                  fill = "right") %>%
    mutate(Type = ifelse(Type=="fake",Temp,Type)) %>% 
    select(-Temp) %>%
    # Change term names--------------------------------------------
    mutate(term = case_when(
      term == "SimCond" ~ "MainEffect",
      term == "Block1" ~ "Block",
      term == "SimCond:Block1" ~ "Interaction",
      TRUE ~ NA_character_)) %>%
    # Decide whether the row is power or type1 ------------------
    mutate(PorT = case_when(
      term == "Interaction" ~ case_when(
        Type %in% c("all", "interaction") ~ "power",
        Type %in% c("main", "none") ~ "type1",
        TRUE ~ NA_character_),
      term == "MainEffect" ~ case_when(
        Type %in% c("all", "main", "power") ~ "power",
        Type %in% c("interaction", "none", "type1") ~ "type1",
        TRUE ~ NA_character_),
      term == "Block" ~ "power", # block is always
      TRUE ~ NA_character_)) %>%
    # Determine what sign things should be ------------------------
  { if ("Block1Code" %in% names(.))
    mutate(., p_sign = case_when(
      term %in% c("MainEffect", "Interaction") ~ 1,
      term == "Block" & Block1Code < 0 ~ -1,
      term == "Block" & Block1Code > 0 ~ 1,
      TRUE ~ NA_real_))
    else 
      mutate(., p_sign = 1)
  }  %>%
    # Gets og N
    group_by(!!!groupers, Analysis, Type, term) %>% 
    mutate(og_n=n(), 
           singfit=ifelse(grepl("singular[^A-Za-z]{0,1} fit", messages), 1, 0)) %>% 
    ungroup()
  
  # Filter out convergence warnings, etc. ------------------------
  df <- df %>% exclude_nonconverged()
  if (no_messages)
    df <- df %>% filter(messages == "")
  if (no_sing_fit)
    df <- df %>% filter(!grepl("singular[^A-Za-z]{0,1} fit", messages))
  
  { if (isTRUE(is_study4)) 
    df %>%
      # Adding the int stuff for better type1 grouping
      mutate(type1grouper = case_when(
        term=="Interaction" ~ case_when(
          Type %in% c("main","all") ~ "alltype1er-int",
          Type %in% c("none","interaction") ~ "inttype1er",
          TRUE ~ NA_character_),
        term=="Interaction" ~ case_when(
          Type %in% c("interaction","all") ~ "alltype1er-main",
          Type %in% c("none","main") ~ "maintype1er",
          TRUE ~ NA_character_),
        TRUE ~ "whocares")
      )
    else
      df
  } %>%
    # Groups without type first
    group_by(!!!group_wo_type) %>% 
    # This gets the new alpha
    # STOP MAKE SURE alpha_correcter can handle null/empty vals cuz of block
    mutate(new_alpha = alpha_correcter(p.value[PorT=="type1"])) %>%
    # Adds the Type back in
    group_by(!!!full_grouping) %>% 
    summarise(all =   mean(ifelse(p.value < 0.05, 1, 0)),
              pos =   mean(ifelse(p.value < 0.05 & estimate > 0, 1, 0)),
              neg =   mean(ifelse(p.value < 0.05 & estimate < 0, 1, 0)),
              mean =  mean(estimate),
              n=n(),
              og_n=first(og_n),
              n_singfits = sum(singfit),
              # Get the more complicated stuf-------------------------------------
              right = mean(case_when(
                p_sign > 0 ~ ifelse(p.value < 0.05 & estimate > 0, 1, 0),
                p_sign < 0 ~ ifelse(p.value < 0.05 & estimate < 0, 1, 0),
                TRUE ~ NA_real_)),
              corrected_right = mean(case_when(
                p_sign > 0 ~ ifelse(p.value < new_alpha & estimate > 0, 1, 0),
                p_sign < 0 ~ ifelse(p.value < new_alpha & estimate < 0, 1, 0),
                TRUE ~ NA_real_)),
              corrected_all = mean(ifelse(p.value < new_alpha, 1, 0))) %>%
    ungroup() %>%
    mutate(Stat = ifelse(PorT == "power", right, all),
           CorrectedStat = ifelse(PorT == "power", corrected_right, corrected_all)) %>%
    select(-p_sign) %>%
    select(term, Analysis, Type, everything())
}


# Cleans up some of the more hard-to-understand values/cols
col_cleaner <- function(df) {
  cols <- names(df)
  
  # Turn 'Bin' into 'Mean RT'
  if ("Bin" %in% cols)
    df <- df %>% 
      mutate(`mean RT` = case_when(
        Bin == "bin1" ~ "high",
        Bin == "bin2" ~ "medium",
        Bin == "bin3" ~ "low",
        TRUE ~ NA_character_) %>%
          factor(levels=c("high", "medium", "low"))) %>%
      select(-Bin)
  
  # Turn Block1Code to 'IntDirection'
  if ("Block1Code" %in% cols)
    df <- df %>% 
      mutate(IntDirection = case_when(
        Block1Code == -1 ~ intgoodtext1,
        Block1Code == 1 ~ intdeptext1,
        TRUE ~ NA_character_) %>%
          factor(levels=c(intgoodtext1, intdeptext1)))
  
  # Turn "Space" into the following:
  if ("Space" %in% cols)
    df <- df %>% 
      mutate(Space = case_when(
        Space == "" ~ "RawRTs",
        Space == "_log_space" ~ "log(RTs)",
        TRUE ~ NA_character_) %>%
          factor(levels=c("RawRTs", "log(RTs)")))
  
  df %>% as_tibble()
}


# Shows the tallies of different raised conditions
# replacement list is named list where names = patterns and vals = replacements
break_down_conditions <- function(df, col, ..., replacement_list=NULL, sep="--------") {
  col <- enquo(col)
  break_down_by <- rlang::enquos(...)
  colname <- paste0(quo_name(col))
  
  df %>% 
    select(k, type, !!col, !!!break_down_by) %>%
    filter(!!col != "") %>%
    mutate(z = seq_along(k)) %>%
    tidyr::separate_rows(!!col, sep=sep) %>%
    group_by(z) %>%
    distinct(!!col, .keep_all=T) %>%
    ungroup() %>%
    {
      temp <- .
      l <- replacement_list
      patterns <- names(l)
      if (length(l) > 0)
        for (i in 1:length(l))
          temp <- mutate(temp, !!colname := gsub(patterns[[i]], l[[i]], !!col))
      temp
    } %>%
    group_by(!!col, type, !!!break_down_by) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    arrange(-n)
}
#------------------------------------------------------------

# For making pretty tables ------------------------------------------------
# For making the pretty tables
z_colorer <- function(x, named_color_list) {
  beginning = "<span style=\"display: block; padding: 0 4px; border-radius: 4px; text-align: center; background-color: "
  end <- "</span>"
  paste0(beginning, named_color_list[x], "\">", x, end)
}
z_barer <- function(x, named_color_list="NA",
                    default_color="lightgreen",
                    fun = formattable::proportion) {
  beginning = "<span style=\"display: inline-block; direction: rtl; border-radius: 4px; padding-right: 2px; background-color: "
  colors <- ifelse(is.na(named_color_list[x]),
                   default_color, named_color_list[x])
  widths <- paste0("width: ", scales::percent(fun(as.numeric(x))))
  end <- "</span>"
  paste0(beginning,
         colors, "; ",
         widths, "\">", x, end)
}
# --------------------------------------------------------------


# Plotters --------------------------------------------------------------
add_results_ornaments <- function(p, top_margin=20) {
  plot_grob <- ggplotGrob(
    p + theme(plot.margin = margin(t=top_margin, 5.5, 5.5, 5.5))
  )
  plot_grob %>%
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 5, r = 7, clip="on")  %>%
    gtable::gtable_add_grob(textGrob("pre-residualization", 0.5, 0.75),
                            t = 1, b = 7, l = 5, r = 7, clip="on") %>%
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 9, r = 11, clip="on")  %>%
    gtable::gtable_add_grob(textGrob("post-residualization", 0.5, 0.75),
                            t = 1, b = 7, l = 9, r = 11, clip="on")  %>%
                            {
                              cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(.))
                            }
}

# For having a 1x4 facet wrap where the top is shared
ornaments_split_1x4 <- function(p, top_margin=20, texts=None) {
  plot_grob <- ggplotGrob(
    p + theme(plot.margin = margin(t=top_margin, 5.5, 5.5, 5.5))
  )
  plot_grob %>%
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 5, r = 9, clip="on")  %>%
    gtable::gtable_add_grob(textGrob(texts[[1]], 0.5, 0.75),
                            t = 1, b = 7, l = 5, r = 9, clip="on") %>%
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 13, r = 17, clip="on")  %>%
    gtable::gtable_add_grob(textGrob(texts[[2]], 0.5, 0.75),
                            t = 1, b = 7, l = 13, r = 17, clip="on")  %>%
                            {
                              cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(.))
                            }
}
# for having a 1x6 grouped into twos:
ornaments_split_1x6 <- function(p, top_margin=20, texts=None) {
  plot_grob <- ggplotGrob(
    p + theme(plot.margin = margin(t=top_margin, 5.5, 5.5, 5.5))
  )
  plot_grob %>%
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 5, r = 9, clip="on")  %>%
    gtable::gtable_add_grob(textGrob(texts[[1]], 0.5, 0.75),
                            t = 1, b = 7, l = 5, r = 9, clip="on") %>%
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 13, r = 17, clip="on")  %>%
    gtable::gtable_add_grob(textGrob(texts[[2]], 0.5, 0.75),
                            t = 1, b = 7, l = 13, r = 17, clip="on")  %>%
    
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 21, r = 25, clip="on")  %>%
    gtable::gtable_add_grob(textGrob(texts[[3]], 0.5, 0.75),
                            t = 1, b = 7, l = 21, r = 25, clip="on")  %>%
                            {
                              cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(.))
                            }
}
# twiddled to get looking good. HACK
spec_fig_ornaments_1x6 <- function(p, top_margin=40, texts=c(intdeptext_full, intgoodtext_full,"scale-independent")) {
  plot_grob <- ggplotGrob(
    p + theme(plot.margin = margin(t=top_margin, 5.5, 5.5, 5.5))
  )
  plot_grob %>%
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 5, r = 9, clip="on")  %>%
    gtable::gtable_add_grob(textGrob(texts[[1]], 0.5, 0.65),
                            t = 1, b = 7, l = 5, r = 9, clip="on") %>%
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 13, r = 17, clip="on")  %>%
    gtable::gtable_add_grob(textGrob(texts[[2]], 0.5, 0.65),
                            t = 1, b = 7, l = 13, r = 17, clip="on")  %>%
    
    gtable::gtable_add_grob(segmentsGrob(0, 0.99, 0.99, 0.99),
                            t = 3, b = 8, l = 21, r = 25, clip="on")  %>%
    gtable::gtable_add_grob(textGrob(texts[[3]], 0.5, 0.60),
                            t = 1, b = 7, l = 21, r = 25, clip="on")  %>%
                            {
                              cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(.))
                            }
}

add_top_facet_title <- function(p, text, n_pan=2, top_margin=20) {
  plot_grob <- ggplotGrob(
    p + theme(plot.margin = margin(t=top_margin, 5.5, 5.5, 5.5))
  )
  r_param <- (n_pan-1)*2+5
  plot_grob %>%
    gtable::gtable_add_grob(textGrob(text, 0.5, 0.75),
                            t = 1, b = 7, l = 5, r = r_param, clip="on") %>%
    {
      cowplot::ggdraw() +
        cowplot::draw_grob(grid::grobTree(.))
    }
}
add_mean_rt_top <- function(p, ...) {
  add_top_facet_title(p=p, text="Mean RT", n_pan=3, ...)
}
add_distro_top <- function(p, ...) {
  add_top_facet_title(p=p, text="Distribution", n_pan=2, ...)
}
add_direction_top <- function(p, ...) {
  add_top_facet_title(p=p, text="Interaction direction", n_pan=2, ...)
}


get_gg_range <- function(g) {
  x_range <- layer_scales(g)$x$range
  y_range <- layer_scales(g)$y$range
  
  x_discrete <- "RangeDiscrete" %in% class(x_range)
  y_discrete <- "RangeDiscrete" %in% class(y_range)
  
  list(x=x_range$range, y=y_range$range, x_discrete=x_discrete, y_discrete=y_discrete)
}
share_ranges <- function(x, ...) {
  l <- list(...)
  if (length(l) == 0) l <- x
  else l <- append(list(x), l)
  
  if (some(l, ~.$x_discrete==TRUE)) {
    xmin <- xmax <- NA
    warning("One of x ranges is discrete, defaulting to NA")
  } else {
    xmin <- reduce(l, ~min(.x, .y$x[[1]]), .init = Inf)
    xmax <- reduce(l, ~max(.x, .y$x[[2]]), .init = -Inf)
  }
  
  if (some(l, ~.$y_discrete==TRUE)) {
    ymin <- ymax <- NA
    warning("One of y ranges is discrete, defaulting to NA")
  } else {
    ymin <- reduce(l, ~min(.x, .y$y[[1]]), .init = Inf)
    ymax <- reduce(l, ~max(.x, .y$y[[2]]), .init = -Inf)
  }
  list(x=c(xmin,xmax), y=c(ymin,ymax))
}
get_shared_range <- function(x, ...) {
  l <- list(...)
  if (length(l) == 0) l <- x
  else l <- append(list(x), l)
  
  map(l, get_gg_range) %>%
    share_ranges()
}




plot_singular <- function(raw_data, ..., shape_col, facet_col,
                          shape_name, facet_name,
                          data_is_prepped = FALSE) {
  
  shape_col <- enquo(shape_col)
  facet_col <- enquo(facet_col)
  grouping_vars <- rlang::enquos(...)
  
  if (data_is_prepped == FALSE) {
    raw_data <- raw_data %>% 
      exclude_nonconverged() %>%
      tidyr::separate(type, c("Analysis","Type"))
    
    raw_data <- raw_data %>%
      mutate(!!shape_col := as.factor(!!shape_col),
             !!facet_col := as.factor(!!facet_col))
  }
  
  raw_data %>%
    group_by(!!!grouping_vars, Analysis, Type) %>%
    summarise(total_n=n(),
              n=sum(ifelse(grepl("singular[^A-Za-z]{0,1} fit", messages), 1, 0))) %>%
    ungroup() %>%
    {cbind(., as_tibble(Hmisc::binconf(.$n, .$total_n)))} %>%
    mutate(Nsubj = as.factor(Nsubj),
           Nitems = as.factor(Nitems)) %>%
    group_by(Nsubj) %>%
    mutate(width = n_distinct(!!!grouping_vars, Analysis, Type)) %>% 
    ungroup() %>%
    mutate(Analysis = analysis_labels[Analysis] %>% set_names(NULL)) %>%
    mutate(Analysis = factor(Analysis, levels=set_names(analysis_labels, NULL))) %>%
    ggplot(aes(x = Nsubj,
               color = Nitems,
               shape = !!shape_col,
               alpha = Type,
               y = PointEst)) +
    geom_pointrange(aes(ymax=Upper, ymin=Lower), position = position_dodge(0.5)) +
    facet_grid(vars(Analysis), vars(!!facet_col), labeller = labeller(.rows = label_both, .cols = label_value)) +
    # ggtitle("Singular fits") +
    ylab(singfit_y_text) +
    xlab(subject_text) +
    scale_alpha_manual("", values = type_alphas, labels = type_labels) +
    scale_shape(shape_name) +
    scale_color_discrete(item_text)
}

plot_convergence <- function(raw_data, ..., shape_col, facet_col,
                             shape_name=NA, facet_name=NA,
                             data_is_prepped = FALSE,
                             show_type1 = TRUE) {
  
  shape_col <- enquo(shape_col)
  facet_col <- enquo(facet_col)
  if (is.na(shape_name)) shape_name <- quo_name(shape_col)
  if (is.na(facet_name)) facet_name <- quo_name(facet_col)
  grouping_vars <- rlang::enquos(...)
  
  if (data_is_prepped == TRUE)
    stop("Data should not be prepped for convergence plots")
  
  d <- raw_data %>%
    mutate(!!shape_col := as.factor(!!shape_col),
           !!facet_col := as.factor(!!facet_col)) %>%
    tidyr::separate(type, c("Analysis","Type")) %>%
    mutate(group_id = paste(!!!grouping_vars, Analysis, Type)) %>%
    group_by(!!!grouping_vars, Analysis, Type, group_id) %>% mutate(total_n = n()) %>%
    exclude_nonconverged() %>%
    summarise(n = n(), total_n=first(total_n)) %>% ungroup() %>%
    {cbind(., as_tibble(Hmisc::binconf((.$total_n-.$n), .$total_n)))} %>%
    mutate(Nsubj = as.factor(Nsubj),
           Nitems = as.factor(Nitems)) %>%
    group_by(Nsubj) %>%
    mutate(width = n_distinct(group_id)) %>% ungroup() %>%
    mutate(Analysis = analysis_labels[Analysis] %>% set_names(NULL)) %>%
    mutate(Analysis = factor(Analysis, levels=set_names(analysis_labels, NULL)))
  
  if (show_type1)
    p <- d %>% 
    ggplot(aes(x = Nsubj,       color = Nitems,   shape = !!shape_col, 
               width = 0.1 * width, y = PointEst, alpha = Type)) +
    scale_alpha_manual("", values = type_alphas, labels = type_labels)
  else
    p <- d %>% filter(Type=="power") %>% 
    ggplot(aes(x = Nsubj,       color = Nitems,   shape = !!shape_col, 
               width = 0.1 * width, y = PointEst))
  
  p +
    geom_pointrange(aes(ymax=Upper,ymin=Lower), position = position_dodge(0.5)) +
    facet_grid(vars(Analysis), vars(!!facet_col), 
               labeller = labeller(.rows = label_both, .cols = label_value)) +
    # ggtitle("Failure of convergence rate") +
    ylab(convergence_y_text) +
    xlab(subject_text) +
    scale_shape(shape_name) +
    scale_color_discrete(item_text)
}

plot_type1 <- function(raw_data, ..., shape_col, facet_col,
                       shape_name, facet_name,
                       data_is_prepped = FALSE) {
  shape_col <- enquo(shape_col)
  facet_col <- enquo(facet_col)
  grouping_vars <- rlang::enquos(...)
  group_without_effectsize <- grouping_vars %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  
  if (data_is_prepped == FALSE) {
    raw_data <- raw_data %>%
      exclude_nonconverged() %>%
      tidyr::separate(type, c("Analysis","Type"))
    
    raw_data <- raw_data %>%
      mutate(!!shape_col := as.factor(!!shape_col),
             !!facet_col := as.factor(!!facet_col))
  }
  
  raw_data %>%
    filter(Type == "type1") %>%
    group_by(!!!group_without_effectsize) %>%
    filter(EffectSize == first(EffectSize)) %>%
    ungroup() %>%
    mutate(Nsubj = factor(Nsubj)) %>%
    mutate(issignificant = ifelse(p.value < 0.05, 1, 0)) %>%
    group_by(Nsubj) %>% mutate(width=0.1*n_distinct(!!!group_without_effectsize)) %>% ungroup() %>%
    mutate(Analysis = analysis_labels[Analysis] %>% set_names(NULL)) %>%
    mutate(Analysis = factor(Analysis, levels=set_names(analysis_labels, NULL))) %>%
    ggplot(aes(x = Nsubj,
               color = as.factor(Nitems),
               width = width,
               # shape = !!shape_col,
               y = issignificant)) +
    zplyr::stat_errorbar(position = position_dodge(0.5)) +
    facet_grid(vars(Analysis), vars(!!facet_col),
               labeller = labeller(.rows = label_both, .cols = label_value)) +
    geom_hline(yintercept=0.05, linetype="dashed") +
    # ggtitle("Type I error rate") +
    ylab("Type I error rate") +
    xlab(subject_text) +
    # scale_shape(shape_name) +
    scale_color_discrete(item_text)
}

# Probably should turn `facet_col` back to `shape_col` and let the original `facet_col` be flexible again
plot_power <- function(data, ..., shape_col, facet_col,
                       shape_name=NA, facet_name=NA,
                       data_is_prepped = FALSE) {
  facet_col <- enquo(facet_col)
  shape_col <- enquo(shape_col)
  if (is.na(shape_name)) shape_name <- NULL
  if (is.na(facet_name)) facet_name <- NULL
  grouping_vars <- rlang::enquos(...)
  
  if (data_is_prepped == FALSE) processed_data <- summarise_mixed_dfs(data, !!!grouping_vars)
  else processed_data <- data
  
  p.power_corrected <- processed_data %>%
    filter(Nitems==Nsubj) %>%
    filter(Type=="power") %>%
    # mutate(EffectSize = ifelse(EffectSize %in% c(15, 56), "small 15/56", "large 35/80")) %>%
    # mutate(EffectSize = factor(EffectSize, levels=c("small 15/56", "large 35/80"))) %>%
    ggplot(aes(x = Nsubj,
               y = CorrectedStat,
               color = Analysis)) +
    geom_line(size=2) +
    facet_grid(vars(!!shape_col), vars(!!facet_col)) +
    scale_y_continuous(corrected_power_text,
                       sec.axis = dup_axis(name=shape_name, 
                                           breaks=c(-Inf, Inf), labels = NULL)) +
    scale_x_continuous(subject_text,
                       sec.axis = dup_axis(name=facet_name, 
                                           breaks=c(-Inf, Inf), labels = NULL)) +
    scale_color_manual(values = analysis_colors,
                       labels = analysis_labels)
  
  p.power <- processed_data %>%
    filter(Nitems==Nsubj) %>%
    filter(Type=="power") %>%
    # mutate(EffectSize = ifelse(EffectSize %in% c(15, 56), "small 15/56", "large 35/80")) %>%
    # mutate(EffectSize = factor(EffectSize, levels=c("small 15/56", "large 35/80"))) %>%
    ggplot(aes(x = Nsubj,
               y = Stat,
               color = Analysis)) +
    geom_line(size=2) +
    facet_grid(vars(!!shape_col), vars(!!facet_col)) +
    scale_y_continuous(power_text,
                       sec.axis = dup_axis(name=shape_name, 
                                           breaks=c(-Inf, Inf), labels = NULL)) +
    scale_x_continuous(subject_text,
                       sec.axis = dup_axis(name=facet_name, 
                                           breaks=c(-Inf, Inf), labels = NULL)) +
    scale_color_manual(values = analysis_colors,
                       labels = analysis_labels)
  
  return(list(p.power = p.power, p.power_corrected = p.power_corrected))
  
}


data_report <- function(raw_data, ..., shape_col, facet_col,
                        shape_name, facet_name,
                        preprocessed_data = NULL) {
  
  shape_col <- enquo(shape_col)
  facet_col <- enquo(facet_col)
  grouping_vars <- rlang::enquos(...)
  
  group_without_effectsize = grouping_vars %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  
  if (is.null(preprocessed_data)) processed_data <- summarise_mixed_dfs(raw_data, !!!grouping_vars)
  else processed_data <- preprocessed_data
  
  kept_raw_data <- raw_data
  
  raw_data <- raw_data %>%
    exclude_nonconverged() %>%
    tidyr::separate(type, c("Analysis","Type"))
  
  # For convenience:
  processed_data <- processed_data %>%
    mutate(!!shape_col := as.factor(!!shape_col),
           !!facet_col := as.factor(!!facet_col))
  raw_data <- raw_data %>% 
    mutate(!!shape_col := as.factor(!!shape_col),
           !!facet_col := as.factor(!!facet_col))
  # Needs to be below here
  type1_diffs <- type1_checker(processed_data, !!!grouping_vars)
  
  p.all_returned <- kept_raw_data %>%
    group_by(!!!grouping_vars, type) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    mutate(Nsubj = as.factor(Nsubj),
           Nitems = as.factor(Nitems),
           perc_tot = format(n/max(n), digits=3),
           i = seq_along(Nsubj)) %>%
    ggplot(aes(x = Nsubj,
               fill = Nitems,
               color = !!shape_col,
               label = perc_tot,
               y = n)) +
    geom_bar(stat="identity", position = "dodge") +
    geom_text(aes(y=n/2, group=i), position=position_dodge(1),
              color="black") +
    facet_grid(vars(type), vars(!!facet_col),
               labeller = labeller(.rows = label_both, .cols = label_value)) +
    ggtitle("# of models") +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank())
  
  p.convergence <- kept_raw_data %>%
    plot_convergence(!!!grouping_vars, 
                     shape_col = !!shape_col, shape_name = shape_name, 
                     facet_col = !!facet_col, facet_name = facet_name,
                     data_is_prepped = FALSE) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank()) +
    smaller_legends_theme
  
  p.sing_fit <- raw_data %>%
    plot_singular(!!!grouping_vars,  
                  shape_col = !!shape_col, shape_name = shape_name, 
                  facet_col = !!facet_col, facet_name = facet_name,
                  data_is_prepped = TRUE) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank()) +
    smaller_legends_theme
  
  p.typeI <- raw_data %>% 
    plot_type1(!!!grouping_vars,  
               shape_col = !!shape_col, shape_name = shape_name, 
               facet_col = !!facet_col, facet_name = facet_name,
               data_is_prepped = TRUE) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank())
  
  p.typeI.comparison <- processed_data %>%
    filter(Type=="type1") %>%
    zplyr::left_join(type1_diffs) %>%
    mutate(logit = qlogis(Stat)) %>%
    rowwise() %>%
    mutate(differs_from_05 = differs_from_05[[Analysis]]) %>%
    ungroup() %>%
    ggplot(aes(x=Type1sDiffer.p.v, y=logit, shape = differs_from_05,
               color = LowerType1,
               group=paste0(!!!grouping_vars))) +
    geom_point(size=3) +
    geom_line() +
    geom_hline(yintercept = qlogis(0.05), linetype="dashed") +
    geom_vline(xintercept = 0.05) +
    xlab("p-value of whether linear vs. log differ") +
    scale_shape_discrete("sig. diff from .05") +
    scale_color_discrete("lower Type I") +
    ylab("Type I error rate (log-odds)") +
    ggtitle("Comparing Type I error rates between analyses\n and whether each significantly differs from 0.05")
  
  power_graphs <- plot_power(processed_data, !!!grouping_vars, 
                             facet_col = !!facet_col, shape_col = !!shape_col,
                             facet_name = facet_name, shape_name=shape_name,
                             data_is_prepped = TRUE)
  
  p.power_corrected <- power_graphs$p.power_corrected 
  p.power           <- power_graphs$p.power
  
  list(
    p.all_returned = p.all_returned,
    p.convergence = p.convergence,
    p.sing_fit = p.sing_fit,
    p.typeI = p.typeI,
    p.typeI.comparison = p.typeI.comparison,
    p.power = p.power,
    p.power_corrected = p.power_corrected
  )
}
# --------------------------------------------------------------

# New model/plotting functions ---------------------------------

# Takes the potentially raw data and makes it ready for plotting type 1 data
build_type1_plot_data <- function(.data, grouping_quos, data_is_prepped = FALSE) {
  if (data_is_prepped == FALSE) 
    .data <- .data %>%
      exclude_nonconverged() %>%
      tidyr::separate(type, c("Analysis","Type"))
  
  group_without_effectsize <- grouping_quos %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  
  .data %>%
    filter(Type == "type1") %>%
    group_by(!!!group_without_effectsize) %>%
    filter(EffectSize == first(EffectSize)) %>%
    ungroup() %>%
    mutate(Nsubj = factor(Nsubj), Nitems = as.factor(Nitems)) %>%
    mutate(issignificant = ifelse(p.value < 0.05, 1, 0)) %>%
    group_by(!!!grouping_quos, Analysis, Type) %>%
    summarise(n=n(), are_sig=sum(issignificant)) %>% ungroup() %>%
    mutate(Stat=are_sig/n)
}

# Prepare summarized data
prepare_summed_type1_data <- function(.data, grouping_quos) {
  group_without_effectsize <- grouping_quos %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  .data %>%
    filter(Type == "type1") %>%
    group_by(!!!group_without_effectsize) %>%
    filter(EffectSize == first(EffectSize)) %>%
    ungroup() %>%
    mutate(Analysis = analysis_labels[Analysis] %>% set_names(NULL)) %>%
    mutate(Analysis = factor(Analysis, levels=set_names(analysis_labels, NULL))) %>%
    mutate(y=Stat*n) %>%
    {cbind(., as_tibble(Hmisc::binconf(.$y, .$n)))} %>%
    mutate(Nsubj = as.factor(Nsubj), Nitems = as.factor(Nitems)) %>%
    mutate(sigDiffFrom05 = factor(as.character(Upper < 0.05 | Lower > 0.05),
                                     levels=c("FALSE","TRUE")))
}

# Plots base plot
base_plot_type1 <- function(.data, grouping_quos) {
  prepare_summed_type1_data(.data, grouping_quos) %>% 
    ggplot(aes(x = Nsubj, y = PointEst)) +
    # geom_line() +
    geom_pointrange(aes(ymax=Upper,ymin=Lower), 
                    position = position_dodge(def_dodge_width)) +
    geom_hline(yintercept=0.05, linetype="dashed") +
    ylab("Type I error rate") +
    xlab(subject_text)
}

# p <- summarized_reg_data %>%
#   filter(Space=="", Nsubj==Nitems) %>%
#   base_plot_type1(giant_grouping_vars) +
#   aes(color=Analysis) + 
#   facet_grid(vars(NULL), vars(SampleName),
#              labeller = label_value) +
#   scale_color_manual(values = analysis_colors,
#                      labels = analysis_labels)

# Takes the potentially raw data and makes it ready for plotting singfit data
build_singfit_plot_data <- function(.data, grouping_quos, data_is_prepped = FALSE) {
  if (data_is_prepped == FALSE) 
    .data <- .data %>%
      exclude_nonconverged() %>%
      tidyr::separate(type, c("Analysis","Type"))
  
  group_without_effectsize <- grouping_quos %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  
  .data %>%
    mutate(Nsubj = factor(Nsubj), Nitems = as.factor(Nitems)) %>%
    mutate(issignificant = ifelse(p.value < 0.05, 1, 0)) %>%
    group_by(!!!grouping_quos, Analysis, Type) %>%
    summarise(n=n(), 
              n_singfits=sum(ifelse(grepl("singular[^A-Za-z]{0,1} fit", messages), 1, 0))) %>%
    ungroup()
}

# Prepare the summarized data for singular fit plots and models
prepare_summed_singfit_data <- function(.data, grouping_quos) { 
  group_without_effectsize <- grouping_quos %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  .data %>%
    filter(!(EffectSizeName=="large" & Type=="type1")) %>%
    mutate(EffectSizeName = ifelse(Type=="type1", "null", 
                                   as.character(EffectSizeName)) %>%
             factor(levels=c("null", "medium", "large"))) %>%
    mutate(Analysis = analysis_labels[Analysis] %>% set_names(NULL)) %>%
    mutate(Analysis = factor(Analysis, levels=set_names(analysis_labels, NULL))) %>%
    mutate(y=n_singfits) %>%
    {cbind(., as_tibble(Hmisc::binconf(.$n_singfits, .$n)))} %>%
    mutate(Nsubj = as.factor(Nsubj), Nitems = as.factor(Nitems))
}

# Plots base plot
base_plot_singfit <- function(.data, grouping_quos) {
  .data %>%
    prepare_summed_singfit_data(grouping_quos) %>%
    ggplot(aes(x = Nsubj, y = PointEst)) +
    geom_pointrange(aes(ymax=Upper,ymin=Lower), 
                    position = position_dodge(0.5)) +
    ylab(singfit_y_text) +
    xlab(subject_text)
}

# p <- summarized_fetal_data %>%
#   filter(Space=="", Nsubj==Nitems) %>%
#   base_plot_singfit(giant_grouping_vars) +
#   aes(color=Nitems, shape=EffectSizeName, alpha=Type) + 
#   facet_grid(vars(Analysis), vars(SampleName),
#              labeller = label_value) +
#   scale_y_continuous(sec.axis = dup_axis(name="Analysis approach", 
#                                          breaks=c(-Inf, Inf), labels=NULL)) +
#   scale_shape("Effect size") +
#   scale_alpha_manual("", values = type_alphas, labels = type_labels) +
#   scale_color_discrete(guide=FALSE)

# Takes the potentially raw data and makes it ready for plotting convergence data
build_conv_plot_data <- function(.data, grouping_quos) {
  if ("type" %in% names(.data))
    .data <- .data %>%
      tidyr::separate(type, c("Analysis","Type"))
  
  group_without_effectsize <- grouping_quos %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  
  .data %>%
    group_by(!!!grouping_quos, Analysis, Type) %>%
    mutate(og_n=n()) %>%
    exclude_nonconverged() %>%
    group_by(!!!grouping_quos, Analysis, Type) %>%
    summarise(n=n()) %>% ungroup()
}

prepare_summed_conv_data <- function(.data, grouping_quos) {
  # Gets the grouping quos without the effect size
  group_without_effectsize <- grouping_quos %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  
  .data %>%
    # removes the LARGE effect size TYPE1 data  (don't need type 1 twice)
    filter(!(EffectSizeName=="large" & Type=="type1")) %>%
    # changes effect size names to "null", "medium", and "large"
    mutate(EffectSizeName = ifelse(Type=="type1", "null", 
                                   as.character(EffectSizeName)) %>%
             factor(levels=c("null", "medium", "large"))) %>%
    # changes the analsis label names
    mutate(Analysis = analysis_labels[Analysis] %>% set_names(NULL)) %>%
    mutate(Analysis = factor(Analysis, levels=set_names(analysis_labels, NULL))) %>%
    mutate(y=og_n-n, n=og_n) %>%
    
    {cbind(., as_tibble(Hmisc::binconf(.$y, .$n)))} %>%
    mutate(Nsubj = as.factor(Nsubj), Nitems = as.factor(Nitems))
}
# ZACHEDIT 2022
# for the data with the int stuff
prepare_summed_conv_data_int <- 



base_plot_convergence <- function(.data, grouping_quos, show_type1 = TRUE) {
  .data <- .data %>% prepare_summed_conv_data(grouping_quos)
  
  if (show_type1==FALSE)
    .data <- .data %>% filter(Type=="power")
  
  .data %>%
    ggplot(aes(x = Nsubj, y = PointEst)) +
    geom_pointrange(aes(ymax=Upper,ymin=Lower), 
                    position = position_dodge(def_dodge_width)) +
    ylab(convergence_y_text_si) +
    xlab(subject_text)
}

# Prepare the summarized power data for modeling and plots
prepare_summed_power_data <- function(.data, grouping_quos, 
                                      is_corrected=FALSE,
                                      is_type1=FALSE
                                      ) {
  group_without_effectsize <- grouping_quos %>% 
    purrr::keep(~ quo_name(.) != "EffectSize" & quo_name(.) !=  "EffectSizeName")
  
  .data <- .data %>%
    mutate(Analysis = analysis_labels[Analysis] %>% set_names(NULL)) %>%
    mutate(Analysis = Analysis %>%
             factor(levels=set_names(analysis_labels, NULL))) %>%
    filter(Nitems==Nsubj)
  
  power_or_type1 <- "power"
  
  if (is_corrected) {
    .data <- .data %>% mutate(y=round(CorrectedStat*n))
    if (is_type1) stop("No such thing as 'corrected' Type I data")
  } else {
    .data <- .data %>% mutate(y=round(Stat*n))
    if (is_type1) power_or_type1 <- "type1"
  }
  
  .data %>%
    filter(Type==power_or_type1) %>%
    {cbind(., as_tibble(Hmisc::binconf(.$y, .$n)))}
}

base_plot_power <- function(.data, grouping_quos, is_corrected = FALSE) {
  if (is_corrected) y_text <- corrected_power_text
  else y_text <- power_text
  .data %>% 
    prepare_summed_power_data(grouping_quos, is_corrected) %>%
    mutate(Nsubj=as.factor(Nsubj)) %>%
    ggplot(aes(x = Nsubj, y = PointEst)) +
    geom_line() +
    geom_pointrange(aes(ymax=Upper, ymin=Lower)) +
    # geom_linerange(aes(ymax=Upper, ymin=Lower)) +
    # geom_ribbon()
    ylab(y_text) + xlab(subject_text)
}

base_plot_advantage <- function(.data, grouping_quos, 
                                is_corrected = FALSE, 
                                in_log_odds = TRUE,
                                is_type1 = FALSE) {
  y_var <- quo(`log-odds difference`)
  if (is_corrected) {
    if (in_log_odds) y_text <- corrected_log_odds_diff_text
    else {
      y_text <- corrected_odds_ratio_text
      y_var <- quo(`odds ratio`)
    }
  } else {
    if (in_log_odds) y_text <- odds_ratio_text
    else {
      y_text <- log_odds_ratio_text 
      y_var <- quo(`odds ratio`)
    }
  }
  
  if (is_type1) {
    # Ugh just deal with it. HACK
    y_text <- "XXXXXXXXXXXXXXX"
    group_without_effectsize <- grouping_quos %>% 
      purrr::keep(~ quo_name(.) != "EffectSize" & 
                    quo_name(.) != "EffectSizeName")
    .data <- .data %>% 
      group_by(!!!group_without_effectsize, Analysis) %>%
      filter(EffectSize==first(EffectSize)) %>%
      mutate(EffectSizeName="null") %>%
      ungroup()
  }
  
  .data <- .data %>% 
    prepare_summed_power_data(grouping_quos = grouping_quos, 
                              is_corrected = is_corrected,
                              is_type1 = is_type1) %>%
    mutate(Nsubj=as.factor(Nsubj)) %>%
    group_by(!!!grouping_quos, Analysis) %>%
    filter(!any(PointEst %in% c(1,0))) %>%
    ungroup() %>%
    select(!!!grouping_quos, Analysis, PointEst, Lower, Upper) %>%
    mutate(Odds = (PointEst)/(1-PointEst)) %>%
    select(-PointEst, -Lower, -Upper)
  
  if (all(unique(.data$Analysis) %in% analysis_labels)) {
    v <- names(analysis_labels) %>% set_names(analysis_labels)
    .data$Analysis <- v[.data$Analysis] %>% set_names(NULL)
  }
  
  if ("lgshift" %in% .data$Analysis) {
    .data <- 
      bind_rows(
        "lgshift" = .data %>% 
          filter(Analysis!="log") %>% 
          mutate(Analysis=ifelse(Analysis=="linear","linear","log")),
        "log" = .data %>% filter(Analysis!="lgshift"),
        .id="log_type"
      )
  }
  .data %>% 
    tidyr::pivot_wider(names_from = Analysis,
                       values_from = Odds) %>%
    mutate(`odds ratio` = log/linear,
           # The additional minus is so that it makes it 'advantage' 
           #  (ie small = good)
           `log-odds difference` = log(log) - log(linear)) %>%
           {
             if (is_type1) {
               mutate(., `odds ratio` = 1/`odds ratio`,
                      `log-odds difference` = -1 * `log-odds difference`)
             } else {
               .
             }
           }%>%
    ggplot(aes(x = Nsubj, y = !!y_var)) +
    geom_line() +
    xlab(subject_text) +
    { if (in_log_odds == FALSE)
      list(
        geom_hline(yintercept=1, linetype="dashed"),
        scale_y_log10(y_text,
                      # breaks=c(1, 1.5, 2, 4),
                      limits =c(0.75, NA))
      )
      else
        geom_hline(yintercept=0, linetype="dashed")
    }
}



# summarized_fetal_data %>%
#   base_plot_power(giant_grouping_vars) +
#   aes(color=Analysis) +
#   facet_grid(vars(EffectSizeName), vars(SampleName),
#              labeller = label_value) +
#   scale_y_continuous(sec.axis = dup_axis(
#     name="Effect size", 
#     breaks=c(-Inf, Inf), labels = NULL)) +
#   scale_color_manual(values = analysis_colors,
#                      labels = analysis_labels)

# --------------------------------------------------------------

# Anova functions ----------------------------------------------
make_formulae <- function(e, ...) {
  v <- enexprs(...)
  e <- enexpr(e)
  full_thing = reduce(v, ~expr(!!.x*!!.y))
  
  if (length(v) > 1) 
    formulae <- (length(v)-1):1 %>%
    map(
      function(n) {
        combn(v, n, simplify = FALSE) %>%
          map(~reduce(.x, function(x,y) expr(!!x*!!y))) %>%
          reduce(~expr(!!.x+!!.y*!!e), .init=expr(1+!!full_thing))
      })
  else
    formulae <- list()
  
  list(list(expr(1+!!full_thing*!!e)), 
       formulae, 
       list(expr(1+!!full_thing+!!e)),
       list(expr(1+!!full_thing))) %>%
    flatten() %>%
    map(~new_formula(expr(cbind(y, n-y)), .x))
}


contrast_namer <- function(fac, contr_fn=contr.sum) {
  row_names <- contrasts(fac) %>% 
    row.names()
  contr <- contr_fn(n_distinct(fac))
  attr(contr, "dimnames") <- list(row_names, 
                                  row_names[1:(length(row_names)-1)]) 
  fac %>% `contrasts<-`(value=contr)
}


good_sigfiger <- function(x,digits=3) {
  formatC(signif(x, digits=digits), digits=digits,format="fg", flag="#")
}

full_regression_summarizer <- function(m, as_char=TRUE) {
  l <- c(
    Nsubj = "SampSize-",
    SampleName = "DataPrep-",
    EffectSizeName = "Effect-",
    EffectSize = "Effect-", 
    Type = "Effects-",
    Direction = "Dir-",
    Space = "Space-",
    Analysis = "Anlys-",
    `RT ~ ...` = "stndrd",
    `mean RT` = "MeanRT-",
    `log(RT) ~ ...` = "logged",
    unproblematic = "unprbl",
    `scale-dependent` = "scaledep"
    )
  res <- broom::tidy(m)
  res$term <- reduce2(l, names(l),
          ~gsub(..3, ..2, ..1),
          .init=res$term)
  res %>%
  {
    if (as_char==FALSE) {
      mutate_at(., vars(estimate:statistic), ~signif(.x, digits=3))
    } else {
      mutate_at(., vars(estimate:statistic), good_sigfiger)
    }
  } %>%
    mutate(` ` = map_chr(p.value, ~get_asterisks(
      .x, mapping = c("0.05"="*", "0.01"="**", "0.001"="***", "")))) %>%
    mutate(p.value = format.pval(p.value, digits=2))
}


auto_anova <- function(df, critvar, ...,
                       return_everything = FALSE,
                       plus_one_smoothing = TRUE) {
  
  critvar <- enexpr(critvar)
  critname1 <- paste0(quo_name(critvar),"1")
  preds <- enexprs(...)
  lm_df <- df %>%
    select(!!!preds, !!critvar, y, n) %>%
    mutate_at(vars(!!!preds, !!critvar), factor) %>%
    mutate_at(vars(!!!preds, !!critvar), 
              ~`contrasts<-`(.x, value=contr.sum(n_distinct(.x))))
  
  if (plus_one_smoothing)
    lm_df <- mutate(lm_df, y=y+1, n=n+2)
  
  pretty_lm_df <- lm_df %>%
    mutate_at(vars(!!!preds, !!critvar), contrast_namer)
  
  formulae <- make_formulae(!!critvar, !!!preds)
  models <- map(formulae, ~glm(.x, data=lm_df, family = "binomial"))
  full_m <- models[[1]]
  models_lagged <- models[2:length(models)]
  models_lead <- models[1:(length(models)-1)]
  models <- models[2:length(models)]
  chi_sqs <- map2(models_lead, models_lagged, ~(anova(.x, .y, test="Chisq")))
  omni_tests <- imap_dfr(chi_sqs, function(x,i) {
    as.data.frame(x)[2,] %>%
      transmute(i = i,
                df = abs(Df),
                deviance = abs(Deviance),
                p = `Pr(>Chi)`,
                sig = map_chr(p, get_asterisks))
  })
  pretty_regsum <- glm(formulae[[1]], data=pretty_lm_df, family="binomial") %>%
    full_regression_summarizer()
    
  
  p_vals <- map_dbl(chi_sqs, ~.$`Pr(>Chi)`[[2]])
  
  if (return_everything==FALSE) {

  intercept <- coefficients(full_m)[["(Intercept)"]]
  diff_in_logodds <-  coefficients(full_m)[[critname1]]*2
  diff_at_05 <- plogis(qlogis(0.05)+diff_in_logodds/2) - 
    plogis(qlogis(0.05)-diff_in_logodds/2)
  diff_at_intercept <- plogis(intercept+diff_in_logodds/2) - 
    plogis(intercept-diff_in_logodds/2)
  
  main_effect_data <- summary(full_m)$coefficients %>% 
    as.data.frame() %>% 
    rownames_to_column("term") %>% 
    filter(term %in% c(critname1, "(Intercept)"))
  
  return(list(p_vals=p_vals, diff_in_logodds=diff_in_logodds, 
              diff_at_05=diff_at_05, 
              diff_at_intercept=diff_at_intercept,
              main_effect_data=main_effect_data,
              omni_tests = omni_tests))
  } else {
    return(list(
      omni_tests = omni_tests,
      full_model = full_m,
      models = models,
      p_vals=p_vals,
      data=lm_df,
      regsum=pretty_regsum
    ))
  }
}

type1_sig_rater <- function(df, preds, offset=0.05) {
  preds <- enexpr(preds)
  real_preds <- expr(1 + offset(o) + !!preds)
  formula <- new_formula(expr(cbind(y, n-y)), real_preds)
  glm(formula, family="binomial",
      data=mutate(df, o=qlogis(offset))) %>%
    broom::tidy() %>%
    filter(term=="(Intercept)")
}

format_omnibus <- function(df) {
  df %>%
    mutate(sig=ifelse(sig==" n.s.","",sig)) %>%
    mutate_at(vars(p, deviance, df),
              ~format(., digits=3, big.mark = ",")) %>%
    rename(`$\\emph{p}$`     = p,
           `$\\Delta_{dev}$` = deviance,
           `$\\Delta_{df}$` = df,
           ` ` = sig)
}

get_descriptive_stats <- function(df, ...,
                                  say_greater=FALSE) {
  grouping_quos <- rlang::enquos(...)
  `%comp%` <- if (say_greater) `>` else `<`
  comparator <- if (say_greater) "greater" else "less than" 
  
  df <- df %>%
    mutate(SA = names(analysis_labels[match(Analysis, analysis_labels)])) %>%
    group_by(!!!grouping_quos) %>%
    mutate(Stat = y/n,
           SmoothStat = (y+1)/(n+2),
           qStat=qlogis(Stat),
           qSmoothStat=qlogis(SmoothStat))
  
  number_at_edge <- df %>% 
    filter(all(Stat==1) | all(Stat==0)) %>% 
    {nrow(.)/2}
  
  number_equal <-   df %>% 
    filter(!(all(Stat==1) | all(Stat==0))) %>% 
    filter(Stat[SA=="log"] == Stat[SA=="linear"]) %>% {nrow(.)/2}
  
  df %>% 
    filter(!(all(Stat==1) | all(Stat==0))) %>% 
    summarise(log_odds_diff = qStat[SA=="log"] - qStat[SA=="linear"],
              is_less       =  Stat[SA=="log"] %comp%  Stat[SA=="linear"],
              smooth_odds_diff = qSmoothStat[SA=="log"] - qSmoothStat[SA=="linear"],
              diff = Stat[SA=="log"] - Stat[SA=="linear"]) %>%
    ungroup() %>%
    summarise(log_odds = mean(log_odds_diff),
      smooth_odds = mean(smooth_odds_diff),
      diff = mean(diff),
      y = sum(ifelse(is_less,1,0)),
      n = n(), 
      is_less_ratio = mean(ifelse(is_less,1,0))) %>%
    mutate(ceiling_floor = number_at_edge,
           equal = number_equal,
           comp = comparator) %>%
    select(comp, y, n, everything())
}

get_stat_closer_to_05_df <- function(df, ...,
                                    say_greater=FALSE,
                                    closer_val=0.05
                                     
                                      ) {
  grouping_quos <- rlang::enquos(...)
  `%comp%` <- if (say_greater) `>` else `<`
  comparator <- if (say_greater) "greater" else "less than" 
  
  df <- df %>%
    mutate(SA = names(analysis_labels[match(Analysis, analysis_labels)])) %>%
    group_by(!!!grouping_quos) %>%
    mutate(Stat = y/n,
           SmoothStat = (y+1)/(n+2),
           qStat=qlogis(Stat),
           qSmoothStat=qlogis(SmoothStat)) %>%
    filter(n_distinct(Stat) != 1)
  
  df %>% 
    summarise(closer_diff =    abs(qStat[SA=="log"] - qlogis(closer_val)) -
                               abs(qStat[SA=="linear"] - qlogis(closer_val)),
              is_closer = abs(qStat[SA=="log"] - qlogis(closer_val)) %comp%
                abs(qStat[SA=="linear"] - qlogis(closer_val))) %>%
    ungroup() %>%
    mutate(total_is_closer = sum(ifelse(is_closer==TRUE,1,0)),
           total = n(),
           comp = comparator) %>%
    select(closer_diff, is_closer, total_is_closer, total, everything())
}


output_descriptive_stats <- function(df, 
                                     intro_text,
                                     stat_text,
                                     say_greater=FALSE,
                                     at_val = NA,
                                     used_smoothing=TRUE) {
  
  temp_ceiling_text <- paste0(" (of which there were ", 
                              df$ceiling_floor, ")")
  temp_equal_text <- paste0(english::Words(df$equal), 
                            " simulations had equal ",
                            stat_text, ". ")
  if (is.na(at_val))
    at_val_text <- ". "
  else {
    diff <- abs(plogis(qlogis(at_val) + df$smooth_odds/2) - 
                  plogis(qlogis(at_val) - df$smooth_odds/2))
    diff <- format(diff, digits=2)
    at_val_text <-  paste0(", which corresponds to a difference of ",
                           diff, " when centered around ", 
                           as.character(at_val), " in proportion space. ")
  }
  if (used_smoothing) {
    log_odds_text <- paste0(
      format(abs(df$smooth_odds), digits=3),
      " in log-odds (using plus-one smoothing to avoid any infinite values)"
    )
  } else {
    log_odds_text <- paste0(
      format(abs(df$log_odds), digits=3),
      " in log-odds"
    )
  }
  
  
  comp_word <- if (say_greater) "higher " else "lower "
  equal_text <- if (df$equal > 0) temp_equal_text else ""
  ceiling_text <- if (df$ceiling_floor > 0) temp_ceiling_text else ""
  is_higher_text <- if (df$smooth_odds > 0) "higher " else "lower "
  
  with(
    df,
    {
      paste0(
        intro_text,
        "the log-transformed analysis approach had ",
        comp_word, stat_text, " than the standard approach for ", 
        y, " of the ", n, " comparisons that were not both at ceiling or floor",
        ceiling_text, ". ",
        equal_text,
        "The log-transformed analysis approach had ", is_higher_text, stat_text,
        " on average, with a mean difference of ", log_odds_text,
        at_val_text
      )
    }
  )
}
# --------------------------------------------------------------


# Fancy Plotter --------------------------------------------------------------
fancy_plots <- function(df, ..., main_plot_title,
                        grouping_quosures, 
                        col_col, col_name,
                        color_col, color_name,
                        row_col=NULL, row_name=NA,
                        is_corrected=FALSE,
                        return_components = FALSE) {
  
  main_plot_filterers <- rlang::enquos(...)
  
  col_col   <- enquo(col_col)
  color_col <- enquo(color_col)
  row_col   <- enquo(row_col)
  
  if (is.na(col_name)) col_name <- NULL
  if (is.na(row_name)) row_name <- NULL
  
  
  if (quo_is_null(row_col)) no_row <- TRUE
  else no_row <- FALSE
  
  if (is_corrected) {
    df <- df %>% mutate(Stat=CorrectedStat)
    main_y_text <- corrected_power_text
    satellite_y_text <- corrected_odds_ratio_text
  } else {
    main_y_text <- power_text
    satellite_y_text <- odds_ratio_text
  }
  
  # Line parameters
  line_width <- 1
  v_lineend <- "butt"
  v_linejoin <- "mitre"
  v_linecolor <- "grey60"
  v_linetype <- "solid" # "dotted"
  
  # Getting the line-stuff -------------------------------
  equal_items_subjects <- df %>%
    filter(Type == "power", Nsubj == Nitems)
  
  main_df1 <- equal_items_subjects %>%
    filter(!!!main_plot_filterers)
  
  end_points <- filter(main_df1, Nsubj==64, Nitems==64) %>%
    select(!!!grouping_quosures, Analysis, Stat) %>%
    tidyr::spread(Analysis, Stat)
  
  main_df1_odds_difference <- (end_points$log/(1 - end_points$log)) / (end_points$linear/(1 - end_points$linear))
  
  end_points$odds_diff <- main_df1_odds_difference
  end_points$text <- paste0("", round(main_df1_odds_difference, 1), " greater odds")
  
  # The annotation for the main plot
  main_df1_annotation <- c(
    geom_segment(aes(x=Nsubj, xend=Nsubj, 
                     yend=log, y=linear),
                 data = end_points, inherit.aes = FALSE, 
                 linetype = v_linetype,
                 size = line_width,
                 lineend = v_lineend,
                 linejoin = v_linejoin,
                 color = v_linecolor,
                 arrow = arrow(ends="both", length = unit(0.35, "cm"))
    ),
    geom_label_repel(data = end_points, inherit.aes = FALSE,
                     aes(x=Nsubj, y=mean(c(linear, log)), label=text),
                     min.segment.length = unit(0, "lines"),
                     force=10,
                     nudge_x = -60, nudge_y = 0.05))
  
  satellite_annotation <- c(
    geom_segment(aes(x=Nsubj, xend=Nsubj, 
                     yend = odds_diff, y = 1),
                 data = end_points, 
                 inherit.aes = FALSE, 
                 linetype = v_linetype,
                 size = line_width,
                 lineend = v_lineend,
                 linejoin = v_linejoin,
                 color = v_linecolor,
                 arrow = arrow(ends="last", length = unit(0.35, "cm")))
  )
  
  # the satellite plots
  main_df1_satellites <- equal_items_subjects %>% 
    select(!!!grouping_quosures, Analysis, Stat) %>%
    mutate(Odds = (Stat)/(1-Stat)) %>%
    select(-Stat) %>%
    tidyr::spread(Analysis, Odds) %>% 
    mutate(`odds ratio` = log/linear) %>%
    ggplot(aes(x=Nsubj, y = `odds ratio`, group := !!color_col, color:= !!color_col)) +
    geom_line(size=2) +
    geom_hline(yintercept=1, linetype="dashed") +
    satellite_annotation +
    xlab(subject_text) +
    ylab(satellite_y_text) +
    facet_grid(vars(!!row_col), vars(!!col_col)) +
    scale_y_log10(breaks=c(1, 1.5, 2, 4),
                  limits =c(0.75, NA),
                  sec.axis = dup_axis(name=row_name, breaks=c(-Inf,Inf), labels=NULL)) +
    scale_x_continuous(breaks = c(8,16,32,64),
                       limits = c(6,70),
                       sec.axis = dup_axis(name=col_name, breaks=c(-Inf,Inf), labels=NULL))+
    theme( panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           legend.background = element_blank(),
           legend.position = "right") +
    scale_color_discrete(color_name)
  
  # The main plot
  main_df1_facet <- main_df1 %>%
    ggplot(aes(x = Nsubj, 
               y = Stat, 
               color = Analysis)) +
    geom_line(size=2) +
    main_df1_annotation + 
    scale_x_continuous(subject_text,
                       breaks = c(8,16,32,64),
                       limits = c(6,70)) +
    ylab(main_y_text) + 
    theme(legend.position = c(.72, 0.20),
          legend.background = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank() ) +
    scale_color_manual(values = analysis_colors,
                       labels = analysis_labels)
  
  if (return_components) {
    list(main_facet = main_df1_facet, satellites=main_df1_satellites)
  } else {
    plot_grid(
      rel_widths = c(1, 1.5),
      main_df1_facet + 
        # coord_cartesian(ylim=shared_range) +
        theme(plot.margin = margin(b=5.5, l=5.5, r=5.5, t=36)),
      main_df1_satellites)
  }
}

# --------------------------------------------------------------



# Distribution plotting functions -------------------------------------------
# Constants and functions
p.distr.theme <-  theme(
  text = element_text(size=9),
  # axis.title.y = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black")
)

make_RTdensity_panel <- function(df, colname, title, ..., geom="abs_text", size=3) {
  ggplot(df, aes_string(x=colname)) +
    geom_density() +
    zplyr::stat_moments(aes(xpos=0.6, ypos=0.80), geom=geom, ...,
                        size=size,
                        sig=T, moment="both", excess_kurtosis = T) +
    scale_x_continuous(title) +
    scale_y_continuous("") +
    p.distr.theme
}
# --------------------------------------------------------------


# 'Jaeger Plot' functions --------------------------------------------------------------
# Gets appropriate asterisks of significance for p value
get_asterisks <- function(p.value, 
                          mapping = c("0.05"="*", "0.01"="**", "0.001"="***", " n.s.")) {
  stopifnot(length(p.value)==1)
  
  named <- rlang::have_name(mapping)
  named_vals <- mapping[named]
  unnamed_val <- mapping[!named]
  threshs <- as.numeric(names(named_vals))
  named_vals <- named_vals[order(threshs)]
  
  if (length(threshs) != length(mapping) - 1 || 
      length(threshs) == 0 || any(threshs >= 1) || 
      any(is.na(threshs)))
    stop("mapping variable not properly specified")
  
  purrr::reduce2(
    sort(threshs), named_vals,
    function(x, v, nom) if (p.value < v) rlang::done(nom) else unnamed_val,
    .init = -Inf) %>% paste0()
}

get_mu_sd <- function(data) {
  data %>%
    group_by(Group, Trial.Order, Structure, Subject) %>%
    dplyr::summarise(RT = mean(RT)) %>% 
    group_by(Group, Trial.Order, Structure) %>%
    dplyr::summarise(mean = mean(RT),
                     sd = sd(RT, na.rm = T)) %>%
    ungroup()
}

# Plots the relationship between mu and sd for a subject
plot_muvsd_per_subj <- function(data, subject_num = NULL,
                                bin_by = c("item","quant_trial","quant_RT")) {
  bin_by <- match.arg(bin_by)
  if (bin_by=="item") grouper=quo(Item)
  else if (bin_by=="quant_trial") grouper=quo(Trial.bin)
  else if (bin_by=="quant_RT") grouper=quo(RT.bin)
  
  
  if (is.null(subject_num))
    subject_num <- sample(unique(data$Subject), 1)
  else 
    subject_num <- as.character(subject_num)
  
  d <- data %>%
    filter(Subject==subject_num) %>%
    mutate(RT.bin = cut(RT, breaks = unique(quantile(RT, probs = seq(0, 1, 0.05))), 
                        include.lowest = T),
           Trial.bin = cut(Trial.Order, breaks = unique(quantile(Trial.Order, probs = seq(0, 1, 0.05))), 
                           include.lowest = T)) %>% 
    group_by(Subject, !!grouper) %>%
    filter(!is.na(RT)) %>%
    dplyr::summarise(
      n = length(RT),
      mean = mean(RT),
      sd = sd(RT)
    )
  
  p <- ggplot(d, aes(x = mean, y = sd)) + 
    geom_point(alpha = .5, size=1) +
    geom_smooth(size = .5, method="gam") +
    ggtitle(paste("Participant ", first(d$Subject))) +
    scale_x_continuous() +
    scale_y_continuous() + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()) 
}

# Mean and SD of RTs at each trial.order (across subjects in a group)
plot_muvsd_for_all <- function(data) {
  d.sd = get_mu_sd(data) 
  r_test = with(d.sd, cor.test(mean, sd))
  sig <- get_asterisks(r_test$p.value)
  
  p <- d.sd %>% 
    ggplot(aes(x = mean, y = sd)) + 
    geom_point(alpha = .5,
               aes(color = Structure,
                   shape = Structure)) +
    geom_smooth(method = "lm", color = "black", se = F, size = .5) +
    zplyr::annotate_abs_text(
      label = paste0("r = ", round(r_test$estimate, 2), sig),
      xpos = 0.25, ypos = 0.8) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_color_discrete("Structure") +
    scale_shape_discrete("Structure")
  # scale_color_manual("Structure",
  #                    breaks = c("Filler", "RC", "MV"),
  #                    values = c(col.FIL, col.RC, col.MV)
  # ) +
  # scale_shape_manual("Structure",
  #                    breaks = c("Filler", "RC", "MV"),
  #                    values = c(shape.FIL, shape.RC, shape.MV)
  # ) 
}

make_jaegerfig_14.bottom <- function(data, bottom_text, has_legend,
                                     return_components = FALSE) {
  bottom_plot <- plot_muvsd_for_all(data)
  if (!has_legend)
    bottom_plot <- bottom_plot + theme(legend.position = "none")
  
  x.grob <- textGrob(bottom_text, gp = gpar(fontface="bold"))
  y.grob <- textGrob("SD", gp = gpar(fontface="bold"), rot=90)
  
  if (return_components) 
    list(bottom_plot = bottom_plot, x.grob = x.grob, y.grob = y.grob)
  else
    grid.arrange(arrangeGrob(bottom_plot, left = y.grob, bottom = x.grob))
}

make_jaegerfig_14.top <- function(data, n_panels, bottom_text, seed=4,
                                  bin_by = c("item","quant_trial","quant_RT"),
                                  return_components = FALSE,
                                  ncols = NA, subjects=NA) {
  if (!is.na(subjects) & !is.null(subjects)) {
    subjs <- subjects
    if (length(subjects) != n_panels) 
      warning("Number of subjects given differs from number of panels supplied")
    n_panels <- length(subjects)
  } else {
    subjs <- sample(unique(data$Subject), n_panels, replace=FALSE)
  }
  if (is.na(ncols)) ncols <- n_panels
  bin_by <- match.arg(bin_by)
  
  set.seed(seed)
  facet_plots <- map(subjs, ~plot_muvsd_per_subj(data, ., bin_by = bin_by))
  
  top_plot <- plot_grid(plotlist = facet_plots, labels=NA, align="h", ncol = ncols)
  x.grob <- textGrob(bottom_text, gp = gpar(fontface = "bold"))
  y.grob <- textGrob("SD", gp = gpar(fontface = "bold"), rot = 90)
  
  if (return_components)
    list(facet_plots = facet_plots, x.grob = x.grob, y.grob = y.grob)
  else
    grid.arrange(arrangeGrob(top_plot, left = y.grob, bottom = x.grob))
}

combine_jaegerfig_parts <- function(p.top, p.bot, rel_heights = c(.5, 0.08, 1)) {
  p <- plot_grid(plotlist = list(p.top, NULL, p.bot), 
                 ncol = 1, labels = c("A", "", "B"), 
                 rel_heights = rel_heights)
  p
}

make_jaegerfig_14 <- function(data, n_panels = 6, seed=4, ncols=NA, return_components = FALSE,
                              top.text = "Mean of length-corrected RTs (by trial)",
                              bottom.text = "Mean of length-corrected RTs (by trial)",
                              bin_by = c("item","quant_trial", "quant_RT"),
                              bottom.legend = TRUE,
                              vspace = 0.08) {
  bin_by <- match.arg(bin_by)
  if (is.na(ncols)) ncols <- n_panels
  
  p.top <- make_jaegerfig_14.top(
    data=data, n_panels=n_panels, seed=seed, bottom_text=top.text, return_components=return_components, ncols=ncols)
  p.bot <- make_jaegerfig_14.bottom(
    data=data, bottom_text=bottom.text, has_legend=bottom.legend, return_components=return_components)
  
  if (return_components) {
    p.facets <- plot_grid(plotlist = p.top$facet_plots, labels=NA, align="h", ncol = n_panels)
    p.final_top_plot <- grid.arrange(arrangeGrob(p.facets, left = p.top$y.grob, bottom = p.top$x.grob))
    p.final_bot_plot <- grid.arrange(arrangeGrob(p.bot$bottom_plot, left = p.bot$y.grob, bottom = p.bot$x.grob))
    p.final <- plot_grid(plotlist = list(p.final_top_plot, NULL, p.final_bot_plot), 
                         ncol = 1, labels = c("A", "", "B"), 
                         rel_heights = c(.5, vspace, 1))
    list(final = p.final, top = p.top, bottom = p.bot)
  } else {
    plot_grid(plotlist = list(p.top, NULL, p.bot), 
              ncol = 1, labels = c("A", "", "B"), 
              rel_heights = c(.5, vspace, 1))
  }
}
# --------------------------------------------------------------

