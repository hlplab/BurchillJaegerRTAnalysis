library(dplyr)

mark_progress <- function(text, val="a", extra_text=NULL) {
  saveRDS(val,
          paste("~/tmp/prog", Sys.info()[["nodename"]], text, extra_text, sep="_") %>%
            paste0(".RDS"))
}

iterate_over_df <- function(df, pf, f) {
  stopifnot(length(formals(f)) == 1)
  formals(f) <- alist(".dontuse" = )
  orig_bod <- body(f)
  # gives a variable .orig_row if you want to use the original row
  new_bod <- substitute({ .orig_row <- df[.dontuse, ]; with(data = df[.dontuse, ], orig_bod ) })
  body(f) <- new_bod

  pf(1:nrow(df), f)
}

remaining_files <- function(df, bare_predicate) {
  bare_predicate <- rlang::enquo(bare_predicate)
  res <- cs::futureBeep({
    iterate_over_df(
      df,
      purrr::map_dfr,
      function(i) {
        .orig_row %>% mutate(there = !! bare_predicate)
      }) %>%
      filter(there==FALSE)
  })
  value(res)
}

nonempty_remaining_files <- function(df, filecol) {
  filecol <- rlang::enquo(filecol)
  bare_predicate <- quo(file.exists(!!filecol) & file.size(!!filecol) !=0)
  remaining_files(df, !!bare_predicate)
}

file_is_readable <- function(file) {
  res <- NULL
  suppressWarnings(tryCatch({l<-readRDS(file); res<-TRUE}, error=function(x) res<<-FALSE))
  return(res)
}


