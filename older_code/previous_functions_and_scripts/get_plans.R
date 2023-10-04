library(dplyr)
library(cs)

path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/functions/"
source(paste0(path, "amlap_constants.R"))
source(paste0(path, "utils.R"))


gateway_node <- backup_login


# good_nodes <- test_nodes(paste0("node",72:92), gateway_node)
# Good nodes:
good_nodes <- c("node74", "node75", "node76", "node77", "node79", "node80", 
  "node81", "node83", "node84", "node85", "node86", "node87", "node88", 
  "node89", "node91", "node92")

# Make the elite_squad
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = gateway_node),
  tweak(cluster, workers = good_nodes),
  tweak(multiprocess, workers = function() { max(round(future::availableCores()/2), 4) })
)); good_node_plan <- plan("list")
stop("CURRENTLY THIS IS ALL YOU CAN GET")


node_df <- gateway_node %>%
  get_nodes_info(timeout_sec = 30,
                 check_node = c("node74")) %>%
  default_cleanup() %>%
  default_filter() %>%
  arrange(-availmem) %>%
  select(number, percent_free, availmem, ncpus) %>%
  mutate(nodename = paste0("node",number)) %>%
  filter(!(number %in% janky_nodes))

heavy_hitters <- node_df %>%
  filter(ncpus > 8)
decent_nodes <- node_df %>%
  filter(ncpus >= 8)
weak_links <- node_df %>%
  filter(ncpus <= 8)
super_weakos <- node_df %>%
  filter(ncpus < 8)


# Make the elite_squad
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = gateway_node),
  tweak(cluster, workers = heavy_hitters$nodename),
  tweak(multiprocess, workers = function() { max(round(future::availableCores()/2), 4) })
)); elite_squad_plan <- plan("list")

# Make the losers
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = gateway_node),
  tweak(cluster, workers = weak_links$nodename),
  tweak(multiprocess, workers = function() { max(round(future::availableCores()/2), 1) })
)); b_squad_plan <- plan("list")

# Make the super-weaklings
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = gateway_node),
  tweak(cluster, workers = super_weakos$nodename),
  tweak(multiprocess, workers = function() { max(round(future::availableCores()/2), 1) })
)); c_squad_plan <- plan("list")

# Make the node92 plane
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = gateway_node),
  tweak(cluster, workers = "node92"),
  tweak(multiprocess, workers = function() { max(round(future::availableCores()*3/4), 4) })
)); node92_plan <- plan("list")

# Experimental
# Flat, no hierachy, each according to his own
commie_iterator <- function(dfx, pf, f, n_chunks=3) {
  breaks <- floor(seq(0, nrow(dfx), length.out = n_chunks+1))
  breaks <- map2(breaks[-1], lag(breaks)[-1]+1, ~.y:.x)
  message(paste(breaks, collapse=" - "))
  pf(breaks, function(.what) {
    dfz <- dfx[.what,]
    iterate_over_df2(dfz, pf, f)
  }
  )
}
# Subs in the data frame
iterate_over_df2 <- function (df, pf, f) {
  stopifnot(length(formals(f)) == 1)
  formals(f) <- alist(.dontuse = )
  orig_bod <- body(f)
  new_bod <- substitute({
    .orig_row <- df[.dontuse, ]
    with(data = df[.dontuse, ], orig_bod)
  }, list(df=df, orig_bod=orig_bod))
  body(f) <- new_bod
  pf(1:nrow(df), f)
}
commie_workers <- node_df %>%
  arrange(-ncpus) %>%
  # filter(ncpus <= 8) %>%
  mutate(ncpus = ifelse(ncpus>=24, ncpus+7, ncpus)) %>%
  mutate(workers = map2(ncpus, nodename, ~rep(list(.y), times=floor(.x/2)))) %>%
  tidyr::unnest(workers) %>% 
  {
    mutate(., the_server = sample(c(server, backup_server, tertiary_server),
             size=nrow(.), replace=TRUE))
  }
f<-bquote(
  function() {
    dfx <- .(commie_workers)
    temp <- Sys.info()[["nodename"]] 
    if (temp==.(server))
      dplyr::filter(dfx, the_server==.(server))$nodename
    else if (temp==.(backup_server))
      dplyr::filter(dfx, the_server==.(backup_server))$nodename
    else if (temp==.(tertiary_server))
      dplyr::filter(dfx, the_server==.(tertiary_server))$nodename
    else 1 })
ff<-eval(f)
plan(list(
  multisession,
  multisession,
  tweak(remote, workers=c(remote_login, backup_login, tertiary_login)),
  tweak(cluster, workers = ff)
)); communist_plan <- plan("list")

# Make the plan with numbers > 72; essentially "ok_nodes"
plan(list(
  multisession, # for the beep!
  tweak(remote,  workers = gateway_node),
  tweak(cluster, workers = decent_nodes$nodename),
  tweak(multiprocess, workers = function() { max(round(future::availableCores()/2), 4) })
)); decent_plan <- plan("list")




# d_test <- data.frame(x=1:135, y=paste0("a",1:135))
# mm_listsx %beep% {
#   iterate_over_df(
#     d_test,
#     future_map,
#     function(zaza) {
#       Sys.info()[["nodename"]]
#     })
# } %plan% b_squad_plan
# done(mm_listsx)
