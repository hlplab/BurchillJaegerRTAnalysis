library(tidyverse)


orig_dir <- "~/Downloads/NSC_RT_data/"

local_save_path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/"

output_filename <- "nsc_sim_df.RDS"
output_filename2 <- "nsc_sim_df_no_RT_exclusions.RDS"

unix_path <- "~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/PowerSimulations/"
path_on_server <- "/u/zburchil/workspace/launchpad/"


# Approximate where the sentence ends are by looking for ".", "?" and "!"
sentence_breaks <- read.csv(paste0(orig_dir, "naturalstories-master/words.tsv"),
                            sep="\t", header=F) %>%
  as_tibble() %>%
  filter(V2 %in% c(".", "?", "!")) %>%
  distinct(V1, V2) %>%
  tidyr::separate(V1, c("item","zone","blah")) %>%
  select(item, zone) %>%
  mutate_all(~as.numeric(as.character(.))) %>%
  group_by(item) %>%
  mutate(maxWords = zone-lag(zone, default=0),
         lowerrange = lag(zone, default=0)+1,
         upperrange = zone,
         sennum = seq_along(zone)) %>%
  ungroup()

# Get the sentences breaks in such a way that we can left_join them
expanded_breaks <- sentence_breaks %>%
  rowwise() %>%
  mutate(zone2 = paste(lowerrange:upperrange, collapse=",")) %>%
  ungroup() %>%
  tidyr::separate_rows(zone2, sep=",") %>%
  mutate(zone = as.numeric(zone2)) %>%
  select(item, zone, sennum)

# Get the words to join to the RTs
word.df <- read.csv(paste0(orig_dir, 'all_stories.tok'), sep = '\t')

# Get the original RTs
d.orig <- rbind(read.csv(paste0(orig_dir, 'batch2_pro.csv')),
            read.csv(paste0(orig_dir, 'batch1_pro.csv'))) %>% 
  group_by(WorkerId) %>% 
  mutate(item_order = paste(unique(item), collapse=" ")) %>%
  ## subtract 2 from zone to properly align region.
  ## but the RTs seem to line up correctly in plots
  mutate(zone = zone - 2) %>%
  ungroup() %>% 
  #read in story words and region
  #item is story (1-10), zone is RT region
  left_join(word.df, by=c("item","zone")) %>%
  zplyr::left_join(expanded_breaks) %>%
  #remove regions that do not have words
  filter(!is.na(word)) %>%
  
  # Add in the word order stuff
  group_by(WorkerId) %>%
  mutate(Temp = seq_along(item)) %>%
  group_by(WorkerId, item, sennum) %>%
  mutate(Temp = min(Temp, na.rm = TRUE)) %>%
  group_by(WorkerId) %>%
  mutate(Trial.Order = dense_rank(Temp)) %>%
  # Make sentence-level IDs:
  group_by(WorkerId, item) %>%
  mutate(ItemToBe = paste0(item, "-", dense_rank(Trial.Order))) %>%
  
  # Make it in the standardized format
  group_by(WorkerId, item, sennum) %>%
  mutate(Word.Order = zone - min(zone, na.rm = TRUE) + 1) %>%
  ungroup() %>%
  mutate(Word.Length = nchar(as.character(word)),
         item = factor(item),
         ItemToBe = factor(ItemToBe)) %>%
  rename(Subject = WorkerId,
         Story = item,
         Item = ItemToBe,
         Word = word) %>%
  select(Story, Subject, Item, Trial.Order, Word.Order, Word.Length, Word, RT, 
         item_order, correct, zone)
  
d.preexcl <- d.orig %>%
  #exclude stories where subject does not get more than 4/6 correct
  filter(correct > 4) %>%
  # Also exclude any participants with more than 3/4s of their RTs sub 100 or above 2000
  mutate(isbad = ifelse(RT <= 100 | RT >= 2000, 1, 0)) %>%
  group_by(Subject) %>%
  filter(mean(isbad) < 0.75) %>%
  ungroup() %>% select(-isbad) %>%
  select(-correct, -item_order)
  
#exclude data points less than 50 ms, greater than 3000 ms
d.post_RT_excl <- d.preexcl %>%
  filter(RT > 100, RT < 2000)

# Reminder that our definitions of sentences isn't really theirs:
reported_sentences <- c(57,37,55,55,45,64,48,33,48,43)
reported_words <- c(1073,990,1040,1085,1013,1089,999,980,1038,938)

d.orig %>%
  transmute(item = Story, zone = zone) %>%
  distinct(item, zone) %>%
  left_join(mutate(expanded_breaks, item=factor(item))) %>%
  group_by(item) %>%
  summarise(n = n_distinct(sennum),
            numwords = n()) %>%
  arrange(item) %>%
  mutate(reported_sens = reported_sentences,
         reported_words = reported_words) %>%
  filter(n!=reported_sens | numwords!=reported_words)

saveRDS(d.post_RT_excl, paste0(local_save_path, output_filename))
saveRDS(d.preexcl,      paste0(local_save_path, output_filename2))

system(
  paste0(
    "rsync -avz ", paste0(unix_path, output_filename), " zburchil@cycle1.cs.rochester.edu:", path_on_server))

system(
  paste0(
    "rsync -avz ", paste0(unix_path, output_filename2), " zburchil@cycle1.cs.rochester.edu:", path_on_server))

