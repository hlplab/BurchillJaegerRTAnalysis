library(tidyverse)

orig_file <- "/Users/zburchill/Downloads/untitled folder/13428_2012_313_MOESM1_ESM/selfpacedreading.RT.txt"
local_save_path <- "~/Box Sync/Power simulations for RTs/PowerSimulations/"
trial_cutoff <- 144
output_filename <- "franketal_sim_df.RDS"
output_filename2 <- "franketal_sim_df_no_RT_exclusions.RDS"

unix_path <- "~/Box\\ Sync/Power\\ simulations\\ for\\ RTs/PowerSimulations/"
path_on_server <- "/u/zburchil/workspace/launchpad/"



d.preexcl <- read.csv(orig_file, sep = "\t", header = TRUE) %>%
  as_tibble() %>%
  transmute(Subject = factor(subj_nr),
            Item = factor(sent_nr),
            Trial.Order = sent_pos,
            Word.Order = word_pos,
            Correct = correct,
            Word = word,
            RT = RT) %>% 
            {
              d.correctness <- filter(., Correct!="-") %>%
                distinct(Subject, Item, Correct) %>% 
                group_by(Subject) %>%
                summarise(AvgAccuracy = mean(ifelse(Correct=="c", 1, 0)))
              
              zplyr::left_join(., d.correctness, by="Subject")
              
            } %>%
  # Adding word length
  mutate(AlphaWord = gsub("[,\\.$]", "", Word),
         Word.Length = nchar(AlphaWord)) %>%
  select(-AlphaWord)

d.reg_excl <- d.preexcl %>%
  # The accuracy-based subject exclusions (in the original paper)
  filter(AvgAccuracy > 0.75) %>%
  group_by(Subject) %>%
  filter(trial_cutoff %in% Trial.Order) %>%
  filter(Trial.Order <= trial_cutoff) %>%
  ungroup()

d.RT_excl <- d.reg_excl %>%
  # Our way of doing it
  filter(RT > 100 & RT < 2000)
  

saveRDS(d.RT_excl, paste0(local_save_path, output_filename))
saveRDS(d.reg_excl, paste0(local_save_path, output_filename2))

system(
  paste0(
    "rsync -avz ", paste0(unix_path, output_filename), " zburchil@cycle1.cs.rochester.edu:", path_on_server))

system(
  paste0(
    "rsync -avz ", paste0(unix_path, output_filename2), " zburchil@cycle1.cs.rochester.edu:", path_on_server))

