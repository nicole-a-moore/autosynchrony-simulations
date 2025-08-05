select = dplyr::select

## read in avonet to see how many species are resident
avo = read.csv("data-raw/AVONET_Raw.csv")

names(avo)

avo = avo %>% 
  select(Species1, Migration) %>%
  distinct()

avo$Migration
#1 = Sedentary. 
#2 = Partially migratory, i.e. minority of population migrates long distances, or most of population undergoes short-distance migration, nomadic movements, distinct altitudinal migration, etc.
#3 = Migratory, i.e. majority of population undertakes long-distance migration

bb$genus_sp

bb = left_join(bb, avo, by = c("genus_sp" = "Species1"))

bb = bb %>%
  mutate(Migration = ifelse(Migration == 1, "Resident", ifelse(Migration == 2, "Partial migrant",
                                                               ifelse(Migration == 3, "Full migrant", NA))))

bb %>% 
  select(genus_sp, Migration) %>%
  unique() %>%
  ggplot(aes(x = Migration)) + geom_bar()


list = bb %>% 
  select(genus_sp, Migration) %>%
  unique() 


## read in range shift data 
sync = read.csv("outputs/data-processed/bbs-synchrony-stats.csv")
ac = read.csv("outputs/data-processed/bbs-autocorrelation-stats.csv")
rs = read.csv("outputs/data-processed/bbs-range-shift-rates.csv")

rs <- left_join(rs, list) %>%
  left_join(., ac) %>%
  left_join(., sync)


rs %>%
  mutate(Migration = ifelse(Migration == "Partial migrant", "Resident", Migration)) %>%
  ggplot(aes(x = mean_beta, y = slope, colour = mean_mean_pearson_corr)) +
  geom_point() +
  facet_grid(Migration~range_edge)








