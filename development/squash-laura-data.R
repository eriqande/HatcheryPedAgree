

library(tidyverse)

cvsh <- read_rds("~/Downloads/cvsh.rds")

# first identify which rendition of an individual to keep
rendition_genos <- cvsh %>%
  count(NMFS_DNA_ID) %>%
  filter(n != 1) %>%
  semi_join(cvsh, ., by = "NMFS_DNA_ID") %>%
  group_by(NMFS_DNA_ID) %>%
  mutate(rendition = 1:n()) %>%
  select(rendition, everything()) %>%
  select(-(BOX_ID:NMFS_ID_VLOOKUP))

keepers <- rendition_genos %>%
  gather(key = "locus", value = "allele_int", -rendition, -NMFS_DNA_ID) %>%
  group_by(NMFS_DNA_ID, rendition) %>%
  summarise(num_non_miss = sum(allele_int != 0) / 2) %>%
  arrange(NMFS_DNA_ID, desc(num_non_miss)) %>%
  slice(1)

# make a wide genos data frame of those individuals
final_dupie_genos <- rendition_genos %>%
  semi_join(keepers, by = c("NMFS_DNA_ID", "rendition")) %>%
  ungroup() %>%
  select(-rendition)

# now, we want to make a similar wide file of all the rest,
# dropping the ones that were duplicated.
non_duped_genos <- cvsh %>%
  anti_join(rendition_genos, by = "NMFS_DNA_ID") %>%
  select(-(BOX_ID:NMFS_ID_VLOOKUP))

wide_genos <- bind_rows(non_duped_genos, final_dupie_genos)
