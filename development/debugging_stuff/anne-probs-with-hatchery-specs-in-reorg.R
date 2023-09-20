

rm(list = ls())
library(tidyverse)
library(HatcheryPedAgree)


load("~/Downloads/for_eric.rda")

reorganize_matching_samples(genos_na, rrsh_metadata, for_histo$clusters)

# oh my goodness.
ggplot(for_histo$pairs, aes(x = num_match / num_non_miss)) +
  geom_histogram()


# step by step

genotypes <- genos_na
metadata <- rrsh_metadata
clusters <- for_histo$clusters
