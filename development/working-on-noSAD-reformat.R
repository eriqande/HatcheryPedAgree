
library(tidyverse)
library(HatcheryPedAgree)


S <- readr::read_rds("private/rrsh_metadata.rds")
# get a data set to play with:
D <- slurp_snppit(
  DIR = "~/Documents/UnsyncedData/HatcheryPedAgreeOutputs/RussianRiver-July-6-2020/snppit-run-noSAD",
  S = S
  )


D2 <- reformat_no_sex_or_date_results(D)
