#' Read the outpur from snppit into a format useful for comparison
#'
#' More later...
#' @param DIR the directory where the output lives
#' @param S the tibble of sex and date.  Must have columns
#' Indiv,   SEX, and    collection_date
slurp_snppit <- function(DIR, S) {
  P <- read_tsv(file.path(DIR, "snppit_output_ParentageAssignments.txt"), trim_ws = TRUE, na = "---") %>%
    left_join(S %>% rename(kid_sex = SEX, kid_date = collection_date), by = c("Kid" = "Indiv")) %>%
    left_join(S %>% rename(pa_sex = SEX, pa_date = collection_date), by = c("Pa" = "Indiv")) %>%
    left_join(S %>% rename(ma_sex = SEX, ma_date = collection_date), by = c("Ma" = "Indiv")) %>%
    select(ends_with("sex"), ends_with("date"), everything())
}

