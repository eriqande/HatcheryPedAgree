#' Read the output from snppit into a format useful for comparison
#'
#' More later...
#' @param DIR the directory where the output lives
#' @param S the tibble of sex and date.  Must have columns
#' indiv,   sex, spawner_group, and hatchery
#' @export
slurp_snppit <- function(DIR, S) {
  S2 <- S %>%
    select(indiv, hatchery, sex, spawner_group, year)

  P <- read_tsv(file.path(DIR, "snppit_output_ParentageAssignments.txt"), trim_ws = TRUE, na = "---") %>%
    left_join(
      rename(
        S2,
        kid_sex = sex,
        kid_sg = spawner_group,
        kid_year = year,
        kid_hatchery = hatchery
      ),
      by = c("Kid" = "indiv")
    ) %>%
    left_join(
      rename(
        S2,
        pa_sex = sex,
        pa_sg = spawner_group,
        pa_year = year,
        pa_hatchery = hatchery
      ),
      by = c("Pa" = "indiv")
    ) %>%
    left_join(
      rename(
        S2,
        ma_sex = sex,
        ma_sg = spawner_group,
        ma_year = year,
        ma_hatchery = hatchery
      ),
      by = c("Ma" = "indiv")
    ) %>%
    rename(kid = Kid, pa = Pa, ma = Ma) %>%
    select(
      kid,
      pa,
      ma,
      ends_with("sex"),
      ends_with("year"),
      ends_with("sg"),
      ends_with("hatchery"),
      everything()
      )

  P
}

