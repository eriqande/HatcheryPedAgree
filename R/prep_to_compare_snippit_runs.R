#' Prepare SNPPIT runs to compare them
#'
#' More later
#' @param L a list of outputs from different snppit runs.  These should be named by whatever
#' you want to call each run.
#' @export
prep_to_compare_snppit_runs <- function(
  L
) {

  type_levels = names(L)

  # get the number of the types of analyses
  num_types <- length(type_levels)

  # combine the list into a single tibble
  D <- bind_rows(L, .id = "type_of_analysis")

  D$type_f = factor(D$type_of_analysis, levels = type_levels)

  # now, we need to fill in values for the FDR (and p-value) for sorting these
  # things along the x-axis.  Namely we need the values for the "base" case (type_f = 1),
  # if those were NA. I will fill those with the average of the remaining values.
  # We also need to get the kid's spawn year in there for facetting.
  D2 <- D %>%
    arrange(kid_hatchery, kid, type_f) %>%
    group_by(kid_hatchery, kid) %>%
    mutate(
      FDR_comp = ifelse(is.na(FDR), mean(FDR, na.rm = TRUE), FDR),
      Pvalue_comp = ifelse(is.na(Pvalue), mean(Pvalue, na.rm = TRUE), Pvalue)
    ) %>%
    mutate(FDR_comp1 = FDR_comp[1],  # these are the "composite" scores
           Pvalue_comp1 = Pvalue_comp[1]) %>%
    arrange(FDR_comp1, Pvalue_comp1, kid) %>%
    ungroup() %>%
    mutate(idx_all_aggregated = as.integer(factor(kid, levels = unique(kid)))) %>%
    mutate(
      kid_min_year = map_int(
        str_split(kid_year, ","),
        .f = function(x) min(as.integer(x))
        ))



  # Now, give each one an index based on the year that they are in, and also
  # on the year and the hatchery that they are in
  D3 <- D2 %>%
    group_by(kid_min_year) %>%
    mutate(idx_by_year = as.integer(factor(kid, levels = unique(kid)))) %>%
    group_by(kid_hatchery, kid_min_year) %>%
    mutate(idx_by_year_and_hatchery = as.integer(factor(kid, levels = unique(kid)))) %>%
    ungroup() %>%
    select(starts_with("idx"), everything())


  return(D3)

  ### Made it to here so far, and just testing
  g <- ggplot(D3 %>% filter(!is.na(kid_hatchery)), aes(x = idx_by_year_and_hatchery, y = FDR + 0.05 * (as.integer(type_f) - 1), colour = type_f)) +
    geom_hline(yintercept = 0.05, colour = "white") +
    geom_point(shape = 21, stroke = 0.2, size = 1.2) +
    facet_grid(kid_min_year ~ kid_hatchery) +
    scale_color_manual(values = c("blue", "red"))
  ggsave(g, filename = "development/outputs/RR-sad-no-sad-fdr-hathchery-by-year-facet-grid.pdf",
         width = 25, height = 14)

  # now, we want to pick out those indivs that are not assigned to the same
  # parents as the base case (type_f == 1) when both of them are non-missing
  base_casers <- D3 %>%
    select(type_f, kid, pa, ma) %>%
    rename(base_pa = pa,
           base_ma = ma) %>%
    filter(as.integer(type_f) == 1) %>%
    select(-type_f)

  D3_plus <- D3 %>%
    left_join(base_casers, by = "kid")

  # finally, we want to get some things to make colors
  D4 <- D3_plus %>%
    mutate(
      parent_compatibility = case_when(
        (!is.na(ma) & !is.na(pa) & !is.na(base_pa) & !is.na(base_ma)) & !((ma == base_ma & pa == base_pa) | (pa == base_ma & ma == base_pa)) ~ "ma_pa_not_matching_base_case",
        is.na(ma_date) | is.na(pa_date)       ~ "date_not_recorded",
        is.na(ma_sex)  | is.na(pa_sex)        ~ "sex_not_recorded",
        pa_sex == ma_sex & pa_date != ma_date ~ "sex_and_date_mismatch",
        pa_sex == ma_sex                      ~ "sex_mismatch",
        pa_date != ma_date                    ~ "date_mismatch",
        TRUE                                  ~ "pair_compatible"
      ),
      trio_type_group = case_when(
        MaxP.Pr.Relat == "C_Se_Se"  ~ "C_Se_Se",
        is.na(MaxP.Pr.Relat)      ~ NA_character_,
        TRUE                      ~ "NOT C_Se_Se"
      )
    )

  D4
}
