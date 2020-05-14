#' Prepare SNPPIT runs to compare them
#'
#' More later
#'
#' @export
prep_to_compare_snppit_runs <- function(
  D,
  type_levels = c("nSD_full", "nSD_dump10", "wSD_full", "wSD_dump10"),
  #type_levels = c("nSD_full", "nSD_dump10", "nSD_dumpE", "nSD_Hil", "wSD_full", "wSD_dump10", "wSD_dumpE", "wSD_Hil")
  GTseqed = sex_and_date$Indiv,
  GT_hi_miss = gt_miss10p
) {

  # get the number of the types of analyses
  num_types <- length(type_levels)

  D$type_f = factor(D$type_of_analysis, levels = type_levels)

  # now, we need to fill in values for the FDR (and p-value) for sorting these
  # things along the x-axis.  Namely we need the values for the "base" case (type_f = 1),
  # if those were NA. I will fill those with the average of the remaining values.
  # We also need to get the kid's spawn here in there for facetting.
  D2 <- D %>%
    mutate(kid_spawn_year = year(mdy(kid_date))) %>%
    arrange(Kid, type_f) %>%
    group_by(Kid) %>%
    mutate(
      FDR_comp = ifelse(is.na(FDR), mean(FDR, na.rm = TRUE), FDR),
      Pvalue_comp = ifelse(is.na(Pvalue), mean(Pvalue, na.rm = TRUE), Pvalue)
    ) %>%
    group_by(Kid) %>%
    mutate(FDR_comp1 = FDR_comp[1],
           Pvalue_comp1 = Pvalue_comp[1]) %>%
    arrange(FDR_comp1, Pvalue_comp1, Kid) %>%
    ungroup() %>%
    mutate(idx_years_aggregated = as.integer(factor(Kid, levels = unique(Kid))))



  # Now, give each one an index based on the year that they are in
  D3 <- D2 %>%
    group_by(kid_spawn_year) %>%
    mutate(idx_by_year = as.integer(factor(Kid, levels = unique(Kid)))) %>%
    ungroup()

  # now, we want to pick out those indivs that are not assigned to the same
  # parents as the base case (type_f == 1) when both of them are non-missing
  base_casers <- D3 %>%
    select(type_f, Kid, Pa, Ma) %>%
    rename(basePa = Pa,
           baseMa = Ma) %>%
    filter(as.integer(type_f) == 1) %>%
    select(-type_f)

  D3_plus <- D3 %>%
    left_join(base_casers, by = "Kid")

  # finally, we want to get some things to make colors
  D4 <- D3_plus %>%
    mutate(
      parent_compatibility = case_when(
        !is.na(Ma) & !is.na(Pa) & !is.na(basePa) & !is.na(baseMa) & ((Ma != basePa & Ma != baseMa) | (Pa != basePa & Pa != baseMa)) ~ "MaPa_not_matching_base_case",
        (Ma %in% GT_hi_miss) | (Pa %in% GT_hi_miss) ~ "Ma_or_Pa_is_potential_GT_Seq_hi_misser",
        !(Kid %in% GTseqed & Ma %in% GTseqed & Pa %in% GTseqed) ~ "full_trio_not_GTseq-ed",
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
