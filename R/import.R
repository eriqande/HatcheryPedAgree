
#' @importFrom dplyr arrange case_when ends_with everything filter group_by left_join mutate recode rename select summarize ungroup
#' @importFrom ggplot2 aes facet_wrap geom_label geom_line geom_point ggplot guide_legend scale_colour_manual scale_fill_manual
#' @importFrom lubridate mdy
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom stringr str_c
#' @importFrom tibble tibble
#' @importFrom tidyr gather spread unite
#' @importFrom utils write.table
NULL


# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if (getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "1",
      "2",
      "FDR",
      "FDR_comp",
      "FDR_comp1",
      "Indiv",
      "Kid",
      "Locus",
      "Ma",
      "Pa",
      "Pvalue",
      "Pvalue_comp",
      "Pvalue_comp1",
      "SEX",
      "collection_date",
      "gene_copy",
      "geno",
      "gt_miss10p",
      "guides",
      "height",
      "idx_by_year",
      "idx_years_aggregated",
      "kid_date",
      "kid_spawn_year",
      "parent_compatibility",
      "sex",
      "sex_and_date",
      "snp_as_int",
      "spawn_date",
      "trio_type_group",
      "type_f",
      "year"
    ))
}
