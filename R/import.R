
#' @importFrom dplyr anti_join arrange bind_rows case_when count desc ends_with everything filter group_by left_join mutate n n_distinct recode rename select semi_join slice starts_with summarize summarise ungroup
#' @importFrom ggplot2 aes facet_grid facet_wrap geom_hline geom_label geom_line geom_point ggplot ggsave guide_legend scale_color_manual scale_colour_manual scale_fill_manual
#' @importFrom lubridate mdy
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom stringr str_c str_replace str_split str_split_fixed
#' @importFrom tibble enframe tibble
#' @importFrom tidyr complete gather pivot_longer pivot_wider spread unite
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
      "Kid",
      "Ma",
      "Pa",
      "Pvalue",
      "Pvalue_comp",
      "Pvalue_comp1",
      "SpawnYear",
      "allele_int",
      "cluster",
      "collection",
      "dummy__",
      "fract",
      "from",
      "gc_wtd_fract",
      "gene_copy",
      "geno",
      "guides",
      "hatchery",
      "height",
      "id",
      "idx_by_year",
      "idx_by_year_and_hatchery",
      "idx_years_aggregated",
      "indiv",
      "indiv1_allele",
      "indiv2_allele",
      "indiv_1",
      "indiv_2",
      "kid",
      "kid_hatchery",
      "kid_id",
      "kid_min_year",
      "kid_spawn_year",
      "kid_year",
      "locus",
      "ma",
      "ma_hatchery",
      "ma_id",
      "ma_year",
      "map_int",
      "max_year",
      "member",
      "min_year",
      "name",
      "new_id",
      "num_discrepant_gene_copies",
      "num_good_genos",
      "num_non_miss",
      "original_id",
      "pa",
      "pa_id",
      "pa_year",
      "parent_compatibility",
      "repunit",
      "retained_id",
      "sample_type",
      "sex",
      "spawner_group",
      "tmp",
      "to",
      "trio_index",
      "trio_type_group",
      "type_f",
      "value",
      "variable",
      "year",
      "years_list"
    ))
}

