
#' convert long format genotypes to rubias format
#'
#' @param genos the long format genotypes tibble with columns
#' `indiv`, `locus`, `gene_copy` and `allele_int` (or another column
#' holding the alleles).
#' @param allele_column unquoted name of the column holding the alleles to
#' use.  By default this is `allele_int`.
#' @param missing quoted symbol used to denote missing data
#' @export
genos2rubias <- function(
  genos,
  allele_column = allele_int,
  missing = "0"
) {

  # check to make sure genos is in the right format
  # .... need to make a general function for this
  allele_column = rlang::enquo(allele_column)

  genos %>%
    mutate(value = ifelse(!! allele_column == missing, NA, !! allele_column)) %>%
    spread(key = locus, value = value)
}
