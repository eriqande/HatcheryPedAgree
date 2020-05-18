
#' convert long format genotypes to rubias format
#'
#' @param genos the long format genotypes tibble with columns
#' `indiv`, `locus`, `gene_copy` and `allele_int`.  Missing data
#' in the allele_int column must be denoted by NA.
#' @param allele_column unquoted name of the column holding the alleles to
#' use.  By default this is `allele_int`.
#' @param missing quoted symbol used to denote missing data
#' @export
genos2rubias <- function(
  genos,
  allele_column = allele_int
) {

  # check to make sure genos is in the right format
  # .... need to make a general function for this

  # enquo allele_columm for NSE later
  allele_column = rlang::enquo(allele_column)

  # making missing data NA, because that is what rubias expects. Then
  # return just the columns we want, and pivot it wider
  tmp <- genos %>%
    mutate(value = !! allele_column) %>%  # just a simple name change
    select(indiv, locus, gene_copy, value) %>%
    mutate(value = as.character(value)) %>%  # rubias wants these columns to be characters
    tidyr::pivot_wider(
      names_from = c(locus, gene_copy),
      values_from = value
    )

  # make first instance of locus name have now "_1" extension, and make the second
  # instance have an "_1" extension.
  names(tmp)[seq(2, ncol(tmp), by = 2)] <- str_replace(names(tmp)[seq(2, ncol(tmp), by = 2)], "_1$", "")
  names(tmp)[seq(3, ncol(tmp), by = 2)] <- str_replace(names(tmp)[seq(3, ncol(tmp), by = 2)], "_2$", "_1")

  # now, add the columns that rubias needs
  tmp %>%
    mutate(
      sample_type = "reference",
      repunit = "dummy_repunit",
      collection = "dummy_collection"
    ) %>%
    select(sample_type, repunit, collection, indiv, everything())

}
