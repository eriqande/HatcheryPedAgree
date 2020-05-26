#' Find matching samples to a single representative
#'
#' Given a specification of pairs of genotypes that must be from the
#' same individual, this identifies connected components and returns a tibble
#' with one column, `indiv`, and another column, `aliases` which is a list column, that
#' includes, for each `indiv`, all the other names it is known by.
#'
#' @param genotypes A tibble like coho_genotypes that has
#' @param ... parameters to be passed to rubias::close_matching_samples.
#' Intended to be used for `min_frac_non_miss` and `min_frac_matching`
#' @export
#' @examples
#' # There are not actually any matching samples in coho_genotypes
#' # but we will just create some pairs that match by cranking
#' # the min_frac_matching down to 80%
#' find_matching_samples(coho_genotypes, min_frac_matching = 0.80)
find_matching_samples <- function(genotypes, return_aliases = FALSE, ...) {

  aliases <- NULL

  # find the pairs
  pairs <- rubias::close_matching_samples(
    D = genos2rubias(genotypes),
    gen_start_col = 5,
    ...
  ) %>%
    select(num_non_miss:indiv_2)

  # now create a graph and find the connectected components
  if (return_aliases == TRUE) {
    aliases <- pairs %>%
      rename(to = indiv_1, from = indiv_2) %>%
      select(to, from) %>%
      igraph::graph_from_data_frame() %>%
      igraph::components() %>%
      .$membership %>%
      enframe() %>%
      rename(indiv = name, cluster = value) %>%
      group_by(cluster) %>%
      mutate(
        csz = n(),
        tmp = list(indiv)
      ) %>%
      ungroup() %>%
      mutate(
        aliases = purrr::map2(.x = tmp, .y = indiv, function(x, y) setdiff(x, y))
      ) %>%
      select(indiv, aliases)
  }

  list(pairs = pairs, aliases = aliases)

}
