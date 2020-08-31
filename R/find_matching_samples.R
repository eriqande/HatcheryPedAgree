#' Find matching samples to a single representative
#'
#' Given a specification of pairs of genotypes that must be from the
#' same individual, this identifies connected components and returns a tibble
#' with one column, `indiv`, and another column, `aliases` which is a list column, that
#' includes, for each `indiv`, all the other names it is known by.
#'
#' @param genotypes A tibble like coho_genotypes that has
#' @param return_cluster Set to TRUE by default, but you might not want
#' to do this if you have a very permissive cutoff.  It makes a graph of
#' the pairs and finds the connected components.
#' @param ... parameters to be passed to rubias::close_matching_samples.
#' Intended to be used for `min_frac_non_miss` and `min_frac_matching`
#' @return Returns a list with three components as follows:
#' * `pairs`: A tibble holding the matching pairs that were found.  It has
#'   the following columns:
#'     * `num_non_miss`: number of loci missing in neither member of the pair
#'     * `num_match`: number of non-missing loci having the same genotype in
#'       in each member of the pair.
#'     * `indiv_1`: the ID of the first member of the pair.
#'     * `indiv_2`:
#' * `clusters`:
#' * `aliases`:
#' @export
#' @examples
#' # There are not actually any matching samples in coho_genotypes
#' # but we will just create some pairs that match by cranking
#' # the min_frac_matching down to 80%
#' find_matching_samples(coho_genotypes, min_frac_matching = 0.80)
find_matching_samples <- function(
  genotypes,
  return_clusters = TRUE,
  ...
) {

  clusters <- NULL
  aliases <- NULL

  # find the pairs
  pairs <- rubias::close_matching_samples(
    D = genos2rubias(genotypes),
    gen_start_col = 5,
    ...
  ) %>%
    select(num_non_miss:indiv_2)

  # now create a graph and find the connectected components
  if (return_clusters == TRUE) {
    clusters <- pairs %>%
      rename(to = indiv_1, from = indiv_2) %>%
      select(to, from) %>%
      igraph::graph_from_data_frame() %>%
      igraph::components() %>%
      .$membership %>%
      enframe() %>%
      rename(indiv = name, cluster = value) %>%
      arrange(cluster, indiv)

    aliases <- clusters %>%
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

  list(pairs = pairs, clusters = clusters, aliases = aliases)

}
