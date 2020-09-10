
#' Relabel matching genotypes to a single ID and reorganize their years and spawner groups
#'
#' If there are matching samples, we want to make sure that we have just a single
#' ID for each cluster, and that the years and spawner groups get organized
#' so that SNPPIT can read all that information in.
#' @param genotypes the genotypes
#' @param metadata the meta data
#' @param clusters the tibble with the matching sample clusters.
#' @return A list with a variety of components:
#' * `matchers_metadata`:
#' * `snppit_meta`:
#' * `snppit_genos`:
#' * `cross_hatchery_matches`:
#' * `cross_sex_matches`:
#' * `geno_discord`:  a list like that returned in [count_discrepancies()].
#' @export
reorganize_matching_samples <- function(genotypes, metadata, clusters) {

  # make sure that NA spawner groups are coded as ?
  metadata <- metadata %>%
    mutate(spawner_group = ifelse(is.na(spawner_group), "?", spawner_group))

  # figure out which individual to take for each cluster.
  # It will be the one with the most non-missing genotypes and, in case
  # of a tie, the one taken will be pretty much arbitrary.  Note that
  # if there are any IDs in genotypes or clusters that are not in meta data this will
  # cause a problem.  So, we semijoin on the metadata at the end.
  ids_key <- genotypes %>%
    semi_join(clusters, by = "indiv") %>%
    group_by(indiv) %>%
    summarise(num_good_genos = sum(!is.na(allele_int) / 2)) %>%
    left_join(clusters, ., by = "indiv") %>%
    arrange(cluster, desc(num_good_genos)) %>%
    group_by(cluster) %>%
    dplyr::do(tibble(retained_id = .$indiv[1], original_id = .$indiv)) %>%
    ungroup() %>%
    semi_join(metadata, by = c("original_id" = "indiv"))


  # now, we attach that result to the meta data
  matchers_meta <- ids_key %>%
    left_join(metadata, by = c("original_id" = "indiv"))

  # use the matching pairs, thus organized, to compute some reports about
  # genotype discordance, etc.
  geno_discord <- count_discrepancies(
    pairs = matchers_meta,
    genotypes = genotypes
  )

  # check for matchers that are in different hatcheries
  # and bark a message if there are any.  They have to be treated
  # a little differently.
  cross_hatchery_matches <- matchers_meta %>%
    group_by(cluster) %>%
    filter(n_distinct(hatchery) > 1)

  if (nrow(cross_hatchery_matches) > 0) {
    message("Some matching genotypes found in different hatcheries. See cross_hatchery_matches in output.")
  }

  # check for matchers that have different sexes listed
  cross_sex_matches <- matchers_meta %>%
    group_by(cluster) %>%
    filter(n_distinct(sex) > 1)

  if (nrow(cross_hatchery_matches) > 0) {
    message("Some matching genotypes found to have different sexes. See cross_sex_matches in output.")
  }

  # To deal with matchers in multiple hatcheries, we still want to use just a single genotype (with the
  # least missing data) but we will have to use a single canonical ID for each of the hatcheries that
  # the matching samples occur in.  The problem that occurs then is that we could assign an individual
  # to itself.  So, to catch those cases visually, later, we will name those canonical IDs in the other hacheries
  # as `retained-id_hatchery`.
  matchers_meta2 <- matchers_meta %>%
    group_by(cluster, retained_id, hatchery) %>%
    mutate(idx_in_hatchery = 1:n()) %>%
    group_by(cluster, retained_id) %>%
    mutate(new_id = case_when(
      hatchery == hatchery[1] ~ retained_id[hatchery == hatchery[1]],
      hatchery != hatchery[1] ~ str_c(retained_id, hatchery[hatchery != hatchery[1]], sep = "_")
    )) %>%
    select(cluster, new_id, everything()) %>%
    ungroup()

  # below here we are mainly compiling up the genotypes.  There is just one thing that
  # we want to do before that, however: for the retained_id fish that had matching samples
  # that showed a homozygote-to-homozygote discordance (i.e. 11 vs 44), we want to mark those
  # loci as missing in those retained_id fish.  We do that just by modifying genotypes.
  if(nrow(geno_discord$alt_homoz_mismatches) > 0) {
    make_these_na <- geno_discord$alt_homoz_mismatches %>%
      select(retained_id, locus) %>%
      rename(indiv = retained_id) %>%
      mutate(dummy__ = 1)

    genotypes %>%
      left_join(make_these_na, by = c("indiv", "locus")) %>%
      mutate(allele_int = ifelse(!is.na(dummy__), NA, allele_int)) %>%
      select(-dummy__)

    message("Rendered ", nrow(make_these_na), " genotypes NA because of mismatching homozygotes in the matching samples analysis")
  }


  # now, we need to prepare the genotypes that we will use for SNPPIT.  Firt, we want to
  # chuck the genotypes of the original_id's that were not retained
  genos_non_canon_dropped <- matchers_meta2 %>%
    filter(retained_id != original_id) %>% # this lets us keep the canonical genotypes IDs
    anti_join(genotypes, ., by = c("indiv" = "original_id"))

  # But, we also need to copy some genotypes for the cases where the same fish
  # was sampled in multiple hatcheries
  genos_newly_constructed <- matchers_meta2 %>%
    filter(new_id != retained_id) %>% # get the fish that have new genotype names
    count(new_id, retained_id) %>% # get just a single instance of each
    select(-n) %>%
    left_join(., genotypes, by = c("retained_id" = "indiv")) %>%
    select(-retained_id) %>%
    rename(indiv = new_id)

  # combine those to get the data set we need for SNPPIT
  snppit_genos <- bind_rows(
    genos_non_canon_dropped,
    genos_newly_constructed
  )

  # now, we need to condense the meta data that we need
  # into a tibble that has combined the years and the spawn
  # groups that each individual belongs to. We can do this by grouping
  # on new_id and hatchery at this point.
  matchers_meta_snppit <- matchers_meta2 %>%
    rename(indiv = new_id) %>%
    group_by(indiv, hatchery) %>%
    summarise(
      year = paste(sort(unique(year)), collapse = ","),
      spawner_group = paste(unique(spawner_group), collapse = ","),
      sex = sex[1]  # just take the sex of the first (canonical) occurrence of the genotype
    ) %>%
    ungroup()

  # and we will need to combine that with all the fish that were not
  # part of the meta_matchers.  We make a data frame of all the individuals
  # that are not handled by the matchers and then we bind_rows() them and
  # arrange them nicely.
  meta_singletons <- metadata %>%
    anti_join(ids_key, by = c("indiv" = "original_id"))

  snppit_meta <- bind_rows(
    matchers_meta_snppit,
    meta_singletons %>% mutate(year = as.character(year))
  ) %>%
    arrange(hatchery, year, spawner_group)

  # now, down at the end here we return a list of outputs.
  # note that we leave the retained_id column in matchers_meta2, as it
  # can be helpful to know that
  list(
    matchers_metadata = matchers_meta2,
    snppit_meta = snppit_meta,
    snppit_genos = snppit_genos,
    cross_hatchery_matches = cross_hatchery_matches,
    cross_sex_matches = cross_sex_matches,
    geno_discord = geno_discord
  )

}
