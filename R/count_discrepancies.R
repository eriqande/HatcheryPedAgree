
#' from pairs of duplicately genotyped samples, prepare a report of discrepancies
#'
#' After running a matching samples analysis we can use the resulting
#' pairs to investigate genotyping discordance rates at different
#' loci.  That is what this does.
#' @param pairs a tibble with at least two columns: `retained_id` and `original_id`.  The comparisons
#'   are made between a single canonical individual in `retained_id` and any other
#'   matching samples in `original_id`.  If `retained_id == original_id` the row is removed.
#' @param genotypes a tibble with columns `indiv`, `locus`, `gene_copy`, and `allele_int`.
#' @return  Returns a list of three components as follows:
#' * `matching_samples_genos`: a tibble with 7 columns.  The genotypes of two different
#'   individuals occupy two different rows.  The first row is the first gene copy and the
#'   second row is the second gene copy.
#'     * `retained_id`: ID of the fish that is used in downstream analyses.
#'     * `original_id`: ID of the other fish whose genotype is being compared to that
#'        of the retained ID.
#'     * `locus`: the locus name.
#'     * `gene_copy`: the gene copy index (1 or 2)
#'     * `indiv1_allele`: the alleles at the retained_id indiv.  These have been
#'        sorted in ascending order within the locus to make it easy to compare with the indiv2allele
#'     * `indiv2_allele`: same as above, but for original_id
#'     * `num_discrepant_gene_copies`: the number of gene copies that are discrepant.  0 = none;
#'        1 = discrepancy where one is heterozygous and the other homozygous; 2 = one individual
#'        is homozygous and the other is homozygous for the other allele.
#' * `genotype_discrepancies_summary`: counts and fractions of discrepancies across all retained_id vs
#'    original_id pairs at a locus.  All columns should be self-explanatory except for `gc_wtd_fract` which
#'    is the average number of discrepant gene copies per genotype at the locus. This might be
#'    considered a suitable estimate of the per-allele genotyping error rate.
#' * `alt_homoz_mismatches`: a tibble recording all the genotypes amongst the retained_id vs original_id
#'    pairs that are discrepancies in which each member of the pair is homozygous (i.e. they are alternate
#'    homozygotes.)  Each row denotes a mismatching locus.
#' @export
count_discrepancies <- function(
  pairs,
  genotypes
) {

  # get a convenient data structure for comparing genotypes
  pg <- pairs %>%
    select(retained_id, original_id) %>%
    filter(retained_id != original_id) %>%
    left_join(genotypes, by = c("retained_id" = "indiv")) %>%
    rename(indiv1_allele = allele_int) %>%
    left_join(
      genotypes,
      by = c(
        "original_id" = "indiv",
        "locus",
        "gene_copy"
        )
    ) %>%
    rename(indiv2_allele = allele_int) %>%
    group_by(retained_id, original_id, locus) %>%
    mutate(
      indiv1_allele = {
        if(!any(is.na(indiv1_allele)))
          sort(indiv1_allele)
        else
          indiv1_allele
      },
      indiv2_allele = {
        if(!any(is.na(indiv2_allele)))
          sort(indiv2_allele)
        else
          indiv2_allele
      }
    ) %>%
  mutate(num_discrepant_gene_copies = sum(indiv1_allele != indiv2_allele)) %>%
  ungroup()

  # now we can tally up the number and fraction of times each locus has
  # 0, 1, or 2 discrepant alleles in all the comparisons.
  report <- pg %>%
    group_by(locus, num_discrepant_gene_copies) %>%
    summarise(
      n = n() / 2
    ) %>%
    filter(!is.na(num_discrepant_gene_copies)) %>%
    complete(
      num_discrepant_gene_copies = c(0:2),
      fill = list(n = 0)
    ) %>%
    mutate(
      fract = n / sum(n),
      gc_wtd_fract = sum(fract * num_discrepant_gene_copies)
    ) %>%
    arrange(desc(gc_wtd_fract))

  # finally, we record the retained_indiv by locus combos that
  # have mismatches between alternate homozygotes:
  alt_homoz_mismatches <- pg %>%
    filter(num_discrepant_gene_copies == 2) %>%
    group_by(retained_id, locus) %>%
    slice(1) %>%
    ungroup()

  # finally, return a list
  list(
    matching_samples_genos = pg,
    genotype_discrepancies_summary = report,
    alt_homoz_mismatches = alt_homoz_mismatches
  )

}












