#' prepare long-format genotypes for snppit
#'
#' More later...
#' @param G a long format data frame of SNPs.  It must have, at a minimum, the columns
#'   `indiv`, `locus`, `gene_copy`, `allele_int`.  Missing data should be represented as
#'   NA.
#' @param S a data frame of meta data with the columns indiv,   sex, and    spawner_group
#' @export
#' @examples
#' prep_for_snppit(G = coho_genotypes, S = coho_metadata)
prepare_snppit_infile <- function(G,
                            S,
                            use_spawner_group = TRUE,
                            use_sex = TRUE,
                            min_age = 2,
                            max_age = 4,
                            geno_err = 0.005,
                            outf = tempfile()) {

  # these will get changed if we don't use them
  sex_col <- "POPCOLUMN_SEX"
  spawn_group_col <- "POPCOLUMN_SPAWN_GROUP"


  # first, make a genotype matrix. Put spaces between the gene copies in a locus
  M <- G %>%
    select(indiv, locus, gene_copy, allele_int) %>%
    mutate(allele_int = ifelse(is.na(allele_int), 0, allele_int)) %>%
    spread(key = gene_copy, value = allele_int) %>%
    unite(col = "geno", `1`, `2`, sep = " ") %>%
    spread(key = locus, value = geno)

  # get the marker names
  SNP_names <- colnames(M)[-1]

  # prep the sex and spawn date data
  meta <- S %>%
    mutate(
      sex = recode(sex, Female = "F", Male = "M"),
      sex = ifelse(is.na(sex), "?", sex),
    ) %>%
    arrange(year, spawner_group) %>%
    select(indiv, sex, year, spawner_group)

  # break the genos into candidate parents and offspring.
  # Basically if your spawn year is less than the min year + the
  # min age, then you don't get to be an offspring.  Likewise, if
  # your spawn year is greater than max year - min age, you don't get
  # to be a candidate parent
  Offs <- meta %>%
    filter(year >= min(year) + min_age) %>%
    select(-sex, -spawner_group) %>%
    mutate(age_range = str_c(min_age, "-", max_age)) %>%
    left_join(., M, by = "indiv")

  Pars <- meta %>%
    filter(year <= max(year) - min_age) %>%
    left_join(., M, by = "indiv")


  if (use_sex == FALSE) {
    sex_col <- ""
    Pars <- Pars %>% select(-sex)
  }
  if (use_spawner_group == FALSE) {
    spawn_group_col <- ""
    Pars <- Pars %>% select(-spawn_date)
  }


  # now assemble the text ouput
  preamble <- glue::glue(
    "NUMLOCI {length(SNP_names)}
MISSING_ALLELE 0
{sex_col}
POPCOLUMN_REPRO_YEARS
{spawn_group_col}
OFFSPRINGCOLUMN_SAMPLE_YEAR
OFFSPRINGCOLUMN_AGE_AT_SAMPLING\n
"
  )

  cat(preamble, file = outf)
  write.table(cbind(SNP_names, geno_err),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    file = outf, append = TRUE
  )
  cat("POP Parents\n", file = outf, append = TRUE)
  write.table(Pars,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    file = outf, append = TRUE, sep = "\t"
  )
  cat("OFFSPRING Offspring Parents\n", file = outf, append = TRUE)
  write.table(Offs,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    file = outf, append = TRUE, sep = "\t"
  )

  # message to user
  message("snppit input file written to ", outf)
}