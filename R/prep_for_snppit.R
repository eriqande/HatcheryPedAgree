#' prepare long-format genotypes for snppit
#'
#' More later...
#' @param G a long format data frame of SNPs.  It must have, at a minimum, the columns
#'   `Indiv`, `Locus`, `gene_copy`, `snp_as_int`.  Missing data should be represented as
#'   a 0.
#' @param S a data frame of meta data with the columns Indiv,   SEX, and    collection_date
prep_for_snppit <- function(G, S,
                            use_spawn_date = TRUE,
                            use_sex = TRUE,
                            minAge = 2,
                            maxAge = 4,
                            geno_err = 0.005,
                            outf = "snippit_out.txt") {

  # these will get changed if we don't use them
  SEX_COL <- "POPCOLUMN_SEX"
  SPAWN_GROUP_COL <- "POPCOLUMN_SPAWN_GROUP"


  # first, make a genotype matrix. Put spaces between the gene copies in a locus
  M <- G %>%
    select(Indiv, Locus, gene_copy, snp_as_int) %>%
    spread(key = gene_copy, value = snp_as_int) %>%
    unite(col = "geno", `1`, `2`, sep = " ") %>%
    spread(key = Locus, value = geno)

  # get the marker names
  SNP_names <- colnames(M)[-1]

  # prep the sex and spawn date data
  meta <- sex_and_date %>%
    mutate(
      spawn_date = mdy(collection_date),
      sex = recode(SEX, Female = "F", Male = "M"),
      sex = ifelse(is.na(sex), "?", sex),
      year = year(spawn_date)
    ) %>%
    arrange(spawn_date) %>%
    select(Indiv, sex, year, spawn_date)

  # break the genos into candidate parents and offspring.
  # Basically if your spawn year is less than the min year + the
  # min age, then you don't get to be an offspring.  Likewise, if
  # your spawn year is greater than max year - min age, you don't get
  # to be a candidate parent
  Offs <- meta %>%
    filter(year >= min(year) + minAge) %>%
    select(-sex, -spawn_date) %>%
    mutate(age_range = str_c(minAge, "-", maxAge)) %>%
    left_join(., M, by = "Indiv")

  Pars <- meta %>%
    filter(year <= max(year) - minAge) %>%
    left_join(., M, by = "Indiv")


  if (use_sex == FALSE) {
    SEX_COL <- ""
    Pars <- Pars %>% select(-sex)
  }
  if (use_spawn_date == FALSE) {
    SPAWN_GROUP_COL <- ""
    Pars <- Pars %>% select(-spawn_date)
  }


  # now assemble the text ouput
  preamble <- glue::glue(
    "NUMLOCI {length(SNP_names)}
MISSING_ALLELE 0
{SEX_COL}
POPCOLUMN_REPRO_YEARS
{SPAWN_GROUP_COL}
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
}
