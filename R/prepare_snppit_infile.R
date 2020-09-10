#' prepare long-format genotypes for snppit
#'
#' The meta data, S, determines who gets put into the snppit infile via a left
#' join.  So, if it turns out that we are requsting the addition of individuals
#' that are not in `G`, this will throw an error.
#'
#' This function is called within [run_snippit()].
#'
#' @param G a long format data frame of SNPs.  It must have, at a minimum, the columns
#'   `indiv`, `locus`, `gene_copy`, `allele_int`.  Missing data should be represented as
#'   NA.
#' @param S a data frame of meta data with the columns `indiv`, `year`, `sex`, `spawner_group`,
#' and `hatchery`
#' @return Function does not return a value of any importance.
#' The side effect of the function is to write out a snppit input file!
#' @export
#' @examples
#' prepare_snppit_infile(G = coho_genotypes, S = coho_metadata)
prepare_snppit_infile <- function(G,
                                  S,
                                  use_spawner_group = TRUE,
                                  use_sex = TRUE,
                                  min_age = 1,
                                  max_age = 6,
                                  geno_err = 0.005,
                                  outf = tempfile()) {

  # these will get changed if we don't use them
  sex_col <- "POPCOLUMN_SEX"
  spawn_group_col <- "POPCOLUMN_SPAWN_GROUP"


  if(any(is.na(S$year))) {
    stop("Missing data not allowed in the year column of S")
  }

  # first, make a genotype matrix. Put spaces between the gene copies in a locus
  M <- G %>%
    select(indiv, locus, gene_copy, allele_int) %>%
    mutate(allele_int = ifelse(is.na(allele_int), 0, allele_int)) %>%
    spread(key = gene_copy, value = allele_int) %>%
    unite(col = "geno", `1`, `2`, sep = " ") %>%
    spread(key = locus, value = geno)

  # get the marker names
  SNP_names <- colnames(M)[-1]

  # prep the sex and spawn date data, and also deal with years, which could
  # be a character string.
  meta <- S %>%
    mutate(
      sex = recode(sex, Female = "F", Male = "M"),
      sex = ifelse(is.na(sex), "?", sex),
    ) %>%
    arrange(year, spawner_group) %>%
    select(hatchery, indiv, sex, year, spawner_group) %>%
    mutate(
      years_list = str_split(year, ","),
      min_year = map_int(years_list, function(x) min(as.integer(x))),
      max_year = map_int(years_list, function(x) max(as.integer(x)))
      )

  # make sure genotype data are available for all the individuals in the metadata
  these_have_no_genos <- meta %>%
    anti_join(M, by = "indiv")
  if (nrow(these_have_no_genos) > 0) {
    stop(
      "Error.  S in prepare_snppit_file is requesting addition of fish with no genotypes in G: ",
      paste(these_have_no_genos$indiv, collapse = ", ")
    )
  }

  # break the genos into candidate parents and offspring.
  # Basically if your latest spawn year is less than the min year + the
  # min age, then you don't get to be an offspring.  Likewise, if
  # your earliest spawn year is greater than max year - min age, you don't get
  # to be a candidate parent
  Offs <- meta %>%
    filter(max_year >= min(min_year) + min_age) %>%
    select(-sex, -spawner_group) %>%
    mutate(age_range = str_c(min_age, "-", max_age)) %>%
    left_join(., M, by = "indiv")

  Pars <- meta %>%
    filter(min_year <= max(max_year) - min_age) %>%
    left_join(., M, by = "indiv")


  if (use_sex == FALSE) {
    sex_col <- ""
    Pars <- Pars %>% select(-sex)
  }
  if (use_spawner_group == FALSE) {
    spawn_group_col <- ""
    Pars <- Pars %>% select(-spawner_group)
  }


  # now assemble the text ouput
  preamble <- glue::glue("
    NUMLOCI {length(SNP_names)}
    MISSING_ALLELE 0
    {sex_col}
    POPCOLUMN_REPRO_YEARS
    {spawn_group_col}
    OFFSPRINGCOLUMN_SAMPLE_YEAR
    OFFSPRINGCOLUMN_AGE_AT_SAMPLING\n
    ")

  cat(preamble, file = outf)
  write.table(
    cbind(SNP_names, geno_err),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    file = outf,
    append = TRUE
  )

  # cycle over the different hatcheries that might be involved
  # and put in a block of parents for each
  dump <- lapply(split(Pars, Pars$hatchery), function(x) {
    cat("POP ", x$hatchery[1], "\n", file = outf, append = TRUE)
    write.table(select(x, -hatchery, -years_list, -min_year, -max_year),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      file = outf, append = TRUE, sep = "\t"
    )
  })



  # then do the same thing for offspring from the different hatcheries.  And,
  # for now, assume that the offspring could have come from any of the hatcheries
  # (hence the ? after the Offspring name)
  dump <- lapply(split(Offs, Offs$hatchery), function(x) {
    cat("OFFSPRING Offspring-", x$hatchery[1], " ?\n", file = outf, append = TRUE, sep = "")
    write.table(select(x, -hatchery, -years_list, -min_year, -max_year),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      file = outf, append = TRUE, sep = "\t"
    )
  })

  # message to user
  message("snppit input file written to ", outf)
}
