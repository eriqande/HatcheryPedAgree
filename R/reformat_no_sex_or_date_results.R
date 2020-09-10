

#' After slurping in a SNPPIT run done without sex or date constraints, harmonize the declared sex and dates with the output
#'
#' The main thing going on here is to make sure that the pa column are males and the ma column are females
#' and to declare a reported sex mismatch if they are both males or females.
#' @param D the output of  [slurp_snppit()].
#' @export
reformat_no_sex_or_date_results <- function(D) {

  # create a column that has the index for each trio
  D2 <- D %>%
    ungroup() %>%
    mutate(trio_index = 1:n()) %>%
    select(trio_index, everything())

  # first break D into reported trios and non-reported trios
  trios <- D2 %>%
    filter(!is.na(ma) & !is.na(pa))

  NA_trios <- D2 %>%
    filter(is.na(ma) | is.na(pa))

  # now, we break trios into the "active" columns (which we are going to be changing)
  # and the "inactive" ones, which we will join back on via the trio index later.
  t_active <- trios %>%
    select(trio_index:ma_hatchery) %>%
    rename(
      kid_id = kid,
      pa_id = pa,
      ma_id = ma
    ) %>%
    mutate(
      kid_year = as.character(kid_year),
      ma_year = as.character(ma_year),
      pa_year = as.character(pa_year)
    )
  t_inactive <- trios %>%
    select(-(kid:ma_hatchery))

  # now, we do some pivot operations so that each trio is stored as three
  # rows in order: kid, ma, pa; and the sex, year, spawner group, and hatchery of
  # each member of the trio is in a separate column.
  # Perhaps there is a way to do this as a single pivot_longer, but I
  # ended up doing a pivot_longer followed by a pivot wider
  ta2 <- t_active %>%
    pivot_longer(
      cols = -c(trio_index, SpawnYear),
      names_to = c("member", "variable"),
      names_sep = "_"
    ) %>%
    pivot_wider(
      names_from = variable,
      values_from = value
    ) %>%
    arrange(trio_index, member)

  # now, if we group by trio index, we have three rows: kid, ma, pa.  And now
  # we look at the sex column to see if the ma and pa labels need to be reassigned.
  ta3 <- ta2 %>%
    group_by(trio_index) %>%
    mutate(
      member = {
        ret <- c("kid", "ma", "pa")
        if ( (!is.na(sex[2]) && !is.na(sex[3])) &&    # if both are not NA and sexes are reversed
             (sex[2] == "Male" && sex[3] == "Female")
        ) {
          ret <- c("kid", "pa", "ma")
        } else if (!is.na(sex[3]) && sex[3] == "Female" && is.na(sex[2])) {  # if pa's sex is not NA and Female, and ma's sex is NA,
          ret <- c("kid", "pa", "ma")
        } else if (!is.na(sex[2]) && sex[2] == "Male" && is.na(sex[3])) {  # if ma's sex is not NA and Male, and pa's sex is NA
          ret <- c("kid", "pa", "ma")
        }

        ret
      }
    ) %>%
    ungroup()

  # now, we pivot this wider
  ta4 <- ta3 %>%
    pivot_wider(
      names_from = member,
      values_from = id:hatchery,
      names_sep = "_"
    )

  # and deal with the names.  This is a bit ugly and will break
  # when there are more extra columns in the front there...
  tmpmat <- names(ta4)[-(1:2)] %>%
    str_split_fixed(., "_", 2)
  names(ta4)[-(1:2)] <- paste(tmpmat[, 2], tmpmat[, 1], sep = "_")

  # rename a few things and join the inactive columns back on there
  ta_done <- ta4 %>%
    rename(
      kid = kid_id,
      pa = pa_id,
      ma = ma_id
    ) %>%
    select(trio_index, kid, pa, ma, everything()) %>%
    left_join(t_inactive, by = "trio_index")

  # then return those along with the NA columns
  bind_rows(
    NA_trios,
    ta_done
  ) %>%
    arrange(trio_index) %>%
    select(-trio_index)


}
