
# internal function to pick out the meaning of a single bit
pick_bit <- function(zero, one, bit_as_int, V) {
  c(zero, one)[any(as.integer(intToBits(bit_as_int) & intToBits(V))) + 1]
}

# internal function to print out the meanings of the different bit-masks
categories_from_bit_mask_integer <- function(BI, FDR1, FDR2) {
  if(any(as.integer(intToBits(1) & intToBits(BI)))) {
    Run1 <- "Run1: parents absent"
  } else {
    Run1 <- paste0(
      "Run1: parents present; MaxPostTrio ",
      pick_bit("is", "is NOT", 4, BI),
      " C_Se_Se; FDR ",
      pick_bit("<= ", "> ", 16, BI),
      FDR1,
      collapse = ""
    )
  }

  if(any(as.integer(intToBits(2) & intToBits(BI)))) {
    Run2 <- "Run2: parents absent"
  } else {
    Run2 <- paste0(
      "Run2: parents ",
      pick_bit("present", "absent", 2, BI),
      "; MaxPostTrio ",
      pick_bit("is", "is NOT", 8, BI),
      " C_Se_Se; FDR ",
      pick_bit("<= ", "> ", 32, BI),
      FDR2,
      collapse = ""
    )
  }

  # finally, get the number of shared parents.  This will either
  # be 0, 1, 2, or NA (it is NA if either run did not identify any parents.)
  if (any(as.integer((intToBits(1) | intToBits(2)) & intToBits(BI)))) {
    num_shared <- NA
  } else if (any(as.integer(intToBits(64) & intToBits(BI)))) {
    num_shared <- 0
  } else if (any(as.integer(intToBits(128) & intToBits(BI)))) {
    num_shared <- 1
  } else if (!any(as.integer(intToBits(64 + 128) & intToBits(BI)))) {
    num_shared <- 2
  }

  Shared <- paste0(
    "Number of shared parents = ",
    num_shared,
    collapse = ""
  )

  # finally return a string:
  paste0(
    BI, "|",
    Run1,
    "|",
    Run2,
    "|",
    Shared
  )
}


#' Partition the results of two SNPPIT runs by comparisons of their features
#'
#' This function is intended to break down the results of two SNPPIT runs
#' into non-overlapping groups of kids.  The primary purpose for this is
#' to compare runs with sex-and-date to runs without the sex and date.  After
#' doing this you can use a later function to pick out the trio from either
#' run1 or run2 according to which category it is in.
#'
#' This is implemented in terms of bit-masks.  Here are the bits, given by the
#' "1" value.  We set this up so that good trios shared by both runs will be 0s.
#' - 1: Parents not identified (NA) in Run1
#' - 2: Parents not identified (NA) in Run2
#' - 4: Max posterior trio in Run1 is not C_Se_Se
#' - 8: Max posterior trio in Run2 is not C_Se_Se
#' - 16: FDR > FDR1 in Run1
#' - 32: FDR > FDR2 in Run2
#' - 64: Run1 and Run2 parent pairs share exactly 0 individuals
#' - 128: Run1 and Run2 parent pairs share exactly 1 individual
#'
#' If neither of bits 1, 2, 64, or 128 are set, then there were two parents matching between the groups.
#' @param Run1 the output of [slurp_snppit()] for the first (typically with sex and date) snppit run.
#' @param Run2 the output of [slurp_snppit()] or [reformat_no_sex_or_date_results()] for the second
#' (typically _without_ sex and date) snppit run.
#' @param FDR1 the FDR threshold to consider for Run1 (defaults to 0.01).
#' @param FDR2 the FDR threshold to consider for Run2 (defaults to 0.01).
#' @export
partition_two_snppit_results <- function(
  Run1,
  Run2,
  FDR1 = 0.01,
  FDR2 = 0.01
) {
  # do a full join on both snppit results. Join on kid and name all remaining columns with a _1 and _2
  # and only retain the columns that we need to make the bitmasks
  JJ <- full_join(
    Run1,
    Run2,
    by = "kid",
    suffix = c("_1", "_2")
  ) %>%
    select(
      kid,
      pa_1,
      pa_2,
      ma_1,
      ma_2,
      MaxP.Pr.Relat_1,
      MaxP.Pr.Relat_2,
      FDR_1,
      FDR_2
    )

  # Now, we can quite simply form the bitmasks.  This was harder than I
  # thought it would be because the & does not short-circuit evaluate the
  # way I thought it would, so NAs are a hassle
  JJ2 <- JJ %>%
    mutate(
      bit_category = case_when(
        (is.na(pa_1) | is.na(ma_1)) & (is.na(pa_2) | is.na(ma_2)) ~ 3L,
        (is.na(pa_1) | is.na(ma_1)) ~ 1L,
        (is.na(pa_2) | is.na(ma_2)) ~ 2L,
        TRUE ~ 4L * (!is.na(MaxP.Pr.Relat_1) & (MaxP.Pr.Relat_1 != "C_Se_Se")) +
          8L * (!is.na(MaxP.Pr.Relat_2) & (MaxP.Pr.Relat_2 != "C_Se_Se")) +
          16L * (!is.na(FDR1) & (FDR_1 > FDR1)) +
          32L * (!is.na(FDR2) & (FDR_2 > FDR2)) +
          64L * (
            (!is.na(pa_1) & !is.na(ma_1) & !is.na(pa_2) & !is.na(ma_2))
            & (( (pa_1 == pa_2) + (pa_1 == ma_2) + (ma_1 == pa_2) + (ma_1 == ma_2)) == 0)
          ) +
          128L * (
            (!is.na(pa_1) & !is.na(ma_1) & !is.na(pa_2) & !is.na(ma_2))
            & (( (pa_1 == pa_2) + (pa_1 == ma_2) + (ma_1 == pa_2) + (ma_1 == ma_2)) == 1)
          )
      )
    )

  JJ2 %>%
    mutate(
      verbal_category = purrr::map_chr(
        .x = bit_category,
        .f = categories_from_bit_mask_integer,
        FDR1 = FDR1,
        FDR2 = FDR2
      )
    )

}
