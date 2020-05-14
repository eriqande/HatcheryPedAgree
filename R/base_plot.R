#' plot to compare snppit runs
#'
#'  more later
#'  @export
base_plot <- function(D,
                      by_year = TRUE,
                      type_sep_y = 0.05,
                      point_size = 1,
                      legend_point_size = 2,
                      point_stroke = 0.3,
                      line_size = 0.08) {

  # stack the values out by the type_sep_y, and make the
  # pp_comp and trio types factors so we can sort them the
  # way we would like to
  Da <- D %>%
    mutate(FDR = FDR + type_sep_y * (as.integer(type_f) - 1)) %>%
    mutate(
      parent_compatibility = factor(
        parent_compatibility,
        levels = c(
          "MaPa_not_matching_base_case",
          "pair_compatible",
          "date_mismatch",
          "sex_mismatch",
          "sex_and_date_mismatch",
          "full_trio_not_GTseq-ed",
          "Ma_or_Pa_is_potential_GT_Seq_hi_misser",
          "date_not_recorded",
          "sex_not_recorded"
        )
      )
    ) %>%
    mutate(kid_spawn_year = str_c("Kid Spawn Year = ", kid_spawn_year)) # Make more informative for panel titles

  # make a data frame with no missing points for drawing contiguous
  # lines.
  Da_nona <- Da %>%
    filter(!is.na(Pa) & !is.na(Ma) & !is.na(FDR))

  # get a data frame for putting labels on the left hand side
  label_tibble <- Da %>%
    group_by(type_f) %>%
    summarize(height = mean(type_sep_y * (as.integer(type_f) - 1)))

  # declare our color choices
  pp_comp_colors <- c(
    MaPa_not_matching_base_case = "black",
    pair_compatible = "white",
    date_mismatch = "blue",
    sex_mismatch = "orange",
    sex_and_date_mismatch = "turquoise",
    `full_trio_not_GTseq-ed` = "#FF00FF",
    Ma_or_Pa_is_potential_GT_Seq_hi_misser = "#9966FF",
    date_not_recorded = "gray",
    sex_not_recorded = "gray"
  )

  tt_colors <- c(
    C_Se_Se = "black",
    `NOT C_Se_Se` = "red"
  )

  if (by_year == FALSE) {
    g <- ggplot(Da, aes(x = idx_years_aggregated,
                        y = FDR))
  } else {
    g <- ggplot(Da, aes(x = idx_by_year,
                        y = FDR)) +
      facet_wrap(~ kid_spawn_year, ncol = 1)
  }

  g +
    geom_line(data = Da_nona,
              mapping = aes(group = type_f),
              colour = "black",
              size = line_size,
    ) +
    geom_point(
      aes(
        fill = parent_compatibility,
        colour = trio_type_group
      ),
      shape = 21,
      size = point_size,
      stroke = point_stroke
    ) +
    scale_fill_manual(values = pp_comp_colors) +
    scale_colour_manual(values = tt_colors) +
    geom_label(
      data = label_tibble,
      mapping = aes(y = height, label = type_f),
      x = 0,
      hjust = 1
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = legend_point_size)),
      colour = guide_legend(override.aes = list(size = legend_point_size))
    )
}
