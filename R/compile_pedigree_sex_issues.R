
#' Compile the connected components of parent pairs with sex conflicts or missing data
#'
#' After a pedigree is inferred, some of the parents might both be of the same sex.
#' This function looks through the pedigree and first warns the user if a pair is
#' compatible (opposite sexes) but the pa is listed as Male or the ma is listed as Female.
#' Then it finds all parent pairs whose members have the same reported sex and it
#' finds the connected components of all of these, and it makes a plot of those as well.
#' It does the same thing for fish with no reported sex.
#' @param Ped A tibble that has the inferred pedigree in it.
#' **This should already be filtered using whatever criteria you want to use (FDR, etc.)**
#' This will typically be the output from `reformat_no_sex_of_date_results()`. Whatever it is,
#' it must have the columns `kid`, `pa`, `ma`, `pa_sex`, `ma_sex`, `pa_year`, and `ma_year`
#' @param strip_pattern A pattern that should be removed from IDs.  This is
#' used for the case where certain samples appear in multiple hatcheries and the occurrence
#' in one of the hatcheries was denoted by adding a suffix to the ID.  The default is set
#' up for that as an example
#' @param scc_pointsize The number passed into ggraph for the size of node points for the sex-clashing clusters
#' plot.  For the example data, this is set as 7 so that plot works well when printed at 30 inches by 30 inches
#' as a PDF.
#' @param suc_pointsize The number passed into ggraph for the size of node points for the sex-unknown clusters
#' plot.  For the example data, this is set as 7 so that plot works well when printed at 30 inches by 30 inches
#' as a PDF.
#' @return Returns a list with five components:
#' - `parent_pairs_graph`: a tidygraph object that holds the adjencies between parents and some of their meta data
#' - `sex_clash_clusters`: the portion of the tidygraph that has all the information about the connected
#' components that include at least one sex-clashing parent pair.
#' - `sex_unknown_clusters`: the portion of the tidygraph that has all the information about the connected
#' components that include at least one individual with sex unknown (i.e. not recorded in the data).
#' - `sex_clash_clusters_plot`: a faceted ggplot object.  Each panel shows one of the connected components with a
#' sex-clashing pair.
#' - `sex_unknown_clusters_plot`: a faceted ggplot object.  Each panel shows one of the connected components with a
#' sex-unknown member.
#' @examples
#' # the package comes with a pedigree for this example:
#' ped <- ped_with_sex_issues
#' issues <- compile_pedigree_sex_issues(ped)
#'
#' # to look at the graph you would do like:
#' \dontrun{
#' # sex clashing clusters
#' ggsave(
#'     issues$sex_clash_clusters_plot,
#'     filename = "clash.pdf",
#'     width = 30,
#'     height = 30
#' )
#' # sex unknown clusters
#' ggsave(
#'     issues$sex_unknown_clusters_plot,
#'     filename = "unknown.pdf",
#'     width = 15,
#'     height = 15
#' )
#' }
#' @export
compile_pedigree_sex_issues <- function(
  Ped,
  strip_pattern = "_Coyote-Valley$|_Warm-Springs$",
  scc_pointsize = 7,
  suc_pointsize = 7
  ) {

  # first check to make sure we have all the right columns
  miss_col <- setdiff(
    c("kid", "pa", "ma", "pa_sex", "ma_sex", "pa_year", "ma_year"),
    names(Ped)
  )
  if(length(miss_col) > 0) {
    stop("Ped seems to be missing columns: ", paste(miss_col, collapse = ", "))
  }

  # now, give warning if sex labels don't match ma and pa, but otherwise don't conflict
  test <- Ped %>%
    filter(
      (ma_sex == "Male" & pa_sex == "Female") |
        (ma_sex == "Male" & is.na(pa_sex)) |
        (pa_sex == "Female" & is.na(ma_sex))
      )
  if(nrow(test) > 0) {
    warning("Parent sex specifications in Ped do not accord with parent type (ma or pa) in cases that don't have a sex conflict.  Did you forget to run this through reformat_no_sex_or_date_results()?")
  }

  #### now we are ready to find the connected components ####

  # first, get the canonical ordering for the parent pairs
  cano_rr <- Ped %>%
    mutate(
      tmp1 = str_replace_all(pa, strip_pattern, ""),
      tmp2 = str_replace_all(ma, strip_pattern, ""),
    ) %>%
    mutate(
      cano1 = case_when(
        tmp1 < tmp2 ~ tmp1,
        TRUE ~ tmp2
      ),
      cano2 = case_when(
        tmp1 < tmp2 ~ tmp2,
        TRUE ~ tmp1
      ),
      sex1 = case_when(
        tmp1 < tmp2 ~ pa_sex,
        TRUE ~ ma_sex
      ),
      sex2 = case_when(
        tmp1 < tmp2 ~ ma_sex,
        TRUE ~ pa_sex
      )
    ) %>%
    select(cano1, cano2, sex1, sex2, everything())

  # enumerate the parent pairs
  pairs <- cano_rr %>%
    count(cano1, cano2)

  # make a graph and identify every parent's the connected component
  conn_comps <- pairs %>%
    rename(from = cano1, to = cano2) %>%
    select(to, from) %>%
    igraph::graph_from_data_frame() %>%
    igraph::components() %>%
    .$membership %>%
    enframe() %>%
    rename(indiv = name, cluster = value) %>%
    arrange(cluster, indiv) %>%
    group_by(cluster) %>%
    mutate(cluster_size = n()) %>%
    ungroup()

  # make a tidygraph of all the parent-pair connections, and then also join the
  # sex and spawn year for each fish on there
  tgraph <- pairs %>%
    rename(from = cano1, to = cano2) %>%
    select(to, from) %>%
    igraph::graph_from_data_frame() %>%
    as_tbl_graph()

  # get sex info for each fish from the info in the pedigree
  sex_info <- bind_rows(
    Ped %>% select(kid, kid_sex) %>% rename(indiv = kid, sex = kid_sex),
    Ped %>% select(ma, ma_sex)  %>% rename(indiv = ma, sex = ma_sex),
    Ped %>% select(pa, pa_sex) %>% rename(indiv = pa, sex = pa_sex)
  ) %>%
    count(indiv, sex) %>%
    select(-n)

  # check to make sure that no individual has multiple sexes
  multi_sex_check <- sex_info %>%
    count(indiv) %>%
    filter(n > 1)
  if(nrow(multi_sex_check) > 0) {
    stop("Problem: these individuals appear to have multiple listed sexes: ", paste(multi_sex_check$indiv, collapse = ", "))
  }

  # also get the spawn info
  spawn_info <- bind_rows(
    Ped %>% select(ma, ma_year)  %>% rename(indiv = ma, spawn_year = ma_year),
    Ped %>% select(pa, pa_year) %>% rename(indiv = pa, spawn_year = pa_year)
  ) %>%
    count(indiv, spawn_year) %>%
    select(-n)

  # Now, joint that sex and spawn information on the nodes in the tidygraph
  tgraph2 <- tgraph %>%
    left_join(sex_info, by = c("name" = "indiv")) %>%
    left_join(spawn_info, by = c("name" = "indiv")) %>%
    left_join(conn_comps, by  = c("name" = "indiv"))

  # Now, find the IDs of the sex-pair clashers
  tmp <- Ped %>%
    filter(ma_sex == pa_sex | ma_sex == "Male" | pa_sex == "Female")
  sex_pair_clashers <- tibble(
    indiv = unique(c(tmp$pa, tmp$ma)),
    sex_clasher = TRUE
  )

  # join those on
  tgraph3 <- tgraph2 %>%
    left_join(sex_pair_clashers, by = c("name" = "indiv")) %>%
    mutate(sex_clasher = ifelse(is.na(sex_clasher), FALSE, sex_clasher))

  # now, filter down to connected components that include at least one sex-clasher
  sex_clash_clusters <- tgraph3 %>%
    group_by(cluster) %>%
    filter(any(sex_clasher == TRUE)) %>%
    ungroup()

  # and also filter down to connected components that involve some sex-unkown indivs
  sex_unknown_clusters <- tgraph3 %>%
    group_by(cluster) %>%
    filter(any(is.na(sex))) %>%
    ungroup()


  # now, make some plots
  sex_clash_clusters_plot <- ggraph(sex_clash_clusters, layout = "kk") +
    geom_edge_link() +
    geom_node_point(aes(shape = sex, fill = sex_clasher, colour = sex_clasher), size = scc_pointsize) +
    geom_node_point(aes(shape = sex), fill = NA, colour = "black", size = scc_pointsize) +
    geom_node_text(aes(label = name), repel = TRUE) +
    facet_nodes(~ cluster, scales = "free") +
    scale_shape_manual(values = c(Male = 22, Female = 21), na.value = 23)

  sex_unknown_clusters_plot <- ggraph(sex_unknown_clusters, layout = "kk") +
    geom_edge_link() +
    geom_node_point(aes(shape = sex, fill = sex), size = suc_pointsize) +
    geom_node_point(aes(shape = sex), fill = NA, colour = "black", size = suc_pointsize) +
    geom_node_text(aes(label = name), repel = TRUE) +
    facet_nodes(~ cluster, scales = "free") +
    scale_shape_manual(values = c(Male = 22, Female = 21), na.value = 23) +
    scale_fill_manual(values = c(Male = "skyblue", Female = "pink"), na.value = "white")


  # return a number of things in a list:
  list(
    parent_pairs_graph = tgraph3,
    sex_clash_clusters = sex_clash_clusters,
    sex_unknown_clusters = sex_unknown_clusters,
    sex_clash_clusters_plot = sex_clash_clusters_plot,
    sex_unknown_clusters_plot = sex_unknown_clusters_plot
  )

}
