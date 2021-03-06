% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reorganize_matching_samples.R
\name{reorganize_matching_samples}
\alias{reorganize_matching_samples}
\title{Relabel matching genotypes to a single ID and reorganize their years and spawner groups}
\usage{
reorganize_matching_samples(genotypes, metadata, clusters)
}
\arguments{
\item{genotypes}{the genotypes}

\item{metadata}{the metadata}

\item{clusters}{the tibble with the matching sample clusters.}
}
\value{
A list with a variety of components:
\itemize{
\item \code{matchers_metadata}: a tibble of all the matching sample clusters and their metadata
\item \code{snppit_meta}: a tibble of metadata in long format that is preped for use in \code{\link[=prepare_snppit_infile]{prepare_snppit_infile()}}
\item \code{snppit_genos}: a tibble of genotype data in long format that is preped for use in \code{\link[=prepare_snppit_infile]{prepare_snppit_infile()}}
\item \code{cross_hatchery_matches}: a tibble of matching samples that were spawned at more than one hatchery
\item \code{cross_sex_matches}: a tibble of matching samples that are marked as both sexes
\item \code{geno_discord}:  a list like that returned in \code{\link[=count_discrepancies]{count_discrepancies()}}.
}
}
\description{
If there are matching samples, we want to make sure that we have just a single
ID for each cluster, and that the years and spawner groups get organized
so that SNPPIT can read all that information in.
}
