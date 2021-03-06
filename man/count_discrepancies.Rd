% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_discrepancies.R
\name{count_discrepancies}
\alias{count_discrepancies}
\title{from pairs of duplicately genotyped samples, prepare a report of discrepancies}
\usage{
count_discrepancies(pairs, genotypes)
}
\arguments{
\item{pairs}{a tibble with at least two columns: \code{retained_id} and \code{original_id}.  The comparisons
are made between a single canonical individual in \code{retained_id} and any other
matching samples in \code{original_id}.  If \code{retained_id == original_id} the row is removed.}

\item{genotypes}{a tibble with columns \code{indiv}, \code{locus}, \code{gene_copy}, and \code{allele_int}.}
}
\value{
Returns a list of three components as follows:
\itemize{
\item \code{matching_samples_genos}: a tibble with 7 columns.  The genotypes of two different
individuals occupy two different rows.  The first row is the first gene copy and the
second row is the second gene copy.
\itemize{
\item \code{retained_id}: ID of the fish that is used in downstream analyses.
\item \code{original_id}: ID of the other fish whose genotype is being compared to that
of the retained ID.
\item \code{locus}: the locus name.
\item \code{gene_copy}: the gene copy index (1 or 2)
\item \code{indiv1_allele}: the alleles at the retained_id indiv.  These have been
sorted in ascending order within the locus to make it easy to compare with the indiv2allele
\item \code{indiv2_allele}: same as above, but for original_id
\item \code{num_discrepant_gene_copies}: the number of gene copies that are discrepant.  0 = none;
1 = discrepancy where one is heterozygous and the other homozygous; 2 = one individual
is homozygous and the other is homozygous for the other allele.
}
\item \code{genotype_discrepancies_summary}: counts and fractions of discrepancies across all retained_id vs
original_id pairs at a locus.  All columns should be self-explanatory except for \code{gc_wtd_fract} which
is the average number of discrepant gene copies per genotype at the locus. This might be
considered a suitable estimate of the per-allele genotyping error rate.
\item \code{alt_homoz_mismatches}: a tibble recording all the genotypes amongst the retained_id vs original_id
pairs that are discrepancies in which each member of the pair is homozygous (i.e. they are alternate
homozygotes.)  Each row denotes a mismatching locus.
}
}
\description{
After running a matching samples analysis we can use the resulting
pairs to investigate genotyping discordance rates at different
loci.  That is what this does.
}
