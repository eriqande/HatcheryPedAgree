

#' prepare data for SNPPIT and run it
#'
#' This returns the directory where it ran, so that you can slurp stuff out of there.
#'
#' @param genotypes A tibble of genotypes.  Should have been reoranized with [reorganize_matching_samples()] if there were matching samples.
#' @param metadata A tibble of metadata.  Should have been reoranized with [reorganize_matching_samples()] if there were matching samples.
#' @param additional_args  A quoted string of additional (i.e. beyond the -f option) options and arguments
#' to pass to snppit.
#' @param outdir The directory to run SNPPIT in and in which to write the output.
#' @param ... Additional arguments to pass to `prepare_snppit_infile()`.  Don't pass `outf` here, though!
#' @export
run_snppit <- function(
  genotypes,
  metadata,
  additional_args = "",
  outdir = tempfile(),
  ...
) {

  dir.create(outdir, recursive = TRUE)

  outf <- file.path(outdir, "snppit_input.txt")
  message("Preparing SNPPIT data into file ", outf)

  dump <- prepare_snppit_infile(
    G = genotypes,
    S = metadata,
    outf = outf,
    ...
  )

  message("Running snppit in directory ", outdir)

  snppit_path <- system.file("bin/snppit-Darwin", package = "HatcheryPedAgree")

  call <- paste(
    "cd", outdir, ";",
    snppit_path, "-f", outf, additional_args
    )

  system(call)

  return(outdir)

}
