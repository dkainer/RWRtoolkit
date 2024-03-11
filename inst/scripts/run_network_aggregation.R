########################################################################
# Parse arguments
########################################################################

parse_arguments <- function() {
  # [TODO] Need to verify 'help' messages (I ripped several of these from
  # other scripts followed by liberal copy/pasta).
  suppressPackageStartupMessages(library(optparse))
  option_list = list(
    make_option(c("-d", "--data"),
      action = "store",
      default = NULL,
      type = "character",
      help = "path to the .Rdata file for your combo of
        underlying functional networks. This file is produced
        by RWR_make_MHobject.R"
    ),
    make_option(
      c("-o", "--outdir"),
      action = "store",
      default = "netstats",
      type = "character",
      help = 'Output file name (default "netstats")'
    ),
    make_option(
      c("--merged_with_all_edges"),
      action = "store_true",
      default = FALSE,
      help = "A boolean denoting a return of a merged down multiplex network along with  network edge counts and vertex counts."
    ),
    make_option(
      c("--merged_with_edgecounts"),
      action = "store_true",
      default = FALSE,
      help = "A boolean denoting a return of a merged down multiplex, but simplified with edge weights denoting the total number of layers in which that edge existed."
    ),
    make_option(
      c("-v", "--verbose"),
      action = "store_true",
      default = FALSE,
      help = "A boolean denoting the verbosity of output."
    )
  )

  opt <- optparse::parse_args(
    optparse::OptionParser(
      option_list = option_list
    )
    )

  # Verify opt.
  errors <- FALSE

  # Require that the user passes either --data or --flist (not both).
  if (!is.null(opt$data) & !is.null(opt$flist)) {
  errors <- TRUE
  message(
    "\n[ERROR] You may provide either --data or --flist but not both.\n"
  )
  }

  # Require that the user provided at least one input.
  if (is.null(opt$data) ) {

  errors <- TRUE
  message(
    "\n[ERROR] You must supply one of the following arguments:\n",
    "  --data\n"
  )}

  if (errors) {
    stop("Exiting due to errors (see above).")
  }

  return(opt)
}

########################################################################
# Main
########################################################################

main <- function() {
  # Get args.
  
  opt <- parse_arguments()
  print("RUNNING NETSTATS WITH ")
  print(opt)
  RWRtoolkit::RWR_network_aggregation()(
    data = opt$data,
    merged_with_edgecounts = opt$merged_with_edgecounts,
    merged_with_all_edges = opt$merged_with_all_edges,
    verbose = opt$verbose
  )
}

status  <-  main()
quit(save = "no", status = status)
