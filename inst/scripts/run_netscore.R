parse_arguments <- function() {
  suppressPackageStartupMessages(require(optparse))
  option_list <- list(
    make_option(
      c("-g", "--gold"),
      action = "store",
      default = NULL,
      type = "character",
      help = "full path to the gold network file. This is the network that
            contains known truth edges.  Input file is a table with first 3 cols
            <GeneA> <GeneB> <score> column names don't matter"
    ),
    make_option(
      c("-n", "--network"),
      action = "store",
      default = NULL,
      type = "character",
      help = "full path to your network to be evaluated. This is a table
            with first 3 cols <GeneA> <GeneB> <score> column names don't matter"
    ),
    make_option(
      c("-r", "--refgenes"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Only required if using --perm. Full path to a file containing
            all gene/node labels that were used in the process of generating the
            network provided in --network."
    ),
    make_option(
      c("-p", "--perm"),
      action = "store",
      default = 0,
      type = "numeric",
      help = "number of random permutations to use for estimating significance.
            If you set this, you must also set --refgenes.
            default = %default"
    ),
    make_option(
      c("-t", "--threads"),
      action = "store",
      default = 1,
      type = "numeric",
      help = "only set this is you want to do permutations.
            default = %default"
    ),
    make_option(
      c("-o", "--outdir"),
      action = "store",
      default = NULL,
      type = "character",
      help = "output directory"
    ),
    make_option(
      c("-v", "--verbose"),
      action = "store_true",
      default = FALSE,
      type = "logical",
      help = "log more stuff"
    )
  )

  desc <- "This tool performs a network intersect between your input network
    (--network) and a gold truth network (--gold). It scores the strength of
    that intersect with multiple metrics"
  opt <- parse_args(
    OptionParser(option_list = option_list, description = desc),
    convert_hyphens_to_underscores = TRUE
  )

  # Check for required arguments.
  errors <- 0
  if (is.null(opt$gold)) {
    message("ERROR:: You must provide a gold network to score against (--gold)")
    errors <- errors + 1
  }
  if (is.null(opt$network)) {
    message("ERROR:: You didn't specify a network to be scored (--network)")
    errors <- errors + 1
  }
  if (is.null(opt$outdir)) {
    message("ERROR:: You didn't specify an output directory (--outdir)")
    errors <- errors + 1
  }
  if (opt$perm > 0 & is.null(opt$refgenes)) {
    message("ERROR:: You provided --perm but did not provide --refgenes.")
    errors <- errors + 1
  }

  # Print the arguments for the user.
  if (opt$verbose) {
    message("You passed the following arguments:")
    for (name in names(opt)) {
      message(sprintf("%s = %s \n", name, opt[[name]]))
    }
  }

  # Exit if there are errors in the arguments.
  if (errors > 0) {
    message("Found errors in your arguments (see above). Use --verbose to inspect your arguments.")
    quit(status = 1)
  }

  return(opt)
}


main <- function() {
  opt <- parse_arguments()
  RWRtoolkit::RWR_netscore(
    gold = opt$gold,
    network = opt$network,
    reference_geneset = opt$refgenes,
    perm = opt$perm,
    threads = opt$threads,
    outdir = opt$outdir,
    verbose = opt$verbose,
    write_to_file = TRUE
  )
  return(0)
}


status <- main()
quit(status = status)
