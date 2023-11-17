########################################################################
# Perform K-fold Cross Validation on a gene set using RWR to find the RWR rank of the left-out genes
# - Input: Pre-computed multiplex network and a geneset
# - Output: Table with the ranking of each gene in the gene set when left out, along with AUPRC and AUROC curves
# Author: David Kainer
########################################################################

parse_arguments <- function() {
  suppressPackageStartupMessages(require(optparse))

  option_list <- list(
    make_option(c("-d", "--data"),
      action = "store",
      default = NULL,
      type = "character",
      help = "The path to the .Rdata file for your combo of underlying
              functional networks. This file is produced by RWR_make_multiplex."
    ),
    make_option(c("-g", "--geneset"),
      action = "store",
      default = NULL,
      type = "character",
      help = "The path to the gene set file. It must have the following
                    first two columns with no headers tab-delimited:
                    <setid> <gene> <weight>."
    ),
    make_option(c("--method"),
      action = "store",
      default = "kfold",
      type = "character",
      help = "Cross-validation method. `kfold`, `loo`, or `singletons`.
                    [default %default]"
    ),
    make_option(c("-f", "--folds"),
      action = "store",
      default = 5,
      type = "numeric",
      help = "Number (k) of folds to use in k-fold CV.
                    [default %default]"
    ),
    make_option(c("-r", "--restart"),
      action = "store",
      default = 0.7,
      type = "numeric",
      help = "Set the restart parameter [0,1). Higher value means the
                    walker will jump back to a seed node more often.
                    [default %default]"
    ),
    make_option("--tau",
      action = "store",
      default = "1.0",
      help = "comma-separated list of values between that MUST add
              up to the number of network layers in the .Rdata file.
              One value per network layer that determines the probability
              that the random walker will restart in that layer. e.g.
              if there are three layers (A,B,C) in your multiplex network,
              then --tau '0.2,1.3,1.5' will mean that layer A is less likely
              to be walked on after a restart than layers B or C.
              [default %default]"
),
    make_option(c("-n", "--numranked"),
      action = "store",
      default = 1.0,
      type = "numeric",
      help = "proportion of ranked genes to return [0,1]. e.g. 0.1
                    will return the top 10%.
                    [default %default]"
    ),
    make_option(c("-o", "--outdir"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Path to the output directory. Both 'fullranks' and
              'meanranks' will be saved here with auto-generated filenames.
              (--out-fullranks and --out-meanranks override this.)"
    ),
    make_option(c("-m", "--modname"),
      action = "store",
      default = "default",
      type = "character",
      help = "String to include in output filename. (--out-fullranks
                    and --out-meanranks override this.)"
    ),
    make_option(c("-p", "--plot"),
      action = "store_true",
      default = FALSE,
      help = "Output plots of ROC, PRC,  etc.
                    [default %default]"
    ),
    make_option(c("--out-fullranks"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Specify the full path for full results.
                    Ignore --outdir and --modname and use this instead."
    ),
    make_option(c("--out-meanranks"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Specify the full path for mean results.
                    Ignore --outdir and --modname and use this instead."
    ),
    make_option(c("-t", "--threads"),
      action = "store",
      default = parallel::detectCores() - 1,
      type = "numeric",
      help = "Number of threads to use. Default for your system is
                    all cores - 1. [default %default]"
    ),
    make_option(c("-v", "--verbose"),
      action = "store_true",
      default = FALSE,
      help = "Verbose mode. [default %default]"
    )
  )

  opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

  errors <- 0
  # Check whether all necessary args have been set by the user.
  if (is.null(opt$data)) {
    message("ERROR:: --data is required.")
    errors <- errors + 1
  }
  if (is.null(opt$geneset)) {
    message("ERROR:: --geneset is required.")
    errors <- errors + 1
  }
  if (is.null(opt$outdir) & is.null(opt$out_fullranks)) {
    message("ERROR:: You must provide either --outdir or --out-fullranks.")
    errors <- errors + 1
  }
  if (is.null(opt$outdir) & is.null(opt$out_meanranks)) {
    message("ERROR:: You must provide either --outdir or --out-meanranks")
    errors <- errors + 1
  }

  if (opt$verbose) {
    # Do this before quitting for errors so the user can inspect args.
    print(opt)
  }

  if (errors > 0) {
    quit()
  }

  return(opt)
}


main <- function(opt) {
  # # Try to use PATH_TO_RWRtoolkit; warn the user if it's set but doesn't exist.
  # path_to_RWRtoolkit = Sys.getenv('PATH_TO_RWRtoolkit', unset=NA)
  # if (is.na(path_to_RWRtoolkit)) {
  #     path_to_RWRtoolkit = '.'
  # } else if (!dir.exists(path_to_RWRtoolkit)) {
  #     message(sprintf('WARNING: PATH_TO_RWRtoolkit is set but the directory does not exist: %s', path_to_RWRtoolkit))
  #     message('WARNING: Defaulting to local directory.')
  #     path_to_RWRtoolkit = '.'
  # }
  #
  # # Source the utils.R file; exit with error if it does not exist.
  # utils_path = file.path(path_to_RWRtoolkit, 'utils.R')
  # if (file.exists(utils_path)) {
  #     source(utils_path)
  # } else {
  #     message(sprintf('ERROR: Cannot find utils.R at path: %s', utils_path))
  #     return(1)
  # }

  opt <- parse_arguments()

  # suppressPackageStartupMessages({
  #     library(RandomWalkRestartMH)
  #     library(igraph)
  #     library(patchwork)
  #     library(dplyr)
  #     library(ggplot2)
  #     library(foreach)
  # })
  # options(dplyr.summarise.inform = FALSE)

  RWRtoolkit::RWR_CV(
    data = opt$data,
    geneset_path = opt$geneset,
    method = opt$method,
    folds = opt$folds,
    restart = opt$restart,
    tau = opt$tau,
    numranked = opt$numranked,
    outdir = opt$outdir,
    modname = opt$modname,
    plot = opt$plot,
    out_full_ranks = opt$out_fullranks,
    out_mean_ranks = opt$out_meanranks,
    threads = opt$threads,
    verbose = opt$verbose,
    write_to_file = TRUE
  )

  return(0)
}

status <- main()
quit(save = "no", status = status)
