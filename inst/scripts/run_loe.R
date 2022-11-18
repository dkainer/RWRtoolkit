library(RWRtoolkit)

## Parse args function
parse_arguments <- function() {
  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option(c("-d", "--data"),
      action = "store",
      default = NULL,
      type = "character",
      help = "path to the .Rdata file for your combo of
                   underlying functional networks. This file is produced
                   by RWR_make_MHobject.R"
    ),
    make_option(c("-g", "--seed_geneset"),
      action = "store",
      default = NULL,
      type = "character",
      help = "path to the first gene set file. This file must
                   have the following first two cols without heading:
                   <setid> <gene> If you wish to include weights for the
                   genes then the third col should be numeric: <setid>
                   <gene> <weight> Note, the weights can be on any
                   numeric scale (they will be normalised) but should
                   all be > 0"
    ),
    make_option(c("-p", "--query_geneset"),
      action = "store",
      default = NULL,
      type = "character",
      help = "<optional> path to the second (target) geneset
                   file.   If it is not provided then the seeds from
                   geneset1 will just be used in a context analysis"
    ),
    make_option(c("-r", "--restart"),
      action = "store",
      default = 0.7,
      type = "numeric",
      help = "set the restart parameter. default [default %default]"
    ),
    make_option("--tau",
      action = "store",
      default = "1.0",
      help = "comma-separated list of values between that
                   MUST add up to the number of network layers in the
                   .Rdata file.  One value per network layer that
                   determines the probability that the random walker
                   will restart in that layer.  e.g. if there are three
                   layers (A,B,C) in your multiplex network, then --tau
                   '0.2,1.3,1.5' will mean that layer A is less likely
                   to be walked on after a restart than layers B or
                   C."
    ),
    make_option(c("-o", "--outdir"),
      action = "store",
      default = NULL,
      type = "character",
      help = "full path to the output file directory. Two
                   output files will be generated with different
                   suffixes"
    ),
    make_option(c("-n", "--numranked"),
      action = "store",
      default = 1.0,
      type = "numeric",
      help = "proportion of ranked genes to return. e.g. 0.1
                   will return the top 10%. default [default %default]"
    ),
    make_option(c("-e", "--eval"),
      action = "store_true",
      default = FALSE,
      help = "include this parameter if you want to output a
                   PNG plot of ROC and PRC. [default %default]"
    ),
    make_option(c("-m", "--modname"),
      action = "store",
      default = "default",
      type = "character",
      help = "alias for this run. Useful for output."
    ),
    make_option(c("-c", "--cyto"),
      action = "store",
      default = 0,
      type = "numeric",
      help = "Specify a number N > 0 if you wish to see a
                   pretty network of the seeds and top N ranked genes.
                   Cytoscape must already be running first! [default %default]"
    ),
    make_option(c("-v", "--verbose"),
      action = "store_true",
      default = FALSE,
      type = "logical",
      help = "Print additional logs."
    )
  )

  desc <- paste(
    "There are two main ways to use RWR_LOE: 1) provide",
    "one geneset (g) and you will get all genes in the",
    "network ranked according to connectivity to genes in",
    "g. 2) provide two genesets (g and p) and you will get",
    "the genes in p ranked according to connectivity to",
    "genes in g, including optional precision/recall evaluation"
  )
  opt <- parse_args(
    OptionParser(
      option_list = option_list,
      description = desc
    ),
    convert_hyphens_to_underscores = TRUE
  )
  return(opt)
}

## Processing arguments
opt <- parse_arguments()
print(opt)

## Call to LOE
RWRtoolkit::RWR_LOE(data = opt$data, seed_geneset = opt$seed_geneset, query_geneset = opt$query_geneset, restart = opt$restart, tau = opt$tau, outdir = opt$outdir, numranked = opt$numranked, eval = opt$eval, modname = opt$modname, cyto = opt$cyto, verbose = opt$verbose)
