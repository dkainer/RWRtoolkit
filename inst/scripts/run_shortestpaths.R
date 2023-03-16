########################################################################
# RWR_shortestpaths.
# Find shortest paths between genes in the given gene sets.
########################################################################

library(RWRtoolkit)

## Parse args function
parse_arguments <- function() {
  suppressPackageStartupMessages(require(optparse))
  option_list <- list(
    make_option(
      c("-d", "--data"),
      action = "store",
      default = NULL,
      type = "character",
      help = "path to the .Rdata file for your combo of underlying functional
            networks. This file is produced by RWR_make_MHobject.R"
    ),
    make_option(
      c("-g", "--source_geneset"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Path to the gene set file. If given with --target_geneset,
            find shortest paths between genes in --source_geneset and --target_geneset.
            Otherwise, find shortest paths among genes in --source_geneset only.
            It must have the following cols without heading: <setid> <gene>"
    ),
    make_option(
      c("-p", "--target_geneset"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Path to the second geneset file. If it is not provided then
            the shortest paths will be calculated between genes in source_geneset"
    ),
    make_option(
      c("-o", "--outdir"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Full path to the output file directory. Two output files will
            be generated with different suffixes"
    ),
    make_option(
      c("--out-path"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Specify the full path for output. Ignore --outdir and --modname
            and use this instead."
    ),
    make_option(
      c("-t", "--threads"),
      action = "store",
      default = 1,
      type = "numeric",
      help = "Number of threads to use. default [default %default]"
    ),
    make_option(
      c("-c", "--cyto"),
      action = "store_true",
      default = FALSE,
      help = "Include this parameter if you wish to see the shortest paths in
                    Cytoscape. Cytoscape must already be running first!
                    If your geneset(s) are large (e.g. 100+) then loading all the
                    shortest paths fully into cytoscape can take a while."
    ),
    make_option(
      c("-v", "--verbose"),
      action = "store_true",
      default = FALSE,
      help = "Log more stuff."
    )
  )
  desc <- paste(
    "Find shortest paths between genes. These will be output as a",
    "table of edges. There are two ways to use this:",
    "\n  1) provide one geneset (g) and you will get all shortest paths between genes in g.",
    "\n  2) provide two genesets (g and p) and you will get the shortest paths between genes from g to p"
  )

  opt <- parse_args(
    OptionParser(
      option_list = option_list,
      description = desc
    ),
    convert_hyphens_to_underscores = TRUE
  )

  if (opt$verbose) {
    print(opt)
  }

  # Check required args.
  errors <- 0
  if (is.null(opt$data)) {
    message("ERROR:: You must provide an Rdata file for your networks.")
    errors <- errors + 1
  }
  if (is.null(opt$source_geneset)) {
    message("ERROR:: You must provide --source_geneset.")
    errors <- errors + 1
  }
  if (is.null(opt$outdir) && is.null(opt$out_path)) {
    message("ERROR:: You must provide either --outdir or --out-path.")
    errors <- errors + 1
  }

  if (errors > 0) {
    quit()
  }

  if (opt$cyto == TRUE) {
    ping <- NULL
    ping <- RCy3::cytoscapePing()
    if (is.null(ping)) {
      message("ERROR:: You set --cyto=TRUE but Cytoscape is not running. Open cytoscape then try again.")
    }
  }
  return(opt)
}

main <- function() {
  ## Processing arguments
print("parse options")
  opt <- parse_arguments()

    print("OPTINOS")
    print(opt)
  ## Call to ShortestPaths
  RWRtoolkit::RWR_ShortestPaths(
    data = opt$data,
    source_geneset = opt$source_geneset,
    target_geneset = opt$target_geneset,
    outdir = opt$outdir,
    out_path = opt$out_path,
    threads = opt$threads,
    cyto = opt$cyto,
    verbose = opt$verbose,
    write_to_file = TRUE
  )
  return(0)
}

status <- main()
quit(save = "no", status = status)
