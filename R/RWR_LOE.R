###############################################################################
# Perform RWR from one geneset and returns either (1) ranking of the rest of
# the genes in the network or (2) ranking of genes in a provided second geneset.
# - Input: Pre-computed multiplex network and at least one geneset of seeds.
# - Output: Table with the ranking of genes from (1) network or (2) second
#           geneset
# Copyright (C) 2022  David Kainer
#
# This file is part of RWRtoolkit.
#
# RWRtoolkit is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# RWRtoolkit is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# RWRtoolkit. If not, see <https://www.gnu.org/licenses/>.
###############################################################################

#' @importFrom dplyr %>%
########################################################################
# Internal Functions
########################################################################

#' @keywords internal
calc_metrics_loe <- function(res, seed_geneset, query_geneset) {
  # PR example: 5 genes left out need to be detected from 10 rankings
  # ranking:  1    2    3    4   5   6   7    8    9     10
  # truth:    1    0    1    0   0   1   0    0    1     1
  # recall:   0.2  0.2  0.4  0.4 0.4 0.6 0.6  0.6  0.8   1.0
  # Prec:     1    0.5  0.67 0.5 0.4 0.5 0.43 0.38 0.44  0.5
  # Avg Prec: 1 +  0  + 0.67 + 0 +0 +0.5 + 0 + 0 + 0.44 +0.5) / 5 = 0.62

  # if there are genes in geneset1 that are also in geneset2...
  # then manually append them to the top of RWR results with the best
  # possible rank. This is because RWR_LOE does not provide ranks to
  # genes in the seed geneset, but in this scenario we need them to have
  # a top rank because technically they represent exact overlap between
  # geneset1 and geneset2

  dupes <- dplyr::slice(
    query_geneset,
    which(query_geneset$gene %in% seed_geneset$gene)
  )

  if (nrow(dupes) > 0) {
    warning(sprintf("WARNING: %i genes in the query geneset are also found in the seed geneset. They will be removed from the query geneset before running eval\n", length(dupes)), appendLF = T) # nolint
    res <- rbind(
      data.frame(
        NodeNames = dupes$gene,
        Score = 1,
        rank = 1,
        num_in_network = dplyr::first(res$num_in_network),
        num_seeds = dplyr::first(res$num_seeds),
        networks = dplyr::first(res$networks),
        modname = dplyr::first(res$modname),
        seed_geneset = dplyr::first(res$seed_geneset),
        query_geneset = dplyr::first(res$query_geneset),
        InQueryGeneset = 1
      ),
      res
    )
  }

  res <- calc_ROCPRC(
    df = res,
    scorescol = "rank",
    labelscol = "InQueryGeneset"
  )

  filtervalue <- dplyr::filter(res, rank == nrow(!!query_geneset))
  outputvalue <- filtervalue %>% dplyr::pull(PREC) # nolint PREC is column name

  # 1. precision @ size(query_geneset) (aka R-PREC)
  output <- data.frame(
    value = ifelse(length(outputvalue) > 0,
      outputvalue,
      NA
    ),
    measure = "P@SizeQueryGeneSet"
  )

  # 2. Avg PRC (uses Average Precision, not interpolated precision)
  output <- rbind(
    output,
    dplyr::filter(res, TP == 1) %>%
      dplyr::summarise(
        value = sum(PREC) / sum(InQueryGeneset),
        measure = "AvgPrec"
      )
  )

  # 3. AUPRC
  output <- rbind(
    output,
    dplyr::summarise(
      res,
      value = area_under_curve(
        REC,
        PREC,
        method = "trapezoid",
        ties = "max"
      ),
      measure = "AUPRC"
    )
  )

  # 4. AUROC
  output <- rbind(
    output,
    dplyr::summarise(
      res,
      value = sum(REC) / dplyr::n(), measure = "AUROC"
    )
  )

  # # 5. NDCG
  # output <- rbind(
  #   output,
  #   dplyr::filter(
  #     res, rank == nrow(!!query_geneset)
  #   ) %>%
  #     dplyr::summarise(
  #       value = ndcg,
  #       measure = "NDCG@SizeQueryGeneSet"
  #     )
  # )

  # output <- rbind(
  #   output,
  #   dplyr::summarise(
  #     res,
  #     value = area_under_curve(
  #       REC,
  #       ndcg,
  #       method = "trapezoid",
  #       ties = "max"
  #     ),
  #     measure = "AUNDCG"
  #   )
  # )

  return(list(summary = output, results = res))
}

#' @keywords internal
view_top_network_loe <- function(results,
                                 nw_mpo,
                                 seed_geneset,
                                 query_geneset = NULL,
                                 ntop = 199,
                                 cyto = 1,
                                 modname = "") {
  if (cyto) {
    topresults <- RandomWalkRestartMH::create.multiplexNetwork.topResults(
      results,
      nw_mpo,
      k = ntop
    )

    igraph::V(topresults)$label.cex <- 0.6 # nolint: igraph method
    if (!is.null(query_geneset)) {
      type_data <- dplyr::case_when(
        igraph::V(topresults)$name %in% seed_geneset$gene ~ "InSet1",
        igraph::V(topresults)$name %in% query_geneset$gene ~ "InSet2",
        TRUE ~ "Other"
      )


      igraph::V(topresults)$type <- type_data
    } else {
      type_data <- dplyr::case_when(
        igraph::V(topresults)$name %in% seed_geneset$gene ~ "InSet1",
        TRUE ~ "Other"
      )
      igraph::V(topresults)$type <- type_data
    }

    igraph::V(topresults)$Rank <- NA # nolint: igraph methods
    tmp <- match(igraph::V(topresults)$name, results$RWRM_Results$NodeNames)
    igraph::V(topresults)$Rank <- results$RWRM_Results$rank[tmp] # nolint: igraph methods

    RCy3::cytoscapePing()
    RCy3::createNetworkFromIgraph(
      topresults,
      title = paste0("seeds_to_top", ntop),
      collection = modname
    )
    RCy3::layoutNetwork("kamada-kawai")
    RCy3::setNodeBorderColorDefault(new.color = "#666666")
    RCy3::setNodeBorderWidthDefault(new.width = 4)
  }
}

#' @keywords internal
save_plots_loe <- function(metrics,
                           seed_geneset,
                           query_geneset,
                           outdir,
                           modname) {
  # this uses custom scoring and avoids the use of Precrec library
  if (is.null(outdir)) {
    message("You must provide a path to --outdir to save the metrics and plots.") # nolint message
    return(1)
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  p1 <- ggplot2::ggplot(metrics$results) +
    ggplot2::geom_path(
      ggplot2::aes(x = REC, y = PREC),
      col = "darkorange",
      alpha = 0.5
    ) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::xlab("Recall") +
    ggplot2::ylab("Precision")

  ### ROC
  #  plot ROC curve including all folds and average
  p2 <- ggplot2::ggplot(metrics$results) +
    ggplot2::geom_line(
      ggplot2::aes(x = FPR, y = REC),
      col = "darkorange",
      alpha = 1
    ) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1),
      linetype = "dashed",
      col = "darkgrey"
    ) +
    ggplot2::xlab("FPR") +
    ggplot2::ylab("TPR")

 
  # plot ranking distribution of hits in bins of 100 for each fold
  p4 <- ggplot2::ggplot(
    metrics$results %>% dplyr::filter(InQueryGeneset == 1)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(x = rank),
      alpha = 1,
      binwidth = 100
    )

  out_path <- get_file_path(seed_geneset$setid[1],
    query_geneset$setid[1],
    modname,
    outdir = outdir,
    ext = "metrics.png"
  )

  png(out_path, width = 1200, height = 1000)

  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
  print(p1, vp = vplayout(1, 1)) # Top left
  print(p2, vp = vplayout(1, 2)) # Top right
  print(p4, vp = vplayout(2, 1:2)) # Bottom Right
  dev.off()
  message(paste("File saved:", out_path))
}

########################################################################
# Main Function
########################################################################
#' RWR LOE
#'
#' `RWR_LOE` perform RWR lines-of-evidence.  RWR_LOE has two possible functions.
#' Given one geneset of seeds, rankings for all other genes in the network will
#' be returned.  Given a second geneset of genes to be queried, ranking for just
#' the genes in that geneset will be returned.  This can be used to build
#' multiple lines of evidence from the various input networks to relate the two
#'  gene sets.
#'
#' @param data The path to the Rdata object generated from RWR_make_multiplex
#' @param seed_geneset The path to the gene set file. It must have the following
#'                     first two columns with no headers tab-delimited:
#'                      {<}setid{>} {<}gene{>} {<}weight{>}
#' @param query_geneset The path to the optional second gene set file. It must
#'                      have the following first two columns with no headers
#'                      tab-delimited:
#'                      {<}setid{>} {<}gene{>} {<}weight{>}
#' @param restart Sets the restart parameter [0,1). Higher value means the
#'                walker will jump back to a seed node more often.
#'                Defaulted to 0.7
#' @param tau Comma-separated list of values between that MUST add up to the
#'            number of network layers in the .Rdata file. One value per network
#'            layer that determines the probability that the random walker will
#'            restart in that layer. e.g. if there are three layers (A,B,C) in
#'            your multiplex network, then --tau '0.2,1.3,1.5' will mean that
#'            layer A is less likely to be walked on after a restart than layers
#'            B or C.
#'            Default 1.0
#' @param outdir Full path to the output file directory
#' @param numranked Proprtion of ranked genes to return \[0,1\].  e.g. 0.1 will
#'                  return the top 10%. Default 1.0
#' @param plot Output PNG plot of ROC and PRC of the query geneset with respect the
#'              seed geneset. 
#' @param modname String to include in output file name. Default "default"
#' @param cyto Specify a number N > 0 if you wish to see a network of the seeds
#'             and top N ranked genes (cytoscape must already be running)
#' @param verbose Verbose mode (Default False)
#' @return Returns a list of:
#'                  (1) a dataframe of results (scores for genes) and
#'                  (2) the seed genes.
#' @examples
#' # An example of a default RWR LOE with one geneset which will
#' # return a ranked list
#'
#' extdata.dir <- system.file("example_data", package = "RWRtoolkit")
#' multiplex_object_filepath <- paste(
#'   extdata.dir,
#'   "/string_interactions.Rdata",
#'   sep = ""
#' )
#' geneset1_filepath <- paste(extdata.dir, "/geneset1.tsv", sep = "")
#' geneset2_filepath <- paste(extdata.dir, "/geneset2.tsv", sep = "")
#' outdir <- "./rwr_loe"
#'
#' # for all other genes in the network
#' loe_gene_seed_list <- RWR_LOE(
#'   data = multiplex_object_filepath,
#'   seed_geneset = geneset1_filepath,
#'   tau = "1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0",
#'   outdir = outdir
#' )
#'
#' # An example of an RWR LOE with two genesets, the first defines the seeds
#' # and the output will be rankings for the genes in the second geneset
#' loe_gene_seed_list2 <- RWR_LOE(
#'   data = multiplex_object_filepath,
#'   seed_geneset = geneset1_filepath,
#'   query_geneset = geneset2_filepath,
#'   tau = "1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0",
#'   outdir = outdir
#' )
#'
#' @export
RWR_LOE <- function(data = NULL, # nolint PACKAGE FUNCTION NAME
                    seed_geneset = NULL,
                    query_geneset = NULL,
                    restart = 0.7,
                    tau = "1.0",
                    outdir = "./",
                    numranked = 1.0,
                    plot = FALSE,
                    modname = "default",
                    cyto = 0,
                    verbose = FALSE) {
  # Check parameters and print key parameters on verbose mode
  if (is.null(seed_geneset)) {
    stop("ERROR: Mandatory arguement seed_geneset is missing.")
  }

  if (verbose) {
    cat(
      paste(
        "\nRdata file:     ",
        data,
        "\nSeed geneset: ",
        seed_geneset,
        "\nQuery geneset: ",
        query_geneset,
        "\noutput: ",
        outdir,
        "\nproportion of ranked nodes to be saved: ",
        numranked,
        "\n\n"
      )
    )
  }

  # load the data:
  data_list <- load_multiplex_data(data)
  nw_mpo <- data_list$nw.mpo
  nw_adjnorm <- data_list$nw.adjnorm


  # Set tau which must match the network layers (utility function)
  tau <- get_or_set_tau(nw_mpo, optTau = tau)

  # Load the gene sets, or set to NULL (utility function)
  seed_geneset <- load_geneset(seed_geneset, nw_mpo, verbose)$geneset
  query_geneset <- load_geneset(
    query_geneset,
    nw_mpo,
    verbose
  )$geneset # will return NULL if not set by user

  # Core of method
  message("\nBeginning RWR_LOE ...")

  # Set network_names to be the concatenation of individual network names
  network_names <- paste(
    names(nw_mpo)[1:nw_mpo$Number_of_Layers],
    collapse = "_"
  )

  # Run RWR on multiplex network starting at seed genes defined by
  # the seed_geneset variable.
  results <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(
    x = nw_adjnorm,
    MultiplexObject = nw_mpo,
    Seeds = seed_geneset$gene,
    r = restart
  )

  # Add a rank column to the results (min_rank gives the minimum rank on a
  # tie, skips next ranking)
  # Any gene that scored 0, give it worst possible rank
  results$RWRM_Results <- results$RWRM_Results %>% # nolint RWRMH package
    dplyr::mutate(rank = dplyr::min_rank(-Score)) %>%
    dplyr::mutate(
      rank = dplyr::if_else(
        Score == 0,
        true = nw_mpo$Number_of_Nodes_Multiplex,
        false = rank
      )
    )
  # Add in additional columns to results
  results$RWRM_Results <- results$RWRM_Results %>% # nolint RWRMH package
    dplyr::mutate(
      num_in_network = nrow(seed_geneset),
      num_seeds = length(results$Seed_Nodes),
      networks = network_names,
      modname = modname,
      seed_geneset = seed_geneset$setid[1]
    )

  # If a second geneset was input, flag genes that are in query_geneset
  if (!is.null(query_geneset)) {
    query_message <- "A query geneset was provided.  Results will be filtered for only genes from the query geneset." # nolint message
    message(query_message)
    results$RWRM_Results <- results$RWRM_Results %>% # nolint RWRMH package
      tibble::add_column(query_geneset = query_geneset$setid[1]) %>%
      tibble::add_column(
        InQueryGeneset = as.numeric(
          results$RWRM_Results$NodeNames %in% query_geneset$gene
        )
      )
  } else {
    message("A query geneset was not provided.")
  }
  message("... RWR_LOE complete.")
  print(head(results$RWRM_Results))

  # Define output path
  if (!is.null(query_geneset)) {
    out_path <- get_file_path(
      "RWR-LOE",
      seed_geneset$setid[1],
      "vs",
      query_geneset$setid[1],
      modname,
      outdir = outdir,
      ext = ".ranks.tsv"
    )
  } else {
    out_path <- get_file_path("RWR-LOE",
      seed_geneset$setid[1],
      modname,
      outdir = outdir,
      ext = ".ranks.tsv"
    )
  }
  message(paste("Output path:", out_path))

  # Save the table of results, ie, scores and ranks
  # for each gene in the network.
  result_table <- results$RWRM_Results %>% dplyr::slice_head(prop = numranked)
  write_table(result_table, out_path, verbose = TRUE)

  # If evaluation mode is turned on and both a query and seed geneset are present, plot ROC and PRC.
  if (!is.null(query_geneset) && plot) {
    message("evaluating metrics for finding query genes from seeds")
    metrics <- calc_metrics_loe(
      results$RWRM_Results,
      seed_geneset,
      query_geneset
    )

    out_path <- get_file_path(
      seed_geneset$setid[1],
      query_geneset$setid[1],
      modname,
      outdir = outdir,
      ext = "metrics.tsv"
    )

    write_table(metrics$summary, out_path)
    message(paste("Saved metrics summary to file:", out_path))
    save_plots_loe(metrics, seed_geneset, query_geneset, outdir, modname)
  } else if (is.null(query_geneset) && plot) {
    warning(sprintf("Plot mode is turned ON but there is no query geneset, evaluation will not occur.\n")) # nolint message
  }

  # If cytoscape mode is turned on, open an interactive cytoscape session.
  if (cyto > 0) {
    message("opening top n network in cytoscape")
    view_top_network_loe(
      results,
      nw_mpo,
      seed_geneset,
      query_geneset,
      ntop = cyto,
      cyto = cyto,
      modname = modname
    )
  }

  return(results)
}
