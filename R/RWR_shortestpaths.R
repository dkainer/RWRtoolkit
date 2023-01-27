########################################################################
# Find shortest paths between genes in gene sets. Given a single gene set,
#   find the shortest paths between the genes in that gene set. Given two
#   gene sets, find the shortest paths for pairs of genes between gene sets.
# - Input: Pre-computed multiplex network and one or two genesets
# - Output: Edge list table
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
########################################################################

#' @importFrom dplyr %>%
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%

########################################################################
# Internal Functions
########################################################################

merge_networks <- function(nw_mpo) {
  # Merge the subnetworks into one big network.
  message("Merging network layers down...")
  nl <- nw_mpo$Number_of_Layers
  nw_dflist <- lapply(nw_mpo[1:nl], igraph::as_data_frame)
  nw_df <- do.call(rbind, nw_dflist)
  nw_df <- nw_df %>%
    dplyr::group_by(type) %>%
    dplyr::mutate(weightnorm = weight / sum(weight))
  nw_merged <- igraph::graph_from_data_frame(nw_df, directed = FALSE)
  message("Done.")
  return(nw_merged)
}

get_shortest_paths <- function(nw_merged,
                               source_geneset,
                               target_geneset,
                               threads = 1) {
  # For each gene in source_geneset, get the
  # shortest path to each gene in target_geneset.

  targets <- which(igraph::V(nw_merged)$name %in% target_geneset$gene)
  wt <- 1 - igraph::E(nw_merged)$weightnorm

  message(sprintf(
    "Calculating Shortest paths for %d x %d = %d gene pairs",
    nrow(source_geneset),
    length(targets),
    nrow(source_geneset) * length(targets)
  ))

  doParallel::registerDoParallel(cores = threads)
  res <- foreach::foreach(g = source_geneset$gene, .combine = rbind) %:%
    foreach::foreach(t = targets, .combine = rbind) %dopar% {
      # Only get the shortest path if the start and end
      # nodes are not the same!
      if (igraph::V(nw_merged)$name[t] != g) {
        sp <- igraph::shortest_paths(
          nw_merged,
          from = g,
          to = t,
          output = "vpath",
          weights = wt
        )


        sg_df <- igraph::as_data_frame(
          igraph::induced_subgraph(
            nw_merged,
            vids = sp$vpath[[1]],
            impl = "create_from_scratch"
          )
        )


        sg_df$pathname <- paste(g,
          igraph::V(nw_merged)$name[t],
          sep = "_"
        )

        sg_df$pathlength <- length(sp$vpath[[1]])

        path_string <- paste(
          lapply(
            sp$vpath[[1]],
            function(x) igraph::V(nw_merged)$name[x]
          ),
          collapse = "->"
        )

        sg_df$pathelements <- path_string
        sg_df
      }
    }
  doParallel::stopImplicitCluster()

  message("Finished calculating Shortest paths.")
  return(res)
}

save_results <- function(rwr_result,
                         source_geneset = NULL,
                         target_geneset = NULL,
                         outdir = NULL,
                         out_path = NULL) {
  # First, generate the out_path.
  if (is.null(outdir) && is.null(out_path)) {
    warning("You must provide either `outdir` or `out_path`.")
  } else if (!is.null(out_path)) {
    out_path <- out_path
  } else {
    if (is.null(target_geneset)) {
      out_path <- get_file_path(
        "RWR-SPATHS_",
        source_geneset$setid[1],
        outdir = outdir
      )
    } else {
      out_path <- get_file_path(
        "RWR-SPATHS_",
        source_geneset$setid[1],
        target_geneset$setid[1],
        outdir = outdir
      )
    }
  }

  # If out_path is still NULL, there's nothing to do.
  if (!is.null(out_path)) {
    if (!dir.exists(dirname(out_path))) {
      dir.create(dirname(out_path), recursive = TRUE)
    }
    write_table(rwr_result, out_path)
    message(sprintf("Saved results to file:\n  %s", out_path))
  }
}

open_cytoscape <- function(res, source_geneset, target_geneset) {
  # Visualize in Cyto.
  ig <- igraph::graph_from_data_frame(res, directed = F)
  igraph::V(ig)$endpoint <- igraph::V(ig)$name %in%
    c(source_geneset$gene, target_geneset$gene)
  RCy3::cytoscapePing()
  RCy3::createNetworkFromIgraph(ig, "myIgraph")
  RCy3::setNodeBorderColorDefault(new.color = "#666666")
  RCy3::setNodeBorderWidthDefault(new.width = 4)
}



extract_node_from_row <- function(path_df,
                                  node_name,
                                  known_path = NULL,
                                  is_penultimate = FALSE,
                                  ultimate = NULL) {
  rows <- path_df[(path_df$to == node_name & !path_df$from %in% known_path) |
    (path_df$from == node_name & !path_df$to %in% known_path), ]

  node_in_to_data <- if (is_penultimate) {
    node_name %in% rows$to && ultimate %in% rows$from
  } else {
    node_name %in% rows$to
  }

  from_column <- if (node_in_to_data) "to" else "from"
  from_node <- unique(rows[[from_column]])

  to_column <- if (node_in_to_data) "from" else "to"
  to_node <- unique(rows[[to_column]])

  to_node <- unlist(
    lapply(
      to_node,
      function(x) if (x %in% known_path) NULL else x
    )
  )

  from_node <- unlist(
    lapply(
      from_node,
      function(x) if (x %in% known_path) NULL else x
    )
  )

  # TODO: Refactor to be less convoluted.
  if (is_penultimate && !is.null(ultimate)) {
    # the penultimate node has nodes within the from and to as it

    if (ultimate %in% from_node) {
      from_node <- unlist(
        lapply(
          from_node,
          function(x) if (x == ultimate) x else NULL
        )
      )
    }

    if (ultimate %in% to_node) {
      to_node <- unlist(
        lapply(
          to_node,
          function(x) if (x == ultimate) x else NULL
        )
      )
    }
  }

  list(
    from_column = from_column,
    from = from_node,
    to_column = to_column,
    to = to_node
  )
}



########################################################################
# Main Function
########################################################################

#'  Extract Highest Scoring Path
#'
#' `extract_highest_scoring_path`  parses through the dataframe output of
#' RWR_ShortestPaths and extracts the highest weighted edge per layer of a
#' user defined path. The method with which these edges are extracted can
#' be either by the normalized weight (default) or by the non-normalized
#' weight. If the edges are unweighted, this will not yield any meaningful
#' results
#'
#' @param shortest_paths_df A dataframe denoting the shortest paths of
#'                          between a series of predefined genes. The output
#'                          of RWR_ShortestPaths
#' @param desired_path      A string denoting the desired path that the user
#'                          wishes to extract (must match pathname of an
#'                          existing path within the shorest_paths_df).
#' @param weight_type       A string determining the weight used to determine
#'                          the highest weighted path. Options:
#'                          "weightnorm" (Default)  uses normalized edges by
#'                                      edgeweight_i / N_vertices_in_layer
#'
#'                          "weight" uses edge weight alone.
#' @return Returns a data frame following the path with the highest edge weights
#' @examples
#'
#' # An example of Running `extract_highest_scoring_path`
#'
#' extdata.dir <- system.file("example_data", package = "RWRtoolkit")
#' multiplex_object_filepath <- paste(extdata.dir,
#'   "/string_interactions.Rdata",
#'   sep = ""
#' )
#' geneset1_filepath <- paste(extdata.dir, "/geneset1.tsv", sep = "")
#' geneset2_filepath <- paste(extdata.dir, "/geneset2.tsv", sep = "")
#' outdir <- paste(extdata.dir, "/out/rwr_shortestpath", sep = "")
#'
#' rwr_shortest_path_output <- RWR_ShortestPaths(
#'   data = multiplex_object_filepath,
#'   source_geneset = geneset1_filepath,
#'   target_geneset = geneset2_filepath,
#'   write_to_file = TRUE,
#'   outdir = outdir
#' )
#'
#' optimal_path <- extract_highest_scoring_path(
#'   rwr_shortest_path_output,
#'   desired_path = "TPI1_PMM2"
#' )
#'
#' @export
extract_highest_scoring_path <- function(shortest_paths_df,
                                         desired_path,
                                         weight_type = "normalized") {
  weight_column <- if (weight_type == "weightnorm") "weightnorm" else "weight"
  path_rows_df <- shortest_paths_df[
    shortest_paths_df$pathname == desired_path,
  ]

  total_path <- unique(path_rows_df$pathelements)

  split_path <- unlist(stringr::str_split(total_path, "->"))

  known_path <- c()
  path_rows <- NULL

  for (node_idx in seq(1, length(split_path) - 1)) {
    is_penultimate <- if (node_idx == length(split_path) - 1) T else F
    ultimate <- if (is_penultimate) split_path[length(split_path)] else NULL

    node <- split_path[node_idx]
    node_information <- extract_node_from_row(
      path_rows_df,
      node,
      known_path,
      is_penultimate,
      ultimate
    )

    from_column <- node_information$from_column
    to_column <- node_information$to_column

    desired_rows <- path_rows_df[
      path_rows_df[[from_column]] == node_information$from &
        path_rows_df[[to_column]] == node_information$to,
    ]

    top_row <- desired_rows[which.max(desired_rows[[weight_column]]), ]

    top_row$from <- node_information$from
    top_row$to <- node_information$to

    known_path <- append(known_path, node)
    path_rows <- rbind(path_rows, top_row)
  }
  path_rows
}

#' RWR Shortest Paths
#'
#' `RWR_ShortestPaths` Find shortest paths between genes in the given gene
#' sets. These will be output as a table of edges. There are two ways to use
#' this:
#' 1) provide one geneset (g) and you will get all shortest paths between genes
#' in g.
#' 2) provide two genesets (g and p) and you will get the shortest paths
#' between genes in g and genes in p.
#'
#' @param data Path to the .Rdata file for your combo of
#' underlying functional networks. This file is produced by
#' RWR_make_MHobject.R
#' @param source_geneset Path to the gene set file. If given with
#' --target_geneset, find shortest paths between genes in --source_geneset and
#' --target_geneset. Otherwise, find shortest paths among genes in
#' --source_geneset only. It must have the following cols without heading:
#' {<}setid{>} {<}gene{>}.
#' @param target_geneset Path to the second geneset file. If it is not provided
#' then the shortest paths will be calculated between genes in source_geneset.
#' @param outdir Full path to the output file directory. Two output files will
#' be generated with different suffixes.  Default "."
#' @param out_path Specify the full path for output. Ignore --outdir and
#' --modname and use this instead.
#' @param threads Number of threads to use.
#' @param cyto Include this parameter if you wish to see the shortest paths in
#' Cytoscape. Cytoscape must already be running first!  If your geneset(s) are
#' large (e.g. 100+) then loading all the shortest paths fully into cytoscape
#' can take a while. Default FALSE
#' @param verbose Verbose mode.  Default FALSE
#' @param write_to_file Also write the results to a file.  Default FALSE
#' @return Returns a data frame.
#' @examples
#'
#' # An example of Running RWR_ShortestPaths
#'
#' extdata.dir <- system.file("example_data", package = "RWRtoolkit")
#' multiplex_object_filepath <- paste(extdata.dir,
#'   "/string_interactions.Rdata",
#'   sep = ""
#' )
#' geneset1_filepath <- paste(extdata.dir, "/geneset1.tsv", sep = "")
#' geneset2_filepath <- paste(extdata.dir, "/geneset2.tsv", sep = "")
#' outdir <- paste(extdata.dir, "/out/rwr_shortestpath", sep = "")
#'
#' rwr_shortest_path_output <- RWR_ShortestPaths(
#'   data = multiplex_object_filepath,
#'   source_geneset = geneset1_filepath,
#'   target_geneset = geneset2_filepath,
#'   write_to_file = TRUE,
#'   outdir = outdir
#' )
#'
#' @export
RWR_ShortestPaths <- function( # nolint
                              data = NULL,
                              source_geneset = NULL,
                              target_geneset = NULL,
                              outdir = ".",
                              out_path = NULL,
                              threads = 1,
                              cyto = FALSE,
                              verbose = FALSE,
                              write_to_file = FALSE) {
  nw_mpo <- NULL
  # This is a saved version of the multiplex network and adj matrix.
  load(data)
  nw_mpo <- nw.mpo # nolint - loaded from filepath
  if (verbose) {
    message("Loaded data:")
    print(ls())
  }

  # If geneset is NULL, get shortest paths
  # from source_geneset to source_geneset.
  if (is.null(target_geneset)) {
    target_geneset <- source_geneset
  }
  # Gene set 1.
  message("Getting source gene set ...")
  source_geneset_plus_extras <- load_geneset(source_geneset, nw_mpo, verbose)
  source_genes <- source_geneset_plus_extras[[1]]
  message(head(source_genes))

  # Gene set 2.
  message("Getting target gene set ...")
  target_geneset_plus_extras <- load_geneset(target_geneset, nw_mpo, verbose)
  target_genes <- target_geneset_plus_extras[[1]]
  message(head(target_genes))

  nw_merged <- merge_networks(nw_mpo)

  res <- get_shortest_paths(nw_merged, source_genes, target_genes, threads)

  if (write_to_file) {
    save_results(
      res,
      source_geneset = source_genes,
      target_geneset = target_genes,
      outdir = outdir,
      out_path = out_path
    )
  }

  if (cyto) {
    open_cytoscape(res, source_genes, target_genes)
  }

  return(res)
}
