###############################################################################
# Perform aggregation on a multiplex objects from RWR_make_multiplex 
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



#' Merged with All Edges
#'
#' Merge down all layers in a multiplex object, but don't
#' aggregate the edges (i.e. keep all edges). This function can take a dummy
#' multiplex or a real multiplex.
#'
#' @param mpo A multiplex network object.
#' @param verbose Print progress to console.
#'
#' @return An igraph network object.
#' @importFrom dplyr %>%
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#'
merged_with_all_edges <- function(mpo, verbose=FALSE) {
  message(sprintf("merging %d network layers ...\n", mpo$Number_of_Layers))
  nl    <- mpo$Number_of_Layers
  nw_dflist <- lapply(mpo[1:nl], igraph::as_data_frame)
  nw_df   <- dplyr::bind_rows(nw_dflist)
  nw_df   <- nw_df %>%
          dplyr::group_by(type) %>%           # nolint
          dplyr::mutate(
            weightnorm = weight / sum(weight) # nolint
          )

  nw_merged <- igraph::graph_from_data_frame(nw_df, directed = FALSE)
  message("merging network layers DONE.")
  message(sprintf("Merged network has %d edges and %d rows.\n",
      igraph::ecount(nw_merged), igraph::vcount(nw_merged)
  ))

  return(
    list(merged_network = nw_merged,
        edge_count = igraph::ecount(nw_merged),
        vertex_count = igraph::vcount(nw_merged))
  )
}

#' Merged With Edgecounts
#'
#' Merge down all layers and aggregate multi-edges into one edge
#' where edge weight is the number of layers in which the two nodes are
#' connected.
#' 
#' @param mpo A multiplex network object. Edge weights are the sum of the
#'       weights of the edges in the multiplex network object.
#' @param inv Set the edge weight to the reciprocal of the sum of the weights.
#' @param verbose Print progress to console.
#'
#' @return An igraph network object.
#' @importFrom dplyr %>%
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
merged_with_edgecounts <- function(mpo, inv=FALSE, verbose=FALSE) {

  n_nodes <- mpo$Number_of_Nodes_Multiplex

  unioned_networks <- NULL
  for (i in 1:mpo$Number_of_Layers) {
    igraph::E(mpo[[i]])$weight <- 1
    unioned_networks <-  igraph::union(mpo[[i]], unioned_networks)
  }

  weight_key_mask <- grepl(
    ".*weight.*",
    names(igraph::edge.attributes(unioned_networks))
  )

  layer_names <- names(
    igraph::edge.attributes(unioned_networks)
    )[weight_key_mask]

  summed_edge_counts <- rep(0, length(igraph::E(unioned_networks)))
  for (layer in layer_names) {
    layer_weights <- igraph::get.edge.attribute(unioned_networks, layer)
    weights_no_na <- replace(layer_weights, is.na(layer_weights), 0)

    summed_edge_counts <- summed_edge_counts + weights_no_na
  }

  igraph::E(unioned_networks)$weight <- summed_edge_counts

  aggr <- igraph::simplify(unioned_networks)


  list(
    merged_network = aggr,
    edge_count = igraph::ecount(aggr),
    vertex_count = igraph::vcount(aggr)
  )
}


#' RWR_network_aggregation
#' 
#' This function acts as an aggregator for RWR multiplex objects. 
#' @param merged_with_all_edges     A boolean denoting a return of a merged
#'                    down multiplex network along with 
#'                    network edge counts and vertex counts.
#'                    Default False
#' @param merged_with_edgecounts    A boolean denoting a return of a merged
#'                    down multiplex, but simplified with
#'                    edge weights denoting the total number
#'                    of layers in which that edge existed.
#'                    Default False
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#' 
RWR_network_aggregation <- function(
    data = NULL,
    flist= NULL,
    outdir = NULL,
    merged_with_all_edges = F,
    merged_with_edgecounts = F,
    verbose = F) {

  utils_output_list <- list()
  `%notin%` <- Negate(`%in%`)


  # Load the multiplex.
  if ( !is.null(data) ) {
    # Get 'nw_mpo' from the .Rdata file.
    nw_mpo <- load_multiplex_data(data)$nw.mpo
  } else if (!is.null(flist)) {
    # Create 'nw_mpo' from the flist.
    nw_mpo <- make_dummy_multiplex(flist, verbose = verbose)
  } else {
    nw_mpo <- NULL
  }

  if (merged_with_all_edges &&
    parameters_exist(
      mpo = nw_mpo,
      required_net = "mpo",
      function_name = "merged_with_all_edges")) {
    utils_output_list$merged_with_all_edges <- merged_with_all_edges(
      nw_mpo,
      verbose = verbose)

    write_networks_to_file_if_fp(
      outdir,
      "merged_with_all_edges.tsv",
      utils_output_list$merged_with_all_edges$merged_network,
      verbose)
  }

  if (merged_with_edgecounts  &&
    parameters_exist(
      mpo = nw_mpo,
      required_net = "mpo",
      function_name = "merged_with_edgecounts")) {
    utils_output_list$merged_with_edgecounts <- merged_with_edgecounts(
      nw_mpo,
      verbose = verbose)

    write_networks_to_file_if_fp(
      outdir,
      "merged_with_edgecounts.tsv",
      utils_output_list$merged_with_edgecounts$merged_network,
      verbose)
  }

  return(utils_output_list)
}
