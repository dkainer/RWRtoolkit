###############################################################################
# Generates a homogenous or heterogeneous multiplex network for downstream
# functions
# - Input: A file list of a set of network
# - Output: A multiplex object
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
#' @importFrom foreach %do%

########################################################################
# Internal Functions
########################################################################
get_neighboring_edge_ids <- function(g, node){
    neighbors <- igraph::neighbors(g, node)$name
    connected_nodes <- c()
    for (neighbor in neighbors){
        connected_nodes <- c(connected_nodes, node, neighbor)
    }
    edge_ids <- igraph::get.edge.ids(g, connected_nodes)
    return(edge_ids)
}


remove_edges_connected_to_node <- function(g, node){
    if( ! (node %in% igraph::V(g)$name )) return(g)
    edge_ids_to_remove <- get_neighboring_edge_ids(g, node)
    g <- igraph::delete_edges(g, edge_ids_to_remove)
    g
}


drop_edges_from_nodes_in_graph <- function(g, elements_to_remove){
    for (element in elements_to_remove){
        g <- remove_edges_connected_to_node(g, element)
    }
    g
}



get_class_of_nw_groups <- function(nw_groups){
  groups_class <- class(nw_groups)
  if (length(groups_class) > 1){
    if ('data.frame' %in% groups_class) return('data.frame')  
  }

  unique_classes <- unique(purrr::map(nw_groups, class))
  if (length(unique_classes) > 1 ) {
    stop('List of graph elements ought to be only igraph objects')
  }
  unique_classes[[1]]
}

get_network_name <- function(g, index, previous_names){
  graph_name_does_not_exist <- is.null(igraph::graph_attr(g, "name"))
  matches_other_name_in_list <- if (!graph_name_does_not_exist)  
    igraph::graph_attr(g, "name") %in% previous_names 
    else FALSE
    
  if (graph_name_does_not_exist || matches_other_name_in_list) {
    return(paste("graph", index, sep='_'))
  }
  return( igraph::graph_attr(g, "name") )
}

# get_graph_list_names <- function(graph_list){
#   name_list <- c()
#   for (g_idx in 1:length(graph_list)){
#     g <- graph_list[[g_idx]]
#     name <- get_network_name(g, g_idx, name_list)
#     name_list <- c(name_list, name)
#   }
#   name_list
# }

make_multiplex <- function(nwdf, knockout_nodes=c()) {
  # Make_multiplex is a wrapper around iGraph and RandomWalkRestartMH

  # Preparing data for create.multiplex function call
  # For each network, read as a datatable, convert to igraph, attribute 
  # edges with nwnames

 
  nwlist <- if (get_class_of_nw_groups(nwdf) != "igraph") foreach::foreach(d = iterators::iter(nwdf, by = "row")) %do% {
    nw_g <- load_network(
      d$nwfile,
      type = d$nwname,
      name= d$nwname,
      col_names = c("from", "to", "weight"),
      select = 1:3, 
    )
    nw_g
  } else nwdf

  name_vector <- c()
  for (g_idx in 1:length(nwlist)){
    g <- nwlist[[g_idx]]
    name <- get_network_name(g, g_idx, name_vector)
    g <- igraph::set_graph_attr(g, 'name', name)
    name_vector <- c(name_vector, name)


    if (length(knockout_nodes) > 0){
      g <- drop_edges_from_nodes_in_graph(g, knockout_nodes)

    }

    g <- igraph::set_edge_attr(g, "type", seq(1, igraph::ecount(g)), name)
    nwlist[[g_idx]] <- g
  }
  
    


  names(nwlist) <- name_vector
  # nwlist 


  # Make the multiplex network object (return as mpo)
  print("constructing a multiplex network...be patient if there are lots of layers or layers are big")
  mpo <-  RandomWalkRestartMH::create.multiplex(LayersList = nwlist)
  print("constructing a multiplex network...DONE")
  return(mpo)
}

make_homogenous_network <- function(nw.groups, delta, out=NULL, knockout_nodes = c(), verbose=F) {
  # Takes in group data from flist, converts to multiplex network, and saves to a file

  # Constructs Multiplex Network
  # cat("constructing multiplex network from 1 group of layers: ", nw.groups[[1]]$nwname, "\n")
  nw.mpo <- make_multiplex(nw.groups, knockout_nodes)
  
  # Create the adjacency matrix and normalize the data for the network
  print("constructing the adjacency matrix...be VERY patient if there are lots of layers or layers are big")
  nw.adj <- RandomWalkRestartMH::compute.adjacency.matrix(nw.mpo, delta = delta)
  nw.adjnorm <- RandomWalkRestartMH::normalize.multiplex.adjacency(nw.adj)

  if (is.null(out)){
    return(list(
      "nw.mpo" = nw.mpo, 
      "nw.adj" = nw.adj, 
      "nw.adjnorm" = nw.adjnorm
    ))
  }
  # Save data to file with presupplied filename or default: network.Rdata
  if (!dir.exists(dirname(out))) {
    dir.create(dirname(out), recursive = TRUE)
  }
  save(nw.mpo, nw.adj, nw.adjnorm, list = c("nw.mpo", "nw.adj", "nw.adjnorm"), file = out)
  message("\nDONE - Homogenous Multiplex object saved to RData file for use in further functions.")
  message(paste("File path: ", out))
}

make_heterogeneous_multiplex <- function(nw.groups, delta, lambda, out, verbose) {
  # - Implementation functional - no downstream analysis available
  # - Takes in group data from flist, converts to multiplex heterogenous network, and saves to file.
  # TODO: Impelement multi-delta if desired
  # TODO: complete

  cat(
    "constructing multiplex heterogeneous network from 2 groups of layers plus bipartite links:\n",
    nw.groups[[1]]$nwname, "\n",
    nw.groups[[2]]$nwname, "\n",
    nw.groups[[3]]$nwname, "\n"
  )

  # Create 2 groups of layers
  nw.mpo1 <- make_multiplex(nw.groups[[1]])
  nw.mpo2 <- make_multiplex(nw.groups[[2]])

  # Link layers one and two with the bipartite group in nw.groups[[3]]
  cat("prepping bipartite links between groups 1 and 2\n")
  bipartite_links <- foreach::foreach(d = iterators::iter(nw.groups[[3]], by = "row"), .combine = rbind) %do% {
    print(d$nwname)
    nw <- data.table::fread(d$nwfile, col.names = c("from", "to", "weight"))
    nw
  }

  # Retrieve links between graph networks.
  bipartite_links <- bipartite_links %>%
    dplyr::filter(from %in% nw.mpo1$Pool_of_Nodes) %>%
    dplyr::filter(to %in% nw.mpo2$Pool_of_Nodes)
  cat("putting it all together...\n")

  # Combine layers and bipartite links into mutiplex heterogeneous network
  nw.mph <- create.multiplexHet(Multiplex_object_1 = nw.mpo1, Multiplex_object_2 = nw.mpo2, Nodes_relations = bipartite_links)
  cat("constructing full supra-adjacency matrix...be VERY patient if there are lots of layers\n")
  nw.adj <- compute.transition.matrix(nw.mph, delta1 = delta, delta2 = delta, lambda = lambda)

  # Save data to file with presupplied filename or default: network.Rdata
  if (!dir.exists(dirname(out))) {
    dir.create(dirname(out), recursive = TRUE)
  }
  save(nw.mph, nw.mpo1, nw.mpo2, bipartite_links, nw.adj, file = out)
  message("\nDONE - Heterogeneous Multiplex object saved to RData file for use in further functions.")
  message(paste("File path: ", out))
}

read_flist <- function(flist) {
  flist_table <- tryCatch(
      {
        data.table::fread(flist, header = F)
      },
      error = function(cond) {
        # TODO: Depending on verbosity, update error message?
        cat("Error in reading in file list:\n", cond$message)
        stop(cond)
      }
    )

    if (ncol(flist_table) < 2 ) {
      stop("flist files must have at least a file path and a layer name.")
    }

    # add extra column if only 2 supplied. 
    if (ncol(flist_table) == 2) {
      flist_table$V3 <- 1
    }

    # Extract only the first threee columns
    col_sliced_flist <- flist_table[ , c('V1', 'V2', 'V3')]
    colnames(col_sliced_flist) <- c("nwfile", "nwname", "nwgroup")

    return(col_sliced_flist)

}

########################################################################
# Main Function
########################################################################

#' RWR Make Multiplex
#'
#' `RWR_make_multiplex` creates a multiplex network.
#'
#' @param flist  Table describing network files to use.  File columns:
#'               {<}path to file{>} {<}short name of network{>}.  {<}group{>}.
#'               'groups' are either 1, 2, or 3.  All 1's will form one
#'               multiplex network (e.g. gene-to-gene), All 2's will form a
#'               separate multiplex network (e.g. disease-to-disease), And all
#'               3's will be used to join the 1's and 2's together (e.g.
#'               gene-to-disease) You don't have to have both 1's and 2's.
#'               But if you do have 1's and 2's, you SHOULD have at least one
#'               3 to join them up.  Can be delimited by comma, tab, space,
#'               pipe, or semicolon.
#' @param delta  Probability to change between layers at the next step \[0,1\].
#'               If delta = 0, the particle will always remain in the same layer
#'               after a non-restart iteration.  If delta = 1, the particle will
#'               always change between layers, therefore not following the
#'               specific edges of each layer. Default is 0.5.
#' @param output Output file name (default "network.Rdata")
#' @param test   Runs an example. Default FALSE
#' @param verbose Verbose mode. Default FALSE
#' @return Mutliplex object is saved to a file (.rdata) to load into subsequent
#'          functions.
#' @examples
#' #
#' # An example of a default RWR Make Multiplex with an output "network.Rdata"
#' extdata.dir <- system.file("example_data", package = "RWRtoolkit")
#' outdir <- "./rwr_make_multiplex"
#'
#' layers.path <- paste(extdata.dir, "/layers/", sep = "")
#' layers <- list.files(layers.path)
#' layer_with_paths <- paste(layers.path, layers, sep = "")
#' layer_names <- sub(pattern = "(.*)\\..*$",
#'                    replacement = "\\1", basename(layers))
#' groups <- rep(1, length(layer_names))
#' flistdatatable <- data.table::data.table(layer_with_paths,
#'                                          layer_names,
#'                                          groups)
#'
#' outfile <- paste(outdir, "/multiplex.Rdata", sep = "")
#' write.table(flistdatatable,
#'   row.names = FALSE,
#'   col.names = FALSE,
#'   sep = "\t",
#'   file = "example.flist", quote = FALSE
#' )
#'
#' RWR_make_multiplex(
#'   flist = "example.flist",
#'   output = outfile
#' )
#'
#'
#'
#' # An example of an RWR Make Multiplex with a non-default delta and
#' # with a specified output filename.
#' outfile <- paste(outdir, "/multiplex_d25_l75.Rdata", sep = "")
#' RWR_make_multiplex(
#'   flist = "example.flist",
#'   delta = 0.25,
#'   output = outfile
#' )
#'
#' system("rm example.flist")
#'
#' @export
RWR_make_multiplex <- function(flist = "", delta = 0.5, output = "network.Rdata",  knockout_nodes = c(), verbose = FALSE) {
  if (flist == "") {
    stop("Please provide a path to your flist, or pass test=TRUE to view an example")
  }

  # Read flist into datatable, fails on file read err
  inDF <- read_flist(flist)

  # Split dataframe into groups
  nw.groups <- inDF %>% dplyr::group_split(nwgroup)

  # Call appropriate network creation function based on groups
  if (length(nw.groups) == 1 && all(nw.groups[[1]]$nwgroup == 1)) {
    make_homogenous_network(nw.groups[[1]], delta, output, knockout_nodes, verbose)
  } else if (length(nw.groups) == 3 && nw.groups[[1]]$nwgroup == 1 && nw.groups[[2]]$nwgroup == 2 && nw.groups[[3]]$nwgroup == 3) {
    warning(
      paste("Hetergeneous Multiplexes are capable of being made",
      " however, the reaminder of the methods in RWRtoolkit have yet",
      "to be validated with respect to the networks.",
      "\n\nThis is planned for V2 of RWRtoolkit."))
    lambda = 0.5
    make_heterogeneous_multiplex(nw.groups, delta, lambda, output, verbose)
  } else {
    # TODO: Add example of flist for homogenous and heterogeneous networks
    stop("Error: Please ensure your fList file is properly formatted", call = F)
  }

  return(0)
}
