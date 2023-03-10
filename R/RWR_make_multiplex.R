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

make_multiplex <- function(nwdf) {
  # Make_multiplex is a wrapper around iGraph and RandomWalkRestartMH

  # Preparing data for create.multiplex function call
  # For each network, read as a datatable, convert to igraph, attribute edges with nwnames
  nwlist <- foreach::foreach(d = iterators::iter(nwdf, by = "row")) %do% {
    nw <- data.table::fread(d$nwfile, col.names = c("from", "to", "weight"), select = 1:3)
    nw.g <- igraph::graph_from_data_frame(nw, directed = F)
    igraph::edge_attr(nw.g, "type") <- d$nwname
    nw.g
  }
  names(nwlist) <- nwdf$nwname

  # Calculate the pairwise intersection of edges between each network layer
  # This is useful to understand if any layers seem redundant or not
  # overlaps <- combn(nwlist, 2, simplify=T, FUN = function(x) {
  #       ig <- igraph::graph.intersection(x[[1]],x[[2]],keep.all.vertices=F)
  #       ug <- igraph::union(x[[1]],x[[2]])
  #
  #       return(100*ecount(ig)/ecount(ug))
  # })

  # Make the multiplex network object (return as mpo)
  print("constructing a multiplex network...be patient if there are lots of layers or layers are big")
  mpo <- suppressWarnings(RandomWalkRestartMH::create.multiplex(LayersList = nwlist))
  print("constructing a multiplex network...DONE")
  return(mpo)
}

make_homogenous_network <- function(nw.groups, delta, out, verbose) {
  # Takes in group data from flist, converts to multiplex network, and saves to a file

  # Constructs Multiplex Network
  cat("constructing multiplex network from 1 group of layers: ", nw.groups[[1]]$nwname, "\n")
  nw.mpo <- make_multiplex(nw.groups[[1]])

  # Create the adjacency matrix and normalize the data for the network
  print("constructing the adjacency matrix...be VERY patient if there are lots of layers or layers are big")
  nw.adj <- RandomWalkRestartMH::compute.adjacency.matrix(nw.mpo, delta = delta)
  nw.adjnorm <- RandomWalkRestartMH::normalize.multiplex.adjacency(nw.adj)

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
#' @param flist Table describing network files to use.  File columns: {<}path to file{>} {<}short name of network{>} {<}group{>}.
#' 'groups' are either 1, 2, or 3.  All 1's will form one multiplex network (e.g. gene-to-gene), All 2's will form a
#' separate multiplex network (e.g. disease-to-disease), And all 3's will be used to join the 1's and 2's together
#' (e.g. gene-to-disease) You don't have to have both 1's and 2's.  But if you do have 1's and 2's, you SHOULD have
#' at least one 3 to join them up.  Can be delimited by comma, tab, space, pipe, or semicolon.
#' @param delta Probability to change between layers at the next step \[0,1\]. If delta = 0, the particle will always
#' remain in the same layer after a non-restart iteration.  If delta = 1, the particle will always change between
#' layers, therefore not following the specific edges of each layer. Default is 0.5.
#' @param lambda For heterogeneous networks only.  Probability \[0,1\] the walker can jump between layer groups when it is
#' at a node with a bipartite link. If lambda=1 then walker will oscillate between groups every time it is at a node
#' with a bipartite link.  Default is 0.5.
#' @param output Output file name (default "network.Rdata")
#' @param test Runs an example. Default FALSE
#' @param verbose Verbose mode. Default FALSE
#' @return Mutliplex object is saved to a file (.rdata) to load into subsequent functions.
#' @examples
#' #
#' # An example of a default RWR Make Multiplex with an output "network.Rdata"
#' extdata.dir <- system.file("example_data", package = "RWRtoolkit")
#' outdir <- "./rwr_make_multiplex"
#'
#' layers.path <- paste(extdata.dir, "/layers/", sep = "")
#' layers <- list.files(layers.path)
#' layer_with_paths <- paste(layers.path, layers, sep = "")
#' layer_names <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(layers))
#' groups <- rep(1, length(layer_names))
#' flistdatatable <- data.table::data.table(layer_with_paths, layer_names, groups)
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
#' # lambda with a specified output filename.
#' outfile <- paste(outdir, "/multiplex_d25_l75.Rdata", sep = "")
#' RWR_make_multiplex(
#'   flist = "example.flist",
#'   delta = 0.25,
#'   lambda = 0.75,
#'   output = outfile
#' )
#'
#' system("rm example.flist")
#'
#' @export
RWR_make_multiplex <- function(flist = "", delta = 0.5, lambda = 0.5, output = "network.Rdata",  verbose = FALSE) {
  if (flist == "") {
    stop("Please provide a path to your flist, or pass test=TRUE to view an example")
  }

  # Read flist into datatable, fails on file read err
  inDF <- read_flist(flist)

  # Split dataframe into groups
  nw.groups <- inDF %>% dplyr::group_split(nwgroup)

  # Call appropriate network creation function based on groups
  if (length(nw.groups) == 1 && all(nw.groups[[1]]$nwgroup == 1)) {
    make_homogenous_network(nw.groups, delta, output, verbose)
  } else if (length(nw.groups) == 3 && nw.groups[[1]]$nwgroup == 1 && nw.groups[[2]]$nwgroup == 2 && nw.groups[[3]]$nwgroup == 3) {
    make_heterogeneous_multiplex(nw.groups, delta, lambda, output, verbose)
  } else {
    # TODO: Add example of flist for homogenous and heterogeneous networks
    stop("Error: Please ensure your fList file is properly formatted", call = F)
  }

  return(0)
}
