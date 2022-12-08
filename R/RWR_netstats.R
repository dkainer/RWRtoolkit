########################################################################
# [DONE] Split this script into:
#   1. R/RWR_netstats.R
#   2. inst/scripts/run_netstats.R
# [TODO] I stripped out a lot of the verbosity from the original script. Need to add that back in.
# [TODO] What is the expected output? Print to stdout? Write tables of metrics to file(s)?
# [DONE] Docstrings.
########################################################################

########################################################################
# Helper functions.
########################################################################
## Template Roxygen documentation:
# @title Get the number of edges in a network.
#
# @description Get the number of edges in a network.
#
# @param type The edge type.
# @param verbose Print progress to console.
#
# @return An igraph network object.
#
# @export

#' Load Network
#' 
#' Loads an igraph network from a filepath to an edge list.
#'
#' @param path_to_edgelist A path to a network file in edgelist format. Must
#' have three columns: 'from', 'to', and 'weight'.
#' @param type The edge type.
#' @param name The name of the network.
#' @param col.names The column names to use for the edgelist.
#' @param select The columns to select from the edgelist.
#' @param directed Whether the network is directed.
#' @param verbose Print progress to console.
#'
#' @return An igraph network object.
load_network <- function(
    path_to_edgelist,
    type=NULL,
    name=NULL,
    col.names=c("from","to","weight"),
    select=1:3,
    header='auto',
    directed=FALSE,
    verbose=FALSE
) {
    edgelist <- data.table::fread(path_to_edgelist, col.names=col.names, select=select, header=header)

    g <- igraph::graph_from_data_frame(edgelist, directed=directed)

    if (!is.null(name)) {
        igraph::graph_attr(g, 'name') <- name
    }

    if (!is.null(type)) {
        igraph::edge_attr(g, "type") <- type
    }

    return(g)
}

#' @title Load an flist file.
#'
#' @description An flist file is a tab-delimited file that describes a
#' multiplex network. It should have the following three columns:
#'
#' <path to file> <short name of network> <group>
#'
#' 'groups' are either 1, 2, or 3.  All 1's will form one multiplex network
#' (e.g. gene-to-gene), All 2's will form a separate multiplex network (e.g.
#' disease-to-disease), And all 3's will be used to join the 1's and 2's
#' together (e.g. gene-to-disease) You don't have to have both 1's and 2's. But
#' if you do have 1's and 2's, you SHOULD have at least one 3 to join them up.
#'
#' @param path_to_flist A path to an flist file.
#' @param verbose Print progress to console.
#'
#' @return A data.table
#'
#' @export
load_flist <- function(path_to_flist, verbose=FALSE) {
    flist = data.table::fread(
        path_to_flist,
        header=F,
        col.names=c("nwfile","nwname","nwgroup"),
        select = 1:3
    )
    if (verbose) {
        print(flist)
    }
    return(flist)
}

#' @title Make a *dummy* multiplex network object.
#'
#' @description Create a multiplex network object from the provided flist. This
#' is a '*dummy*' multiplex network object because it is simply a list of
#' network layers with a single extra attribute `Number_of_Layers`.
#'
#' @param flist_or_path A data.table or path to an flist file.
#' @param verbose Print progress to console.
#'
#' @return An multiplex network object.
#'
#' @export
make_dummy_multiplex <- function(
    flist_or_path,
    verbose=FALSE
) {
    # Create a dummy multiplex object (it's just a list of igraph objects).
    # Also, add Number_of_Layers attribute, but other attributes are not added
    # (see RWR_make_multiplex.R for that).
    require(foreach)

    if ( is.data.frame(flist_or_path) ) {
        flist <- flist_or_path
    } else {
        flist <- load_flist(flist_or_path)
    }

    mpo <- foreach::foreach(row_=iterators::iter(flist, by='row')) %do%
    {
        g = load_network(
            row_$nwfile,
            type=row_$nwname,
            name=row_$nwname,
            col.names=col.names,
            select=select,
            header=header,
            directed=directed,
            verbose=verbose
        )
        g
    }
    names(mpo) <- flist$nwname

    # For the comparison functions (below), we need to be able to access the
    # network names and Number_of_Layers. The other helper functions
    # (merge_mpo) require additional attributes.
    mpo$Number_of_Layers = length(mpo)

    # Attributes that are missing from this dummy multiplex:
    # mpo$Number_of_Nodes_Multiplex  # RWR_make_multiplex.R
    # mpo$Pool_of_Nodes              # RWR_make_multiplex.R

    return(mpo)
}


#' @title Merge a multiplex network object and keep all edges.
#'
#' @description Merge down all layers in a multiplex object, but don't
#' aggregate the edges (i.e. keep all edges). This function can take a dummy
#' multiplex or a real multiplex.
#'
#' @param mpo A multiplex network object.
#' @param verbose Print progress to console.
#'
#' @return An igraph network object.
#'
#' @export
merged_with_all_edges <- function(mpo, verbose=FALSE) {
    message(sprintf("merging %d network layers ...", mpo$Number_of_Layers))
    nl        <- mpo$Number_of_Layers
    nw.dflist <- lapply(mpo[1:nl], igraph::as_data_frame)
    nw.df     <- dplyr::bind_rows(nw.dflist)
    nw.df     <- nw.df %>% dplyr::group_by(type) %>% dplyr::mutate(weightnorm = weight/sum(weight))
    nw.merged <- igraph::graph_from_data_frame(nw.df , directed = FALSE)
    message("merging network layers DONE.")
    message(sprintf("Merged network has %d edges and %d rows.", igraph::ecount(nw.merged), igraph::vcount(nw.merged) ))
    return(nw.merged)
}

#' @title Merge a multiplex network object and aggregate edges.
#'
#' @description Merge down all layers and aggregate multi-edges into one edge
#' where edge weight is the number of layers in which the two nodes are
#' connected.
#'
#' @param mpo A multiplex network object. Edge weights are the sum of the
#' weights of the edges in the multiplex network object. This CANNOT be a
#' '*dummy*' multiplex (see RWR_make_multiplex.R to create the multiplex).
#' @param inv Set the edge weight to the reciprocal of the sum of the weights.
#' @param verbose Print progress to console.
#'
#' @return An igraph network object.
#'
#' @export
merged_with_edgecounts <- function(mpo, inv=FALSE, verbose=FALSE) {
    # IMPORTANT: This one has to be a REAL multiplex (ie, use
    # RWR_make_multiplex.R, not the dummy function in this script)!
    require(Matrix)
    n <- mpo$Number_of_Nodes_Multiplex
    A <- as( matrix(0, ncol=n, nrow=n), "dgCMatrix")
    for(i in 1:mpo$Number_of_Layers)
    {
        A <- A + igraph::as_adjacency_matrix(mpo[[i]], sparse=TRUE)
    }
    aggr <- igraph::graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE, diag=FALSE)
    aggr <- igraph::set.vertex.attribute(
        aggr,
        name="name",
        index=igraph::V(aggr),
        value=mpo$Pool_of_Nodes
    )
    if (inv) {
        aggr <- igraph::set.edge.attribute(
            aggr,
            "weight",
            index = E(aggr),
            value=( 1 / get.edge.attribute(aggr, "weight", index = E(aggr)) )
        )
    }
    return(aggr)
}

#' @title Get the network name.
#'
#' @param g An igraph network object.
#' @param default The default network name.
#'
#' @return A string.
#'
#' @export
get_name <- function(g, default='<G>') {
    if (is.null(igraph::get.graph.attribute(g, 'name'))) {
        g_name = default
    } else {
        g_name = igraph::get.graph.attribute(g, 'name')
    }
    return(g_name)
}

########################################################################
# Metrics for single networks. In it's current form, this is just a report that
# prints to stderr.
########################################################################

# [TODO] number of connected components.
# [TODO] number of isolates.
# [TODO] average path length.
# [TODO] degree distribution.

#' @title Calculate basic network statistics for a network and print to screen.
#'
#' @param g An igraph object.
#' @param name A name for the network (optional).
#' @param directed Logical, whether directed or undirected paths are to be
#' considered. This is ignored for undirected graphs.
#' @param unconnected Logical, what to do if the graph is unconnected. If
#' FALSE, the function will return a number that is one larger the largest
#' possible diameter, which is always the number of vertices. If TRUE, the
#' diameters of the connected components will be calculated and the largest one
#' will be returned.
#' @param weights Optional positive weight vector for calculating weighted
#' distances. If the graph has a weight edge attribute, then this is used by
#' default.
#' @param verbose Print progress to console.
#'
#' @return NULL
#'
#' @seealso \code{\link{igraph::vcount}}, \code{\link{igraph::ecount}},
#' \code{\link{igraph::diameter}}.
#'
#' @export
basic_statistics <- function(
    g,
    name=NULL,
    directed=FALSE,
    unconnected=TRUE,
    weights=NULL,
    verbose=FALSE
) {
    # Calculate basic network statistics for monoplex network.
    # igraph::diameter
    # graph
    #     The graph to analyze.
    # directed
    #     Logical, whether directed or undirected paths are to be considered.
    #     This is ignored for undirected graphs.
    # unconnected
    #     Logical, what to do if the graph is unconnected. If FALSE, the
    #     function will return a number that is one larger the largest possible
    #     diameter, which is always the number of vertices. If TRUE, the
    #     diameters of the connected components will be calculated and the
    #     largest one will be returned.
    # weights
    #     Optional positive weight vector for calculating weighted distances.
    #     If the graph has a weight edge attribute, then this is used by
    #     default.
    if (is.null(name)) {
        name = get_name(g, default="<G>")
    }
    title = sprintf("Network stats for network %s", name)
    hrule = paste0(rep("=", nchar(title)))
    message(title)
    message(hrule)
    message(sprintf("Number of nodes : %d", igraph::vcount(g)))
    message(sprintf("Number of edges : %d", igraph::ecount(g)))
    message(sprintf("Diameter        : %.2f", igraph::diameter(g, directed=directed, unconnected=unconnected, weights=weights)))
    return(NULL)
}

#' @title Calculate basic network statistics for each layer of a multiplex
#' network and print to screen.
#'
#' @param mpo A multiplex object.
#' @param verbose Print progress to console.
#'
#' @return NULL
#'
#' @export
basic_statistics_multiplex <- function(mpo, verbose=FALSE) {
    # Calculate basic network statistics on each layer of a multiplex network.
    for (i in 1:mpo$Number_of_Layers) {
        basic_statistics(mpo[[i]], name=names(mpo)[[i]], verbose=verbose)
    }
    return(NULL)
}

########################################################################
# Metrics for comparing one network to another. To add another metric, first
# create function for your metric. That function should take at least three
# inputs: g, h, verbose. Then, add your function to the list of functions
# called by compare_networks().
########################################################################

#' @title Jaccard similarity coefficient of edges for two networks.
#'
#' @param g Reference network.
#' @param h Network of interest.
#' @param verbose Print progress to console.
#'
#' @return Jaccard similarity coefficient.
#'
#' @export
jaccard_score_edges <- function(g, h, verbose=FALSE) {
    i <- igraph::graph.intersection(g, h)
    u <- igraph::graph.union(g, h)
    score <- igraph::ecount(i) / igraph::ecount(u)
    if (verbose) {
        message(sprintf("Jaccard score for edges (%s vs %s): %.2f", get_name(g, default="<G>"), get_name(h, default="<H>"), score))
    }
    return(score)
}

#' @title Sum of edge weights of the intersection divided by the number of
#' edges in the network of interest.
#'
#' @param g Reference network.
#' @param h Network of interest.
#' @param verbose Print progress to console.
#'
#' @return Overlap score.
#'
#' @export
overlap_score <- function(g, h, verbose=FALSE) {
    # Sum of edges weights from the intersection divided by
    # the number of edges in the network of interest.
    # g : reference network, e.g., a "gold-standard"
    # h : network of interest
    # [TODO] If h is very large relative to g, then the score will be small
    # (the intersect cannot be larger than the smaller network). Is this desirable?
    i <- igraph::graph.intersection(g, h)
    sum_of_i_edge_weights <- sum( igraph::E(i)$weight_1 )
    h_n_edges <- igraph::ecount(h)
    score = sum_of_i_edge_weights / h_n_edges
    if (verbose) {
        message(sprintf("Overlap score (%s vs %s): %.2f", get_name(g, default="<G>"), get_name(h, default="<H>"), score))
    }
    return(score)
}

#' @title Compare two networks with the given metric.
#'
#' @param g Reference network.
#' @param h Network of interest.
#' @param metric Name of the metric to use for comparison.
#' @param verbose Print progress to console.
#'
#' @return Score
#'
#' @export
compare_networks <- function(g, h, metric='overlap', verbose=FALSE) {
    if (metric == 'overlap') {
        score <- overlap_score(g, h, verbose=verbose)
    } else if (metric == 'jaccard') {
        score <- jaccard_score_edges(g, h, verbose=verbose)
    } else {
        message(sprintf("Unknown metric: %s", metric))
        score = NULL
    }

    return(score)
}

########################################################################
# Core functions.
########################################################################

#' @title Compare two networks with the given metric.
#'
#' @param g Reference network.
#' @param h Network of interest.
#' @param metric Metric to use for comparison.
#' @param verbose Print progress to console.
#'
#' @return Score
#'
#' @export
overlap_pair <- function(...) {
    # Wrapper for compare_networks
    compare_networks(...)
}

#' @title Calculate overlap scores between all layers of a multiplex network.
#'
#' @param mpo A multiplex object.
#' @param metric Metric to use for comparison.
#' @param verbose Print progress to console.
#'
#' @return numeric vector
#'
#' @export
overlap_many_pairwise <- function(mpo, metric='overlap', verbose=FALSE) {
    n <- mpo$Number_of_Layers
    mat <- matrix(0,nrow=n, ncol=n)
    if (verbose>0) {
        message("Calculating pairwise Jaccard Indices ...")
    }
    for(i in 1:n){
        for(j in i:n){
            if (verbose>0) {
                message(sprintf("%s vs %s ...", names(mpo)[i], names(mpo)[j]))
            }
            mat[i,j] <- compare_networks(mpo[[i]], mpo[[j]], metric=metric, verbose=verbose)
            mat[j,i] <- mat[i,j]
        }
    }
    rownames(mat) = names(mpo)[1:n]
    colnames(mat) = names(mpo)[1:n]
    if (verbose>0) {
        message("Calculating pairwise Jaccard Indices DONE")
    }
    return(mat)
}

#' @title Calculate overlap scores between a multiplex network and a reference network.
#'
#' @param mpo multiplex network
#' @param reference_network reference network
#' @param metric metric to use for calculating overlap score
#' @param verbose Print progress to console.
#'
#' @return numeric vector
#'
#' @export
overlap_many_vs_reference <- function(mpo, reference_network, metric='overlap', verbose=FALSE) {
    # Initialize vectors.
    scores_ <- vector(length = mpo$Number_of_Layers)
    names_ <- vector(length = mpo$Number_of_Layers)
    for (i in 1:mpo$Number_of_Layers) {
        scores_[i] <- compare_networks(reference_network, mpo[[i]], metric=metric, verbose=verbose)
        names_[i] <- names(mpo)[i]
    }
    names(scores_) <- names_
    return(scores_)
}

#' @title Calculate tau vector for the given multiplex network based on the reference network.
#'
#' @param mpo Multiplex network
#' @param reference_network Reference network
#' @param verbose Print progress to console.
#'
#' @return Tau vector (numeric)
#'
#' @export
get_tau <- function(mpo, reference_network, verbose=FALSE) {
    scores <- overlap_many_vs_reference(mpo, reference_network, metric='overlap', verbose=verbose)
    if ( !is.null(scores) ) {
        tau = mpo$Number_of_Layers * (scores/sum(scores))
    } else {
        message('[WARNING] scores is NULL. Cannot calculate tau.')
        tau = NULL
    }

    if ( verbose & !is.null(tau) ) {
        message(sprintf("Tau vector: %s", paste(round(tau, digits = 3), collapse = ", ")))
    }

    return(tau)
}

#' @title Exclusivity
#'
#' @description Exclusivity is proportion of all edges in the multiplex that
#' are found in only one layer. Here we show how many edges are found in 1
#' layer, 2 layers, 3 layers etc...
#'
#' @param mpo
#' @param verbose Print progress to console.
#'
#' @return NULL
#'
#' @export
exclusivity <- function(mpo, verbose=FALSE) {
    # IMPORTANT: This one has to be a REAL multiplex (ie, use
    # RWR_make_multiplex.R, not the dummy function in this script)!
    # Get the number of edges between each node pair.
    merged <- merged_with_edgecounts(mpo)
    for (i in 1:mpo$Number_of_Layers)
    {
        exc <- round(
            sum( igraph::E(merged)$weight==i ) / igraph::ecount(merged),
            4
        )
        message(sprintf("proportion of all edges found in %d layers: %.4f", i, exc))
    }
    return(NULL)
}

########################################################################
# Command-line interface
########################################################################

#' @title Command-line interface for RWR_netstats.R
#'
#' @param opt List of options
#'
#' @return Integer.
#'
#' @export
RWR_netstats <- function(opt) {

    # Load the multiplex.
    if ( !is.null(opt$data) ) {
        # Get 'nw.mpo' from the .Rdata file.
        load(opt$data)
    } else if ( !is.null(opt$flist) ) {
        # Create 'nw.mpo' from the flist.
        nw.mpo = make_dummy_multiplex(opt$flist, verbose=opt$verbose)
    } else {
        nw.mpo = NULL
    }

    # Load reference networks(s).
    if (!is.null(opt$network_1)) {
        network_1 = load_network(
            opt$network_1,
            name='network_1',
            verbose=opt$verbose
        )
    } else {
        network_1 = NULL
    }

    if (!is.null(opt$network_2)) {
        network_2 = load_network(
            opt$network_2,
            name='network_2',
            verbose=opt$verbose
        )
    } else {
        network_2 = NULL
    }

    # Calculate the requested metrics.
    if (opt$basic_statistics) {
        if ( !is.null(network_1) ) {
            basic_statistics(network_1, verbose=opt$verbose)
        }
        if ( !is.null(network_2) ) {
            basic_statistics(network_2, verbose=opt$verbose)
        }
        if ( !is.null(nw.mpo) ) {
            basic_statistics_multiplex(nw.mpo, verbose=opt$verbose)
        }
    }

    if (opt$overlapSimilarityMultiplex) {
        overlap_many_pairwise(nw.mpo, metric='jaccard', verbose=opt$verbose)
    }

    if (opt$overlapSimilarityMultiplexLayer) {
        overlap_many_vs_reference(nw.mpo, network_1, metric='overlap', verbose=opt$verbose)
    }

    if (opt$overlapSimilarityLayerLayer) {
        overlap_pair(network_1, network_2, metric='overlap', verbose=opt$verbose)
        overlap_pair(network_1, network_2, metric='jaccard', verbose=opt$verbose)
    }

    if (opt$overlapScore) {
        overlap_many_pairwise(nw.mpo, metric='overlap', verbose=opt$verbose)
    }

    if (opt$getTau) {
        get_tau(nw.mpo, network_1, verbose=opt$verbose)
    }

    if (opt$merged_with_all_edges) {
        merged_with_all_edges(nw.mpo, verbose=opt$verbose)
    }

    if (opt$merged_with_edgecounts) {
        merged_with_edgecounts(nw.mpo, verbose=opt$verbose)
    }

    if (opt$exclusivity) {
        exclusivity(nw.mpo, verbose=opt$verbose)
    }

    return(0)
}

# END