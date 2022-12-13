########################################################################
# [DONE] Split this script into:
#   1. R/RWR_netstats.R
#   2. inst/scripts/run_netstats.R
# [TODO] I stripped out a lot of the verbosity from the original script.
#        Need to add that back in.
# [TODO] What is the expected output? Print to stdout? Write tables of
#        metrics to file(s)?
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

#' @title Load Network
#'
#' @description Loads an igraph network from a filepath to an edge list.
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
    type = NULL,
    name = NULL,
    col.names = c("from", "to", "weight"),
    select = 1:3,
    header = "auto",
    directed = FALSE,
    verbose = FALSE
) {
    edgelist <- data.table::fread(
        path_to_edgelist,
        col.names = col.names,
        select = select,
        header = header)

    g <- igraph::graph_from_data_frame(edgelist, directed = directed)

    if (!is.null(name)) {
        igraph::graph_attr(g, "name") <- name
    }

    if (!is.null(type)) {
        igraph::edge_attr(g, "type") <- type
    }

    return(g)
}

#' @title Load flist
#'
#' @description `load_flist` loads an flist file from a given path.
#' An flist file is a tab-delimited file that describes a
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
    flist <- data.table::fread(
        path_to_flist,
        header = F,
        col.names = c("nwfile", "nwname", "nwgroup"),
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
    # TODO: Ask @izaakm why we are not making a real multiplex here:
    # Create a dummy multiplex object (it's just a list of igraph objects).
    # Also, add Number_of_Layers attribute, but other attributes are not added
    # (see RWR_make_multiplex.R for that).
    require(foreach)

    if (is.data.frame(flist_or_path)) {
        flist <- flist_or_path
    } else if (file.exists(flist_or_path)) {
        flist <- load_flist(flist_or_path)
    } else {
        stop("Input must be either a path to an flist or a flist dataframe.")
    }

    mpo <- foreach::foreach(row_ = iterators::iter(flist, by = "row")) %do% {
        g <- load_network(
            row_$nwfile,
            type = row_$nwname,
            name = row_$nwname
            # col.names = col.names, # TODO: ask @izaakm how these
            # select = select,       #       were intended to work
            # header = header,
            # directed = directed,
            # verbose = verbose
        )
        g
    }


    names(mpo) <- flist$nwname

    unioned_layers <- NULL
    for (layer in seq(1, length(mpo))) {
        unioned_layers <- igraph::union(mpo[[layer]], unioned_layers)
    }

    # For the comparison functions (below), we need to be able to access the
    # network names and Number_of_Layers. The other helper functions
    # (merge_mpo) require additional attributes.

    mpo$Number_of_Layers = length(mpo)                                 # nolint
    mpo$Pool_of_Nodes <- igraph::V(unioned_layers)                     # nolint
    mpo$Number_of_Nodes_Multiplex <- length(igraph::V(unioned_layers)) # nolint

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
    nw_dflist <- lapply(mpo[1:nl], igraph::as_data_frame)
    nw_df     <- dplyr::bind_rows(nw_dflist)
    nw_df     <- nw_df %>%
                    dplyr::group_by(type) %>%
                    dplyr::mutate(
                        weightnorm = weight / sum(weight)
                    )

    nw_merged <- igraph::graph_from_data_frame(nw_df, directed = FALSE)
    message("merging network layers DONE.")
    message(sprintf("Merged network has %d edges and %d rows.",
            igraph::ecount(nw_merged), igraph::vcount(nw_merged)
    ))

    return(
        list(merged_network = nw_merged,
                edge_count = igraph::ecount(nw_merged),
                vertex_count = igraph::vcount(nw_merged))
    )
}

#' @title Merge a multiplex network object and aggregate edges.
#'
#' @description Merge down all layers and aggregate multi-edges into one edge
#' where edge weight is the number of layers in which the two nodes are
#' connected.
#'
#' @param mpo A multiplex network object. Edge weights are the sum of the
#'             weights of the edges in the multiplex network object.
#' @param inv Set the edge weight to the reciprocal of the sum of the weights.
#' @param verbose Print progress to console.
#'
#' @return An igraph network object.
#'
#' @export
merged_with_edgecounts <- function(mpo, inv=FALSE, verbose=FALSE) {

    n_nodes <- mpo$Number_of_Nodes_Multiplex
    adj <- Matrix::Matrix(
                0,
                ncol = n_nodes,
                nrow = n_nodes,
                sparse = T
            )

    unioned_networks <- NULL
    for (i in 1:mpo$Number_of_Layers) {
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

    aggr
}

#' @title Get the network name.
#'
#' @param g An igraph network object.
#' @param default The default network name.
#'
#' @return A string.
#'
#' @export
get_name <- function(g, default="<G>") {
    if (is.null(igraph::get.graph.attribute(g, "name"))) {
        g_name <- default
    } else {
        g_name <- igraph::get.graph.attribute(g, "name")
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
        name <-  get_name(g, default = "<G>")
    }
    title <- sprintf("Network stats for network %s", name)
    hrule <- paste0(rep("=", nchar(title)))
    vertex_count <- igraph::vcount(g)
    edge_count <- igraph::ecount(g)
    network_diameter <- igraph::diameter(
        g,
        directed = directed,
        unconnected = unconnected,
        weights = weights
    )
    if (verbose) {
        message(title)
        message(hrule)
        message(sprintf("Number of nodes : %d",   vertex_count))
        message(sprintf("Number of edges : %d",   edge_count))
        message(sprintf("Diameter        : %.2f", network_diameter))
    }

    return(list(
        network_name = name,
        number_of_nodes = vertex_count,
        number_of_edges = edge_count,
        diameter = network_diameter
    ))
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
        basic_statistics(mpo[[i]], name = names(mpo)[[i]], verbose = verbose)
    }
    return(NULL)
}

########################################################################
# Metrics for comparing one network to another. To add another metric, first
# create function for your metric. That function should take at least three
# inputs: g, h, verbose. Then, add your function to the list of functions
# called by compare_networks().
########################################################################

check_weighted_edges <- function(network, network_name) {
    `%notin%` <- Negate(`%in%`)

    if ("weight" %notin% names(igraph::edge.attributes(network))) {
        warning(
            paste(
                "Network",
                network_name,
                "has no weighted edges. All scores will be zero."
            )
        )
    }
}

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
        message(sprintf("Jaccard score for edges (%s vs %s): %.2f",
                get_name(g, default = "<G>"),
                get_name(h, default = "<H>"),
                score))
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
    # (the intersect cannot be larger than the smaller network).
    # Is this desirable?

    check_weighted_edges(g, get_name(g, "<G>"))
    check_weighted_edges(h, get_name(h, "<H>"))

    i <- igraph::graph.intersection(g, h)
    sum_of_i_edge_weights <- sum(igraph::E(i)$weight_1)
    h_n_edges <- igraph::ecount(h)
    score <- sum_of_i_edge_weights / h_n_edges
    if (verbose) {
        message(sprintf("Overlap score (%s vs %s): %.2f",
        get_name(g, default = "<G>"),
        get_name(h, default = "<H>"),
        score))
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
compare_networks <- function(g, h, metric="overlap", verbose=FALSE) {
    if (metric == "overlap") {
        score <- overlap_score(g, h, verbose = verbose)
    } else if (metric == "jaccard") {
        score <- jaccard_score_edges(g, h, verbose = verbose)
    } else {
        message(sprintf("Unknown metric: %s", metric))
        score <- NULL
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
overlap_many_pairwise <- function(mpo, metric="overlap", verbose=FALSE) {
    n <- mpo$Number_of_Layers
    mat <- matrix(0, nrow = n, ncol = n)

    #TODO: @izaakm - is this verbose correct?
    # @verbose was verbose > 0, still curious verbosity - given that metric
    # is optional
    if (verbose) {
        message("Calculating pairwise Jaccard Indices ...")
    }

    for (i in seq(1, n)) {
        for (j in seq(i, n)) {
            if (verbose > 0) {
                message(sprintf("%s vs %s ...", names(mpo)[i], names(mpo)[j]))
            }

            network_a <- mpo[[i]]
            network_b <- mpo[[j]]

            mat[i, j] <- compare_networks(
                            network_a,
                            network_b,
                            metric = metric,
                            verbose = verbose)

            mat[j, i] <- mat[i, j]

        }
    }
    rownames(mat) <- names(mpo)[1:n]
    colnames(mat) <- names(mpo)[1:n]
    if (verbose > 0) {
        message("Calculating pairwise Jaccard Indices DONE")
    }

    return(mat)
}

#' @title Calculate Many Vs. Reference
#'
#' Calculate overlap scores between a multiplex network and a reference network.
#'
#' @param mpo multiplex network
#' @param reference_network reference network
#' @param metric metric to use for calculating overlap score
#' @param verbose Print progress to console.
#'
#' @return numeric vector
#'
#' @export
overlap_many_vs_reference <- function(
                                mpo,
                                reference_network,
                                metric = "overlap",
                                verbose = FALSE) {
    # Initialize vectors.
    scores_ <- vector(length = mpo$Number_of_Layers)
    names_ <- vector(length = mpo$Number_of_Layers)
    for (i in 1:mpo$Number_of_Layers) {
        scores_[i] <- compare_networks(
            reference_network,
            mpo[[i]],
            metric = metric,
            verbose = verbose)
        names_[i] <- names(mpo)[i]
    }
    names(scores_) <- names_
    return(scores_)
}

#' @title Calculate Tau
#'
#' Calculate tau vector for the given multiplex
#' network based on the reference network.
#'
#' @param mpo Multiplex network
#' @param reference_network Reference network
#' @param verbose Print progress to console.
#'
#' @return Tau vector (numeric)
#'
#' @export
calculate_tau <- function(mpo, reference_network, verbose=FALSE) {
    scores <- overlap_many_vs_reference(mpo,
                    reference_network,
                    metric = "overlap",
                    verbose = verbose)

    if (!is.null(scores)) {
        tau <- mpo$Number_of_Layers * (scores / sum(scores))
    } else {
        message("[WARNING] scores is NULL. Cannot calculate tau.")
        tau <- NULL
    }

    if (verbose & !is.null(tau)) {
        message(sprintf(
            "Tau vector: %s",
            paste(
                round(tau, digits = 3),
                collapse = ", ")))
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
    # Get the number of edges between each node pair.
    merged <- merged_with_edgecounts(mpo)
    edge_count_of_merged <- igraph::ecount(merged)

    n_found_matrix <- matrix(0, nrow = mpo$Number_of_Layers, ncol = 2)
    for (i in 1:mpo$Number_of_Layers) {
        sum_of_edges_in_layer <- sum(igraph::E(merged)$weight == i)
        exc <- round(
                sum_of_edges_in_layer / edge_count_of_merged,
            4
        )

        n_found_matrix[i, ] <- c(i, exc)
        if (verbose) {
            message(
                sprintf(
                    "proportion of all edges found in %d layers: %.4f",
                    i,
                    exc
                )
            )
        }
    }

    pct_found_df <- data.frame(n_found_matrix)
    colnames(pct_found_df) <- c("n_layers", "pct_found")
    return(pct_found_df)
}


parameters_exist <- function(mpo = NULL,
                            flist = NULL,
                            network_1 = NULL,
                            network_2 = NULL,
                            required_net,
                            function_name) {

    return_val <- FALSE
    if (required_net == "mpo") return_val <- !is.null(mpo)
    if (required_net == "flist") return_val <- !is.null(flist)
    if (required_net == "network_1") return_val <- !is.null(network_1)
    if (required_net == "network_2") return_val <- !is.null(network_2)

    if (!return_val)
        warning(
            paste(
                required_net,
                "required for ",
                function_name
            )
        )
    return_val
}


########################################################################
# Main Function
########################################################################

#' @title Command-line interface for RWR_netstats.R
#'
#' @param data                          The filepath to an mpo object.
#' @param flist                         An flist. Currently creates "faux mpo"
#'                                      object
#'                                      for fast statistical inferences.
#' @param network_1                     A path to an edgelist. Used for basic
#'                                      statistics, overlap_sim_multiplex_layer,
#'                                      overlap_pair, and calculate tau
#' @param network_2                     A path to an edgelist. Used for basic
#'                                      statistics and overlap_pair.
#' @param basic_statistics              A boolean denoting a return for basic
#'                                      statistics concerning supplied networks,
#'                                      or flists.
#' @param overlap_sim_multiplex         A boolean denoting a return of jaccard
#'                                      similarity metrics for the supplied
#'                                      multiplex
#' @param overlap_sim_multiplex_layer   A boolean denoting a return of the
#'                                      calculated edge weight overlap between
#'                                      a multiplex network and a reference
#'                                      network (supplied as "network_1")
#' @param overlap_sim_layer_layer       A boolean denoting a return of jaccard
#'                                      and edge weight overlap between two
#'                                      supplied networks (network_1 and
#'                                      network_2)
#' @param overlap_score                 A boolean denoting a return of a matrix
#'                                      of overlap scores between all layers
#'                                      of the multiplex.
#' @param calculate_tau                 A boolean denoting a return of the
#'                                      distribution of "tau", with respect
#'                                      to the network layers,  calculated via
#'                                      edge overlap weight / total edgeweight
#'                                      multipled by the total number of layers
#' @param merged_with_all_edges         A boolean denoting a return of a merged
#'                                      down multiplex network along with 
#'                                      network edge counts and vertex counts.
#' @param merged_with_edgecounts        A boolean denoting a return of a merged
#'                                      down multiplex, but simplified with 
#'                                      edge weights denoting the total number
#'                                      of layers in which that edge existed.
#' @param exclusivity                   A boolean denoting a return of total
#'                                      percentage of edges that exist within
#'                                      all n layers of the multiplex.
#' @param verbose                       A boolean denoting the verbosity of
#'                                      output
#' @return A list of data frames, igraph networks, and vectors denoting various
#'         parameter defined statistics.
#'
#' @export
RWR_netstats <- function(
    data = NULL,
    flist  = NULL,
    network_1 = NULL,
    network_2 = NULL,
    basic_statistics = F,
    overlap_sim_multiplex = F,
    overlap_sim_multiplex_layer = F,
    overlap_sim_layer_layer = F,
    overlap_score = F,
    calculate_tau = F,
    merged_with_all_edges = F,
    merged_with_edgecounts = F,
    exclusivity = F,
    verbose = F) {

    netstat_output <- list()

    # Load the multiplex.
    if ( !is.null(data) ) {
        # Get 'nw_mpo' from the .Rdata file.
        nw_mpo <- RWRtoolkit::load_multiplex_data(data)
    } else if (!is.null(flist)) {
        # Create 'nw_mpo' from the flist.
        nw_mpo <- make_dummy_multiplex(flist, verbose = verbose)
    } else {
        nw_mpo <- NULL
    }

    # Load reference networks(s).
    print("Network 1")
    print(network_1)
    if (!is.null(network_1)) {
        network_1 <- load_network(
            network_1,
            name = "network_1",
            verbose = verbose
        )
    } else {
        network_1 <- NULL
    }
    print(network_1)

    if (!is.null(network_2)) {
        network_2 <- load_network(
            network_2,
            name = "network_2",
            verbose = verbose
        )
    } else {
        network_2 <- NULL
    }


    # ensure that something is available:
    if (is.null(data) &
        is.null(flist) &
        is.null(network_1) &
        is.null(network_2)
        ) {

        stop( paste(
            "[ERROR] You must supply one of the following arguments:",
            "       data",
            "       flist",
            "       (data or flist) and network_1",
            "       network_1 and/or network_2",
            sep = "\n"
        ))
    }
    # Calculate the requested metrics.
    if (basic_statistics) {
        if (!is.null(network_1)) {
            netstat_output$base_stats_net1 <- basic_statistics(
                                                network_1,
                                                verbose = verbose)
        }
        if (!is.null(network_2)) {
            netstat_output$base_stats_net2 <- basic_statistics(
                                                network_2,
                                                verbose = verbose)
        }
        if (!is.null(nw_mpo)) {
            netstat_output$base_stats_mpo <- basic_statistics_multiplex(
                                                nw_mpo,
                                                verbose = verbose)
        }
    }

    if (overlap_sim_multiplex &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "overlap_sim_multiplex")
        ) {
        netstat_output$overlap_sim_multiplex_jaccard <- overlap_many_pairwise(
                                            nw_mpo,
                                            metric = "jaccard",
                                            verbose = verbose)
    }

    if (overlap_sim_multiplex_layer &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "overlap_sim_multiplex_layer") &&
        parameters_exist(
            network_1 = network_1,
            required_net = "network_1",
            function_name = "overlap_sim_multiplex_layer")
        ) {

        netstat_output$overlap_sim_multiplex_layer <- overlap_many_vs_reference(
            nw_mpo,
            network_1,
            metric = "overlap",
            verbose = verbose
        )
    }

    if (overlap_sim_layer_layer &&
        parameters_exist(
            network_1 = network_1,
            required_net = "network_1",
            function_name = "overlap_sim_layer_layer")
        && parameters_exist(
            network_2 = network_2,
            required_net = "network_2",
            function_name = "overlap_sim_layer_layer")
        ) {
        netstat_output$overlap_pair_overlap_weight <- overlap_pair(network_1,
                     network_2,
                     metric = "overlap",
                     verbose = verbose)
        netstat_output$overlap_pair_jaccard <- overlap_pair(network_1,
                     network_2,
                     metric = "jaccard",
                     verbose = verbose)
    }

    if (overlap_score  &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "mpo_overlap_score")) {
        netstat_output$mpo_overlap_score <- overlap_many_pairwise(
                                                        nw_mpo,
                                                        metric = "overlap",
                                                        verbose = verbose)
    }

    if (calculate_tau &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "calculate_tau")) {
        netstat_output$calculated_tau <- calculate_tau(
                                                nw_mpo,
                                                network_1,
                                                verbose = verbose)
    }

    if (merged_with_all_edges &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "merged_with_all_edges")) {
        netstat_output$merged_with_all_edges <- merged_with_all_edges(
                                                        nw_mpo,
                                                        verbose = verbose)
    }

    if (merged_with_edgecounts  &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "merged_with_edgecounts")) {
        netstat_output$merged_with_edgecounts <- merged_with_edgecounts(
                                                nw_mpo,
                                                verbose = verbose)
    }

    if (exclusivity &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "exclusivity")) {
        netstat_output$exclusivity <- exclusivity(
                                                nw_mpo,
                                                verbose = verbose)
    }

    return(netstat_output)
}

# END
