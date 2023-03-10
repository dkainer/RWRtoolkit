########################################################################
# Helper functions.
########################################################################

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

#' Load flist
#'
#' `load_flist` loads an flist file from a given path.
#' An flist file is a tab-delimited file that describes a
#' multiplex network. It should have the following three columns:
#'
#' "path to file" "short name of network" 'group"
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

#' Make a *dummy* multiplex network object.
#'
#' Create a multiplex network object from the provided flist. This
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


#' Merge a multiplex network object and keep all edges.
#'
#' Merge down all layers in a multiplex object, but don't
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
    message(sprintf("merging %d network layers ...\n", mpo$Number_of_Layers))
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
    message(sprintf("Merged network has %d edges and %d rows.\n",
            igraph::ecount(nw_merged), igraph::vcount(nw_merged)
    ))

    return(
        list(merged_network = nw_merged,
            edge_count = igraph::ecount(nw_merged),
            vertex_count = igraph::vcount(nw_merged))
    )
}

#' Merge a multiplex network object and aggregate edges.
#'
#' Merge down all layers and aggregate multi-edges into one edge
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

#' Get the network name.
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

#' Calculate basic network statistics for a network and print to screen.
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
calculate_basic_statistics <- function(
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
    title <- sprintf("Network stats for network %s\n", name)
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
        message(sprintf("Number of nodes : %d\n",   vertex_count))
        message(sprintf("Number of edges : %d\n",   edge_count))
        message(sprintf("Diameter        : %.2f\n", network_diameter))
    }

    return(data.frame(list(
        network_name = name,
        number_of_nodes = vertex_count,
        number_of_edges = edge_count,
        diameter = network_diameter
    )))
}

#' Calculate basic network statistics for each layer of a multiplex
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

    network_name_list <- c()
    number_of_nodes_list <- c()
    number_of_edges_list <- c()
    diameter_list <- c()
    for (i in seq(1, mpo$Number_of_Layers)) {
        outlist <- calculate_basic_statistics(
            mpo[[i]],
            name = names(mpo)[[i]],
            verbose = verbose)

        network_name_list <- c(network_name_list, outlist$network_name)
        number_of_nodes_list <- c(number_of_nodes_list, outlist$number_of_nodes)
        number_of_edges_list <- c(number_of_edges_list, outlist$number_of_edges)
        diameter_list <- c(diameter_list, outlist$diameter)

    }
    return(data.frame(list(
        network_name = network_name_list,
        number_of_nodes = number_of_nodes_list,
        number_of_edges = number_of_edges_list,
        diameter = diameter_list
    )))
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

#' Jaccard similarity coefficient of edges for two networks.
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
        message(sprintf("Jaccard score for edges (%s vs %s): %.2f\n",
                get_name(g, default = "<G>"),
                get_name(h, default = "<H>"),
                score))
    }
    return(score)
}

#' Sum of edge weights of the intersection divided by the number of
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
        message(sprintf("Overlap score (%s vs %s): %.2f\n",
        get_name(g, default = "<G>"),
        get_name(h, default = "<H>"),
        score))
    }
    return(score)
}

#' Compare two networks with the given metric.
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
        message(sprintf("Unknown metric: %s\n", metric))
        score <- NULL
    }

    return(score)
}

########################################################################
# Core functions.
########################################################################

#' Compare two networks with the given metric.
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

#' Calculate overlap scores between all layers of a multiplex network.
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
        for (j in seq(1, n)) {
            if (verbose > 0) {
                message(sprintf("%s vs %s ...\n", names(mpo)[i], names(mpo)[j]))
            }

            network_a <- mpo[[i]]
            network_b <- mpo[[j]]

            mat[i, j] <- compare_networks(
                            network_a,
                            network_b,
                            metric = metric,
                            verbose = verbose)

        }
    }
    rownames(mat) <- names(mpo)[1:n]
    colnames(mat) <- names(mpo)[1:n]
    if (verbose > 0) {
        message("Calculating pairwise Jaccard Indices DONE")
    }

    return(mat)
}

#' Calculate Many Vs. Reference
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

#' Calculate Tau
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
            "Tau vector: %s\n",
            paste(
                round(tau, digits = 3),
                collapse = ", ")))
    }

    return(tau)
}

#' Exclusivity
#'
#' Exclusivity is proportion of all edges in the multiplex that
#' are found in only one layer. Here we show how many edges are found in 1
#' layer, 2 layers, 3 layers etc...
#'
#' @param mpo     A multiplex object 
#' @param verbose Print progress to console.
#'
#' @return NULL
#'
#' @export
exclusivity <- function(mpo, verbose=FALSE) {
    # Get the number of edges between each node pair.
    merged <- merged_with_edgecounts(mpo)
    merged_net <- merged$merged_network
    edge_count_of_merged <- merged$edge_count

    n_found_matrix <- matrix(0, nrow = mpo$Number_of_Layers, ncol = 2)
    for (i in 1:mpo$Number_of_Layers) {
        num_edges_found <- sum(igraph::E(merged_net)$weight == i)
        exc <- round(
                num_edges_found / edge_count_of_merged,
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

write_stats_to_file_if_fp <- function(
    outdir_path = NULL,
    filename = NULL,
    netstats = NULL,
    write_rows = F,
    verbose = FALSE ){
    
    if(is.null(outdir_path)) return()

    outfile_path <- paste(outdir_path, filename,  sep = "/")

    write_table(netstats, outfile_path, write_rows, verbose = verbose)
}

write_networks_to_file_if_fp <- function(
    outdir_path = NULL,
    filename = NULL,
    netstats_network = NULL,
    verbose = FALSE ){
    
    if(is.null(outdir_path)) return()

    outfile_path <- paste(outdir_path, filename,  sep = "/")

    edgelist <- igraph::get.edgelist(netstats_network, names=T)

    from <- edgelist[,1]
    to <- edgelist[,2]
    weight <- igraph::E(netstats_network)$weight

    edge_attributes <- names(igraph::edge.attributes(netstats_network))
    
    if (length(edge_attributes) > 1){
        weight_norm <- igraph::E(netstats_network)$weightnorm
        type <- igraph::E(netstats_network)$type

        write_table(
            data.frame(
                from = from,
                to = to,
                weight = weight,
                weightnorm = weight_norm,
                type = type
            ),
            outfile_path,
            verbose
        )
        return()
    }

    write_table(
        data.frame(
            from = from,
            to = to,
            weight = weight
        ),
        outfile_path,
        verbose
    )
}

########################################################################
# Main Function
########################################################################

#' Command-line interface for RWR_netstats.R
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
#' @param outdir_path                   If a directory is supplied, all output
#'                                      data are saved to tsv within that dir.
#' @param basic_statistics              A boolean denoting a return for basic
#'                                      statistics concerning supplied networks,
#'                                      or flists. Default True.
#' @param scoring_metric                A string denoting either "jaccard", 
#'                                      "overlap" or "both" to describe. 
#'                                      Default is "jaccard". 
#' @param pairwise_between_mpo_layer    A boolean denoting a return of the 
#'                                      pairwise score (defined by 
#'                                      `scoring_metric`) between each 
#'                                      layer of the supplied
#'                                      multiplex. Default False
#' @param multiplex_layers_to_refnet    A boolean denoting a return of the
#'                                      calculated score (defined by 
#'                                      `scoring_metric`: default jaccard) 
#'                                      between a multiplex network and a 
#'                                      reference network (supplied as 
#'                                      "network_1"). Default False
#' @param net_to_net_similarity         A boolean denoting a scoring between two
#'                                      supplied networks (network_1 and
#'                                      network_2) (scoring defined by 
#'                                      `scoring_metric`: default jaccard). 
#'                                      Default False
#' @param calculate_tau                 A boolean denoting a return of the
#'                                      distribution of "tau", with respect
#'                                      to the network layers,  calculated via
#'                                      edge overlap weight / total edgeweight
#'                                      multipled by the total number of layers
#'                                      Default False
#' @param merged_with_all_edges         A boolean denoting a return of a merged
#'                                      down multiplex network along with 
#'                                      network edge counts and vertex counts.
#'                                      Default False
#' @param merged_with_edgecounts        A boolean denoting a return of a merged
#'                                      down multiplex, but simplified with
#'                                      edge weights denoting the total number
#'                                      of layers in which that edge existed.
#'                                      Default False
#' @param calculate_exclusivity_for_mpo A boolean denoting a return of total
#'                                      percentage of edges that exist within
#'                                      all n layers of the multiplex.
#'                                      Default False
#' @param verbose                       A boolean denoting the verbosity of
#'                                      output
#' @return A list of data frames, igraph networks, and vectors denoting various
#'         parameter defined statistics.
#' @examples
#' 
#' # An example of running netstats with all statistics: 
#' extdata.dir <- system.file("example_data", package="RWRtoolkit")
#' mpo_path <- paste(extdata.dir, "string_interactions.Rdata", sep = "/")
#' gold <- paste(extdata.dir, "netstat/combined_score-random-gold.tsv", sep="/")
#' test <- paste(extdata.dir, "netstat/combined_score-random-test.tsv", sep="/")
#' 
#' output_nstats <- RWR_netstats(
#'      data = mpo_path,
#'      network_1 = gold,
#'      network_2 = test,
#'      basic_statistics = T,
#'      overlap_sim_multiplex = T,
#'      overlap_sim_multiplex_layer = T,
#'      overlap_sim_multiplex_layer_jaccard = T,
#'      overlap_sim_layer_layer = T,
#'      overlap_score = T,
#'      calculate_tau_for_mpo = T,
#'      merged_with_all_edges = T,
#'      merged_with_edgecounts = T,
#'      calculate_exclusivity_for_mpo = T,
#'      outdir_path = "./",
#'      verbose = T
#' )
#' 
#' print(output_nstats)
#' 
#' @export
RWR_netstats <- function(
    data = NULL,
    flist  = NULL,
    network_1 = NULL,
    network_2 = NULL,
    outdir_path = NULL,
    basic_statistics = T,
    scoring_metric = "jaccard",
    pairwise_between_mpo_layer = F,
    multiplex_layers_to_refnet = F,
    net_to_net_similarity = F,
    calculate_tau_for_mpo = F,
    merged_with_all_edges = F,
    merged_with_edgecounts = F,
    calculate_exclusivity_for_mpo = F,
    verbose = F) {
    netstat_output <- list()

    `%notin%` <- Negate(`%in%`)
    scoring_metric <- tolower(scoring_metric)
    if (scoring_metric %notin% c("jaccard", "overlap", "both")){
        stop("Incorrect usage: scoring_metric must be one of the following: ['jaccard', 'overlap', 'both']") #nolint
    }

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




    # Load reference networks(s).
    if (!is.null(network_1)) {
        network_1 <- load_network(
            network_1,
            name = "network_1",
            verbose = verbose
        )
    } else {
        network_1 <- NULL
    }


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
            netstat_output$base_stats_net1 <- calculate_basic_statistics(
                                                network_1,
                                                verbose = verbose)
        }
        if (!is.null(network_2)) {
            netstat_output$base_stats_net2 <- calculate_basic_statistics(
                                                network_2,
                                                verbose = verbose)
        }
        if (!is.null(nw_mpo)) {
            netstat_output$base_stats_mpo <- basic_statistics_multiplex(
                                                nw_mpo,
                                                verbose = verbose)
        }
        
        #write base stats to file: 
        filestats <- rbind(
            netstat_output$base_stats_net1,
            netstat_output$base_stats_net2,
            netstat_output$base_stats_mpo
        )

        write_stats_to_file_if_fp(
                outdir_path = outdir_path,
                filename = "base_stats.tsv",
                netstats = filestats,
                write_rows = F,
                verbose = verbose)
    }


    if (pairwise_between_mpo_layer &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "pairwise_between_mpo_layer")
        ) {
        
        if (scoring_metric == "both" || scoring_metric == "jaccard"){
            netstat_output$pairwise_between_mpo_layer_jaccard <- overlap_many_pairwise(
                                                nw_mpo,
                                                metric = "jaccard",
                                                verbose = verbose)
            write_stats_to_file_if_fp(
                outdir_path = outdir_path,
                filename = "pairwise_between_mpo_layer_jaccard.tsv",
                netstats = netstat_output$pairwise_between_mpo_layer_jaccard,
                write_rows = T,
                verbose = verbose)
        }
        if (scoring_metric == "both" || scoring_metric == "overlap"){
            netstat_output$pairwise_between_mpo_layer_overlap <- overlap_many_pairwise(
                                                nw_mpo,
                                                metric = "overlap",
                                                verbose = verbose)
            write_stats_to_file_if_fp(
                outdir_path = outdir_path,
                filename = "pairwise_between_mpo_layer_overlap.tsv",
                netstats = netstat_output$pairwise_between_mpo_layer_overlap,
                write_rows = T,
                verbose = verbose)
        }
    }

    if (multiplex_layers_to_refnet &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "multiplex_layers_to_refnet") &&
        parameters_exist(
            network_1 = network_1,
            required_net = "network_1",
            function_name = "multiplex_layers_to_refnet")
        ) {
        
        if (scoring_metric == "both" || scoring_metric == "jaccard"){
            netstat_output$multiplex_layers_to_refnet_jaccard <- overlap_many_vs_reference(
                nw_mpo,
                network_1,
                metric = "overlap",
                verbose = verbose
            )

            write_stats_to_file_if_fp(
                outdir_path = outdir_path,
                filename = "multiplex_layers_to_refnet_jaccard.tsv",
                netstats = netstat_output$multiplex_layers_to_refnet_jaccard,
                write_rows = T,
                verbose = verbose)
        }
        
        if (scoring_metric == "both" || scoring_metric == "overlap"){
            netstat_output$multiplex_layers_to_refnet_overlap <- overlap_many_vs_reference(
                nw_mpo,
                network_1,
                metric = "overlap",
                verbose = verbose
            )

            write_stats_to_file_if_fp(
                outdir_path = outdir_path,
                filename = "multiplex_layers_to_refnet_overlap.tsv",
                netstats = netstat_output$multiplex_layers_to_refnet_overlap,
                write_rows = T,
                verbose = verbose)
        }
    }

    if (net_to_net_similarity &&
        parameters_exist(
            network_1 = network_1,
            required_net = "network_1",
            function_name = "net_to_net_similarity")
        && parameters_exist(
            network_2 = network_2,
            required_net = "network_2",
            function_name = "net_to_net_similarity")
        ) {

        measure_list <- list()
        if (scoring_metric == "both" || scoring_metric == "jaccard"){
            overlap_pair_jaccard <- overlap_pair(network_1,
                        network_2,
                        metric = "jaccard",
                        verbose = verbose)
            measure_list$jaccard = overlap_pair_jaccard

        }
        if (scoring_metric == "both" || scoring_metric == "overlap"){
            overlap_pair_overlap <- overlap_pair(network_1,
                        network_2,
                        metric = "overlap",
                        verbose = verbose)
            measure_list$overlap = overlap_pair_overlap
        }

        netstat_output$net_to_net_similarity <- data.frame(measure_list)


        write_stats_to_file_if_fp(
            outdir_path = outdir_path,
            filename = "net_to_net_similarity.tsv",
            netstats = netstat_output$net_to_net_similarity,
            write_rows = F,
            verbose = verbose)

    }

    if (calculate_tau_for_mpo &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "calculate_tau")) {
        netstat_output$calculated_tau <- calculate_tau(
                                                nw_mpo,
                                                network_1,
                                                verbose = verbose)

        write_stats_to_file_if_fp(
            outdir_path = outdir_path,
            filename = "calculated_tau.tsv",
            netstats = netstat_output$calculated_tau,
            write_rows = T,
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

        write_networks_to_file_if_fp(
            outdir_path,
            "merged_with_all_edges.tsv",
            netstat_output$merged_with_all_edges$merged_network,
            verbose )
    }

    if (merged_with_edgecounts  &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "merged_with_edgecounts")) {
        netstat_output$merged_with_edgecounts <- merged_with_edgecounts(
                                                nw_mpo,
                                                verbose = verbose)

        write_networks_to_file_if_fp(
            outdir_path,
            "merged_with_edgecounts.tsv",
            netstat_output$merged_with_edgecounts$merged_network,
            verbose )
    }

    if (calculate_exclusivity_for_mpo &&
        parameters_exist(
            mpo = nw_mpo,
            required_net = "mpo",
            function_name = "exclusivity")) {
        netstat_output$exclusivity <- exclusivity(
                                                nw_mpo,
                                                verbose = verbose)


        write_stats_to_file_if_fp(
            outdir_path = outdir_path,
            filename = "exclusivity.tsv",
            netstats = netstat_output$exclusivity,
            write_rows = F,
            verbose = verbose)
    }

    return(netstat_output)
}

# END
