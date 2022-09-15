########################################################################
# Find shortest paths between genes in gene sets. Given a single gene set, 
#   find the shortest paths between the genes in that gene set. Given two gene sets, 
#   find the shortest paths for pairs of genes between gene sets.
# - Input: Pre-computed multiplex network and one or two genesets
# - Output: Edge list table
# Copyright (C) 2022  David Kainer
# 
# This file is part of RWRtoolkit.
# 
# RWRtoolkit is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
# 
# RWRtoolkit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with RWRtoolkit. 
# If not, see <https://www.gnu.org/licenses/>.########################################################################

#' @importFrom dplyr %>%
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%

########################################################################
# Internal Functions
########################################################################

merge_networks = function(nw.mpo) {
    # Merge the subnetworks into one big network.
    message("Merging network layers down...")
    nl          <- nw.mpo$Number_of_Layers
    nw_dflist   <- lapply(nw.mpo[1:nl], igraph::as_data_frame)
    nw_df       <- do.call(rbind, nw_dflist)
    nw_df       <- nw_df %>% dplyr::group_by(type) %>% dplyr::mutate(weightnorm = weight/sum(weight))
    nw_merged   <- igraph::graph_from_data_frame(nw_df , directed = FALSE)
    message("Done.")
    return(nw_merged)
}

get_shortest_paths = function(nw_merged, source_geneset, target_geneset, threads=NULL) {
    # For each gene in source_geneset, get the shortest path to each gene in target_geneset.
    targets <- which(igraph::V(nw_merged)$name %in% target_geneset$gene)
    wt <- 1-igraph::E(nw_merged)$weightnorm
    # message('Beginning shortest paths ...')
    # message(class(source_geneset))
    # message(head(source_geneset))
    message(sprintf(
        "Calculating Shortest paths for %d x %d = %d gene pairs",
        nrow(source_geneset),
        length(targets),
        nrow(source_geneset)*length(targets)
    ))
    doParallel::registerDoParallel(cores=threads)
    res <- foreach::foreach(g = source_geneset$gene, .combine=rbind) %:%
        foreach::foreach(t = targets, .combine=rbind) %dopar% {
                # Only get the shortest path if the start and end nodes are not the same!
                if( igraph::V(nw_merged)$name[t] != g ) {
                    sp <- igraph::shortest_paths(nw_merged, from=g, to=t, output="vpath", weights=wt)
                    sg <- igraph::as_data_frame(igraph::induced_subgraph(nw_merged, vids = sp$vpath[[1]], impl="create_from_scratch"))
                    sg$pathname <- paste(g,igraph::V(nw_merged)$name[t],sep="_")

                    #### TODO: Discuss shortest path length with @dkainer
                    sg$pathlength <- length(sp$vpath[[1]])
                    sg
                }
            }
    doParallel::stopImplicitCluster()

    message("Finished calculating Shortest paths.")
    return(res)
}

save_results = function(rwr_result, source_geneset=NULL, target_geneset=NULL, outdir=NULL, out_path=NULL) {
    # First, generate the out_path.
    if (is.null(outdir) && is.null(out_path)) {
        warning('You must provide either `outdir` or `out_path`.')
    } else if (!is.null(out_path)) {
        out_path = out_path
    } else {
        out_path = get_file_path(
            "RWR-SPATHS",
            slug="edges",
            modname=NULL,
            outdir=outdir
        )
    }

    # If out_path is still NULL, there's nothing to do.
    if (!is.null(out_path)) {
        if (!dir.exists(dirname(out_path))) {
            dir.create(dirname(out_path), recursive=TRUE)
        }
        write_table(rwr_result, out_path)
        message(sprintf('Saved results to file:\n  %s', out_path))
    }
}

open_cytoscape = function(res, source_geneset, target_geneset) {
    # Visualize in Cyto.
    ig <- igraph::graph_from_data_frame(res, directed=F)
    igraph::V(ig)$endpoint <- igraph::V(ig)$name %in% c(source_geneset$gene,target_geneset$gene)
    RCy3::cytoscapePing()
    RCy3::createNetworkFromIgraph(ig,"myIgraph")
    RCy3::setNodeBorderColorDefault(new.color = "#666666")
    RCy3::setNodeBorderWidthDefault(new.width = 4)
}

########################################################################
# Main Function
########################################################################

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
#' @param data Path to the .Rdata file for your combo of underlying functional
#' networks. This file is produced by RWR_make_MHobject.R
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
#' extdata.dir <- system.file("example_data", package="RWRtoolkit")
#' multiplex_object_filepath <- paste(extdata.dir, '/string_interactions.Rdata', sep='')
#' geneset1_filepath <- paste(extdata.dir, '/geneset1.tsv', sep='')
#' geneset2_filepath <- paste(extdata.dir, '/geneset2.tsv', sep='')
#' outdir <- paste(extdata.dir, '/out/rwr_shortestpath', sep='')
#' 
#' rwr_shortest_path_output <- RWR_ShortestPaths(data=multiplex_object_filepath,
#'                                                  source_geneset=geneset1_filepath,
#'                                                  target_geneset=geneset2_filepath,
#'                                                  write_to_file=TRUE,
#'                                                  outdir=outdir)
#'
#'
#'
#' @export
RWR_ShortestPaths <- function(
        data=NULL,
        source_geneset=NULL,
        target_geneset=NULL,
        outdir='.',
        out_path=NULL,
        threads=parallel::detectCores-1,
        cyto=FALSE,
        verbose=FALSE,
        write_to_file=FALSE
    ) {

    # This is a saved version of the multiplex network and adj matrix.
    load(data)

    if (verbose) {
        message('Loaded data:')
        print(ls())
    }

    # If geneset is NULL, get shortest paths from source_geneset to source_geneset.
    if (is.null(target_geneset)) {
        target_geneset = source_geneset
    }
    # Gene set 1.
    message('Getting source gene set ...')
    source_geneset_plus_extras = load_geneset(source_geneset, nw.mpo, verbose)
    source_genes = source_geneset_plus_extras[[1]]
    message(head(source_genes))

    # Gene set 2.
    message('Getting target gene set ...')
    target_geneset_plus_extras = load_geneset(target_geneset, nw.mpo, verbose)
    target_genes = target_geneset_plus_extras[[1]]
    message(head(target_genes))

    nw_merged = merge_networks(nw.mpo)

    res = get_shortest_paths(nw_merged, source_genes, target_genes, threads)

    if (write_to_file) {
        save_results(res, source_geneset=source_genes, target_geneset=target_genes, outdir=outdir, out_path=out_path)
    }

    if (cyto) {
        open_cytoscape(res, source_geneset, target_geneset)
    }

    return(res)
}
