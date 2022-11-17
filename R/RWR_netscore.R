########################################################################
# Performs a network intersect between an input network and a gold truth network. It scores the strength of that intersect with multiple metrics.
# - Input: A gold standard network as reference and another network to compare to the gold standard.
# - Output: Table of metrics as a comparison between networks
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
# Main Function
########################################################################

#' RWR netscore
#'
#' `RWR_netscore` This tool performs a network intersect between your input
#' network (--network) and a gold truth network (--gold). It scores the
#' strength of that intersect with multiple metrics.
#'
#' @param gold Full path to the gold network file. This is the network that
#' contains known truth edges.  Input file is a table with first 3 cols {<}GeneA{>}
#' {<}GeneB{>} {<}score{>} column names don't matter.
#' @param network Full path to your network to be evaluated. This is a table
#' with first 3 cols {<}GeneA{>} {<}GeneB{>} {<}score{>} column names don't matter.
#' @param reference_geneset Only required if using `perm`. Full path to a file
#' containing all gene/node labels that were used in the process of generating
#' the network provided in `network`.
#' @param perm Number of random permutations to use for estimating
#' significance. If you set this, you must also set `refgenes`.
#' @param threads Only set this is you want to do permutations.
#' @param outdir Output directory. Default "."
#' @param write_to_file Also write the result to a file. Default FALSE
#' @param verbose Verbose mode. Default FALSE
#' @return Returns a table of metrics as a comparison between networks with parameters:
#'              network:  Original Network input.
#'              gold:   Gold Network comparison. 
#'              network_nodes : Number of connected nodes within network (i.e. singleton nodes will not be counted)
#'              gold_nodes:     Number of nodes within the gold network
#'              network_edges:  Total number of edges within the network. 
#'              gold_edges:     Total number of edges within the gold network
#'              nodes_in_common:  Number of nodes existing within both input and gold networks.
#'              intersecting_edges: Number of edges existing within both input and gold networks.
#'              intersect_perc: Percentage of intersecting nodes 
#'              intersect_score:  @dkainer
#'              score_per_network_edge: @dkainer
#'
#' @examples
#' 
#' extdata.dir <- system.file("example_data", package="RWRtoolkit")
#' multiplex_object_filepath <- paste(extdata.dir, '/string_interactions.Rdata', sep='')
#' geneset1_filepath <- paste(extdata.dir, '/geneset1.tsv', sep='')
#' geneset2_filepath <- paste(extdata.dir, '/geneset2.tsv', sep='')
#' outdir <- paste(extdata.dir, '/out/rwr_netscore', sep='') 
#' 
#' # An example of Running RWR_ShortestPaths
#' netscore_table <- RWRtoolkit::RWR_ShortestPaths(data=multiplex_object_filepath,
#'                                               source_geneset=geneset1_filepath,
#'                                               target_geneset=geneset2_filepath,
#'                                               write_to_file=TRUE,
#'                                               outdir=outdir)
#'
#'
#'
#' @export
RWR_netscore <- function(
    gold=NULL,
    network=NULL,
    reference_geneset=NULL,
    perm=0,
    threads=NULL,
    outdir='.',
    write_to_file=FALSE,
    verbose=FALSE
    ) {
    
    # Load the GOLD network data and make it into an iGraph object.
    message("loading GOLD network...")
    GOLD        <- data.table::fread(gold, select = 1:3, col.names = c("GeneA", "GeneB", "score"))
    GOLD        <- dplyr::filter(GOLD, !is.na(score))
    GOLD.ig <- igraph::graph_from_data_frame(d=GOLD, directed=F)
    
    # Load the user network and make it into an iGraph object.
    message("loading user network...")
    NW      <- data.table::fread(network, select = 1:3, col.names = c("GeneA", "GeneB", "score"))
    NW      <- dplyr::filter(NW, !is.na(score))
    NW.ig <- igraph::graph_from_data_frame(d=NW, directed=F)
    
    # Intersect the GOLD network and user network.
    message("getting network intersection...")
    intersect <- igraph::intersection(GOLD.ig, NW.ig, keep.all.vertices=F)

    # Ignore the first element 'score':
    vec = igraph::E(intersect)$score_1
    mask = vec != 'score'
    score = sum(as.numeric(vec[mask]))
    
    message("Intersected gold network and user network")
    message("=======================================")
    message(sprintf("GOLD nodes                  : %d", igraph::vcount(GOLD.ig)))
    message(sprintf("Network nodes               : %d", igraph::vcount(NW.ig)))
    message(sprintf("GOLD edges                  : %d", igraph::ecount(GOLD.ig)))
    message(sprintf("Network edges               : %d", igraph::ecount(NW.ig)))
    message(sprintf("Nodes in common             : %d", igraph::vcount(intersect)))
    message(sprintf("Intersecting edges          : %d", igraph::ecount(intersect)))
    message(sprintf("Percent intersected         : %.1f%%", 100*igraph::ecount(intersect) / igraph::ecount(NW.ig)))
    message(sprintf("Intersect score             : %.2f", score))
    message(sprintf("Score per network edge      : %.2f", score / igraph::ecount(NW.ig)))
    message(sprintf("Score per intersecting edge : %.2f", score / igraph::ecount(intersect)))
    
    out <- data.frame(
        network = basename(network),
        gold = basename(gold),
        network_nodes = igraph::vcount(NW.ig),
        gold_nodes = igraph::vcount(GOLD.ig),
        network_edges = igraph::ecount(NW.ig),
        gold_edges = igraph::ecount(GOLD.ig),
        nodes_in_common = igraph::vcount(intersect),
        intersecting_edges = igraph::ecount(intersect),
        intersect_perc = 100*igraph::ecount(intersect) / igraph::ecount(NW.ig),
        intersect_score = score,
        score_per_network_edge = score / igraph::ecount(NW.ig),
        score_per_intersect_edge = score / igraph::ecount(intersect),
        type = "test"
    )
    
    if (perm > 0) {
        message(sprintf("Running %d random permutations ...", perm))
        
        # Load the reference genes list
        refgenes = load_geneset(reference_geneset)$genes$gene
        
        doParallel::registerDoParallel(cores=threads)
        res <- foreach::foreach(i = 1:perm, .combine=rbind) %dopar% {
            gperm <- NW.ig
            igraph::V(gperm)$name <- sample(refgenes, igraph::vcount(gperm))
            gperm.isec  <- igraph::intersection(GOLD.ig, gperm )
            gperm.score <- sum(igraph::E(gperm.isec)$score_1)
            data.frame(
                network = basename(network),
                GOLD = basename(gold),
                network_edges = igraph::ecount(NW.ig),
                intersecting_edges = igraph::ecount(gperm.isec),
                intersect_perc = 100*igraph::ecount(gperm.isec) / igraph::ecount(gperm),
                intersect_score = gperm.score,
                score_per_network_edge = gperm.score / igraph::ecount(gperm),
                score_per_intersect_edge = gperm.score / igraph::ecount(gperm.isec),
                type = "permutation"
            )
        }
        doParallel::stopImplicitCluster()
        
        message('Permutation results')
        message("=======================================")
        message(sprintf("average intersect score = %.2f", mean(res$intersect_score)))
        message(sprintf("average intersect pecentage = %.2f", mean(res$intersect_perc)))
        message(sprintf("average score per network edge = %.2f", mean(res$score_per_network_edge)))
        message(sprintf("average score per intersected edge = %.2f", mean(res$score_per_intersect_edge)))
        message(sprintf("test network fc(score) compared to permutations = %.2f", out$intersect_score/mean(res$intersect_score)))

        out = dplyr::bind_rows(out, res)
    }
    
    if (write_to_file) {
        out_path = get_file_path(
             "NETSCORE_",
             tools::file_path_sans_ext(basename(network)),
             outdir=outdir
        )
        write_table(out, out_path, verbose=verbose)
    }

    return(out)
}
