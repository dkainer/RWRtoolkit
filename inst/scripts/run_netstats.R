########################################################################
# Parse arguments
########################################################################

parse_arguments <- function() {
    # [TODO] Need to verify 'help' messages (I ripped several of these from
    # other scripts followed by liberal copy/pasta).
    suppressPackageStartupMessages(library(optparse))
    option_list = list(
        make_option(
            c("-f", "--flist"),
            action = "store",
            default = NULL,
            type = "character",
            help = "Table describing network files to use. File columns:

                <path to file> <short name of network> <group>

                'groups' are either 1, 2, or 3.  All 1's will form one
                multiplex network (e.g. gene-to-gene), All 2's will form a
                separate multiplex network (e.g. disease-to-disease), And all
                3's will be used to join the 1's and 2's together (e.g.
                gene-to-disease) You don't have to have both 1's and 2's. But
                if you do have 1's and 2's, you SHOULD have at least one 3 to
                join them up."
        ),
        make_option(c("-d", "--data"),
                action = "store",
                default = NULL,
                type = "character",
                help = "path to the .Rdata file for your combo of
                    underlying functional networks. This file is produced
                    by RWR_make_MHobject.R"
        ),
        make_option(
            c("-o", "--outdir_path"),
            action = "store",
            default = "netstats",
            type = "character",
            help = 'Output file name (default "netstats")'
        ),
        make_option(
            c("--network_1"),
            action = "store",
            default = NULL,
            type = "character",
            help = "A path to an edgelist. Used for basic statistics, overlap_sim_multiplex_layer, overlap_pair, and calculate tau"
        ),
        make_option(
            c("--network_2"),
            action = "store",
            default = NULL,
            type = "character",
            help = "A path to an edgelist. Used for basic statistics and overlap_pair."
        ),
        make_option(
            c("--basic_statistics"),
            action = "store_true",
            default = FALSE,
            help = "A flag denoting a return for basic statistics concerning supplied networks, or flists."
        ),
          make_option(
            c("-s", "--scoring_metric"),
            action = "store",
            default = "jaccard",
            type = "character",
            help = "A metric used for scoring \"jaccard\" similarity, edge weight \"overlap\" or \"both\". Default \"jaccard\"."
        ),
        make_option(
            c("--pairwise_between_mpo_layer"),
            action = "store_true",
            default = FALSE,
            help = "A boolean denoting a return of the pairwise score (defined by `scoring_metric`) between each layer of the suppliedmultiplex. Default False"
        ),
        make_option(
            c("--multiplex_layers_to_refnet"),
            action = "store_true",
            default = FALSE,
            help = " A boolean denoting a return of the calculated score (defined by `scoring_metric`: default jaccard) between a multiplex network and a reference network (supplied as \"network_1\"). Default False"
        ),
        make_option(
            c("--net_to_net_similarity"),
            action = "store_true",
            default = FALSE,
            help = "A boolean denoting a scoring between two supplied networks (network_1 and network_2) (scoring defined by `scoring_metric`: default jaccard). Default False"
        ),
        make_option(
            c("--calculate_tau_for_mpo"),
            action = "store_true",
            default = FALSE,
            help = "A boolean denoting a return of the distribution of 'tau', with respect to the network layers,  calculated via edge overlap weight / total edgeweight multipled by the total number of layers"
        ),
        make_option(
            c("--merged_with_all_edges"),
            action = "store_true",
            default = FALSE,
            help = "A boolean denoting a return of a merged down multiplex network along with  network edge counts and vertex counts."
        ),
        make_option(
            c("--merged_with_edgecounts"),
            action = "store_true",
            default = FALSE,
            help = "A boolean denoting a return of a merged down multiplex, but simplified with edge weights denoting the total number of layers in which that edge existed."
        ),
        make_option(
            c("--calculate_exclusivity_for_mpo"),
            action = "store_true",
            default = FALSE,
            help = "A boolean denoting a return of total percentage of edges that exist within all n layers of the multiplex."
        ),
        make_option(
            c("-v", "--verbose"),
            action = "store_true",
            default = FALSE,
            help = "A boolean denoting the verbosity of output."
        )
    )

    opt <- optparse::parse_args(
                optparse::OptionParser(
                    option_list = option_list
                )
            )

    # Verify opt.
    errors <- FALSE

    # Require that the user passes either --data or --flist (not both).
    if (!is.null(opt$data) & !is.null(opt$flist)) {
        errors <- TRUE
        message(
            "\n[ERROR] You may provide either --data or --flist but not both.\n"
        )
    }

    # Require that the user provided at least one input.
    if (is.null(opt$data) &
        is.null(opt$flist) &
        is.null(opt$network_1) &
        is.null(opt$network_2)
        ) {

        errors <- TRUE
        message(
            "\n[ERROR] You must supply one of the following arguments:\n",
            "  --data\n",
            "  --flist\n",
            "  (--data or --flist) and network_1\n",
            "  --network_1 and/or --network_2"
        )
    }

    # Require that --network_1 is used with --data or --flist.
    if ((!is.null(opt$data) | !is.null(opt$flist)) &
         (is.null(opt$network_1) & !is.null(opt$network_2))) {
        errors <- TRUE
        message(
            "\n[ERROR] You must provide the reference network as --network_1.\n"
        )
    }

    if (errors) {
        stop("Exiting due to errors (see above).")
    }

    return(opt)
}

########################################################################
# Main
########################################################################

main <- function() {
    # Get args.
    
    opt <- parse_arguments()
    print("RUNNING NETSTATS WITH ")
    print(opt)
    RWRtoolkit::RWR_netstats(
        data = opt$data,
        flist = opt$flist,
        network_1 = opt$network_1,
        network_2 = opt$network_2,
        outdir_path = opt$outdir_path,
        basic_statistics = opt$basic_statistics,
        scoring_metric = opt$scoring_metric,
        pairwise_between_mpo_layer = opt$pairwise_between_mpo_layer,
        multiplex_layers_to_refnet = opt$multiplex_layers_to_refnet,
        net_to_net_similarity = opt$net_to_net_similarity,
        calculate_tau_for_mpo = opt$calculate_tau_for_mpo,
        merged_with_edgecounts = opt$merged_with_edgecounts,
        merged_with_all_edges = opt$merged_with_all_edges,
        calculate_exclusivity_for_mpo = opt$calculate_exclusivity_for_mpo,
        verbose = opt$verbose
    )
}

status  <-  main()
quit(save = "no", status = status)
