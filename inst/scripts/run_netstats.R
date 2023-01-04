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
            c("-o", "--outdir"),
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
            help = "Reference network to use for calculating Tau."
        ),
        make_option(
            c("--network_2"),
            action = "store",
            default = NULL,
            type = "character",
            help = "network_2"
        ),
        make_option(
            c("--basic_statistics"),
            action = "store_true",
            default = FALSE,
            help = "basic_statistics"
        ),
        make_option(
            c("--overlap_sim_multiplex"),
            action = "store_true",
            default = FALSE,
            help = "overlap_sim_multiplex"
        ),
        make_option(
            c("--overlap_sim_multiplex_layer"),
            action = "store_true",
            default = FALSE,
            help = "overlap_sim_multiplex_layer"
        ),
        make_option(
            c("--overlap_sim_layer_layer"),
            action = "store_true",
            default = FALSE,
            help = "overlap_sim_layer_layer"
        ),
        make_option(
            c("--overlapScore"),
            action = "store_true",
            default = FALSE,
            help = "overlapScore"
        ),
        make_option(
            c("--getTau"),
            action = "store_true",
            default = FALSE,
            help = "getTau"
        ),
        make_option(
            c("--merged_with_all_edges"),
            action = "store_true",
            default = FALSE,
            help = "merged_with_all_edges"
        ),
        make_option(
            c("--merged_with_edgecounts"),
            action = "store_true",
            default = FALSE,
            help = "merged_with_edgecounts"
        ),
        make_option(
            c("--exclusivity"),
            action = "store_true",
            default = FALSE,
            help = "exclusivity"
        ),
        make_option(
            c("-v", "--verbose"),
            action = "store_true",
            default = TRUE,
            help = "Print extra output [default]"
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

    RWRtoolkit::RWR_netstats(
        data = opt$data,
        flist = opt$flist,
        network_1 = opt$network_1,
        network_2 = opt$network_2,
        basic_statistics = opt$basic_statistics,
        overlap_sim_multiplex = opt$overlap_sim_multiplex,
        overlap_sim_multiplex_layer = opt$overlap_sim_multiplex_layer,
        overlap_sim_layer_layer = opt$overlap_sim_layer_layer,
        overlap_score = opt$overlap_score,
        calculate_tau = opt$calculate_tau,
        merged_with_all_edges = opt$merged_with_all_edges,
        merged_with_edgecounts = opt$merged_with_edgecounts,
        exclusivity = opt$exclusivity,
        outdir_path = opt$outdir,
        verbose = opt$verbose
    )
}

status  <-  main()
quit(save = "no", status = status)
