library(RWRtoolkit)

## Parse args function
parse_arguments <- function() {
    suppressPackageStartupMessages(library(optparse))
    option_list = list(
        make_option(
            c("-f", "--flist"),
            action="store",
            default=NULL,
            type='character',
            help="Table describing network files to use.  File columns: <path
            to file> <short name of network> <group>.  'groups' are either 1,
            2, or 3.  All 1's will form one multiplex network (e.g.
            gene-to-gene), All 2's will form a separate multiplex network (e.g.
            disease-to-disease), And all 3's will be used to join the 1's and 2's
            together (e.g. gene-to-disease) You don't have to have both 1's and
            2's.  But if you do have 1's and 2's, you SHOULD have at least one
            3 to join them up."
        ),
        make_option(
            c("-d","--delta"),
            action="store",
            default=0.5,
            type="numeric",
            help="The parameter delta sets the probability to change between
            layers at the next step. If delta = 0, the particle will always
            remain in the same layer after a non-restart iteration.  If delta = 1, 
            the particle will always change between layers, therefore not following 
            the specific edges of each layer. Default is 0.5.  Note delta must 
            be greater than 0 and less than or equal to 1."
        ),
        make_option(
            c("-l","--lambda"),
            action="store",
            default=0.5,
            type="numeric",
            help="When building a heterogeneous network (i.e. multiple layer
            groups connected with bipartite links), the walker can jump between
            layer groups with probability = lambda when it is at a node with a
            bipartite link. If lambda=1 then walker will oscillate between groups
            every time it is at a node with a bipartite link.  Default is 0.5."
        ),
        make_option(
            c("-t","--test"),
            action="store_true",
            default=FALSE,
            type="character",
            help="Run example to test script."
        ),
        make_option(
            c("-o","--out"),
            action="store",
            default="network.Rdata",
            type="character",
            help='Output file name (default "network.Rdata")'
        ),
        make_option(
            c("-v", "--verbose"),
            action="store_true",
            default=TRUE,
            help="Print extra output [default]"
        )
    )

    opt = parse_args(OptionParser(option_list=option_list))
    return(opt)
}

## Processing arguments
opt = parse_arguments()
print(opt)

## Check whether all necessary args have been set by the user
if (opt$test) {
    opt$flist=""
}

if(is.null(opt$flist)){
    stop("Error. \n - Input file \"flist\" parameter not included.\n", file=stderr())
}

if (opt$verbose) {
    # you can use either the long or short name. opt$a and opt$avar are the same.
    cat("Network files table:           ")
    cat(opt$flist)
}


RWRtoolkit::RWR_make_multiplex(flist=opt$flist,  delta=opt$delta, lambda=opt$lambda, output=opt$out, test=opt$test, verbose=opt$verbose)
