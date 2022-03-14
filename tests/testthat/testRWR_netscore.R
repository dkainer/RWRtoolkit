# context("RWR_netscore Tests")
library(RWRtools)
library(RandomWalkRestartMH)
library(vctrs)
library(igraph)
library(mockery)

########################################################################
# Global Variables/Functions
########################################################################
gold_standard_network = 'tests/testSTRINGDB/rwr_netscore/combined_score-random-gold.tsv'
test_network = 'tests/testSTRINGDB/rwr_netscore/combined_score-random-test.tsv'
reference_geneset = 'tests/testSTRINGDB/rwr_netscore/refgenes.tsv'
outdir = 'tmp'


# path = getwd()
# while (TRUE) {
#     current_dir = basename(path)
#     if (current_dir == 'RWRtools') {
#         setwd(path)
#         message(sprintf('Found RWRtools; setting working directory to %s', path))
#         break
#     } else if (path == '/') {
#         stop('Cannot find RWRtools')
#     } else {
#         path = dirname(path)
#     }
# }
path_to_rwrtools = stringr::str_extract(getwd(), '.*RWRtools')
setwd(path_to_rwrtools)
# message(getwd())

########################################################################
# Tests
########################################################################
describe('basic test', {
    it('runs normally', {
        expect_message(
            RWRtools::RWR_netscore(
                gold=gold_standard_network,
                network=test_network,
                outdir = 'tmp'
            ),
            '.*Intersected gold network and user network.*'
        ) 
    })    
})


describe('permutations test', {
    it('runs normally', {
        expect_message(
            RWRtools::RWR_netscore(
                gold=gold_standard_network,
                network=test_network,
                outdir='tmp',
                reference_geneset=reference_geneset,
                perm=3
            ),
            '.*Permutation results.*'
        ) 
    })    
})

