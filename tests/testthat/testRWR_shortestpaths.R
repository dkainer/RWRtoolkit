# context("RWR_shortestpaths Tests")
library(RWRtools)
library(RandomWalkRestartMH)
library(vctrs)
library(igraph)
library(mockery)

########################################################################
# Global Variables/Functions
########################################################################
data = 'tests/testSTRINGDB/string_interactions.Rdata'
geneset1 = 'tests/testSTRINGDB/geneset1.tsv'
geneset2 = 'tests/testSTRINGDB/geneset2.tsv'
outdir = 'tmp'

# setwd('~/src/RWRtools')  # [TODO] Remove this.
path_to_rwrtools = stringr::str_extract(getwd(), '.*RWRtools')
setwd(path_to_rwrtools)

########################################################################
# Tests
########################################################################
describe('merged multiplex', {
    it('has all expected nodes', {
        expected_nodes = c(
            'CDC5L', 'ENO1', 'ENO3', 'G6PD', 'GAPDH', 'GCK', 'GFPT1', 'GFPT2',
            'GNPNAT1', 'GPI', 'H6PD', 'HK1', 'HK2', 'HK3', 'HKDC1', 'IDNK',
            'MPI', 'PFKL', 'PFKM', 'PFKP', 'PGD', 'PGLS', 'PGM1', 'PKLR',
            'PKM', 'PMM1', 'TALDO1', 'TKT', 'TPI1', 'KHK', 'PMM2'
        )
        load(data)
        geneset1_plus_extras = RWRtools::load_geneset(geneset1, nw.mpo)
        source_genes = geneset1_plus_extras[[1]]
        geneset2_plus_extras = RWRtools::load_geneset(geneset2, nw.mpo)
        target_genes = geneset2_plus_extras[[1]]
        nw_merged = RWRtools::merge_networks(nw.mpo)
        returned_nodes = names(V(nw_merged))

        expect_equal(
            expected_nodes,
            returned_nodes
        ) 
    })    
})


describe('shortest paths result', {
    it('has correct dimensions', {
        #### of dimensionality 
        load(data)
        geneset1_plus_extras = RWRtools::load_geneset(geneset1, nw.mpo)
        source_genes = geneset1_plus_extras[[1]]
        geneset2_plus_extras = RWRtools::load_geneset(geneset2, nw.mpo)
        target_genes = geneset2_plus_extras[[1]]
        nw_merged = RWRtools::merge_networks(nw.mpo)
        res = RWRtools::get_shortest_paths(nw_merged, source_genes, target_genes)

        expect_equal(
            dim(res),
            c(216, 7)
        ) 
    })    
})


describe('basic test', {
    it('runs normally', {
        expect_message(
            RWRtools::RWR_ShortestPaths(
                data=data,
                source_geneset=geneset1,
                target_geneset=geneset2,
                outdir=outdir,
                write_to_file=TRUE
            ),
            '.*Saved results to file.*'
        ) 
    })    
})


# This doesn't do anything ...
teardown(
    {
        system('rm -rf ./tmp')
    },
    env=parent.frame()
)

# END
