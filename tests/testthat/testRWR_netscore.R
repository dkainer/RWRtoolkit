# context("RWR_netscore Tests")
library(RWRtoolkit)
library(RandomWalkRestartMH)
library(vctrs)
library(igraph)
library(mockery)

########################################################################
# Global Variables/Functions
########################################################################
gold_standard_network <- "../testSTRINGDB/rwr_netscore/combined_score-random-gold.tsv"
test_network <- "../testSTRINGDB/rwr_netscore/combined_score-random-test.tsv"
reference_geneset <- "../testSTRINGDB/rwr_netscore/refgenes.tsv"
outdir <- "tmp"

########################################################################
# Tests
########################################################################
describe("basic test", {
  it("runs normally", {
    expect_message(
      RWRtoolkit::RWR_netscore(
        gold = gold_standard_network,
        network = test_network,
        outdir = "tmp"
      ),
      ".*Intersected gold network and user network.*"
    )
  })
})


describe("permutations test", {
  it("runs normally", {
    expect_message(
      RWRtoolkit::RWR_netscore(
        gold = gold_standard_network,
        network = test_network,
        outdir = "tmp",
        reference_geneset = reference_geneset,
        perm = 3
      ),
      ".*Permutation results.*"
    )
  })
})

teardown(
  {
    # setwd(path_to_RWRtoolkit)
    system("rm -rf ./tmp")
  },
  env = parent.frame()
)
