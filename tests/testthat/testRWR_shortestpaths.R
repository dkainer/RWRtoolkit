context("RWR_shortestpaths Tests")
library(RandomWalkRestartMH)
library(vctrs)
library(igraph)
library(mockery)

########################################################################
# Global Variables/Functions
########################################################################
data <- "../testSTRINGDB/string_interactions.Rdata"
geneset1 <- "../testSTRINGDB/geneset1.tsv"
geneset2 <- "../testSTRINGDB/geneset2.tsv"
outdir <- "tmp"

########################################################################
# Tests
########################################################################


describe("shortest paths result", {
  it("has correct dimensions", {
    #### of dimensionality
    load(data)
    geneset1_plus_extras <- load_geneset(geneset1, nw.mpo)
    source_genes <- geneset1_plus_extras[[1]]
    geneset2_plus_extras <- load_geneset(geneset2, nw.mpo)
    target_genes <- geneset2_plus_extras[[1]]
    nw_merged <- merged_with_all_edges(nw.mpo)$merged_network
    res <- get_shortest_paths(
      nw_merged,
      source_genes,
      target_genes
    )

    expect_equal(
      dim(res),
      c(216, 8)
    )
  })
})


describe("basic test", {
  it("runs normally and saves to file", {
    expect_message(
      RWRtoolkit::RWR_ShortestPaths(
        data = data,
        source_geneset = geneset1,
        target_geneset = geneset2,
        outdir = outdir,
        write_to_file = TRUE
      ),
      ".*Saved results to file.*"
    )
  })
})

# ABC Layer 1      # ABC Layer 2            # ABC Layer 3
# A◄───B◄──────C   # G        H             #    H◄─┐
#      ▲       ▲   # │        │             #       │
#      │       │   # └──►D    └─►E       F  #  D◄───E◄┐        F
# D◄───E◄─┬──F─┤   #     │       │       │  #         │        │
#      ▲  │  ▲ │   #     └───►B◄─┘──►C◄──┘  #         ├───►C◄──┘
#      │  │  │ │   #          │      │      #         │
#  G───┘  H ─┴─┘              └──►A◄─┘      #         A

# GET SHORTEST PATH FROM G TO C
describe("test RWR_shortestpaths", {
  it("extracts the shortest paths between two elements", {
    abc_multiplex_filepath <- "../testMultiplex/abc_multiplex.Rdata"
    source_set_filepath <- "../testGenesets/test_abc_genesetG.tsv"
    sink_set_filepath <- "../testGenesets/test_abc_genesetA.tsv"

    actual <- RWRtoolkit::RWR_ShortestPaths(
      data = abc_multiplex_filepath,
      source_geneset = source_set_filepath,
      target_geneset = sink_set_filepath
    )

    expected_pathelements <- "G->E->A"
    expected_pathname <- "G_A"
    expected_pathlength <- 3
    expected_from <- c("E", "E")
    expected_to <- c("G", "A")
    expected_weight <- c(1, 1)
    expected_weightnorm <- c(0.1, 0.2)
    expected_type <- c("layer1", "layer3")

    expected <- data.frame(list(
      from = expected_from,
      to = expected_to,
      weight = expected_weight,
      type = expected_type,
      weightnorm = expected_weightnorm,
      pathname = expected_pathname,
      pathlength = expected_pathlength,
      pathelements = expected_pathelements
    ))

    expect_equal(actual, expected)
  })

  it('extracts multiple paths if they exist', {
    circular_multiplex <- "../testMultiplex/circular_network.Rdata"
    source_set_filepath <- "../testGenesets/test_circular_geneset_1.tsv"
    target_set_filepath <- "../testGenesets/test_circular_geneset_2.tsv"

    actual <- RWRtoolkit::RWR_ShortestPaths(
      data = circular_multiplex,
      source_geneset = source_set_filepath,
      target_geneset = target_set_filepath
    )

    expected_pathelements <- "A->B->D"
    expected_pathname <- "A_D"
    expected_pathlength <- 3
    expected_from <- c("A", "B")
    expected_to <- c("B", "D")
    expected_weight <- c(1, 1)
    expected_weightnorm <- c(0.25, 0.25)
    expected_type <- "layer1"

    expected_1 <- data.frame(list(
      from = expected_from,
      to = expected_to,
      weight = expected_weight,
      type = expected_type,
      weightnorm = expected_weightnorm,
      pathname = expected_pathname,
      pathlength = expected_pathlength,
      pathelements = expected_pathelements
    ))
    
    expected_pathelements <- "A->C->D"
    expected_pathname <- "A_D"
    expected_pathlength <- 3
    expected_from <- c("A", "C")
    expected_to <- c("C", "D")
    expected_weight <- c(1, 1)
    expected_weightnorm <- c(0.25, 0.25)
    expected_type <- "layer1"

    expected_2 <- data.frame(list(
      from = expected_from,
      to = expected_to,
      weight = expected_weight,
      type = expected_type,
      weightnorm = expected_weightnorm,
      pathname = expected_pathname,
      pathlength = expected_pathlength,
      pathelements = expected_pathelements
    ))

    expected <- rbind(expected_1, expected_2)

    expect_equal(actual, expected)
  })
})



describe("tests extract_node_from_row", {
  path_elements <- "G->E->A"
  pathname <- "G_A"
  pathlength <- 3
  from_nodes <- c("E", "E")
  to_nodes <- c("G", "A")
  edgeweight <- c(1, 1)
  edgeweightnorm <- c(0.1, 0.2)
  type <- c("layer1", "layer3")

  shortest_path_df <- data.frame(list(
    from = from_nodes,
    to = to_nodes,
    weight = edgeweight,
    type = type,
    weightnorm = edgeweightnorm,
    pathname = pathname,
    pathlength = pathlength,
    pathelements = path_elements
  ))

  it("extracts the 'to' and 'from' nodes from a shortest path df", {
    starting_node <- "G"

    actual <- extract_node_from_row(shortest_path_df, starting_node)

    expected <- list(
      from_column = "to",
      from = starting_node,
      to_column = "from",
      to = "E"
    )

    expect_equal(actual, expected)
  })



  # it("extracts the 'to' and 'from' nodes from a shortest path df", {
  #     starting_node <- "E"

  #     actual <- extract_node_from_row(shortest_path_df, starting_node)

  #     expected <- list(
  #         from_colum = "from",
  #         from = starting_node,
  #         to_column = "to",
  #         to = "A"
  #     )

  #     expect_equal(actual, expected)
  # })
})


teardown(
  {
    system("rm -rf ./tmp")
  },
  env = parent.frame()
)

# END
