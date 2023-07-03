context("Random Walk Restart Tests")
library(RWRtoolkit)
library(RandomWalkRestartMH)
library(vctrs)
library(igraph)
library(mockery)

# Homogenous Network tibble
nw_tibble <- tibble::tibble(
  "nwfile" = c("../testNetworks/m1.txt", "../testNetworks/m2.txt"),
  "nwname" = c("m1", "m2"),
  "nwgroup" = c(1, 1)
)

nw_tibble_bad_path <- tibble::tibble(
  "nwfile" = c("../testNetworks/iDontExist_m1.txt", "../whereever/m2.txt"),
  "nwname" = c("m1", "m2"),
  "nwgroup" = c(1, 1)
)

nw_groups <- list_of(nw_tibble)

# Heterogeneous Network Tibble
nw_tibble_het1 <- tibble::tibble(
  "nwfile" = c("../testNetworks/m1.txt"),
  "nwname" = c("m1"),
  "nwgroup" = c(1)
)
nw_tibble_het2 <- tibble::tibble(
  "nwfile" = c("../testNetworks/n1.txt"),
  "nwname" = c("n1"),
  "nwgroup" = c(2)
)
nw_tibble_het3 <- tibble::tibble(
  "nwfile" = c("../testNetworks/i1.txt"),
  "nwname" = c("i1"),
  "nwgroup" = c(3)
)


# generates expected matrix for test data from graphs m1 and m2
generate_expected_supraadj <- function(delta) {
  one_minus_delta <- 1 - delta
  # Because of how RandomWalkRestartMH Calculates the edge weights
  # our graph, m1, originally has weights 1 and 2. Those get normalized
  # my the create multiplex to then have the weights 0.5 and 1, respectively.
  w1 <- 1 * one_minus_delta
  w1half <- 0.5 * one_minus_delta

  colnames <- c("0_1", "1_1", "2_1", "3_1", "0_2", "1_2", "2_2", "3_2")
  list_mat <- c(
    c(0, w1half, w1, 0, delta, 0, 0, 0),
    c(w1half, 0, w1half, 0, 0, delta, 0, 0),
    c(w1, w1half, 0, 0, 0, 0, delta, 0),
    c(0, 0, 0, 0, 0, 0, 0, delta),
    c(delta, 0, 0, 0, 0, 0, w1, w1),
    c(0, delta, 0, 0, 0, 0, w1, 0),
    c(0, 0, delta, 0, w1, w1, 0, w1),
    c(0, 0, 0, delta, w1, 0, w1, 0)
  )
  mat <- matrix(list_mat,
                nrow = 8,
                ncol = 8,
                dimnames = list(colnames, colnames))

  # Normalize the matrix by column
  normalized_mat <- mat %*% diag(1 / colSums(mat))
  # non normalized matrix retains the original delta as the intra-layer weight
  expected_nonnormalized_mat <- as(mat, "dgCMatrix")
  # non matrix are normalized along with the [0,1]
  #   normalized edge weights in each layer.
  expected_normalized_mat <- as(normalized_mat, "dgCMatrix")
  colnames(expected_normalized_mat) <- colnames

  return(c(expected_nonnormalized_mat, expected_normalized_mat))
}

run_test_for_diff_graph_data <- function(
  nw_groups,
  delta,
  output_filename,
  verbose) {
  supra_adj_mat <- generate_expected_supraadj(delta)
  expected_nonnormalized_mat <- supra_adj_mat[[1]]
  expected_normalized_mat <- supra_adj_mat[[2]]

  invisible(
    make_homogenous_network(
      nw_groups,
      delta,
      output_filename,
      verbose)
    )
  load(output_filename)

  expect_equal(nw.adjnorm, expected_normalized_mat)   #nolint loaded from file
  expect_equal(nw.adj, expected_nonnormalized_mat)  #nolint loaded from file
}

describe("make_multiplex", {
  it("throws an error when fed an flist tibble with non-existant files", {
    nw_groups <- list_of(nw_tibble_bad_path)
    expect_error(make_multiplex(nw_groups[[1]]))
  })

  it("also makes a multiplex", {
    # expected output as multiplex object:
    # Because of how RandomWalkRestartMH Calculates the edge weights
    # our graph, m1, originally has weights 1 and 2. Those get normalized
    # my the create multiplex to then have the weights 0.5 and 1, respectively.
    graph1 <- make_graph(c("0", "1", "0", "2", "1", "2"), directed = FALSE)
    E(graph1)$weight <- c(0.5, 1.0, 0.5)
    E(graph1)$type <- c("m1", "m1", "m1")
    # Add vertex 3 because RandomWalkRestartMH
    # ensures all layers share common vertices
    graph1 <- add_vertices(graph1, 1, name = c("3"))
    graph2 <- make_graph(
        c("0", "2", "0", "3", "1", "2", "2", "3"),
        directed = FALSE
    )
    E(graph2)$weight <- c(1, 1, 1, 1)
    E(graph2)$type <- c("m2", "m2", "m2", "m2")
    expected_num_layers <- 2
    expected_num_nodes <- 4

    invisible(
      capture.output(
        actual_output <- make_multiplex(nw_groups[[1]])
      )
    )

    expect_equal(actual_output$Number_of_Layers, expected_num_layers)
    expect_equal(actual_output$Number_of_Nodes_Multiplex, expected_num_nodes)
    # Issues with sorted nature of vertices and edges.
    # Check node equality for each layer:
    expect_setequal(V(actual_output$m1)$name, V(graph1)$name)
    expect_setequal(V(actual_output$m2)$name, V(graph2)$name)
    # Check edge equality for each layer:
    expect_setequal(E(actual_output$m1), E(graph1))
    expect_setequal(E(actual_output$m2), E(graph2))
    # Check edge weight equality for each layer:

    expect_equal(E(actual_output$m1)$weight, E(graph1)$weight)
    expect_equal(E(actual_output$m2)$weight, E(graph2)$weight)
    # Check Layer descriptions for each layer:
    expect_equal(E(actual_output$m1)$type, E(graph1)$type)
    expect_equal(E(actual_output$m2)$type, E(graph2)$type)
  })
})

describe("make_homogenous_network", {
  it("takes in an nw_groups object and saves multiplex network data to file", {
    # delta is the non-normalized edge weight of the connecting edges
    delta <- 0.5
    output_filename <- "testthatOutput.txt"
    verbose <- FALSE

    invisible(
      make_homogenous_network(
        nw_groups,
        delta,
        output_filename,
        verbose
      )
    )

    expect_true(output_filename %in% list.files())
  })


  it("ensures multiplex save file matches expected data", {
    delta <- 0.5
    output_filename <- "testthatOutputP5.txt"
    verbose <- FALSE

    run_test_for_diff_graph_data(nw_groups, delta, output_filename, verbose)
  })

  it("ensures multiplex save file matches expected data with 0.7 as delta", {
    delta <- 0.7
    output_filename <- "testthatOutputP7.txt"
    verbose <- FALSE

    run_test_for_diff_graph_data(nw_groups, delta, output_filename, verbose)
  })

  it("ensures multiplex save file matches expected data with 1.0 as delta", {
    output_filename <- "testthatOutput1P0.txt"
    ## By calling 1 as our delta, the network layers m1 and m2 are entirely
    ## ignored, and only our delta intra-layer edges exist.
    delta <- 1
    colnames <- c("0_1", "1_1", "2_1", "3_1", "0_2", "1_2", "2_2", "3_2")
    list_mat <- c(
      c(0, 0, 0, 0, delta, 0, 0, 0),
      c(0, 0, 0, 0, 0, delta, 0, 0),
      c(0, 0, 0, 0, 0, 0, delta, 0),
      c(0, 0, 0, 0, 0, 0, 0, delta),
      c(delta, 0, 0, 0, 0, 0, 0, 0),
      c(0, delta, 0, 0, 0, 0, 0, 0),
      c(0, 0, delta, 0, 0, 0, 0, 0),
      c(0, 0, 0, delta, 0, 0, 0, 0)
    )
    mat <- matrix(
      list_mat,
      nrow = 8,
      ncol = 8,
      dimnames = list(colnames, colnames)
    )
    # Normalize the matrix by column
    normalized_mat <- mat %*% diag(1 / colSums(mat))
    # non normalized matrix retains the original delta as the intra-layer weight
    expected_nonnormalized_mat <- as(mat, "dgCMatrix")
    # non matrix are normalized along with the [0,1]
    #    normalized edge weights in each layer.
    expected_normalized_mat <- as(normalized_mat, "dgCMatrix")
    colnames(expected_normalized_mat) <- colnames

    invisible(
      make_homogenous_network(
        nw_groups,
        delta,
        output_filename,
        verbose
      )
    )
    load(output_filename)

    expect_equal(nw.adj, expected_nonnormalized_mat)
    expect_equal(nw.adjnorm, expected_normalized_mat)
  })

  it("creates multiplexes with tab delimited data", {
    nw_tibble <- tibble::tibble(
      "nwfile" = c(
          "../testNetworks/m1_tabDelim.txt",
          "../testNetworks/m2_tabDelim.txt"
        ),
      "nwname" = c("m1", "m2"),
      "nwgroup" = c(1, 1)
    )

    nw_groups_tab_delim <- list_of(nw_tibble)
    delta <- 0.5
    output_filename <- "testthatOutputP5_fromTabDelmited.txt"
    verbose <- FALSE

    run_test_for_diff_graph_data(
      nw_groups_tab_delim,
      delta,
      output_filename,
      verbose)
  })

  it("creates multiplexes with non-header data", {
    nw_tibble <- tibble::tibble(
      "nwfile" = c(
        "../testNetworks/m1_noHeader.txt",
        "../testNetworks/m2_noHeader.txt"
      ),
      "nwname" = c("m1", "m2"),
      "nwgroup" = c(1, 1)
    )

    nw_groups_nohead <- list_of(nw_tibble)
    delta <- 0.5
    output_filename <- "testthatOutputP5_fromNoHeader.txt"
    verbose <- FALSE

    run_test_for_diff_graph_data(
      nw_groups_nohead,
      delta,
      output_filename,
      verbose
    )
  })

  it("creates multiplexes with mixed delimited data", {
    nw_tibble <- tibble::tibble(
      "nwfile" = c("../testNetworks/m1.txt", "../testNetworks/m2_tabDelim.txt"),
      "nwname" = c("m1", "m2"),
      "nwgroup" = c(1, 1)
    )

    nw_groups_mixed_delim <- list_of(nw_tibble)
    delta <- 0.5
    output_filename <- "testthatOutputP5_fromMixedDelim.txt"
    verbose <- FALSE

    run_test_for_diff_graph_data(
      nw_groups_mixed_delim,
      delta,
      output_filename,
      verbose
    )
  })
  it("creates multiplexes with mixed header data", {
    nw_tibble <- tibble::tibble(
      "nwfile" = c(
        "../testNetworks/m1_noHeader.txt",
        "../testNetworks/m2.txt"
      ),
      "nwname" = c("m1", "m2"),
      "nwgroup" = c(1, 1)
    )

    nw_groups_mixed_header_ <- list_of(nw_tibble)
    delta <- 0.5
    output_filename <- "testthatOutputP5_fromMixedHeader.txt"
    verbose <- FALSE

    run_test_for_diff_graph_data(
      nw_groups_mixed_header_,
      delta,
      output_filename,
      verbose
    )
  })
  it("creates multiplexes with mixed header and mixed delimiter data", {
    nw_tibble <- tibble::tibble(
      "nwfile" = c(
        "../testNetworks/m1_noHeader.txt",
        "../testNetworks/m2_tabDelim.txt"
      ),
      "nwname" = c("m1", "m2"),
      "nwgroup" = c(1, 1)
    )

    nw_groups_mixed_header_delim <- list_of(nw_tibble)
    delta <- 0.5
    output_filename <- "testthatOutputP5_fromMixedHeaderDelim.txt"
    verbose <- FALSE

    run_test_for_diff_graph_data(
      nw_groups_mixed_header_delim,
      delta,
      output_filename,
      verbose
    )
  })

  it("throws an error when it supplied with a delta value of 0.0", {
    # make_homogenous_network <- function(nw_groups, delta, out, verbose) {
    delta <- 0.0
    one_minus_delta <- 1 - delta
    output_filename <- "testthatOutputP7.txt"
    verbose <- FALSE

    expect_error(
      make_homogenous_network(
        nw_groups,
        delta,
        output_filename,
        verbose)
      )
  })
})

describe("make_heterogeneous_multiplex", {
  it("creates a heterogeneous multiplex and saves it to file", {
    ############ QUESTIONS FOR DAVID WITH MAKE MULTIPLEX HET

    #         0_1  1_1       2_1 3_1       4_1 5_1 6_1 7_1 8_1
    # 0_1 .         0.25 0.3333333 0.5 .         0.5 .     .   .
    # 1_1 0.1666667 .    0.1666667 .   0.5000000 .   .     .   .
    # 2_1 0.3333333 0.25 .         .   .         .   0.5   .   .
    # 3_1 0.2500000 .    .         .   0.1666667 .   .     .   .
    # 4_1 .         0.50 .         0.5 .         0.5 .     .   1
    # 5_1 0.2500000 .    .         .   0.1666667 .   .     .   .
    # 6_1 .         .    0.5000000 .   .         .   .     1   .
    # 7_1 .         .    .         .   .         .   0.5   .   .
    # 8_1 .         .    .         .   0.1666667 .   .     .   .


    ## Inter layer matrix edge weidhts depend on the labda.
    ## When given a 0.5, if there existss only one interlayer edge
    ## That edge recieves the weight of lambda. The remainder of the
    ## weights that are not lambda are then calculated 1-lambda
    nw_group_input <- list_of(nw_tibble_het1, nw_tibble_het2, nw_tibble_het3)
    delta <- 1

    ## currently, with package formation, delta appears to do nothing?

    lambda <- 0.6
    out <- "network.Rdata"

    invisible(
      make_heterogeneous_multiplex(
        nw_group_input,
        delta,
        lambda,
        out)
      )
  })
})

describe("Read Flist", {
  file_list <- c("../testNetworks/abc_layer1.tsv", "../testNetworks/abc_layer2.tsv", "../testNetworks/abc_layer3.tsv")
  layer_list <- c("layer1", "layer2", "layer3")
  group_list <- c(1, 1, 1)
  
  it("throws an error if flist is empty", {
    flist <- ""

    expect_error(read_flist(flist))
  })

  it("reads an flist with only 2 columns, adding in a group", {
    flist <- "../testFlists/test_flist_2cols.tsv"

    actual_output <- read_flist(flist)

    expected_output <- data.table::data.table(
      nwfile = file_list,
      nwname = layer_list,
      nwgroup = group_list
    )

    expect_equal(actual_output, expected_output)
  })


  it("reads an flist with 5 columns and pares down last two columns", {
    flist <- "../testFlists/test_flist_5cols.tsv"

    actual_output <- read_flist(flist)

    expected_output <- data.table::data.table(
      nwfile = file_list,
      nwname = layer_list,
      nwgroup = group_list
    )

    expect_equal(actual_output, expected_output)
  })

})

describe("RWR_make_multiplex.R:", {

  it("throws an error if flist elements contain bad paths", {
    bad_flist <- "../testFlists/testFlist_badPaths.txt"

    expect_error(RWRtoolkit::RWR_make_multiplex(bad_flist))
  })

  it("takes flist and makes a homogenous multiplex with default parameters", {
    nw_group_input <- list_of(nw_tibble)
    flist_file_path <- "../testFlists/testFlist.txt"

    make_homogenous_stub <- mock()
    make_heterogenous_stub <- mock()
    stub(
      RWRtoolkit::RWR_make_multiplex,
      "make_homogenous_network",
      make_homogenous_stub
    )
    stub(
      RWRtoolkit::RWR_make_multiplex,
      "make_heterogeneous_multiplex",
      make_heterogenous_stub
    )

    invisible(RWRtoolkit::RWR_make_multiplex(flist_file_path))

    expect_called(make_homogenous_stub, 1)
    expect_args(
      make_homogenous_stub,
      1,
      nw_group_input,
      0.5,
      "network.Rdata",
      FALSE
    )
    expect_called(make_heterogenous_stub, 0)
  })

  it("takes flist and makes a homogenous multiplex with parameters", {
    nw_group_input <- list_of(nw_tibble)
    flist_file_path <- "../testFlists/testFlist.txt"

    delta <- 0.9
    verbose <- TRUE
    output_file <- "myOutput.txt"

    make_homogenous_stub <- mock()
    make_heterogenous_stub <- mock()
    stub(
      RWRtoolkit::RWR_make_multiplex,
      "make_homogenous_network",
      make_homogenous_stub
    )
    stub(
      RWRtoolkit::RWR_make_multiplex,
      "make_heterogeneous_multiplex",
      make_heterogenous_stub
    )

    invisible(
      RWRtoolkit::RWR_make_multiplex(
        flist_file_path,
        delta = delta,
        output = output_file,
        verbose = verbose
      )
    )

    expect_called(make_homogenous_stub, 1)
    expect_args(
      make_homogenous_stub,
      1,
      nw_group_input,
      delta,
      output_file,
      verbose)
    expect_called(make_heterogenous_stub, 0)
  })

  it("fails to create a heterogeneous network due to there not being enough files", { #nolint 
    flist_file_path <- "../testFlists/testFlist_heterogeneous_badGrouping.txt"

    make_homogenous_stub <- mock()
    make_heterogenous_stub <- mock()
    stub(
      RWRtoolkit::RWR_make_multiplex,
      "make_homogenous_network",
      make_homogenous_stub
    )
    stub(
      RWRtoolkit::RWR_make_multiplex,
      "make_heterogeneous_multiplex",
      make_heterogenous_stub
    )

    expect_error(RWRtoolkit::RWR_make_multiplex(flist_file_path))
    expect_called(make_homogenous_stub, 0)
    expect_called(make_heterogenous_stub, 0)
  })

  it("takes flist and makes a heterogeneous multiplex with default parameters", { #notlint
    ## Heterogeneous networks
    nw_group_input <- list_of(nw_tibble_het1, nw_tibble_het2, nw_tibble_het3)
    flist_file_path <- "../testFlists/testFlist_heterogeneous.txt"

    make_homogenous_stub <- mock()
    make_heterogenous_stub <- mock()
    stub(
      RWRtoolkit::RWR_make_multiplex,
      "make_homogenous_network",
      make_homogenous_stub
    )
    stub(
      RWRtoolkit::RWR_make_multiplex,
      "make_heterogeneous_multiplex",
      make_heterogenous_stub
    )

    expect_warning(RWRtoolkit::RWR_make_multiplex(flist_file_path))

    expect_called(make_homogenous_stub, 0)
    expect_called(make_heterogenous_stub, 1)
    expect_args(
      make_heterogenous_stub,
      1,
      nw_group_input,
      0.5,
      0.5,
      "network.Rdata",
      FALSE
    )
  })
})


# removes test files so as not to clutter everything up
teardown(
  {
    system("rm network.Rdata")
    system("rm testthatOutput.txt")
    system("rm testthatOutputP5.txt")
    system("rm testthatOutputP7.txt")
    system("rm testthatOutput1P0.txt")
    system("rm testthatOutputP5_fromTabDelmited.txt")
    system("rm testthatOutputP5_fromNoHeader.txt")
    system("rm testthatOutputP5_fromMixedDelim.txt")
    system("rm testthatOutputP5_fromMixedHeader.txt")
    system("rm testthatOutputP5_fromMixedHeaderDelim.txt")
  },
  env = parent.frame()
)
