context("Random Walk Restart Tests")
library(RWRtoolkit)
library(RandomWalkRestartMH)
library(vctrs)
library(igraph)
library(mockery)



setid <- c("setA", "setA", "setA")
gene <- c("1a", "2b", "3")
weight <- c(1, 2, 3)
geneset_3genes <- data.frame(setid, gene, weight)

setid <- c("setA", "setA", "setA", "setA")
gene <- c("1a", "2b", "3", "4")
weight <- c(1, 2, 3, 4)
geneset_4genes <- data.frame(setid, gene, weight)

setid <- c("setA", "setA", "setA", "setA", "setA")
gene <- c("1a", "2b", "3", "4", "5")
weight <- c(1, 2, 3, 4, 5)
geneset_5genes <- data.frame(setid, gene, weight)
# 5 genes, assuming 3 folds.
chunks <- list(c("1a", "2b"), c("3", "4"), c("5"))
names(chunks) <- c("1", "2", "3")



generate_mock_rwr_output <- function(NodeNames, Score, expected_seed_nodes) {
  fold_rwr_result <- data.frame(NodeNames, Score)
  fold_list <- list(fold_rwr_result, expected_seed_nodes)
  names(fold_list) <- list("RWRM_Results", "Seed_Nodes")
  fold_list
}

generate_expected_rwr_cv_layer <- function(NodeNames, Score, InValset, num_in_network, fold, leftout, seed, networks, geneset_names, method, num_seeds, num_leftout, modname) {
  # hard coded to work with
  rank <- sequence(length(NodeNames))

  networks <- networks
  geneset <- geneset_names

  data.frame(NodeNames, Score, rank, InValset, num_in_network, num_seeds, num_leftout, networks, fold, modname, geneset, seed, leftout, method)
}

describe("update_folds_by_method", {
  # Basic 2 genes
  setid <- c("setA", "setA")
  gene <- c("1", "2")
  weight <- c(1, 1)
  geneset <- data.frame(setid, gene, weight)

  it("throws an error if methdo is no recognized", {
    method <- "unsupported_method"
    numFolds <- 5

    expect_error(update_folds_by_method(geneset, method, numFolds))
  })

  it("sets folds to the number of rows within the geneset for LOO", {
    method <- "loo"
    expected_folds <- 2
    numFolds <- NA
    expected_chunks <- NULL
    expected_method <- "loo"

    output <- update_folds_by_method(geneset, method, numFolds)

    expect_equal(output[[1]], expected_folds)
    expect_equal(output[[2]], geneset)
    expect_equal(output[[3]], expected_chunks)
    expect_equal(output[[4]], expected_method)
  })

  it("sets folds to the number of rows within the geneset for LOO", {
    setid <- c("setA", "setA", "setA")
    gene <- c("2", "3", "4")
    weight <- c(0.4, 0.6, 1)
    geneset <- data.frame(setid, gene, weight)
    method <- "loo"
    numFolds <- NA
    expected_chunks <- NULL
    expected_folds <- 3
    expected_method <- "loo"

    output <- update_folds_by_method(geneset, method, numFolds)

    expect_equal(output[[1]], expected_folds)
    expect_equal(output[[2]], geneset)
    expect_equal(output[[3]], expected_chunks)
    expect_equal(output[[4]], expected_method)
  })

  it("sets folds to the number of rows within the geneset for singletons", {
    method <- "singletons"
    numFolds <- NA
    expected_folds <- 2
    expected_chunks <- NULL
    expected_method <- "singletons"

    output <- update_folds_by_method(geneset, method, NA)

    expect_equal(output[[1]], expected_folds)
    expect_equal(output[[2]], geneset)
    expect_equal(output[[3]], expected_chunks)
    expect_equal(output[[4]], expected_method)
  })

  it("warns user if kfold and there exist more folds than genes", {
    method <- "kfold"
    numFolds <- 5
    gene <- c("2", "3")
    weight <- c(0.4, 0.6)
    geneset <- data.frame(setid, gene, weight)
    expected_folds <- 2
    expected_chunks <- NULL
    expected_method <- "loo"

    output <- expect_warning(update_folds_by_method(geneset, method, numFolds))

    expect_equal(output[[1]], expected_folds)
    expect_equal(output[[2]], geneset)
    expect_equal(output[[3]], expected_chunks)
    expect_equal(output[[4]], expected_method)
  })

  it("randomly shuffles geneset and chunks with method of kfold ", {
    setid <- c("setA", "setA", "setA")
    gene <- c("1a", "2b", "3")
    weight <- c(1, 2, 3)
    geneset <- data.frame(setid, gene, weight)
    method <- "kfold"
    numFolds <- 2
    expected_chunks <- list(c("1a", "2b"), c("3"))
    names(expected_chunks) <- c("1", "2")
    expected_folds <- 2
    expected_method <- "kfold"

    mock_sample <- mock(sequence(3), gene)
    stub(update_folds_by_method, "sample", mock_sample)
    output <- update_folds_by_method(geneset, method, numFolds)

    expect_equal(output[[1]], expected_folds)
    expect_equal(output[[2]], geneset)
    expect_equal(output[[3]], expected_chunks)
    expect_equal(output[[4]], expected_method)
  })
})


describe("extract_lo_and_seed_genes_cv", {
  # Unexpected method is checked earlier in pipeline, redundant to check here.


  it("extracts one seed genes and leaves out the remainder for singletons where fold is 1", {
    method <- "singletons"
    fold <- 1

    expected_seed_genes <- c("1a")
    expected_leftout <- c("2b", "3")

    res <- extract_lo_and_seed_genes_cv(geneset_3genes, method, fold)


    expect_equal(res[[1]], expected_leftout)
    expect_equal(res[[2]], expected_seed_genes)
  })

  it("extracts one seed genes and leaves out the remainder for singletons where fold is 3", {
    method <- "singletons"
    fold <- 3

    expected_seed_genes <- c("3")
    expected_leftout <- c("1a", "2b")

    res <- extract_lo_and_seed_genes_cv(geneset_3genes, method, fold)

    expect_equal(res[[1]], expected_leftout)
    expect_equal(res[[2]], expected_seed_genes)
  })


  it("extracts one seed genes and leaves out the remainder for singletons of w/ geneset of 4 genes", {
    method <- "singletons"
    fold <- 1

    expected_seed_genes <- c("1a")
    expected_leftout <- c("2b", "3", "4")

    res <- extract_lo_and_seed_genes_cv(geneset_4genes, method, fold)

    expect_equal(res[[1]], expected_leftout)
    expect_equal(res[[2]], expected_seed_genes)
  })

  it("extracts seed genes and leaves out only a gene for loo where fold is 1", {
    method <- "loo"
    fold <- 1

    expected_leftout <- c("1a")
    expected_seed_genes <- c("2b", "3", "4")

    res <- extract_lo_and_seed_genes_cv(geneset_4genes, method, fold)

    expect_equal(res[[1]], expected_leftout)
    expect_equal(res[[2]], expected_seed_genes)
  })

  it("extracts seed genes and leaves out only a gene for loo where fold is 3", {
    method <- "loo"
    fold <- 3

    expected_leftout <- c("3")
    expected_seed_genes <- c("1a", "2b", "4")

    res <- extract_lo_and_seed_genes_cv(geneset_4genes, method, fold)

    expect_equal(res[[1]], expected_leftout)
    expect_equal(res[[2]], expected_seed_genes)
  })

  it("extracts seed gene and left out chunks for kfold: fold is 1", {
    method <- "kfold"
    fold <- 1

    expected_leftout <- c("1a", "2b")
    expected_seed_genes <- c("3", "4", "5")

    res <- extract_lo_and_seed_genes_cv(geneset_5genes, method, fold, chunks)

    expect_equal(res[[1]], expected_leftout)
    expect_equal(res[[2]], expected_seed_genes)
  })

  it("extracts seed gene and left out chunks for kfold: fold is 3", {
    method <- "kfold"
    fold <- 3
    expected_leftout <- c("5")
    expected_seed_genes <- c("1a", "2b", "3", "4")

    res <- extract_lo_and_seed_genes_cv(geneset_5genes, method, fold, chunks)

    expect_equal(res[[1]], expected_leftout)
    expect_equal(res[[2]], expected_seed_genes)
  })
})



describe("create_rankings_cv", {
  # TODO:  Test invalid genes getting worst rank.
  it("creates a ranking object for a specific geneset using singletons", {
    ## Generate Mocked data and expected Values
    method <- "singletons"
    test_geneset <- geneset_3genes
    number_of_rows_in_geneset <- nrow(test_geneset)

    # The nodes in the split fold
    NodeNames <- c("0", "1a")
    number_of_nodes_in_rwr_fold <- length(NodeNames)


    networks <- rep("m1_m2", number_of_nodes_in_rwr_fold)
    fold <- 1
    name <- "default"
    seed_genes <- c("1a")
    num_seeds <- rep(1, number_of_nodes_in_rwr_fold)
    expected_seed <- rep(seed_genes, number_of_nodes_in_rwr_fold)
    leftout <- c("2b", "3")
    num_leftout <- rep(2, number_of_nodes_in_rwr_fold)
    expected_leftout <- c("AllBut1", "AllBut1")
    num_in_network <- rep(number_of_rows_in_geneset, number_of_nodes_in_rwr_fold)
    total_nodes_in_faux_multiplex <- 5

    Score <- c(0.03, 0.01) # Some arbitrary score
    fold1_list <- generate_mock_rwr_output(NodeNames, Score, seed_genes)

    # Neither node 0 or 1a are in leftout, returns 0
    InValset <- c(0, 0)
    fold <- rep(1, number_of_nodes_in_rwr_fold)
    expected_fold1 <- generate_expected_rwr_cv_layer(NodeNames, Score, InValset, num_in_network, fold, expected_leftout, expected_seed, networks, c("setA", "setA"), rep(method, 2), num_seeds, num_leftout, name)

    res <- create_rankings_cv(fold1_list, networks, fold, name, test_geneset, method, seed_genes, leftout, total_nodes_in_faux_multiplex)
    expect_equal(res, expected_fold1)
  })

  it("creates a ranking object for a specific geneset using singletons", {
    ## Generate Mocked data and expected Values
    method <- "singletons"
    test_geneset <- geneset_3genes
    number_of_rows_in_geneset <- nrow(test_geneset)

    # The nodes in the split fold
    NodeNames <- c("0", "1a")
    number_of_nodes_in_rwr_fold <- length(NodeNames)


    networks <- rep("m1_m2", number_of_nodes_in_rwr_fold)
    fold <- 1
    name <- "default"
    seed_genes <- c("1a")
    num_seeds <- rep(1, number_of_nodes_in_rwr_fold)
    expected_seed <- rep(seed_genes, number_of_nodes_in_rwr_fold)
    leftout <- c("2b", "3")
    num_leftout <- rep(2, number_of_nodes_in_rwr_fold)
    expected_leftout <- c("AllBut1", "AllBut1")
    num_in_network <- rep(number_of_rows_in_geneset, number_of_nodes_in_rwr_fold)
    total_nodes_in_faux_multiplex <- 5

    Score <- c(0.03, 0.01) # Some arbitrary score
    fold1_list <- generate_mock_rwr_output(NodeNames, Score, seed_genes)

    # Neither node 0 or 1a are in leftout, returns 0
    InValset <- c(0, 0)
    fold <- rep(1, number_of_nodes_in_rwr_fold)
    expected_fold1 <- generate_expected_rwr_cv_layer(NodeNames, Score, InValset, num_in_network, fold, expected_leftout, expected_seed, networks, c("setA", "setA"), rep(method, 2), num_seeds, num_leftout, name)

    res <- create_rankings_cv(fold1_list, networks, fold, name, test_geneset, method, seed_genes, leftout, total_nodes_in_faux_multiplex)
    expect_equal(res, expected_fold1)
  })

  it("creates a ranking object for a specific geneset using loo", {
    method <- "loo"
    test_geneset <- geneset_3genes
    number_of_rows_in_geneset <- nrow(test_geneset)

    NodeNames <- c("0", "1a")
    number_of_nodes_in_rwr_fold <- length(NodeNames)

    networks <- rep("m1_m2", number_of_nodes_in_rwr_fold)
    fold <- 1
    name <- "default"
    seed_genes <- c("2b", "3")
    num_seeds <- rep(2, number_of_nodes_in_rwr_fold)
    expected_seed <- rep("AllBut1", number_of_nodes_in_rwr_fold)
    leftout <- c("1a")
    num_leftout <- rep(1, number_of_nodes_in_rwr_fold)
    expected_leftout <- rep(leftout, number_of_nodes_in_rwr_fold)
    num_in_network <- rep(number_of_rows_in_geneset, number_of_nodes_in_rwr_fold)
    total_nodes_in_faux_multiplex <- 5

    Score <- c(0.03, 0.01) # Some arbitrary score
    fold1_list <- generate_mock_rwr_output(NodeNames, Score, seed_genes)

    # Node 0 isn't in lefout, but  1a is, returns 0,1
    InValset <- c(0, 1)
    fold <- rep(1, number_of_nodes_in_rwr_fold)
    expected_fold1 <- generate_expected_rwr_cv_layer(NodeNames, Score, InValset, num_in_network, fold, expected_leftout, expected_seed, networks, c("setA", "setA"), rep(method, 2), num_seeds, num_leftout, name)

    res <- create_rankings_cv(fold1_list, networks, fold, name, test_geneset, method, seed_genes, leftout, total_nodes_in_faux_multiplex)
    expect_equal(res, expected_fold1)
  })

  it("creates a ranking object for a specific geneset using kmeans", {
    method <- "kmeans"
    test_geneset <- geneset_5genes
    number_of_rows_in_geneset <- nrow(test_geneset)

    # RWR nodes tested in fold
    NodeNames <- c("0", "1a", "2b")
    number_of_nodes_in_rwr_fold <- length(NodeNames)

    networks <- rep("m1_m2", number_of_nodes_in_rwr_fold)
    fold <- 1
    name <- "default"
    method <- "kfold"

    seed_genes <- c("3", "4", "5")
    num_seeds <- rep(3, number_of_nodes_in_rwr_fold)
    expected_seed <- rep("many", number_of_nodes_in_rwr_fold)
    leftout <- c("1a", "2b")

    num_leftout <- rep(2, number_of_nodes_in_rwr_fold)
    expected_leftout <- rep("many", number_of_nodes_in_rwr_fold)
    num_in_network <- rep(number_of_rows_in_geneset, number_of_nodes_in_rwr_fold)
    total_nodes_in_faux_multiplex <- 5

    NodeNames <- c("0", "1a", "2b")
    Score <- c(0.03, 0.01, 0.007) # Some arbitrary score
    fold1_list <- generate_mock_rwr_output(NodeNames, Score, seed_genes)

    # Node 0 isn't in lefout, but  1a is, returns 0,1
    InValset <- c(0, 1, 1)
    fold <- rep(1, number_of_nodes_in_rwr_fold)
    expected_fold1 <- generate_expected_rwr_cv_layer(NodeNames, Score, InValset, num_in_network, fold, expected_leftout, expected_seed, networks, c("setA", "setA", "setA"), rep(method, 3), num_seeds, num_leftout, name)

    res <- create_rankings_cv(fold1_list, networks, fold, name, test_geneset, method, seed_genes, leftout, total_nodes_in_faux_multiplex)


    expect_equal(res, expected_fold1)
  })
})

describe("RWR", {
  setid <- c("setA", "setA", "setA")
  gene <- c("1", "2", "3")
  weight <- c(1, 2, 3)
  geneset_3genes <- data.frame(setid, gene, weight)

  chunks <- list(c("1", "2"), c("3"))
  names(chunks) <- c("1", "2")
  load("../testMultiplex/unitTestMultiplex.Rdata")

  it("runs rwr on all folds and returns list of rankings", {
    method <- "singletons"
    numfolds <- 3 ## immiterial as singletons uses nrows of geneset
    tau <- c(1, 1)
    expected_chunks <- NULL
    expected_folds <- 3
    first_fold_seeds <- c("1")
    first_fold_nodes <- c("0", "2", "3")
    first_fold_scores <- c(0.7, 0.25, 0.07)
    first_inValSet <- c(0, 1, 1)
    second_fold_seeds <- c("2")
    second_fold_nodes <- c("0", "1", "3")
    second_fold_scores <- c(0.6, 0.34, 0.07)
    second_inValSet <- c(0, 1, 1)
    third_fold_seeds <- c("3")
    third_fold_nodes <- c("0", "1", "2")
    third_fold_scores <- c(0.8, 0.14, 0.07)
    third_inValSet <- c(0, 1, 1)
    repititions <- 3
    name <- "default"

    # Create mocks for function
    first_mockExtractReturn <- list(first_fold_nodes, first_fold_seeds)
    second_mockExtractReturn <- list(second_fold_nodes, second_fold_seeds)
    third_mockExtractReturn <- list(third_fold_nodes, third_fold_seeds)
    num_in_network <- rep(4, repititions)
    num_leftout <- rep(1, repititions)
    num_seeds <- rep(1, repititions)
    networks <- rep("m1_m2", repititions)
    leftout <- rep("AllBut1", repititions)

    first_rwr_output <- generate_mock_rwr_output(first_fold_nodes, first_fold_scores, first_fold_seeds)
    second_rwr_output <- generate_mock_rwr_output(second_fold_nodes, second_fold_scores, second_fold_seeds)
    third_rwr_output <- generate_mock_rwr_output(third_fold_nodes, third_fold_scores, third_fold_seeds)

    first_mockLayer <- generate_expected_rwr_cv_layer(first_fold_nodes, first_fold_scores, first_inValSet, num_in_network, rep(1, repititions), leftout, first_fold_seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
    second_mockLayer <- generate_expected_rwr_cv_layer(second_fold_nodes, second_fold_scores, second_inValSet, num_in_network, rep(2, repititions), leftout, second_fold_seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
    third_mockLayer <- generate_expected_rwr_cv_layer(third_fold_nodes, third_fold_scores, third_inValSet, num_in_network, rep(3, repititions), leftout, third_fold_seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)

    mock_update_folds_by_method <- mock(list(expected_folds, geneset_3genes, expected_chunks, method))
    mock_extract_leftout_and_seed_genes_cv <- mock(
      first_mockExtractReturn,
      second_mockExtractReturn,
      third_mockExtractReturn
    )
    mock_RWRMH <- mock(
      first_rwr_output,
      second_rwr_output,
      third_rwr_output
    )
    mock_create_rankings_cv <- mock(
      first_mockLayer,
      second_mockLayer,
      third_mockLayer
    )

    # Stub functions with mocks
    stub(RWR, "update_folds_by_method", mock_update_folds_by_method)
    stub(RWR, "extract_lo_and_seed_genes_cv", mock_extract_leftout_and_seed_genes_cv)
    stub(RWR, "RandomWalkRestartMH::Random.Walk.Restart.Multiplex", mock_RWRMH)
    stub(RWR, "create_rankings_cv", mock_create_rankings_cv)
    expected_response <- list(first_mockLayer, second_mockLayer, third_mockLayer)

    response <- RWR(geneset_3genes, nw.adjnorm, nw.mpo, method, numfolds, tau = tau)

    expect_equal(response, expected_response)
    expect_args(mock_update_folds_by_method, 1, geneset_3genes, method, expected_folds)
    expect_called(mock_update_folds_by_method, 1)
    expect_args(mock_extract_leftout_and_seed_genes_cv, 1, geneset_3genes, method, 1, expected_chunks)
    expect_args(mock_extract_leftout_and_seed_genes_cv, 2, geneset_3genes, method, 2, expected_chunks)
    expect_args(mock_extract_leftout_and_seed_genes_cv, 3, geneset_3genes, method, 3, expected_chunks)
    expect_called(mock_extract_leftout_and_seed_genes_cv, 3)

    expect_args(mock_RWRMH, 1, nw.adjnorm, nw.mpo, first_fold_seeds, 0.7, tau)
    expect_args(mock_RWRMH, 2, nw.adjnorm, nw.mpo, second_fold_seeds, 0.7, tau)
    expect_args(mock_RWRMH, 3, nw.adjnorm, nw.mpo, third_fold_seeds, 0.7, tau)
    expect_called(mock_RWRMH, 3)
    expect_args(mock_create_rankings_cv, 1, first_rwr_output, networks[1], 1, name, geneset_3genes, method, first_fold_seeds, first_fold_nodes, nw.mpo$Number_of_Nodes_Multiplex)
    expect_args(mock_create_rankings_cv, 2, second_rwr_output, networks[1], 2, name, geneset_3genes, method, second_fold_seeds, second_fold_nodes, nw.mpo$Number_of_Nodes_Multiplex)
    expect_args(mock_create_rankings_cv, 3, third_rwr_output, networks[1], 3, name, geneset_3genes, method, third_fold_seeds, third_fold_nodes, nw.mpo$Number_of_Nodes_Multiplex)
    expect_called(mock_create_rankings_cv, 3)
  })
})

describe("Post Processing", {
  ## Base data for rwr_response and ultimate
  method <- "singletons"
  numfolds <- 3 ## immiterial as singletons uses nrows of geneset
  first_fold_seeds <- c("1")
  first_fold_nodes <- c("0", "2", "3")
  first_fold_scores <- c(0.7, 0.25, 0.071)
  first_inValSet <- c(0, 1, 1)
  second_fold_seeds <- c("2")
  second_fold_nodes <- c("0", "1", "3")
  second_fold_scores <- c(0.6, 0.34, 0.07)
  second_inValSet <- c(0, 1, 1)
  third_fold_seeds <- c("3")
  third_fold_nodes <- c("0", "1", "2")
  third_fold_scores <- c(0.8, 0.14, 0.075)
  third_inValSet <- c(0, 1, 1)
  repititions <- 3
  num_in_network <- rep(4, repititions)
  num_leftout <- rep(2, repititions)
  num_seeds <- rep(1, repititions)
  networks <- rep("m1_m2", repititions)
  leftout <- rep("AllBut1", repititions)
  name <- "default"
  first_mockLayer <- generate_expected_rwr_cv_layer(first_fold_nodes, first_fold_scores, first_inValSet, num_in_network, rep(1, repititions), leftout, first_fold_seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
  second_mockLayer <- generate_expected_rwr_cv_layer(second_fold_nodes, second_fold_scores, second_inValSet, num_in_network, rep(2, repititions), leftout, second_fold_seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
  third_mockLayer <- generate_expected_rwr_cv_layer(third_fold_nodes, third_fold_scores, third_inValSet, num_in_network, rep(3, repititions), leftout, third_fold_seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)

  rwr_res <- list(first_mockLayer, second_mockLayer, third_mockLayer)

  describe("post_process_rwr_output_cv", {
    it("binds rows and arranges by rank with no extras", {
      extras <- NULL

      # because all faux layers are ranked in order of 1, 2, 3, then the corresponding expected arrangement is the first from each, then the second, etc.
      expected_nodenames <- c(
        first_fold_nodes[1], second_fold_nodes[1], third_fold_nodes[1],
        first_fold_nodes[2], second_fold_nodes[2], third_fold_nodes[2],
        first_fold_nodes[3], second_fold_nodes[3], third_fold_nodes[3]
      )
      expected_scores <- c(
        first_fold_scores[1], second_fold_scores[1], third_fold_scores[1],
        first_fold_scores[2], second_fold_scores[2], third_fold_scores[2],
        first_fold_scores[3], second_fold_scores[3], third_fold_scores[3]
      )
      expected_ranks <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
      expected_invalset <- c(
        first_inValSet[1], second_inValSet[1], third_inValSet[1],
        first_inValSet[2], second_inValSet[2], third_inValSet[2],
        first_inValSet[3], second_inValSet[3], third_inValSet[3]
      )
      expected_numinnetwork <- rep(4, 9)
      expected_numseeds <- rep(1, 9)
      expected_numleftout <- rep(2, 9)
      expected_networks <- rep("m1_m2", 9)
      expected_fold <- rep(c(1, 2, 3), 3)
      expected_geneset <- rep("setA", 9)
      expected_seed <- rep(c("1", "2", "3"), 3)
      expected_lefout <- rep("AllBut1", 9)
      expected_method <- rep(method, 9)
      expected_name <- rep(name, 9)
      expected_output <- data.frame(
        "NodeNames" = expected_nodenames, "Score" = expected_scores,
        "rank" = expected_ranks, "InValset" = expected_invalset, "num_in_network" = expected_numinnetwork,
        "num_seeds" = expected_numseeds, "num_leftout" = expected_numleftout,
        "networks" = expected_networks, "fold" = expected_fold, "modname" = expected_name, "geneset" = expected_geneset,
        "seed" = expected_seed, "leftout" = expected_lefout, "method" = expected_method
      )

      output <- post_process_rwr_output_cv(rwr_res, extras, numfolds, nw.mpo)

      expect_equal(output, expected_output)
    })

    it("binds rows and arranges by rank with extras as a singleton", {
      setid <- c("setA", "setA")
      gene <- c("5", "6")
      weight <- c(0.6, 0.9)
      extras <- data.frame(setid, gene, weight)

      # because all faux layers are ranked in order of 1, 2, 3, then the corresponding expected arrangement is the first from each, then the second, etc.
      expected_nodenames <- c(
        first_fold_nodes[1], second_fold_nodes[1], third_fold_nodes[1],
        first_fold_nodes[2], second_fold_nodes[2], third_fold_nodes[2],
        first_fold_nodes[3], second_fold_nodes[3], third_fold_nodes[3],
        gene[1], gene[2]
      )
      expected_scores <- c(
        first_fold_scores[1], second_fold_scores[1], third_fold_scores[1],
        first_fold_scores[2], second_fold_scores[2], third_fold_scores[2],
        first_fold_scores[3], second_fold_scores[3], third_fold_scores[3],
        0, 0
      )
      expected_ranks <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4)
      expected_invalset <- c(
        first_inValSet[1], second_inValSet[1], third_inValSet[1],
        first_inValSet[2], second_inValSet[2], third_inValSet[2],
        first_inValSet[3], second_inValSet[3], third_inValSet[3],
        1, 1
      )
      expected_numinnetwork <- rep(4, 11)
      expected_numseeds <- rep(1, 11)
      expected_numleftout <- rep(2, 11)
      expected_networks <- rep("m1_m2", 11)
      expected_fold <- c(rep(c(1, 2, 3), 3), c(1, 2)) ## as the two extras are appended as one of each fold
      expected_geneset <- rep("setA", 11)
      expected_seed <- c(rep(c("1", "2", "3"), 3), c("missing", "missing"))
      expected_lefout <- c(rep("AllBut1", 9), c("missing", "missing"))
      modname <- rep("default", 11)
      expected_method <- rep(method, 11)
      expected_output <- data.frame(
        "NodeNames" = expected_nodenames, "Score" = expected_scores,
        "rank" = expected_ranks, "InValset" = expected_invalset, "num_in_network" = expected_numinnetwork,
        "num_seeds" = expected_numseeds, "num_leftout" = expected_numleftout,
        "networks" = expected_networks, "fold" = expected_fold, "modname" = modname, "geneset" = expected_geneset,
        "seed" = expected_seed, "leftout" = expected_lefout, "method" = expected_method
      )


      output <- post_process_rwr_output_cv(rwr_res, extras, numfolds, nw.mpo)
      expect_equal(output, expected_output)
    })


    it("binds rows and arranges by rank with extras from LOO", {
      method <- "loo"
      numfolds <- 3 ## immiterial as singletons uses nrows of geneset
      load("../testMultiplex/unitTestMultiplex.Rdata") # load nw.mpo

      setid <- c("setA", "setA", "setA")
      gene <- c("1", "2", "3")
      weight <- c(1, 2, 3)
      geneset_3genes <- data.frame(setid, gene, weight)

      repititions <- 3
      seeds <- rep("AllBut1", repititions)
      num_in_network <- rep(3, repititions)
      num_seeds <- rep(2, repititions)
      num_leftout <- rep(1, repititions)

      # layer 1
      first_fold_nodes <- c(0, 1, 2)
      first_fold_scores <- c(0.30, 0.25, 0.11)
      first_fold_invalset <- c(1, 0, 1)
      first_fold_leftout <- c(1)

      # layer 2
      second_fold_nodes <- c(1, 3, 2)
      second_fold_scores <- c(0.22, 0.15, 0.07)
      second_fold_invalset <- c(1, 1, 0)
      second_fold_leftout <- c(2)

      # layer 3
      third_fold_nodes <- c(2, 0, 1)
      third_fold_scores <- c(0.32, 0.18, 0.09)
      third_fold_invalset <- c(1, 1, 1)
      third_fold_leftout <- c(3)

      layer1 <- generate_expected_rwr_cv_layer(first_fold_nodes, first_fold_scores, first_fold_invalset, num_in_network, rep(1, repititions), rep(first_fold_leftout, repititions), seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
      layer2 <- generate_expected_rwr_cv_layer(second_fold_nodes, second_fold_scores, second_fold_invalset, num_in_network, rep(2, repititions), rep(second_fold_leftout, repititions), seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
      layer3 <- generate_expected_rwr_cv_layer(third_fold_nodes, third_fold_scores, third_fold_invalset, num_in_network, rep(3, repititions), rep(third_fold_leftout, repititions), seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
      # RWRtoolkit::RWR(geneset_3genes, nw.adjnorm, nw.mpo, method, numfolds)

      res <- list(layer1, layer2, layer3)
      setid <- c("setA", "setA")
      gene <- c("5", "6")
      weight <- c(0.6, 0.9)
      extras <- data.frame(setid, gene, weight)

      expected_nodenames <- c(
        first_fold_nodes[1], second_fold_nodes[1], third_fold_nodes[1],
        first_fold_nodes[2], second_fold_nodes[2], third_fold_nodes[2],
        first_fold_nodes[3], second_fold_nodes[3], third_fold_nodes[3],
        gene[1], gene[2]
      )
      expected_scores <- c(
        first_fold_scores[1], second_fold_scores[1], third_fold_scores[1],
        first_fold_scores[2], second_fold_scores[2], third_fold_scores[2],
        first_fold_scores[3], second_fold_scores[3], third_fold_scores[3],
        0, 0
      )
      expected_ranks <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4)
      expected_invalset <- c(
        first_fold_invalset[1], second_fold_invalset[1], third_fold_invalset[1],
        first_fold_invalset[2], second_fold_invalset[2], third_fold_invalset[2],
        first_fold_invalset[3], second_fold_invalset[3], third_fold_invalset[3],
        1, 1
      )
      expected_num_in_network <- rep(3, 11)
      expected_num_seeds <- rep(2, 11)
      expected_num_leftout <- rep(1, 11)
      expected_networks <- rep("m1_m2", 11)
      expected_fold <- c(rep(c(1, 2, 3), 3), c(4, 5))
      expected_modname <- rep("default", 11)
      expected_geneset <- rep("setA", 11)
      expected_seed <- c(rep("AllBut1", 9), c("", ""))
      expected_leftout <- c(rep(c("1", "2", "3"), 3), c("5", "6"))
      expected_method <- rep("loo", 11)

      expected_output <- data.frame(
        "NodeNames" = expected_nodenames, "Score" = expected_scores,
        "rank" = expected_ranks, "InValset" = expected_invalset, "num_in_network" = expected_num_in_network,
        "num_seeds" = expected_num_seeds, "num_leftout" = expected_num_leftout,
        "networks" = expected_networks, "fold" = expected_fold, "modname" = expected_modname, "geneset" = expected_geneset,
        "seed" = expected_seed, "leftout" = expected_leftout, "method" = expected_method
      )


      output <- post_process_rwr_output_cv(res, extras, folds, nw.mpo)

      expect_equal(output, expected_output)
    })
  })

  describe("calculate_average_rank_across_folds_cv", {
    it("takes the combined folds from singletons/kfold with no extras and calculates the average of them", {
      nodenames <- c(
        first_fold_nodes[1], second_fold_nodes[1], third_fold_nodes[1],
        first_fold_nodes[2], second_fold_nodes[2], third_fold_nodes[2],
        first_fold_nodes[3], second_fold_nodes[3], third_fold_nodes[3]
      )
      scores <- c(
        first_fold_scores[1], second_fold_scores[1], third_fold_scores[1],
        first_fold_scores[2], second_fold_scores[2], third_fold_scores[2],
        first_fold_scores[3], second_fold_scores[3], third_fold_scores[3]
      )
      ranks <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
      invalset <- c(
        first_inValSet[1], second_inValSet[1], third_inValSet[1],
        first_inValSet[2], second_inValSet[2], third_inValSet[2],
        first_inValSet[3], second_inValSet[3], third_inValSet[3]
      )
      numinnetwork <- rep(4, 9)
      numseeds <- rep(1, 9)
      numleftout <- rep(2, 9)
      networks <- rep("m1_m2", 9)
      fold <- rep(c(1, 2, 3), 3) ## as the two extras are appended as one of each fold
      genesets <- rep("setA", 9)
      seeds <- rep(c("1", "2", "3"), 3)
      leftout <- rep("AllBut1", 9)
      modname <- rep("default", 9)
      method <- rep(method, 9)
      res_combined <- data.frame(
        "NodeNames" = nodenames, "Score" = scores,
        "rank" = ranks, "InValset" = invalset, "num_in_network" = numinnetwork,
        "num_seeds" = numseeds, "num_leftout" = numleftout,
        "networks" = networks, "fold" = fold, "modname" = modname, "geneset" = genesets,
        "seed" = seeds, "leftout" = leftout, "method" = method
      )


      expected_output <- tibble::tibble(
        NodeNames = c("0", "1", "2", "3"),
        meanrank = c(1, 2, 2.5, 3),
        rerank = c(1, 2, 3, 4),
        medrank = c(1, 2, 2.5, 3),
        InValset = c(0, 1, 1, 1),
        geneset = c("setA", "setA", "setA", "setA"),
        num_in_network = c(4, 4, 4, 4)
      )

      output <- calculate_average_rank_across_folds_cv(res_combined)

      expect_equal(output, expected_output)
    })

    it("takes the combined folds from singletons/kfold WITH EXTRAS and calculates the average of them", {
      setid <- c("setA", "setA")
      gene <- c("5", "6")
      weight <- c(0.6, 0.9)
      extras <- data.frame(setid, gene, weight)

      # because all faux layers are ranked in order of 1, 2, 3, then the corresponding expected arrangement is the first from each, then the second, etc.
      nodenames <- c(
        first_fold_nodes[1], second_fold_nodes[1], third_fold_nodes[1],
        first_fold_nodes[2], second_fold_nodes[2], third_fold_nodes[2],
        first_fold_nodes[3], second_fold_nodes[3], third_fold_nodes[3],
        gene[1], gene[2]
      )
      scores <- c(
        first_fold_scores[1], second_fold_scores[1], third_fold_scores[1],
        first_fold_scores[2], second_fold_scores[2], third_fold_scores[2],
        first_fold_scores[3], second_fold_scores[3], third_fold_scores[3],
        0, 0
      )
      ranks <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4)
      invalset <- c(
        first_inValSet[1], second_inValSet[1], third_inValSet[1],
        first_inValSet[2], second_inValSet[2], third_inValSet[2],
        first_inValSet[3], second_inValSet[3], third_inValSet[3],
        1, 1
      )
      numinnetwork <- rep(4, 11)
      numseeds <- rep(1, 11)
      numleftout <- rep(2, 11)
      networks <- rep("m1_m2", 11)
      fold <- c(rep(c(1, 2, 3), 3), c(1, 2)) ## as the two extras are appended as one of each fold
      genesets <- rep("setA", 11)
      seeds <- c(rep(c("1", "2", "3"), 3), c("extra", "extra"))
      leftout <- c(rep("AllBut1", 9), c("extra", "extra"))
      modname <- rep("default", 11)
      method <- rep(method, 11)
      res_combined <- data.frame(
        "NodeNames" = nodenames, "Score" = scores,
        "rank" = ranks, "InValset" = invalset, "num_in_network" = numinnetwork,
        "num_seeds" = numseeds, "num_leftout" = numleftout,
        "networks" = networks, "fold" = fold, "modname" = modname, "geneset" = genesets,
        "seed" = seeds, "leftout" = leftout, "method" = method
      )

      expected_output <- tibble::tibble(
        NodeNames = c("0", "1", "2", "3", "5", "6"),
        meanrank = c(1, 2, 2.5, 3, 4, 4),
        rerank = c(1, 2, 3, 4, 5, 5),
        medrank = c(1, 2, 2.5, 3, 4, 4),
        InValset = c(0, 1, 1, 1, 1, 1),
        geneset = c("setA", "setA", "setA", "setA", "setA", "setA"),
        num_in_network = c(4, 4, 4, 4, 4, 4)
      )

      output <- calculate_average_rank_across_folds_cv(res_combined)
      expect_equal(output, expected_output)
    })

    it("takes the combined folds from LOO and calculates the average of them", {
      method <- "loo"
      numfolds <- 3 ## immiterial as singletons uses nrows of geneset

      setid <- c("setA", "setA", "setA")
      gene <- c("1", "2", "3")
      weight <- c(1, 2, 3)
      geneset_3genes <- data.frame(setid, gene, weight)

      repititions <- 3
      seeds <- rep("AllBut1", repititions)
      num_in_network <- rep(3, repititions)
      num_seeds <- rep(2, repititions)
      num_leftout <- rep(1, repititions)

      # layer 1
      first_fold_nodes <- c(0, 1, 2)
      first_fold_scores <- c(0.30, 0.25, 0.11)
      first_fold_invalset <- c(0, 1, 1)
      first_fold_leftout <- c(1)

      # layer 2
      second_fold_nodes <- c(1, 3, 2)
      second_fold_scores <- c(0.22, 0.15, 0.07)
      second_fold_invalset <- c(1, 1, 1)
      second_fold_leftout <- c(2)

      # layer 3
      third_fold_nodes <- c(2, 0, 1)
      third_fold_scores <- c(0.32, 0.18, 0.09)
      third_fold_invalset <- c(1, 0, 1)
      third_fold_leftout <- c(3)



      layer1 <- generate_expected_rwr_cv_layer(first_fold_nodes, first_fold_scores, first_fold_invalset, num_in_network, rep(1, repititions), rep(first_fold_leftout, repititions), seeds,  networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
      layer2 <- generate_expected_rwr_cv_layer(second_fold_nodes, second_fold_scores, second_fold_invalset, num_in_network, rep(2, repititions), rep(second_fold_leftout, repititions), seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
      layer3 <- generate_expected_rwr_cv_layer(third_fold_nodes, third_fold_scores, third_fold_invalset, num_in_network, rep(3, repititions), rep(third_fold_leftout, repititions), seeds, networks, c("setA", "setA", "setA"), rep(method, repititions), num_seeds, num_leftout, name)
      # RWRtoolkit::RWR(geneset_3genes, nw.adjnorm, nw.mpo, method, numfolds)
      
      res <- list(layer1, layer2, layer3)
      setid <- c("setA", "setA")
      gene <- c("5", "6")
      weight <- c(0.6, 0.9)
      extras <- data.frame(setid, gene, weight)

      nodenames <- c(
        first_fold_nodes[1], second_fold_nodes[1], third_fold_nodes[1],
        first_fold_nodes[2], second_fold_nodes[2], third_fold_nodes[2],
        first_fold_nodes[3], second_fold_nodes[3], third_fold_nodes[3],
        gene[1], gene[2]
      )
      scores <- c(
        first_fold_scores[1], second_fold_scores[1], third_fold_scores[1],
        first_fold_scores[2], second_fold_scores[2], third_fold_scores[2],
        first_fold_scores[3], second_fold_scores[3], third_fold_scores[3],
        0, 0
      )
      ranks <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4)
      invalset <- c(
        first_fold_invalset[1], second_fold_invalset[1], third_fold_invalset[1],
        first_fold_invalset[2], second_fold_invalset[2], third_fold_invalset[2],
        first_fold_invalset[3], second_fold_invalset[3], third_fold_invalset[3],
        1, 1
      )
      num_in_network <- rep(3, 11)
      num_seeds <- rep(2, 11)
      num_leftout <- rep(1, 11)
      networks <- rep("m1_m2", 11)
      fold <- c(rep(c(1, 2, 3), 3), c(4, 5))
      modname <- rep("default", 11)
      geneset <- rep("setA", 11)
      seed <- c(rep("AllBut1", 9), c("", ""))
      leftout <- c(rep(c("1", "2", "3"), 3), c("5", "6"))
      method <- rep("loo", 11)

      res_combined_loo <- data.frame(
        "NodeNames" = nodenames, "Score" = scores,
        "rank" = ranks, "InValset" = invalset, "num_in_network" = num_in_network,
        "num_seeds" = num_seeds, "num_leftout" = num_leftout,
        "networks" = networks, "fold" = fold, "modname" = modname, "geneset" = geneset,
        "seed" = seed, "leftout" = leftout, "method" = method
      )

      expected_output <- tibble::tibble(
        NodeNames = c("0", "1", "3", "2", "5", "6"),
        meanrank = c(1.5, 2, 2, 7 / 3, 4, 4), 
        # Note: On tie of non-extra, count increases on rerank so we get 1, 2, 2, -> 4 (instead of 1,2,2,3)
        rerank = c(1, 2, 2, 4, 5, 5),
        medrank = c(1.5, 2, 2, 3, 4, 4),
        InValset = c(0, 1, 1, 1, 1, 1),
        geneset = c("setA", "setA", "setA", "setA", "setA", "setA"),
        num_in_network = c(3, 3, 3, 3, 3, 3)
      )

      output <- calculate_average_rank_across_folds_cv(res_combined_loo)

  print("EXPECTED EQUAL?")
  print(expected_output)
  print(output)
      expect_equal(output, expected_output)
    })
  })
})
