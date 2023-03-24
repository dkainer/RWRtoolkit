context("Utils Tests")
library(RWRtoolkit)
library(vctrs)
library(mockery)
load("../testMultiplex/unitTestMultiplex.Rdata")

describe("load_network", {
    describe("load_network", {
      expected_from <- c("B", "E", "C", "E", "F", "F", "G", "H", "H", "H")
      expected_to <- c("A", "B", "B", "D", "E", "C", "E", "E", "F", "C")
      expected_elements <- c()
      for (i in seq(1, length(expected_from))) {
        expected_elements <- append(expected_elements, expected_from[i])
        expected_elements <- append(expected_elements, expected_to[i])
      }

      it("loads an igraph object from an edgelist", {
        path_to_edgelist <- "../testNetworks/abc_layer1.tsv"

        actual_network <- load_network(path_to_edgelist)

        expected_network <- igraph::make_graph(expected_elements, directed = F)
        igraph::E(expected_network)$weight <- 1

        expect_setequal(
          igraph::V(actual_network)$name,
          igraph::V(expected_network)$name
        )
        expect_setequal(
          igraph::E(actual_network),
          igraph::E(expected_network)
        )
        expect_setequal(
          igraph::E(actual_network)$weight,
          igraph::E(expected_network)$weight
        )
        expect_true(!is.directed(actual_network))
      })

      it("loads a graph from an edgelist with no weights", {
        path_to_edgelist <- "../testNetworks/abc_layer1_noweight.tsv"

        actual_network <- load_network(path_to_edgelist)

        expected_network <- igraph::make_graph(expected_elements, directed = F)

        expect_setequal(
          V(actual_network)$name,
          V(expected_network)$name
        )
        expect_setequal(
          E(actual_network),
          E(expected_network)
        )
        expect_true(!is.directed(actual_network))
      })

      it("loads a graph from an edgelist with weights and additional cols", {
        path_to_edgelist <- "../testNetworks/abc_layer1_extra_columns.tsv"

        actual_network <- load_network(path_to_edgelist)

        expected_network <- igraph::make_graph(expected_elements, directed = F)
        igraph::E(expected_network)$weight <- 0.5


        expect_setequal(
          V(actual_network)$name,
          V(expected_network)$name
        )
        expect_setequal(
          E(actual_network),
          E(expected_network)
        )
        expect_setequal(
          E(actual_network)$weight,
          E(expected_network)$weight
        )
        expect_true(!is.directed(actual_network))
      })


      it("loads a directed igraph object from an edgelist", {
        path_to_edgelist <- "../testNetworks/abc_layer1.tsv"
        directed <- TRUE
        type <- "testnetwork"
        name <- "TESTLAYER"

        actual_network <- load_network(
          path_to_edgelist = path_to_edgelist,
          type = type,
          name = name,
          directed = directed
        )

        expected_name <- name
        expected_network <- igraph::make_graph(expected_elements, directed = T)
        expected_type <- rep(type, length(E(expected_network)))

        ## established it reads graph, test if it sets attributes
        expect_equal(
          igraph::get.graph.attribute(actual_network)$name,
          expected_name
        )
        expect_equal(E(actual_network)$type, expected_type)
        expect_setequal(E(actual_network), E(expected_network))
        expect_true(igraph::is.directed(actual_network))
      })
    })

})

describe("load_geneset", {
  it("loads genes with sharps within the names", {
    # TODO:
    # Add fill*, trunctes the remainder of the line
    # commend character paramter - metioned by @Izaak comment.char='#'
  })


  it("returns null if geneset path is null", {
    output <- load_geneset(NULL)

    expect_equal(output, NULL)
  })

  it("fails to load the geneset due to bad file path", {
    file_path <- "../testGenesets/bad/file/path.tsv"

    expected_err_msg <- paste(
      "ERROR: geneset file does not exist:", 
      file_path)
    expect_error(load_geneset(file_path), expected_err_msg)
  })

  it("loads bad geneset", {
    file_path <- "../testGenesets/testGeneset_Bad_oneColumn.tsv"
    expected_err_msg <- "Your geneset file is incorrectly formatted. Please see documentation." #nolint msg

    expect_error(load_geneset(file_path), expected_err_msg)
  })

  it("loads csv formatted geneset", {
    file_path <- "../testGenesets/testGeneset1.csv"

    setid <- c("setA", "setA")
    gene <- c("1", "2")
    weight <- c(1, 1)
    expected_geneset <- data.table::data.table(setid, gene, weight)
    expected_extras <- NULL
    expected_output <- list(
      "geneset" = expected_geneset,
      "extras" = expected_extras)

    output <- load_geneset(file_path)
    expect_equal(output, expected_output)
  })

  it("loads a geneset with more than 3 columns, but only reads the first 3", {
    file_path <- "../testGenesets/test_geneset_5_cols.tsv"
     setid <- c("setA", "setA")
     gene <- c("1", "2")
     weight <- c(0.5, 0.2)
     expected_geneset <- data.table::data.table(setid, gene, weight)
     expected_extras <- NULL
     expected_output <- list(
      "geneset" = expected_geneset,
      "extras" = expected_extras)

     output <- load_geneset(file_path)

     expect_equal(output, expected_output)
  })
  
  it("loads unfiltered geneset with no multiplex network", {
    file_path <- "../testGenesets/testGeneset1.tsv"
    setid <- c("setA", "setA")
    gene <- c("1", "2")
    weight <- c(1, 1)
    expected_geneset <- data.table::data.table(setid, gene, weight)
    expected_extras <- NULL
    expected_output <- list(
      "geneset" = expected_geneset,
      "extras" = expected_extras)

    output <- load_geneset(file_path)

    expect_equal(output, expected_output)
  })


  it("loads weighted unfiltered geneset  with no multiplex network", {
    file_path <- "../testGenesets/testGeneset2.tsv"
    setid <- c("setB", "setB")
    gene <- c("2", "3")
    weight <- c(0.4, 0.6)
    expected_geneset <- data.table::data.table(setid, gene, weight)
    expected_extras <- NULL
    expected_output <- list(
      "geneset" = expected_geneset,
      "extras" = expected_extras)

    output <- load_geneset(file_path)

    expect_equal(output, expected_output)
  })

  it("loads un-weighted unfiltered geneset with 3 columns with no multiplex network", {  #nolint title
    ## Has a third column, but words, not weights.
    file_path <- "../testGenesets/testGeneset3.tsv"
    setid <- c("setC", "setC")
    gene <- c("4", "5")
    weight <- c(1, 1)
    expected_geneset <- data.table::data.table(setid, gene, weight)
    expected_extras <- NULL
    expected_output <- list(
      "geneset" = expected_geneset,
      "extras" = expected_extras)

    output <- load_geneset(file_path)

    expect_equal(output, expected_output)
  })

  it("loads filtered geneset filtered by the multiplex pool of nodes", {
    # Even though testGeneset2 has two genes: (2, 3) 
    #   - the pool of nodes has nodes 0, 1, 2.
    # The intersection of the two is specifically gene 2,
    #    which is the only output we expect to see.
    Pool_of_Nodes <- c("0", "1", "2")     #nolint
    nw.mpo <- data.frame(Pool_of_Nodes)   #nolint
    file_path <- "../testGenesets/testGeneset2.tsv"
    setid <- c("setB")
    gene <- c("2")
    weight <- c(0.4)
    expected_geneset <- data.table::data.table(setid, gene, weight)
    gene <- c("3")
    weight <- c(0.6)
    expected_extras <- data.table::data.table(setid, gene, weight)
    expected_output <- list(
      "geneset" = expected_geneset,
      "extras" = expected_extras)

    output <- expect_warning(load_geneset(file_path, nw.mpo))
    expect_equal(output, expected_output)
  })

  it("loads filtered geneset filtered by the multiplex pool of nodes, all filtered", {
    # Even though testGeneset2 has two genes: (2, 3)
    #- the pool of nodes has nodes 0, 1, 2.
    # The intersection of the two is specifically gene 2,
    # which is the only output we expect to see.
    Pool_of_Nodes <- c("0", "1") #nolint
    nw.mpo <- data.frame(Pool_of_Nodes) #nolint
    file_path <- "../testGenesets/testGeneset2.tsv"
    expected_geneset <- data.table::data.table(
        setid = character(),
        gene = character(),
        weight = double())
    setid <- c("setB", "setB")
    gene <- c("2", "3")
    weight <- c(0.4, 0.6)
    expected_extras <- data.table::data.table(setid, gene, weight)
    expected_output <- list(
      "geneset" = expected_geneset,
      "extras" = expected_extras)

    output <- expect_warning(load_geneset(file_path, nw.mpo))

    expect_equal(output, expected_output)
  })


})

describe("get_or_set_tau", {
  it("Returns numeric value if tau is numeric", {
    # Where tau values for layer one = 1.25 and layer 2 = 0.75

    opt_tau <- c(1.25, 0.75)
    expected_output <- opt_tau

    output <- get_or_set_tau(nw.mpo, opt_tau)

    expect_equal(output, expected_output)
  })

  it("Returns numeric values if tau is string list obtained via command line", {
    # Where tau values for layer one = 1.25 and layer 2 = 0.75
    opt_tau <- "1.25,0.75"
    expected_output <- c(1.25, 0.75)

    output <- get_or_set_tau(nw.mpo, opt_tau)

    expect_equal(output, expected_output)
  })

  it("Returns one for each layer with warning because of incorrect sum of tau vals", {  #nolint
    # Where tau values for layer one = 0.25 and layer 2 = 0.75, don't sum to 2
    opt_tau <- c(0.25, 0.75)
    expected_output <- c(1, 1)

    output <- expect_warning(get_or_set_tau(nw.mpo, opt_tau))

    expect_equal(output, expected_output)
  })

  it("Returns one for each layer with warning because of incorrect list size", {
    # Where tau values for layer one = 0.25 and layer 2 = 0.75, don't sum to 2
    opt_tau <- "0.75"
    expected_output <- c(1, 1)

    output <- expect_warning(get_or_set_tau(nw.mpo, opt_tau))

    expect_equal(output, expected_output)
  })
})

describe("chunk", {
  it("takes in a dataframe, and separates it by N values (number of folds for kfold cross val)", {
    gene_lists <- c("a", "b", "c", "d", "e", "f")
    n_folds <- 3
    expected_outputlength <- 3

    output <- chunk(gene_lists, n_folds)

    expect_equal(length(names(output)), expected_outputlength)
    expect_equal(output["1"][[1]], c("a", "b"))
    expect_equal(output["2"][[1]], c("c", "d"))
    expect_equal(output["3"][[1]], c("e", "f"))
  })

  it("Warns users if gene set length is smaller than num folds", {
    gene_lists <- c("a", "b", "c")
    n_folds <- 4
    expected_outputlength <- 3

    output <- expect_warning(chunk(gene_lists, n_folds))

    expect_equal(length(names(output)), expected_outputlength)
  })
})

describe("write_table", {
  it("Fails to write table if no table exists", {
    table <- NA
    path <- "./tableOutput"

    expect_warning(write_table(table, path))
  })

  it("Writes table to existing path", {
    faux_list_of_nodes <- c("A", "B", "C", "D")
    path <- "./someFakeFile.txt"

    write_table_stub <- mock()
    stub(write_table, "write.table", write_table_stub)
    write_table(faux_list_of_nodes, path)

    expect_called(write_table_stub, 1)
    expect_args(write_table_stub, 1, faux_list_of_nodes, path, "\t", F, T, F)
  })

  it("Writes table to nonexisting path", {
    faux_list_of_nodes <- c("A", "B", "C", "D")
    bad_path <- "./someOtherPath/recursiveDepth/output.txt"
    expected_dir <- "recursiveDepth"

    write_table_stub <- mock()
    stub(write_table, "write.table", write_table_stub)
    write_table(faux_list_of_nodes, bad_path)
    check_dir <- system("ls ./someOtherPath/", intern = T)

    expect_equal(check_dir[1], expected_dir)
    expect_called(write_table_stub, 1)
    expect_args(write_table_stub, 1, faux_list_of_nodes, bad_path, "\t", F, T, F)
  })
})

describe("get_base_name", {
  it("extracts a basename from a full file path", {
    full_fp <- "/Users/mattlane/projects/RWRtoolkit/inst/example_data/net.Rdata"

    actual_basename <- get_base_name(full_fp)

    expected_basename <- "net"
    expect_equal(actual_basename, expected_basename)
  })

  it("extracts a basename from a full url with query params", {
    full_url <- "http://www.github.com/username/net.Rdata?raw=True"

    actual_basename <- get_base_name(full_url)

    expected_basename <- "net"
    expect_equal(actual_basename, expected_basename)
  })

  it("extracts a basename from only a basename", {
    base_with_ext <- "net.Rdata"

    actual_basename <- get_base_name(base_with_ext)

    expected_basename <- "net"
    expect_equal(actual_basename, expected_basename)
  })

  
  it("extracts a basename from only a basename with no ext", {
    base_with_ext <- "net"

    actual_basename <- get_base_name(base_with_ext)

    expected_basename <- "net"
    expect_equal(actual_basename, expected_basename)
  })

  it("extracts a basename but leaves non-rdata ext", {
    base_with_nonrdata_ext <- "net.something.r"

    actual_basename <- get_base_name(base_with_nonrdata_ext)

    expected_basename <- "net.something.r"
    expect_equal(actual_basename, expected_basename)
  })

})

describe("get_file_path", {
  it("appends default extension to a filename with no specified out directory", { #nolint
    base_filename <- "someFile"
    expected_basefilename <- "someFile.tsv"

    outfilename <- get_file_path(base_filename)

    expect_equal(outfilename, expected_basefilename)
  })


  it("appends extension to a filename with no specified out directory", {
    base_filename <- "someFile"
    expected_basefilename <- "someFile.csv"
    file_extension <- ".csv"

    outfilename <- get_file_path(base_filename, ext = file_extension)

    expect_equal(outfilename, expected_basefilename)
  })
  it("appends extension to a filename with specified out directory", {
    base_filename <- "someFile"
    outdirpath <- "/some/path/to/file"
    expected_basefilename <- "/some/path/to/file/someFile.tsv"

    outfilename <- get_file_path(base_filename, outdir = outdirpath)

    expect_equal(outfilename, expected_basefilename)
  })
})

describe("dump_layers", {
  it("extracts specific layers and writes them to a table", {
    load("../testMultiplex/network.Rdata")
    fakefile_path <- "somefile_path/network.tsv"
    expected_layer1name <- "m1"
    expected_layer2name <- "m2"
    expected_num_calls <- 2
    from <- c("0", "0", "1")
    to <- c("1", "2", "2")
    weight <- c(0.5, 1.0, 0.5)
    type <- c("m1", "m1", "m1")
    expected_layer1_df <- data.frame(from, to, weight, type)
    from <- c("0", "0", "1", "2")
    to <- c("2", "3", "2", "3")
    weight <- c(1, 1, 1, 1)
    type <- c("m2", "m2", "m2", "m2")
    expected_layer2_df <- data.frame(from, to, weight, type)

    get_filepath_stub <- mock(fakefile_path, fakefile_path)
    write_table_stub <- mock()
    stub(dump_layers, "get_file_path", get_filepath_stub)
    stub(dump_layers, "write_table", write_table_stub)
    dump_layers(nw.mpo)

    expect_called(get_filepath_stub, expected_num_calls)
    expect_args(get_filepath_stub, 1, expected_layer1name, NULL)
    expect_args(get_filepath_stub, 2, expected_layer2name, NULL)
    expect_called(write_table_stub, 2)
    expect_args(write_table_stub, 1, expected_layer1_df, fakefile_path)
    expect_args(write_table_stub, 2, expected_layer2_df, fakefile_path)
  })

  it("writes nothing to a table if no layers exist in object", {
    load("../testMultiplex/network.Rdata")
    expected_num_calls <- 0
    fakefile_path <- "somefile_path/network.tsv"
    nw.mpo$m1 <- NULL # Remove layers from network.
    nw.mpo$m2 <- NULL

    get_filepath_stub <- mock(fakefile_path, fakefile_path)
    write_table_stub <- mock()
    stub(dump_layers, "get_file_path", get_filepath_stub)
    stub(dump_layers, "write_table", write_table_stub)
    dump_layers(nw.mpo)

    expect_called(get_filepath_stub, expected_num_calls)
    expect_called(write_table_stub, expected_num_calls)
  })
})

describe("dump_nodes", {
  it("throws an error if a multiplexMH object is not passed in as mpo", {
    nw.mpo <- NULL #nolint
    expected_num_calls <- 0

    write_table_stub <- mock()
    stub(dump_nodes, "write_table", write_table_stub)
    expect_error(dump_nodes(nw.mpo))

    expect_called(write_table_stub, expected_num_calls)
  })
  it("extracts nodes from the multiplexMH object and saves to file", {
    load("../testMultiplex/network.Rdata")
    mpo <- nw.mpo1 # loaded from the network.Rdata
    mpo_nodes <- c("0", "1", "2")
    expected_dataframe <- data.frame(list(mpoNodes=mpo_nodes))
    expected_num_calls <- 1
    fakefile_path <- "./someFile/path/file.txt"

    write_table_stub <- mock()
    get_filepath_stub <- mock(fakefile_path)
    stub(dump_nodes, "write_table", write_table_stub)
    stub(dump_nodes, "get_file_path", get_filepath_stub)


    dump_nodes(nw.mpo1)
    expect_called(get_filepath_stub, expected_num_calls)
    expect_called(write_table_stub, expected_num_calls)
    expect_args(write_table_stub, 1, expected_dataframe, fakefile_path)
  })
})


# removes test files so as not to clutter everything up
teardown({
    system("rm -rf ./someOtherPath") # nolint removes path to written table via write_table
    system("rm ./tableOutput")
  },
  env = parent.frame()
)
