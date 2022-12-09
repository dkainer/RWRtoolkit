describe("RWR_netstats", {
  describe("load_network", {

    expected_from  <- c("B", "E", "C", "E", "F", "F", "G", "H", "H", "H")
    expected_to <- c("A", "B", "B", "D", "E", "C", "E", "E", "F", "C")
    expected_elements <- c()
    for (i in seq(1, length(expected_from))) {
      expected_elements <- append(expected_elements, expected_from[i])
      expected_elements <- append(expected_elements, expected_to[i])
    }

    it("loads an igraph object from an edgelist", {
      path_to_edgelist <- "../testNetworks/abc_layer1.tsv"

      actual_network <- RWRtoolkit::load_network(path_to_edgelist)

      expected_network <- igraph::make_graph(expected_elements, directed = F)

      expect_setequal(V(actual_network)$name,
                   V(expected_network)$name)
      expect_setequal(E(actual_network),
                  E(expected_network))
      expect_true(!is.directed(actual_network))
    })

    it("loads a directed igraph object from an edgelist", {
      path_to_edgelist <- "../testNetworks/abc_layer1.tsv"
      directed <- TRUE
      type <- "testnetwork"
      name <- "TESTLAYER"

      actual_network <- RWRtoolkit::load_network(
        path_to_edgelist = path_to_edgelist,
        type = type,
        name = name,
        directed = directed
      )

      expected_name <- name
      expected_network <- igraph::make_graph(expected_elements, directed = T)
      expected_type <- rep(type, length(E(expected_network)))

      ## established it reads graph, test if it sets attributes
      expect_equal(igraph::get.graph.attribute(actual_network)$name,
                   expected_name)
      expect_equal(E(actual_network)$type, expected_type)
      expect_setequal(E(actual_network), E(expected_network))
      expect_true(igraph::is.directed(actual_network))
    })
  })

  describe("load_flist", {

    expected_filepaths <- c(
      "../testNetworks/abc_layer1.tsv",
      "../testNetworks/abc_layer2.tsv",
      "../testNetworks/abc_layer3.tsv"
    )
    expected_layers <- c("layer1", "layer2", "layer3")
    expected_group <- rep(1, 3)
    expected_flist <- data.table::data.table(
      data.frame(list(
        nwfile = expected_filepaths,
        nwname = expected_layers,
        nwgroup = expected_group
      ))
    )

    it("loads an flist from a file", {
      flist_path <- "../testFlists/abc_flist.txt"
      actual_flist <- RWRtoolkit::load_flist(flist_path)

      expect_equal(actual_flist, expected_flist)
    })

    it("loads the first 3 columns of an flist from file", {
      flist_path <- "../testFlists/test_flist_5cols.tsv"
      actual_flist <- RWRtoolkit::load_flist(flist_path)

      expect_equal(actual_flist, expected_flist)
    })
  })

  describe("make_dummy_multiplex", {

    filepaths <- c(
      "../testNetworks/abc_layer1.tsv",
      "../testNetworks/abc_layer2.tsv",
      "../testNetworks/abc_layer3.tsv"
    )
    layers <- c("layer1", "layer2", "layer3")
    groups <- rep(1, 3)
    flist <- data.table::data.table(
      data.frame(list(
        nwfile = filepaths,
        nwname = layers,
        nwgroup = groups
      ))
    )

    it("creates a dummy multiplex from a filepath", {
      flist_file_path <- "../testFlists/abc_flist.txt"

      load_flist_stub <- mock(flist)
      load_network_stub <- mock(
          igraph::make_graph(c("A", "B")),
          igraph::make_graph(c("B", "C")),
          igraph::make_graph(c("C", "D")),
          )

      stub(
        RWRtoolkit::make_dummy_multiplex,
        "load_flist",
        load_flist_stub
      )
      stub(
        RWRtoolkit::make_dummy_multiplex,
        "load_network",
        load_network_stub
      )

      RWRtoolkit::make_dummy_multiplex(flist_file_path)

      expect_called(load_flist_stub, 1)
      expect_called(load_network_stub, 3)

      expect_args(load_flist_stub, 1, flist_file_path)
      expect_args(load_network_stub, 1, filepaths[1], layers[1], layers[1])
      expect_args(load_network_stub, 2, filepaths[2], layers[2], layers[2])
      expect_args(load_network_stub, 3, filepaths[3], layers[3], layers[3])
    })

    it("creates a dummy multiplex from an flist data frame", {
      flist_file_path <- flist

      load_flist_stub <- mock()
      load_network_stub <- mock(
        igraph::make_graph(c("A", "B")),
        igraph::make_graph(c("B", "C")),
        igraph::make_graph(c("C", "D")),
      )

      stub(
        RWRtoolkit::make_dummy_multiplex,
        "load_flist",
        load_flist_stub
      )
      stub(
        RWRtoolkit::make_dummy_multiplex,
        "load_network",
        load_network_stub
      )

      RWRtoolkit::make_dummy_multiplex(flist_file_path)

      expect_called(load_flist_stub, 0)
      expect_called(load_network_stub, 3)

      expect_args(load_network_stub, 1, filepaths[1], layers[1], layers[1])
      expect_args(load_network_stub, 2, filepaths[2], layers[2], layers[2])
      expect_args(load_network_stub, 3, filepaths[3], layers[3], layers[3])
    })

    it("throws an error when no flist or multiplex", {
      flist_file_path <- "bad_test_path.tsv"

      load_flist_stub <- mock()
      load_network_stub <- mock()

      stub(
        RWRtoolkit::make_dummy_multiplex,
        "load_flist",
        load_flist_stub
      )
      stub(
        RWRtoolkit::make_dummy_multiplex,
        "load_network",
        load_network_stub
      )

      expect_error(
        RWRtoolkit::make_dummy_multiplex(flist_file_path), 
        "Input must be either a path to an flist or a flist dataframe."
      )
    })
  })

  describe("merging methods: ", {
     # The Multiplex:
     # ABC Layer 1      # ABC Layer 2            # ABC Layer 3
     # A◄───B◄──────C   # G        H             #    H◄─┐
     #      ▲       ▲   # │        │             #       │
     #      │       │   # └──►D    └─►E       F  #  D◄───E◄┐        F
     # D◄───E◄─┬──F─┤   #     │       │       │  #         │        │
     #      ▲  │  ▲ │   #     └───►B◄─┘──►C◄──┘  #         ├───►C◄──┘
     #      │  │  │ │   #          │      │      #         │
     #  G───┘  H ─┴─┘              └──►A◄─┘      #         A
     #
     # Merges down to (not illustrating multi-edges):
     # ┌────────┬───────┐
     # G        ├───H──►F
     # │        ▼       │
     # └─►D◄────E───┐   │
     #    │     │   ▼   │
     #    └─►B◄─┴───C◄──┘
     #       ▲      │
     #       └──►A◄─┘
    mpo_filepath <- "../testMultiplex/abc_multiplex.rdata"
    load(mpo_filepath)

    expected_edges <- c(
      "B", "E",  "B", "C",
      "B", "A",  "E", "F",
      "E", "G",  "E", "H",
      "E", "D",  "C", "F",
      "C", "H",  "F", "H",
      "G", "D",  "E", "H",
      "B", "D",  "B", "E",
      "E", "C",  "C", "F",
      "B", "A",  "C", "A",
      "E", "A",  "E", "H",
      "E", "D",  "C", "A",
      "C", "F"
    )
    expected_type1 <- rep("layer1", 10)
    expected_type2 <- rep("layer2", 8)
    expected_type3 <- rep("layer3", 5)
    expected_type <- c(
      expected_type1,
      expected_type2,
      expected_type3
    )

    expected_weightnorm1 <- rep(1 / 10, 10)
    expected_weightnorm2 <- rep(1 / 8, 8)
    expected_weightnorm3 <- rep(1 / 5, 5)
    expected_weightnorm <- c(
      expected_weightnorm1,
      expected_weightnorm2,
      expected_weightnorm3
    )
    expected_merged <- igraph::make_graph(
      edges = expected_edges,
      directed = F
    )
    E(expected_merged)$weight <- 1
    E(expected_merged)$weightnorm <- expected_weightnorm
    E(expected_merged)$type <- expected_type

    describe("merged_with_all_edges", {
      it("merges all layers of an MPO", {
        actual_merged <- RWRtoolkit::merged_with_all_edges(nw.mpo) #nolint

        expect_setequal(E(actual_merged), E(expected_merged))
        expect_setequal(E(actual_merged)$weight, E(expected_merged)$weight)
        expect_setequal(E(actual_merged)$weightnorm,
                        E(expected_merged)$weightnorm)
        expect_setequal(E(actual_merged)$type, E(expected_merged)$type)
        expect_setequal(V(actual_merged)$name, V(expected_merged)$name)
      })
    })

    describe("merged_with_edgecounts", {
      it("merges layers and removes all duplicate edges w/ edgecount as weight", { #nolint
        actual_merged <- RWRtoolkit::merged_with_edgecounts(nw.mpo)

        expected_simple_merged <- igraph::simplify(expected_merged)

        #  B--E 2
        #  B--C 1
        #  B--D 1
        #  B--A 2
        #  E--C 1
        #  E--F 1
        #  E--G 1
        #  E--H 3
        #  E--D 2
        #  E--A 1
        #  C--F 3
        #  C--H 1
        #  C--A 1
        #  F--H 1
        #  G--D 1

        expected_counts <- c(2, 1, 1, 2, 1, 1, 1, 3, 2, 1, 3, 1, 1, 1, 1)
        E(expected_simple_merged)$weight <- expected_counts

        for (i in seq(1, length(E(actual_merged)))) {
          edge <- igraph::E(actual_merged)[i]
        }

        expect_setequal(
          E(actual_merged)$weight, 
          E(expected_simple_merged)$weight)
        expect_setequal(E(actual_merged), E(expected_simple_merged))
      })
    })
  })

  describe("get_name", {
    unnamed_network <- igraph::make_graph(edges = c("A", "B"))

    it("returns the name of a named network", {
      network_name <- "test_network"
      named_network <- igraph::set.graph.attribute(
                          unnamed_network,
                          "name",
                          network_name)

      actual_name <- RWRtoolkit::get_name(named_network)


      expect_equal(actual_name, network_name)
    })

    it("returns the function default name from an unnamed network", {
      actual_name <- RWRtoolkit::get_name(unnamed_network)

      expected_name <- "<G>"
      expect_equal(actual_name, expected_name)
    })
  })

  describe("basic_statistics", {
    it("calculates basic stats for a network e/v-count, diameter", {
      #    ┌──────1────┐
      #    ▼           ▼
      # ┌──3──┐    ┌───2───┐
      # ▼     ▼    ▼       ▼
      # 6     7 ┌──4──┐    5
      #         ▼     ▼    │
      #         8     9    ▼
      #                   10
      set.seed(42)
      tree_network <- make_tree(n = 10)


      get_name_stub <- mock("TREE")
      stub(
        RWRtoolkit::basic_statistics,
        "get_name",
        get_name_stub
      )

      actual_messages <- capture_messages(
        RWRtoolkit::basic_statistics(tree_network, directed = T))

      expected_messages <- c(
        "Network stats for network TREE\n",
        "==============================\n",
        "Number of nodes : 10\n",
        "Number of edges : 9\n",
        "Diameter        : 3.00\n"
      )

      expect_equal(actual_messages, expected_messages)
      expect_called(get_name_stub, 1)
      expect_args(get_name_stub, 1, tree_network, "<G>")
    })
  })

  describe("basic_statistics_multiplex", {
    it("describes basic statistics for each layer of a multiplex", {
      mp_faux <- list()
      mp_faux$layer1 <- igraph::make_graph(edges = c("A", "B"))
      mp_faux$layer2 <- igraph::make_graph(edges = c("A", "C"))
      mp_faux$Number_of_Layers <- 2 # nolint

      basic_stats_stub <- mock()
      stub(
        RWRtoolkit::basic_statistics_multiplex,
        "basic_statistics",
        basic_stats_stub
      )

      RWRtoolkit::basic_statistics_multiplex(mp_faux)
      expect_called(basic_stats_stub, 2)
      expect_args(basic_stats_stub, 1, mp_faux$layer1, "layer1", F)
      expect_args(basic_stats_stub, 2, mp_faux$layer2, "layer2", F)
    })
  })

  describe("jaccard_score_edges", {

  })

  describe("overlap_score", {

  })

  describe("compare_networks", {

  })

  describe("overlap_pair", {

  })

  describe("overlap_many_pairwise", {

  })

  describe("overlap_many_vs_reference", {

  })

  describe("get_tau", {

  })

  describe("exclusivity", {

  })

  describe("RWR_netstats", {

  })


})
