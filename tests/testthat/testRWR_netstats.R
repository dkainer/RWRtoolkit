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

      actual_network <- load_network(path_to_edgelist)

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
      actual_flist <- load_flist(flist_path)

      expect_equal(actual_flist, expected_flist)
    })

    it("loads the first 3 columns of an flist from file", {
      flist_path <- "../testFlists/test_flist_5cols.tsv"
      actual_flist <- load_flist(flist_path)

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

        expect_setequal(E(actual_merged$merged_network),
                        E(expected_merged))
        expect_setequal(E(actual_merged$merged_network)$weight,
                        E(expected_merged)$weight)
        expect_setequal(E(actual_merged$merged_network)$weightnorm,
                        E(expected_merged)$weightnorm)
        expect_setequal(E(actual_merged$merged_network)$type,
                        E(expected_merged)$type)
        expect_setequal(V(actual_merged$merged_network)$name,
                        V(expected_merged)$name)
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

        for (i in seq(1, length(E(actual_merged$merged_network)))) {
          edge <- igraph::E(actual_merged$merged_network)[i]
        }

        expect_setequal(
          E(actual_merged$merged_network)$weight, 
          E(expected_simple_merged)$weight)
        expect_setequal(E(actual_merged$merged_network), 
                        E(expected_simple_merged))
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

  describe("calculate_basic_statistics", {
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
        RWRtoolkit::calculate_basic_statistics,
        "get_name",
        get_name_stub
      )

      actual_messages <- capture_messages(
        RWRtoolkit::calculate_basic_statistics(tree_network,
                                    directed = T,
                                    verbose = T))

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
        RWRtoolkit::calculate_basic_statistics,
        "get_name",
        get_name_stub
      )

      actual_stats <- RWRtoolkit::calculate_basic_statistics(
                        tree_network,
                        directed = T
      )

      expected_stats <- data.frame(list(
        network_name = "TREE",
        number_of_nodes = 10,
        number_of_edges = 9,
        diameter = 3.00
      ))
      expect_equal(actual_stats, expected_stats)
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
        "calculate_basic_statistics",
        basic_stats_stub
      )

      RWRtoolkit::basic_statistics_multiplex(mp_faux)
      expect_called(basic_stats_stub, 2)
      expect_args(basic_stats_stub, 1, mp_faux$layer1, "layer1", F)
      expect_args(basic_stats_stub, 2, mp_faux$layer2, "layer2", F)
    })
  })

  describe("check weighted edges", {
    it("throws a warning when no weighted edges exist", {
          network <- igraph::make_graph(edges = c("A", "B", "B", "C"))

          expected_warning <- paste(
            "Network <X> has no weighted edges.",
            "All scores will be zero."
          )
          expect_warning(check_weighted_edges(network, "<X>"))
    })

    it("does not throw a warning when weighted edges exist", {
          network_x <- igraph::make_graph(edges = c("A", "B", "B", "C"))
          E(network_x)$weight <- 1

          # No warnings thrown: 
          check_weighted_edges(network_x, "<x>")
          
    })

  })

  describe("Comparison of 2 networks", {
    network_a <- igraph::make_graph(edges = c("A", "B", "B", "C", "C", "D"))
    network_b <- igraph::make_graph(edges = c("C", "D", "D", "E", "E", "F"))
    E(network_a)$weight <- 1
    E(network_b)$weight <- 1

    describe("jaccard_score_edges", {
      # |A n B| / (|A| + |B| - |AnB|)

      it("calculates a jaccard score between two networks", {
        # (C--D) / (A--B, B--C, C--D, E--F, D--E)
        actual_score <- RWRtoolkit::jaccard_score_edges(network_a, network_b)

        expected_score <- 0.2
        expect_equal(actual_score, expected_score)
      })
    })

    describe("overlap_score", {
      it("calculates an overlap score between two netwrosk", {
        # Gets intersection
        # Sums weights
        # divides sum of intersected weights
        actual_overlap <- RWRtoolkit::overlap_score(network_a, network_b)

        sum_intersected_edges <- 1
        total_edges_net_b <- length(E(network_b))
        expected_overlap_score <- sum_intersected_edges / total_edges_net_b

        expect_equal(actual_overlap, expected_overlap_score)
      })

      it("throws a warning when one of the networks has no weights", {
        unweighted_net <- network_b
        unweighted_net <- igraph::remove.edge.attribute(
                                    unweighted_net,
                                    "weight")

        expect_warning(
          RWRtoolkit::overlap_score(network_a, unweighted_net),
          "Network <H> has no weighted edges. All scores will be zero."
        )
      })

    })

    describe("compare_networks", {
      it("compares networks with jaccard", {
        metric <- "jaccard"
        jaccard_stub_output <- 0.5

        jaccard_stub <- mock(jaccard_stub_output)
        overlap_stub <- mock()

        stub(
          RWRtoolkit::compare_networks,
          "jaccard_score_edges",
          jaccard_stub
        )
        stub(
          RWRtoolkit::compare_networks,
          "overlap_score",
          overlap_stub
        )

        actual_output <- RWRtoolkit::compare_networks(
                        network_a,
                        network_b,
                        metric)

        expect_equal(actual_output, jaccard_stub_output)
        expect_called(jaccard_stub, 1)
        expect_called(overlap_stub, 0)
        expect_args(jaccard_stub, 1, network_a, network_b, FALSE)
      })

       it("compares networks with overlap", {
         metric <- "overlap"
         overlap_stub_output <- 0.5

         jaccard_stub <- mock()
         overlap_stub <- mock(overlap_stub_output)

         stub(
           RWRtoolkit::compare_networks,
           "jaccard_score_edges",
           jaccard_stub
         )
         stub(
           RWRtoolkit::compare_networks,
           "overlap_score",
           overlap_stub
         )

         actual_output <- RWRtoolkit::compare_networks(
           network_a,
           network_b,
           metric
         )

         expect_equal(actual_output, overlap_stub_output)
         expect_called(jaccard_stub, 0)
         expect_called(overlap_stub, 1)
         expect_args(overlap_stub, 1, network_a, network_b, FALSE)
       })
    })

    describe("overlap_pair", {
      it("runs compare_networks:", {
        expected_output <- "expected_output"
        compare_stub <- mock(expected_output)
        arg1 <- "A"
        arg2 <- "B"

        stub(
          RWRtoolkit::overlap_pair,
          "compare_networks",
          compare_stub
        )

        actual_output <- RWRtoolkit::overlap_pair(arg1, arg2)

        expect_equal(actual_output, expected_output)
        expect_args(compare_stub, 1, arg1, arg2)
      })
    })
  })

  describe("Overlap w/ mpo objects", {
    mpo_filepath <- "../testMultiplex/abc_multiplex.rdata"
    load(mpo_filepath)

    describe("overlap_many_pairwise", {

      it("calls compare networks 6 times for pairwise comparisons w/ jaccard", {
        metric <- "jaccard"
        comp_networks_stub <- mock(
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0)

        stub(
          RWRtoolkit::overlap_many_pairwise,
          "compare_networks",
          comp_networks_stub
        )

        actual_output <- RWRtoolkit::overlap_many_pairwise(
          nw.mpo,
          metric = metric
        )

        expected_output <- matrix(1, nrow = 3, ncol = 3)
        rownames(expected_output) <- c("layer1", "layer2", "layer3")
        colnames(expected_output) <- c("layer1", "layer2", "layer3")

        expect_equal(actual_output, expected_output)
        expect_called(comp_networks_stub, 9)
        expect_args(comp_networks_stub, 1, nw.mpo[[1]], nw.mpo[[1]], metric, F)
        expect_args(comp_networks_stub, 2, nw.mpo[[1]], nw.mpo[[2]], metric, F)
        expect_args(comp_networks_stub, 3, nw.mpo[[1]], nw.mpo[[3]], metric, F)
        expect_args(comp_networks_stub, 4, nw.mpo[[2]], nw.mpo[[1]], metric, F)
        expect_args(comp_networks_stub, 5, nw.mpo[[2]], nw.mpo[[2]], metric, F)
        expect_args(comp_networks_stub, 6, nw.mpo[[2]], nw.mpo[[3]], metric, F)
        expect_args(comp_networks_stub, 7, nw.mpo[[3]], nw.mpo[[1]], metric, F)
        expect_args(comp_networks_stub, 8, nw.mpo[[3]], nw.mpo[[2]], metric, F)
        expect_args(comp_networks_stub, 9, nw.mpo[[3]], nw.mpo[[3]], metric, F)
      })

      it("cals compare networks 6 times for pairwise comparison w/ overlap", {
        metric <- "overlap"
        comp_networks_stub <- mock(
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0, 
          1.0)
        stub(
          RWRtoolkit::overlap_many_pairwise,
          "compare_networks",
          comp_networks_stub
        )

        actual_output <- RWRtoolkit::overlap_many_pairwise(nw.mpo)

        expected_output <- matrix(1, nrow = 3, ncol = 3)
        rownames(expected_output) <- c("layer1", "layer2", "layer3")
        colnames(expected_output) <- c("layer1", "layer2", "layer3")

        expect_equal(actual_output, expected_output)
        expect_called(comp_networks_stub, 9)
        expect_args(comp_networks_stub, 1, nw.mpo[[1]], nw.mpo[[1]], metric, F)
        expect_args(comp_networks_stub, 2, nw.mpo[[1]], nw.mpo[[2]], metric, F)
        expect_args(comp_networks_stub, 3, nw.mpo[[1]], nw.mpo[[3]], metric, F)
        expect_args(comp_networks_stub, 4, nw.mpo[[2]], nw.mpo[[1]], metric, F)
        expect_args(comp_networks_stub, 5, nw.mpo[[2]], nw.mpo[[2]], metric, F)
        expect_args(comp_networks_stub, 6, nw.mpo[[2]], nw.mpo[[3]], metric, F)
        expect_args(comp_networks_stub, 7, nw.mpo[[3]], nw.mpo[[1]], metric, F)
        expect_args(comp_networks_stub, 8, nw.mpo[[3]], nw.mpo[[2]], metric, F)
        expect_args(comp_networks_stub, 9, nw.mpo[[3]], nw.mpo[[3]], metric, F)
      })
    })

    describe("overlap_many_vs_reference", {
      it("compares against a ref networks", {
        metric <- "overlap"
        ref_net <- igraph::make_graph(edges = c("A", "B", "B", "C"))

        # not necessary for stub to run, but good to see
        E(ref_net)$weight <- 1

        comp_networks_stub <- mock(1.0, 0.5, 0.2)
        stub(
          RWRtoolkit::overlap_many_vs_reference,
          "compare_networks",
          comp_networks_stub
        )

        actual_output <- RWRtoolkit::overlap_many_vs_reference(
                                        nw.mpo,
                                        ref_net)

        expected_output <- c(1.0, 0.5, 0.2) # arbitrary scores*
        names(expected_output) <- c("layer1", "layer2", "layer3")

        expect_equal(actual_output, expected_output)

        expect_called(comp_networks_stub, 3)
        expect_args(comp_networks_stub, 1, ref_net, nw.mpo[[1]], metric, F)
        expect_args(comp_networks_stub, 2, ref_net, nw.mpo[[2]], metric, F)
        expect_args(comp_networks_stub, 3, ref_net, nw.mpo[[3]], metric, F)
      })
    })

    describe("get_tau", {
      it("calculates tau for a multiplex based on ref network", {
        ref_net <- igraph::make_graph(
                            edges = c("C", "F",   "G", "D",   "H", "E", 
                                      "B", "A",   "C", "A",   "F", "G"),
                            directed = F)
        E(ref_net)$weight <- 1
        metric <- "overlap"

        many_v_ref_output <- c(1.0, 0.5, 1.0) #arbitrary scores*
        names(many_v_ref_output) <- c("layer1", "layer2", "layer3")

        overlap_many_v_ref_stub <- mock(many_v_ref_output)

        stub(
          RWRtoolkit::calculate_tau,
          "overlap_many_vs_reference",
          overlap_many_v_ref_stub
        )

        actual_tau <- RWRtoolkit::calculate_tau(nw.mpo, ref_net)

        expected_output <- c(
          1.0 / 2.5 * nw.mpo$Number_of_Layers,
          0.5 / 2.5 * nw.mpo$Number_of_Layers,
          1.0 / 2.5 * nw.mpo$Number_of_Layers
        )

        names(expected_output) <- c("layer1", "layer2", "layer3")
        expect_equal(actual_tau, expected_output)
      })
    })
  })



  describe("exclusivity", {
    it("calls exclusivity and returns a message for each layer", {
      mpo_filepath <- "../testMultiplex/abc_multiplex.rdata"
      load(mpo_filepath)
      
      # 1           2
      # E--A      # B--A B--A
      # E--C      # B--E E--B
      # E--F      # C--A A--C
      # E--G      # E--D E--D
      # D--B      
      # F--H        3
      # B--C      # E--H H--E E--H
      # C--H      # F--C F--C C--F
      # G--D
      
      actual_output <- RWRtoolkit::exclusivity(nw.mpo)

      expected_in_3 <-  2 / 15
      expected_in_2 <-  4 / 15
      expected_in_1 <-  9 / 15
      pct_rows <- c(expected_in_1, expected_in_2, expected_in_3)

      expected_output <- data.frame(
                          list(
                            n_layers = seq(1:3),
                            pct_found = pct_rows)
                          )
      print("Actual")
      print(actual_output)
      print("Expected")
      print(expected_output)
      expect_equal(actual_output, expected_output, tolerance = 1e-3)
    })
  })

  describe("RWR_netstats", {
    mpo_filepath <- "../testMultiplex/abc_multiplex.rdata"
    load(mpo_filepath)
    flist_file_path <- "../testFlists/abc_flist.txt"
    net1_file_path <- "../testNetworks/abc_layer1.tsv"
    net1 <- read.table(net1_file_path)
    net1 <- graph.data.frame(net1)

    net2_file_path <- "../testNetworks/abc_layer2.tsv"
    net2 <- read.table(net2_file_path)
    net2 <- graph.data.frame(net2)

    it("throws error if no networks or paths included", {
      expected_error_message <- paste(
        "[ERROR] You must supply one of the following arguments:",
        "       data",
        "       flist",
        "       (data or flist) and network_1",
        "       network_1 and/or network_2",
        sep = "\n"
      )

      expect_error(
        RWRtoolkit::RWR_netstats(),
        expected_error_message,
        fixed = T
      )
    })

    it("calculates basic statistics with just network 1", {
      load_network_stub <- mock(net1)
      basic_stats_stub <- mock("basic_statistics__mock_output")
      load_flist_stub <- mock()
      make_dummy_multiplex_stub <- mock()
      basic_stats_multiplex_stub <- mock()
      overlap_many_pairwise_stub <- mock()
      overlap_many_vs_reference_stub <- mock()
      overlap_pair_stub <- mock()
      calculate_tau_stub <- mock()
      merged_with_all_edges_stub <- mock()
      merged_with_edgecounts_stub <- mock()
      exclusivity_stub <- mock()

      stub(RWRtoolkit::RWR_netstats,
        "calculate_basic_statistics",
        basic_stats_stub
      )
      stub(RWRtoolkit::RWR_netstats, "load_network", load_network_stub)
      stub(RWRtoolkit::RWR_netstats, "load_flist", load_flist_stub)
      stub(RWRtoolkit::RWR_netstats,
          "make_dummy_multiplex",
          make_dummy_multiplex_stub
      )

      stub(RWRtoolkit::RWR_netstats,
          "basic_statistics_multiplex",
          basic_stats_multiplex_stub
      )
      stub(RWRtoolkit::RWR_netstats,
          "overlap_many_pairwise",
          overlap_many_pairwise_stub
      )
      stub(RWRtoolkit::RWR_netstats,
          "overlap_many_vs_reference",
          overlap_many_vs_reference_stub
      )
      stub(RWRtoolkit::RWR_netstats, "overlap_pair", overlap_pair_stub)
      stub(RWRtoolkit::RWR_netstats, "calculate_tau", calculate_tau_stub)
      stub(RWRtoolkit::RWR_netstats,
          "merged_with_all_edges",
          merged_with_all_edges_stub
      )
      stub(RWRtoolkit::RWR_netstats,
          "merged_with_edgecounts",
          merged_with_edgecounts_stub
      )
      stub(RWRtoolkit::RWR_netstats, "exclusivity", exclusivity_stub)

      actual_output <- RWRtoolkit::RWR_netstats(
          data = NULL,
          flist  = NULL,
          network_1 = net1_file_path,
          network_2 = NULL,
          basic_statistics = T
      )

      expected_output <- list(
        base_stats_net1 = "basic_statistics__mock_output"
      )

      expect_equal(actual_output, expected_output)

      expect_called(load_network_stub, 1)
      expect_called(basic_stats_stub, 1)
      expect_called(load_flist_stub, 0)
      expect_called(make_dummy_multiplex_stub, 0)
      expect_called(basic_stats_multiplex_stub, 0)
      expect_called(overlap_many_pairwise_stub, 0)
      expect_called(overlap_many_vs_reference_stub, 0)
      expect_called(overlap_pair_stub, 0)
      expect_called(calculate_tau_stub, 0)
      expect_called(merged_with_all_edges_stub, 0)
      expect_called(merged_with_edgecounts_stub, 0)
      expect_called(exclusivity_stub, 0)

      expect_args(load_network_stub, 1, net1_file_path, "network_1", F)
      expect_args(basic_stats_stub, 1, net1, F)
    })

    it("calculates basic statistics with all available options", {
      load_multiplex_data_stub <- mock(list(nw.mpo = nw.mpo))

      load_network_stub <- mock(net1, net2)
      basic_stats_stub <- mock(
                    "basic_statistics_mock_output_net1",
                    "basic_statistics_mock_output_net2"
                    )
      load_flist_stub <- mock("load_flist_mock_output")
      make_dummy_multiplex_stub <- mock("make_dummy_multiplex_mock_output")
      basic_stats_multiplex_stub <- mock("basic_stats_multiplex_mock_output")
      overlap_many_pairwise_stub <- mock(
                              "overlap_many_pairwise_mock_output_jaccard",
                              "overlap_many_pairwise_mock_output_mpo_overlap"
                                    )
      overlap_many_vs_reference_stub <- mock(
                                      "overlap_many_vs_reference_mock_output"
                                      )
      overlap_pair_stub <- mock(
        "overlap_pair_mock_output_overlap",
        "overlap_pair_mock_output_jaccard"
        )
      calculate_tau_stub <- mock("calculate_tau_mock_output")
      merged_with_all_edges_stub <- mock("merged_with_all_edges_mock_output")
      merged_with_edgecounts_stub <- mock("merged_with_edgecounts_mock_output")
      exclusivity_stub <- mock("exclusivity_mock_output")


      stub(RWRtoolkit::RWR_netstats,
          "load_multiplex_data",
          load_multiplex_data_stub
      )

      stub(RWRtoolkit::RWR_netstats, "calculate_basic_statistics", basic_stats_stub)
      stub(RWRtoolkit::RWR_netstats, "load_network", load_network_stub)
      stub(RWRtoolkit::RWR_netstats, "load_flist", load_flist_stub)
      stub(RWRtoolkit::RWR_netstats,
          "make_dummy_multiplex",
          make_dummy_multiplex_stub
      )

      stub(RWRtoolkit::RWR_netstats,
          "basic_statistics_multiplex",
          basic_stats_multiplex_stub
      )
      stub(RWRtoolkit::RWR_netstats,
          "overlap_many_pairwise",
          overlap_many_pairwise_stub
      )
      stub(RWRtoolkit::RWR_netstats,
          "overlap_many_vs_reference",
          overlap_many_vs_reference_stub
      )
      stub(RWRtoolkit::RWR_netstats, "overlap_pair", overlap_pair_stub)
      stub(RWRtoolkit::RWR_netstats, "calculate_tau", calculate_tau_stub)
      stub(RWRtoolkit::RWR_netstats,
          "merged_with_all_edges",
          merged_with_all_edges_stub
      )
      stub(RWRtoolkit::RWR_netstats,
          "merged_with_edgecounts",
          merged_with_edgecounts_stub
      )
      stub(RWRtoolkit::RWR_netstats, "exclusivity", exclusivity_stub)

      actual_output <- RWRtoolkit::RWR_netstats(
          data = mpo_filepath,
          flist  = NULL,
          network_1 = net1_file_path,
          network_2 = net2_file_path,
          basic_statistics = T,
          overlap_sim_multiplex_jaccard = T,
          overlap_sim_multiplex_layer = T,
          overlap_sim_layer_layer = T,
          overlap_score = T,
          calculate_tau_for_mpo = T,
          merged_with_all_edges = T,
          merged_with_edgecounts = T,
          calculate_exclusivity_for_mpo = T
      )

      expected_output <- list(
        base_stats_net1 = "basic_statistics_mock_output_net1",
        base_stats_net2 = "basic_statistics_mock_output_net2",
        base_stats_mpo = "basic_stats_multiplex_mock_output",
        overlap_sim_multiplex_jaccard =
              "overlap_many_pairwise_mock_output_jaccard",
        overlap_sim_multiplex_layer = "overlap_many_vs_reference_mock_output",
        overlap_sim_layer_layer = data.frame(
          list(
            overlap = "overlap_pair_mock_output_overlap",
            jaccard = "overlap_pair_mock_output_jaccard"
          )
        ),
        mpo_overlap_score = "overlap_many_pairwise_mock_output_mpo_overlap",
        calculated_tau = "calculate_tau_mock_output",
        merged_with_all_edges = "merged_with_all_edges_mock_output",
        merged_with_edgecounts = "merged_with_edgecounts_mock_output",
        exclusivity = "exclusivity_mock_output"
      )


      expect_equal(actual_output, expected_output)

      expect_called(load_multiplex_data_stub, 1)
      expect_called(load_network_stub, 2)
      expect_called(basic_stats_stub, 2)
      expect_called(load_flist_stub, 0)
      expect_called(make_dummy_multiplex_stub, 0)
      expect_called(basic_stats_multiplex_stub, 1)
      expect_called(overlap_many_pairwise_stub, 2)
      expect_called(overlap_many_vs_reference_stub, 1)
      expect_called(overlap_pair_stub, 2)
      expect_called(calculate_tau_stub, 1)
      expect_called(merged_with_all_edges_stub, 1)
      expect_called(merged_with_edgecounts_stub, 1)
      expect_called(exclusivity_stub, 1)

      expect_args(load_multiplex_data_stub, 1, mpo_filepath)
      expect_args(load_network_stub, 1, net1_file_path, "network_1", F)
      expect_args(load_network_stub, 2, net2_file_path, "network_2", F)
      expect_args(basic_stats_stub, 1, net1, F)
      expect_args(basic_stats_stub, 2, net2, F)
      expect_args(basic_stats_multiplex_stub, 1, nw.mpo, F)
      expect_args(overlap_many_pairwise_stub, 1, nw.mpo, "jaccard", F)
      expect_args(overlap_many_pairwise_stub, 2, nw.mpo, "overlap", F)
      expect_args(overlap_many_vs_reference_stub,
                  1,
                  nw.mpo,
                  net1,
                  "overlap",
                  F)
      expect_args(overlap_pair_stub, 1, net1, net2, "overlap", F)
      expect_args(overlap_pair_stub, 2, net1, net2, "jaccard", F)
      expect_args(calculate_tau_stub, 1, nw.mpo, net1, F)
      expect_args(merged_with_all_edges_stub, 1, nw.mpo, F)
      expect_args(merged_with_edgecounts_stub, 1, nw.mpo, F)
      expect_args(exclusivity_stub, 1, nw.mpo, F)
    })
  })
})
