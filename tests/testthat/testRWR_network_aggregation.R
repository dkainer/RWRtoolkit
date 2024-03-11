
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
    mpo_filepath <- "../testMultiplex/abc_multiplex.Rdata"
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
        actual_merged <- merged_with_all_edges(nw.mpo) #nolint

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
        actual_merged <- merged_with_edgecounts(nw.mpo)

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

