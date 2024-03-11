## These tests presume that the RWRtoolkit package has already
## been installed.
## ANY UPDATES to the code require a fresh installation for these
## to take affect. 


describe("RWR inst/scripts", {
  describe("rwr_make_multiplex", {
    it("test run rwr_make_multiplex", {
      make_multiplex_script <- "../../inst/scripts/run_make_multiplex.R"
      flist <- "../testFlists/abc_flist.txt"
      delta <- 0.5
      outpath <- "./test_make_multiplex.Rdata"
      expected_mpo_path <- "../testMultiplex/abc_multiplex.Rdata"

      script <- paste("Rscript",
                    make_multiplex_script,
                    "--flist",
                    flist,
                    "--delta",
                    delta,
                    "--out",
                    outpath)
      system(script)
      
      load(outpath)
      actual_mpo <- nw.mpo
      actual_adj <- nw.adj
      actual_adjnorm <- nw.adjnorm


      load(expected_mpo_path)
      expected_mpo <- nw.mpo
      expected_adj <- nw.adj
      expected_adjnorm <- nw.adjnorm

      expect_equal(actual_mpo$Number_of_Layers, expected_mpo$Number_of_Layers)
      expect_equal(actual_mpo$Number_of_Nodes_Multiplex,
                  expected_mpo$Number_of_Nodes_Multiplex)
      expect_equal(actual_mpo$Pool_of_Nodes, expected_mpo$Pool_of_Nodes)
      expect_equal(names(V(actual_mpo$layer1)),
                  names(V(expected_mpo$layer1)))
      expect_equal(names(V(actual_mpo$layer2)),
                  names(V(expected_mpo$layer2)))
      expect_equal(names(V(actual_mpo$layer3)),
                  names(V(expected_mpo$layer3)))
      expect_equal(names(E(actual_mpo$layer1)),
                  names(E(expected_mpo$layer1)))
      expect_equal(names(E(actual_mpo$layer2)),
                  names(E(expected_mpo$layer2)))
      expect_equal(names(E(actual_mpo$layer3)),
                  names(E(expected_mpo$layer3)))
      expect_equal(actual_adj, expected_adj)
      expect_equal(actual_adjnorm, expected_adjnorm)
    })

    teardown({
        system("rm test_make_multiplex.Rdata")
      }
    )
  })

  describe("run_cv", {
    it("runs rwr cv from the inst", {
      run_cv_filepath <- "../../inst/scripts/run_cv.R"
      network_filepath <- "../testSTRINGDB/string_interactions.Rdata"
      geneset_filepath <- "../testSTRINGDB/genesetLong.tsv"
      output_directory_script <- "./runCVOutput_script/"
      output_directory_rfunc <- "./runCVOutput_rfunc/"
      output_file_base <- "RWR-CV__setA1_string_interactions_default."
      method <- "kfold"
      folds <- 3
      tau <- "1,1,1,1,1,1,1,1,1"

      # Run script
      script <- paste(
        "Rscript",
        run_cv_filepath,
        "--data",
        network_filepath,
        "--geneset",
        geneset_filepath,
        "--method",
        method,
        "--folds",
        folds,
        "--outdir",
        output_directory_script,
        "--tau",
        tau
      )

      system(script)

      #Run regular
      RWRtoolkit::RWR_CV(
        data = network_filepath,
        geneset_path = geneset_filepath,
        method = method,
        folds = folds,
        outdir = output_directory_rfunc,
        tau = tau,
        write_to_file = T
      )

      output_files_script <- list.files(output_directory_script)
      actual_file_count_script <- length(output_files_script)

      actual_metrics_nrow_script <- nrow(
        read.table(
          paste(output_directory_script,
                output_file_base,
                "metrics",
                ".tsv",
                sep = "")
        )
      )
      actual_summary_nrow_script <- nrow(
        read.table(
          paste(output_directory_script,
                output_file_base,
                "summary",
                ".tsv",
                sep = "")
        )
      )
      actual_fullranks_nrow_script <- nrow(
        read.table(
          paste(output_directory_script,
            output_file_base,
            "fullranks",
            ".tsv",
            sep = ""
          )
        )
      )
      actual_mean_ranks_nrow_script <- nrow(
        read.table(
          paste(output_directory_script,
            output_file_base,
            "meanranks",
            ".tsv",
            sep = ""
          )
        )
      )

      output_files_rfunc <- list.files(output_directory_rfunc)
      actual_file_count_rfunc <- length(output_files_rfunc)
      actual_metrics_nrow_rfunc <- nrow(
        read.table(
          paste(output_directory_rfunc,
                output_file_base,
                "metrics",
                ".tsv",
                sep = "")
        )
      )
      actual_summary_nrow_rfunc <- nrow(
        read.table(
          paste(output_directory_rfunc,
                output_file_base,
                "summary",
                ".tsv",
                sep = "")
        )
      )
      actual_fullranks_nrow_rfunc <- nrow(
        read.table(
          paste(output_directory_rfunc,
            output_file_base,
            "fullranks",
            ".tsv",
            sep = ""
          )
        )
      )
      actual_meanranks_nrow_rfunc <- nrow(
        read.table(
          paste(output_directory_rfunc,
            output_file_base,
            "meanranks",
            ".tsv",
            sep = ""
          )
        )
      )

      expect_equal(output_files_script, output_files_rfunc)
      expect_equal(actual_file_count_script, actual_file_count_rfunc)
      expect_equal(actual_metrics_nrow_script, actual_metrics_nrow_rfunc)
      expect_equal(actual_summary_nrow_script, actual_summary_nrow_rfunc)
      expect_equal(actual_fullranks_nrow_script, actual_fullranks_nrow_rfunc)
      expect_equal(actual_mean_ranks_nrow_script, actual_meanranks_nrow_rfunc)
    })

     teardown({
       system("rm -rf ./runCVOutput_*")
     })
  })

  describe("run_loe", {
    it("Runs rwr loe from the command line", {

      run_loe_filepath <- "../../inst/scripts/run_loe.R"
      network_filepath <- "../testSTRINGDB/string_interactions.Rdata"
      geneset_filepath <- "../testSTRINGDB/geneset1.tsv"
      output_directory_script <- "./rwr_loe_output_script/"
      output_file_base <- "RWR-LOE_setA1_default"
      tau <- "1,1,1,1,1,1,1,1,1"
      script <- paste(
        "Rscript",
        run_loe_filepath,
        "--data",
        network_filepath,
        "--seed_geneset",
        geneset_filepath,
        "--outdir",
        output_directory_script,
        "--tau",
        tau
      )

      system(script)

      actual_output_files <- list.files(output_directory_script)

      expected_output <- paste(output_file_base,'ranks.tsv', sep='.')

      expect_equal(actual_output_files, expected_output)
    })

     teardown({
       system("rm -rf ./rwr_loe_output_*")
     })
  })

  describe("run_shortestpaths", {
    it("Runs rwr shortest paths from the command line", {
      run_shortestpaths_filepath <- "../../inst/scripts/run_shortestpaths.R"
      network_filepath <- "../testSTRINGDB/string_interactions.Rdata"
      source_geneset <- "../testSTRINGDB/geneset1.tsv"
      target_geneset <- "../testSTRINGDB/geneset2.tsv"
      output_directory_script <- "./rwr_shortestpaths_output_script/"
      output_directory_func <- "./rwr_shortestpaths_output_func/"
      output_file_base <- "RWR-shortestpaths_setA1_default"

      tau <- "1,1,1,1,1,1,1,1,1"

      script <- paste(
        "Rscript",
        run_shortestpaths_filepath,
        "--data",
        network_filepath,
        "--source_geneset",
        source_geneset,
        "--target_geneset",
        target_geneset,
        "--outdir",
        output_directory_script
      )
      system(script)


      RWRtoolkit::RWR_ShortestPaths(
        data = network_filepath, 
        source_geneset =  source_geneset, 
        target_geneset = target_geneset, 
        outdir = output_directory_func, 
        write_to_file =  T
      )

      output_files_script <- list.files(output_directory_script)
      output_files_func <- list.files(output_directory_func)

      actual_shortpaths_ranks_script <- nrow(
        read.table(
          paste(output_directory_script,
            output_files_func[1],
            sep = ""
          )
        )
      )


      actual_shortpaths_ranks_func <- nrow(
        read.table(
          paste(output_directory_func,
            output_files_func[1],
            sep = ""
          )
        )
      )

      expect_equal(length(output_files_script), length(output_files_func))
      expect_equal(actual_shortpaths_ranks_script,
                  actual_shortpaths_ranks_func)
    })

    teardown({
      system("rm -rf ./rwr_shortestpaths_output_*")
    })
  })

  describe("RWR Netstats command line", {
    it("Runs RWR netstats from the command line", {
      run_netstats_filepath <- "../../inst/scripts/run_netstats.R"
      network_filepath <- "../testSTRINGDB/string_interactions.Rdata"
      outdir_path <- "./test_rwr_netstats_output_directory_inst_script"
      network1_filepath <- "../testSTRINGDB/netstat/combined_score-random-gold.tsv"
      network2_filepath <- "../testSTRINGDB/netstat/combined_score-random-test.tsv"
      basic_statistics <- TRUE
      scoring_metric <- "both"
      pairwise_between_mpo_layer <- TRUE
      multiplex_layers_to_refnet <- TRUE
      net_to_net_similarity <- TRUE
      calculate_tau_for_mpo <- TRUE
      calculate_exclusivity_for_mpo <- TRUE
      verbose <- TRUE

      script <- paste(
        "Rscript",
        run_netstats_filepath,
        "--data",
        network_filepath,
        "--network_1",
        network1_filepath,
        "--network_2",
        network2_filepath,
        "--outdir",
        outdir_path,
        "--basic_statistics",
        "--scoring_metric",
        scoring_metric,
        "--pairwise_between_mpo_layer",
        "--multiplex_layers_to_refnet",
        "--net_to_net_similarity",
        "--calculate_tau_for_mpo",
        "--calculate_exclusivity_for_mpo")

      system(script)
      actual_output_files <- list.files(outdir_path)


      expected_output_files <- c(
        "base_stats.tsv",
        "pairwise_between_mpo_layer_jaccard.tsv",
        "pairwise_between_mpo_layer_overlap.tsv",
        "multiplex_layers_to_refnet_jaccard.tsv",
        "multiplex_layers_to_refnet_overlap.tsv",
        "net_to_net_similarity.tsv",
        "calculated_tau.tsv",
        "exclusivity.tsv"
         )

      expect_setequal(actual_output_files, expected_output_files)
    })
    
    teardown({
      system("rm -rf ./test_rwr_netstats_output_directory_inst_script")
    })
  })

  describe("RWR Network Aggregation command line", {
    it("Runs RWR network aggregator from the command line", {
      run_aggregator_filepath <- "../../inst/scripts/run_network_aggregation.R"
      network_filepath <- "../testSTRINGDB/string_interactions.Rdata"
      outdir_path <- "./test_rwr_aggregator_output_directory_inst_script"

      verbose <- TRUE
 
      script <- paste(
        "Rscript",
        run_aggregator_filepath,
        "--data",
        network_filepath,
        "--merged_with_all_edges",
        "--merged_with_edgecounts",
        "--outdir",
        outdir_path
    )
      system(script)
      actual_output_files <- list.files(outdir_path)

      expected_output_files <- c(
     "merged_with_all_edges.tsv",
     "merged_with_edgecounts.tsv"
    ) 

      expect_setequal(actual_output_files, expected_output_files)
    })
    
    teardown({ 
      system("rm -rf ./test_rwr_aggregator_output_directory_inst_script")
    })
  })



})
