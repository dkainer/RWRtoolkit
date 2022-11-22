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

      print("actual_mpo")
      print(actual_mpo)
      print("expected_mpo")
      print(expected_mpo)

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
      output_file_base <- "RWR-CV__setA1_string_interactions.Rdata_default."
      method <- "kfold"
      folds <- 3
      tau <- "1,1,1,1,1,1,1,1,1"

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
        output_directory,
        "--tau",
        tau
      )
      system(script)
      RWRtoolkit::RWR_CV(
        data = network_filepath,
        geneset = geneset_filepath,
        
      )

      output_files_script <- list.files(output_directory)
      actual_file_count_script <- length(output_files)
      actual_metrics_nrow_script <- nrow(
        read.table(
          paste(output_directory, output_file_base, "metrics", ".tsv", sep = "")
        )
      )
      actual_summary_nrow_script <- nrow(
        read.table(
          paste(output_directory, output_file_base, "summary", ".tsv", sep = "")
        )
      )
      actual_fullranks_nrow_script <- nrow(
        read.table(
          paste(output_directory,
            output_file_base,
            "fullranks",
            ".tsv",
            sep = ""
          )
        )
      )
      actual_medianranks_nrow_script <- nrow(
        read.table(
          paste(output_directory,
            output_file_base,
            "medianranks",
            ".tsv",
            sep = ""
          )
        )
      )



      expect_equal(actual_file_count, expected_file_count)
      expect_equal(actual_metrics_nrow, expected_metrics_nrow)
      expect_equal(actual_summary_nrow, expected_summary_nrow)
      expect_equal(actual_fullranks_nrow, expected_fullranks_nrow)
      expect_equal(actual_medianranks_nrow, expected_medianranks_nrow)
    })

     teardown({
       system("rm -rf ./runCVOutput")
     })
  })
  describe("run_loe", {
    it("Runs rwr loe from the command line", {

      run_cv_filepath <- "../../inst/scripts/run_loe.R"
      network_filepath <- "../testSTRINGDB/string_interactions.Rdata"
      geneset_filepath <- "../testSTRINGDB/geneset1.tsv"
      output_directory <- "./rwr_loe_output/"
      output_file_base <- "RWR-LOE_setA1_default"

      tau <- "1,1,1,1,1,1,1,1,1"

      script <- paste(
        "Rscript",
        run_cv_filepath,
        "--data",
        network_filepath,
        "--geneset",
        geneset_filepath,
        "--outdir",
        output_directory,
        "--tau",
        tau
      )

      system(script)

      output_files <- list.files(output_directory)

      expect_equal(length(output_files), 1)
      expect_true(stringr::str_detect(string = output_files[1],
                  pattern = output_file_base))
    })

     teardown({
       system("rm -rf ./rwr_loe_output")
     })
  })
  # describe("run_netscore")
  # describe("run_shortestpaths")

})
