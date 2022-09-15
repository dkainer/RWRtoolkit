context("RWR_LOE Tests")
library(RWRtoolkit)

########################################################################
# Global Variables/Functions
########################################################################
setid <- c('setA', 'setA')
gene <- c('1', '2')
weight <- c(1, 1)
genesetA_2genes <- data.frame(setid, gene, weight)

setid <- c('setB', 'setB', 'setB')
gene <- c('0', '1', '3')
weight <- c(1, 1, 1)
genesetB_3genes <- data.frame(setid, gene, weight)

setid <- c('setC', 'setC')
gene <- c('0', '3')
weight <- c(1, 1)
genesetC_2genes <- data.frame(setid, gene, weight)

generate_mock_rwr_result <- function(NodeNames, Score, rank, num_in_network, num_seeds, networks, modname, seed_geneset, query_geneset, InQueryGeneset) {
  # hard coded to work with
  rank <- sequence(length(NodeNames))
  
  networks <- networks

  result <- data.frame(NodeNames, Score, rank, num_in_network, num_seeds, networks, modname, seed_geneset, query_geneset, InQueryGeneset)
  return(result)
}

generate_mock_rwr_result_tibble <- function(result, seeds) {
  result_tibble <-tibble::tibble(
    'RWRM_Results'=result, 
    'Seed_Nodes'=seeds) 
  return(result_tibble)
}

generate_mock_metrics_list <- function() {
  #TODO do I really need all this?
  summary <- data.frame("values"=rep(0,6),"measure"=c("P@SizeGeneset2","AvgPrec","AUPRC","AUROC","NDCG@SizeGeneset2","AUNDCG"))
  result <- data.frame(
    "NodeNames"=c(1,2), "Score"=rep(0.3,2), "rank"=sequence(2), "num_in_network"=rep(2,2), 
    "num_seeds"=rep(2,2), "networks"=rep("name",2), "modname"=rep("name",2), 
    "seed_geneset"=rep("name",2), "query_geneset"=rep("name",2), "InQueryGeneset"=c(0,1),
    "TP"=c(0,1),"FP"=c(0,1),"TN"=c(0,1),"FN"=c(0,1),"cum_TP"=c(0,1),"cum_TN"=c(0,1),
    "FPR"=c(0,1),"PREC"=c(0,1),REC=c(0,1),"dcg"=c(0,1),idcg=c(0,1),ndcg=c(0,1))
  
  metrics_list <-list(
    'summary'=summary, 
    'results'=result) 
  return(metrics_list)
}

#Default dummy rwr result outputs
NodeNames <- c('0', '3')
Score <- c(0.03, 0.01) # Some arbitrary score
rank <- sequence(length(NodeNames))
num_in_network <- rep(2,2)
num_seeds <- rep(2,2)
networks <- rep('m1_m2',2)
modname <- rep('default',2)
seed_geneset <- c('setA', 'setA')
query_geneset <- c('setB', 'setB')
InQueryGeneset <- c(0, 1)
Seeds <- c('1', '2')
res_dummy <- generate_mock_rwr_result(NodeNames, Score, rank, num_in_network, num_seeds, networks, modname, seed_geneset, query_geneset, InQueryGeneset)
res_tibble_dummy <- generate_mock_rwr_result_tibble(res_dummy, Seeds)

########################################################################
# Testing RWR_LOE internal functions
########################################################################
describe('RWR_LOE calc_metrics_loe', {
  it('throws a warning if there are query genes in the seed geneset', {
    expect_warning(RWRtoolkit::calc_metrics_loe(res_dummy,genesetA_2genes,genesetB_3genes))
  })
  
  it('returns a list of two objects back', {
    output <- RWRtoolkit::calc_metrics_loe(res_dummy,genesetA_2genes,genesetC_2genes)
    
    expect_equal(2,length(output))
  })
  
  it('results contain the correct amount of genes (no dupes)', {
    output <- RWRtoolkit::calc_metrics_loe(res_dummy,genesetA_2genes,genesetC_2genes)
    
    expect_equal(nrow(output$result),nrow(genesetC_2genes))
  })
  
  it('results contain the correct amount of genes (with dupes)', {
    output <- expect_warning(RWRtoolkit::calc_metrics_loe(res_dummy,genesetA_2genes,genesetB_3genes))
  
    ## Expected output's length with duplicates is because all 3 one of the query set genes is also 
    ## in the seed genes.  
    dupes <- which(genesetB_3genes$gene %in% genesetA_2genes$gene)
    expect_equal(nrow(output$result),nrow(genesetB_3genes))
  })
})

# #describe('RWR_LOE view_top_network', {
#   #TODO - not sure how to test if cytoscape opened
# #})

describe('RWR_LOE save_plots_loe', {
  it('output file of plots at expected path exists after saving', {
    test_metrics <- generate_mock_metrics_list()
    test_seed_geneset <- genesetA_2genes
    test_query_geneset <- genesetB_3genes
    test_modname <- "testmod"
    test_outdir <- "plots"

    # The save_plots_loe fxn internally calls get_file_path, so we can use the
    # same parameters to test the output file name.
    expected_outputFileName <- RWRtoolkit::get_file_path(
        "RWR-LOE",
        slug="metrics",
        modname=test_modname,
        outdir=test_outdir,
        ext='.png'
    )
    
    RWRtoolkit::save_plots_loe(metrics=test_metrics,
                               seed_geneset=test_seed_geneset,
                               query_geneset=test_query_geneset,
                               outdir=test_outdir,
                               modname=test_modname)

    expect_true(expected_outputFileName %in% list.files(recursive = TRUE))
    system(paste('rm -r plots'))
  })
})

########################################################################
# Testing RWR_LOE parameters
########################################################################
describe('RWR_LOE parameter tests', {
  it('throws an error if all parameters are missing', {
    expect_error(RWRtoolkit::RWR_LOE())
  })
  
  it('throws an error if seed_geneset parameter is missing', {
    test_data <- "../testNetworks/network_m1m2.rdata"
    
    expect_error(RWRtoolkit::RWR_LOE(data=test_data))
  })
  
  it('throws an error if Rdata file does not exist', {
    test_data <- "./fake/path/network.Rdata"
    test_seed_geneset <- "../testGenesets/testGeneset1.tsv"
    
    expect_error(invisible(RWRtoolkit::RWR_LOE(data=test_data,seed_geneset=test_seed_geneset)))
  })
  
  it('if eval mode is turned on, metrics are computed and saved', {
    test_data <- "../testNetworks/network_m1m2.rdata"
    test_seed_geneset <- "../testGenesets/testGeneset1.tsv"
    test_query_geneset <- "../testGenesets/testGeneset2.tsv"
    test_tau <- "1,1"
    test_eval <- TRUE

    #stub functions
    stub_calc_metrics_loe = mock()
    stub_save_plots_loe = mock()
    stub_write_table = mock()
    stub(RWRtoolkit::RWR_LOE, 'calc_metrics_loe', stub_calc_metrics_loe)
    stub(RWRtoolkit::RWR_LOE, 'save_plots_loe', stub_save_plots_loe)
    stub(RWRtoolkit::RWR_LOE, 'write_table', stub_write_table)
    
    invisible(RWRtoolkit::RWR_LOE(data=test_data, seed_geneset=test_seed_geneset, query_geneset=test_query_geneset, tau=test_tau, eval=test_eval))
    expect_called(stub_calc_metrics_loe, 1)
    expect_called(stub_save_plots_loe,  1)
  })
  
  it('if eval mode is turned on but there is no query set, metrics are NOT computed and saved and a warning message is displayed', {
    test_data <- "../testNetworks/network_m1m2.rdata"
    test_seed_geneset <- "../testGenesets/testGeneset1.tsv"
    test_tau <- "1,1"
    test_eval <- TRUE
    
    #stub functions
    stub_calc_metrics = mock()
    stub_save_plots2 = mock()
    stub(RWRtoolkit::RWR_LOE, 'calc_metrics', stub_calc_metrics)
    stub(RWRtoolkit::RWR_LOE, 'save_plots2', stub_save_plots2)
    
    expect_warning(RWRtoolkit::RWR_LOE(data=test_data, seed_geneset=test_seed_geneset, tau=test_tau, eval=test_eval))
    expect_called(stub_calc_metrics, 0)
    expect_called(stub_save_plots2,  0)
  })
  
})

########################################################################
# Testing RWR_LOE functionality
########################################################################
describe('RWR_LOE funtional tests', {
  it('basic parameters return expected scores and seeds', {
    expected_RWRM_Results_NodeNames <- c("0","3")
    expected_RWRM_Results_Score <- c(0.033961993,0.006687154)
    expected_RWRM_Results_rank <- c(1,2)
    expected_RWRM_Results_num_in_network <- c(2,2)
    expected_RWRM_Results_num_seeds <- c(2,2)
    expected_RWRM_Results_networks <- c("m1_m2","m1_m2")
    expected_RWRM_Results_modname <- c("default","default")
    expected_RWRM_Results_seed_geneset <- c("setA","setA")
    expected_RWRM_Results <- data.frame("NodeNames"=expected_RWRM_Results_NodeNames, "Score"=expected_RWRM_Results_Score,
                                        "rank"=expected_RWRM_Results_rank, "num_in_network"=expected_RWRM_Results_num_in_network,
                                        "num_seeds"=expected_RWRM_Results_num_seeds, "networks"=expected_RWRM_Results_networks,
                                        "modname"=expected_RWRM_Results_modname, "seed_geneset"=expected_RWRM_Results_seed_geneset)
    
    expected_seeds=c("1","2")
        
    invisible(capture.output(result <- RWRtoolkit::RWR_LOE(data="../testNetworks/network_m1m2.rdata",seed_geneset="../testGenesets/testGeneset1.tsv", tau="1,1")))
    expect_equal(result$RWRM_Results, expected_RWRM_Results, tolerance=1e-3)
    expect_equal(result$Seed_Nodes, expected_seeds)
  })
  
  it('using query geneset returns expected scores and seeds', {
    expected_RWRM_Results_NodeNames <- c("0","3")
    expected_RWRM_Results_Score <- c(0.033961993,0.006687154)
    expected_RWRM_Results_rank <- c(1,2)
    expected_RWRM_Results_num_in_network <- c(2,2)
    expected_RWRM_Results_num_seeds <- c(2,2)
    expected_RWRM_Results_networks <- c("m1_m2","m1_m2")
    expected_RWRM_Results_modname <- c("default","default")
    expected_RWRM_Results_seed_geneset <- c("setA","setA")
    expected_RWRM_Results_query_geneset <- c("setB","setB")
    expected_RWRM_Results_InQueryGeneset <- c(0,1)
    expected_RWRM_Results <- data.frame("NodeNames"=expected_RWRM_Results_NodeNames, "Score"=expected_RWRM_Results_Score,
                                        "rank"=expected_RWRM_Results_rank, "num_in_network"=expected_RWRM_Results_num_in_network,
                                        "num_seeds"=expected_RWRM_Results_num_seeds, "networks"=expected_RWRM_Results_networks,
                                        "modname"=expected_RWRM_Results_modname, "seed_geneset"=expected_RWRM_Results_seed_geneset,
                                        "query_geneset"=expected_RWRM_Results_query_geneset, "InQueryGeneset"=expected_RWRM_Results_InQueryGeneset)

    expected_seeds=c("1","2")
    
    invisible(capture.output(result <- RWRtoolkit::RWR_LOE(data="../testNetworks/network_m1m2.rdata",seed_geneset="../testGenesets/testGeneset1.tsv", query_geneset="../testGenesets/testGeneset2.tsv", tau="1,1")))
    expect_equal(result$RWRM_Results, expected_RWRM_Results, tolerance=1e-3)
    expect_equal(result$Seed_Nodes, expected_seeds)
  })
})


teardown({
  system('rm RWR-LOE_*')
}, env=parent.frame())
