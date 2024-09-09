

# Function to calculate correlation between corresponding columns of two matrices using multiprocessing
correlation_matrices <- function(mat1, mat2, hpc=FALSE, outputfile='correlation_vector.tsv') {
  
  # Define a function to calculate the correlation for a single column pair
  calc_correlation <- function(i) {
        print("IN CORRELATION", flush=T)
        print(i, flush=T)
	print("WEEEEEEE", flush=T)
    	print(mat1[, i])
	print(typeof(mat1[,i]))
	cor_val <- cor(mat1[, i], mat2[, i]) 
	if (is.na(cor_val)) return(0)

	return(cor_val)	
  }
   
  correlations <- NULL
  col_indices <- as.numeric(1:ncol(mat1))
  check_matrices(mat1, mat2)
  if (hpc){
	output <- task.pull(col_indices, calc_correlation)
	if(comm.rank() == 0){
 # Define the number of cores to use
          correlations <- lapply(col_indices, calc_correlation) 
         # print("OUTPUT from correlations", Flush=T) 
	  # print(output)  
	   #correlations  <- unlist(output)
	   #saveRDS(correlations, 'correlations.rds') 
	   #names(correlations) <- colnames(mat1)
	   correlations <- unlist(correlations)
	   names(correlations) <- colnames(mat1) 
           write.table(data.frame(correlations), file=outputfile, sep='\t')
        }	
	  return()
  } else {
	  # Define the number of cores to use
	  num_cores <- detectCores() - 1
	  cl <- makeCluster(num_cores)

	  # Use mclapply to parallelize the correlation calculation across columns
	  correlations <- mclapply(col_indices, calc_correlation, mc.cores = num_cores)
	  stopCluster(cl)
  }
  # Convert the list of correlations to a vector
  correlations <- unlist(correlations)

  names(correlations) <- colnames(mat1)

  write.table(data.frame(correlations), file=outputfile, sep='\t') 
  return(correlations)
}





extract_and_rename_score_column <- function(df_list){
  lapply(df_list, function(df) {
    seed_value <- unique(df$seed)
    df <- df[, c('NodeNames', 'Score')]
    colnames(df)[which(colnames(df) == "Score")] <- paste0(seed_value)
    df
  })
}

combine_ranks_to_embed_matrix <- function(res){
  extracted_score_and_nodes <- extract_and_rename_score_column(res)
  reduced <- Reduce(function(x, y) merge(x, y, by ='NodeNames', all=T), extracted_score_and_nodes)
  reduced
}


create_rwr_embedding_matrix <- function(data, knockouts, output_file, write_to_file, tau, threads, hpc, verbose){
  data_list <- load_multiplex_data(data, knockouts)
  nw_mpo <- data_list$nw.mpo
  nw_adjnorm <- data_list$nw.adjnorm

  tau <- get_or_set_tau(nw_mpo, tau)

  # Geneset is required.
  geneset_list <- load_geneset(nw_mpo$Pool_of_Nodes, nw_mpo, verbose = verbose)
  geneset <- geneset_list$geneset
  extras <- geneset_list$extras

  # check and update data if necessary 
  folds <- nrow(geneset)
  updated_data_list <- update_folds_by_method(geneset, 'singletons', folds)
  folds <- updated_data_list$folds
  geneset <- updated_data_list$geneset
  chunks <- updated_data_list$chunks
  method <- updated_data_list$method


  ############# run RWR  #####################################################
  res <- RWR(geneset,
             nw_adjnorm,
             nw_mpo,
             method,
             folds,
             chunks,
             restart,
             tau,
             'embedding_mat',
             threads,
             hpc,
             verbose)

  ############# post process results  ########################################
  base_embedded_matrix <- combine_ranks_to_embed_matrix(res)
  
  if (write_to_file) write_table(base_embedded_matrix, path=output_file)
  
  return(base_embedded_matrix) 
}


########################################################################
# Main Function
########################################################################

#' RWR Embedding Comparator
#'
#' `RWR_embedding_comparator` performs a pairwise comparison across vector
#' vectors within the vector embedding matrix defined by random walks 
#' between two matrices: 1, a base matrix acting and a perturbed matrix. 
#' The perturbed embedding matrix comes from knocking out nodes within the
#' multiplex network. 
#'
#' @param data The path to the .Rdata file containing your multiplexed
#'                  functional networks. This file is produced by
#'                  RWR_make_multiplex. Default NULL
#' @param restart Set the restart parameter \[0,1). Higher value means the
#'                walker will jump back to seed node more often. Default 0.7
#' @param tau Comma-separated list of values between that MUST add up to the
#'            number of network layers in the .Rdata file. One value per 
#'            network layer that determines the probability that the random
#'            walker will restart in that layer. e.g. if there are three layers
#'            (A,B,C) in your multiplex network, then --tau '0.2,1.3,1.5' will
#'            mean that layer A is less likely to be walked on after a restart
#'            than layers B or C.  Default 1.0
#' @param outdir Path to the output directory. Both 'fullranks' and 
#'                   'medianranks' will be saved with auto-generated filenames.
#'                   Can be overridden by specifically setting 'out_full_ranks'
#'                   and 'out_mean_ranks' parameters. No defined path will 
#'                   output within the same directory from which the original
#'                   code was run.
#'                   Default NULL
#' @param threads Specify the number of threads to use. Default for your system
#'                is all cores - 1.
#' @param verbose Verbose mode. Default FALSE
#' @return Returns a list of four data tables: fullranks, medianranks, metrics,
#'         and summary.
#' @examples
#'
#' # An example of Running RWR CV
#' # Loads a 10 layer multiplex and does not write to file:
#' extdata.dir <- system.file("example_data", package = "RWRtoolkit")
#' multiplex_object_filepath <- paste(extdata.dir,
#'                                   "/string_interactions.Rdata",
#'                                   sep = "")
#' geneset_filepath <- paste(extdata.dir, "/geneset1.tsv", sep = "")
#' outdir <- "./rwr_cv"
#'
#'
#' cv_examples <- RWR_CV(
#'   data = multiplex_object_filepath,
#'   tau = "1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0",
#'   geneset_path = geneset_filepath,
#'   outdir = outdir,
#'   method = "kfold",
#'   folds = 3
#' )
#'
#' # An example of Running RWR CV with non-default method and writing to file
#' # Loads a 10 layer multiplex and does not write to file:
#' cv_examples <- RWR_CV(
#'   data = multiplex_object_filepath,
#'   tau = "1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0",
#'   geneset_path = geneset_filepath,
#'   outdir = outdir,
#'   method = "singletons",
#'   write_to_file = TRUE
#' )
#'
#' @export
RWR_embedding_comparator <- function(data = NULL,
                   knockouts_to_simulate = c(),
                   restart = 0.7,
                   tau = 1.0,
                   outdir = './',
                   threads = 1,
                   hpc = FALSE,
                   verbose = FALSE,
                   write_to_file = FALSE) {
 
    ############# Initialize  ##################################################
    base_output_file <- paste(out_dir, 'base_embedding_matrix.tsv', sep='/')
    perturbed_output_file <- paste(out_dir, 'perturbed_embedding_matrix.tsv', sep='/')

    base_embedding <- create_rwr_embedding_matrix(
        data, 
        c(), 
        base_output_file, 
        write_to_file, 
        tau, 
        threads, 
        hpc,
        verbose)
    perturbed_embedding <- create_rwr_embedding_matrix(
        data, 
        knockouts_to_simulate, 
        perturbed_output_file, 
        write_to_file, 
        tau, 
        threads,
        hpc, 
        verbose)

    comparison_vector <- correlation_matrices(base_embedded, perturbed_embedding, hpc, outputfile)
}
