########################################################################
# Perform K-fold Cross Validation on a gene set using RWR to find the RWR rank
# of the left-out genes:
# - Input: Pre-computed multiplex network and a geneset
# - Output: Table with the ranking of each gene in the gene set when left out,
#           along with AUPRC and AUROC curves
# Copyright (C) 2022  David Kainer
#
# This file is part of RWRtoolkit.
#
# RWRtoolkit is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# RWRtoolkit is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# RWRtoolkit. If not, see <https://www.gnu.org/licenses/>.
########################################################################

#' @importFrom dplyr %>%
#' @importFrom foreach %dopar%

########################################################################
# Internal Functions
########################################################################

update_folds_by_method <- function(
  geneset,
  method,
  num_folds,
  verbose = FALSE) {
  chunks <- NULL
  if (method %in% c("loo", "singletons")) {
    folds <- nrow(geneset)
    cat("\nCross-validation method is ",
      method,
      ", i.e, folds=",
      folds,
      "\n",
      sep = "",
      file = stderr()
    )
  } else if (method == "kfold") {
    if ((nrow(geneset) / num_folds) < 1) {
      folds <- nrow(geneset)
      method <- "loo"
      base_warning <- "\nWARNING: Cross-validation method was set to k-fold, but `k` is too large for this geneset: " # nolint
      warning_message <- paste(base_warning,
        num_folds,
        "\n",
        "\nCross-validation method is now set to LOO, i.e, folds=",
        folds,
        "\n",
        sep = ""
      )
      warning(warning_message)
    } else {
      folds <- num_folds
      # Randomly shuffle the geneset so that k-folds below is unbiased.
      set.seed(42)
      samples <- sample(nrow(geneset))
      geneset <- geneset[samples, ]
      message("Gene set was randomly shuffled")

      if (verbose) {
        message("Gene set:")
        print(geneset)
      }

      # we shuffle the geneset to break up and pre-ordering
      chunks <- chunk(sample(geneset$gene), folds)
      cat("\nCross-validation method is k-fold; folds=",
        folds,
        "\n",
        sep = "",
        file = stderr()
      )
    }
  } else {
    stop("ERROR: CV method not recognised. must be one of loo, singletons or kfold. Stopping") # nolint message
  }
  return(list(folds=folds, geneset=geneset, chunks=chunks, method=method))
}

extract_lo_and_seed_genes_cv <- function(geneset, method, r, chunks = NULL) {
  if (method == "singletons") {
    seed_genes <- geneset %>%
      dplyr::slice(r) %>%
      dplyr::pull(gene) # nolint       Get one seed gene.  
    leftout <- geneset %>%
      dplyr::filter(gene != seed_genes) %>%     #nolint
      dplyr::pull(gene)
  } else if (method == "loo") {
    leftout <- geneset %>%
      dplyr::slice(r) %>%
      dplyr::pull(gene) # nolint        leave out one gene.
    seed_genes <- geneset %>%
      dplyr::filter(!gene %in% leftout) %>%  #nolint
      dplyr::pull(gene)
  } else {
    leftout <- chunks[[r]] # get the r'th fold for CV
    seed_genes <- geneset %>%
      dplyr::filter(!gene %in% leftout) %>%  #nolint
      dplyr::pull(gene)
  }

  list(leftout, seed_genes)
}

create_rankings_cv <- function(rwr,
                               networks,
                               r,
                               name,
                               geneset,
                               method,
                               seed_genes,
                               leftout,
                               num_nodes_in_multiplex) {
  ## turn Scores into Rankings
  # for any gene that scored 0, give it worst possible rank
  # highest score gets lowest rank
  rwr$RWRM_Results <- rwr$RWRM_Results %>% # nolint output name from rwrmh
    dplyr::mutate(
      rank = dplyr::row_number(-Score),
      InValset = as.numeric(rwr$RWRM_Results$NodeNames %in% leftout)
    ) %>%
    dplyr::mutate(rank = dplyr::if_else(Score == 0,
      true = as.integer(num_nodes_in_multiplex),
      false = rank
    ))

  num_in_network <- nrow(geneset)
  num_seeds <- length(seed_genes)
  num_leftout <- length(leftout)
  networks <- networks
  fold <- r
  modname <- name
  geneset <- geneset$setid[1]
  seed <- if (method == "singletons") {
    seed_genes
  } else if (method == "loo") {
    "AllBut1"
  } else {
    "many"
  }
  leftout <- if (method == "singletons") {
    "AllBut1"
  } else if (method == "loo") {
    leftout
  } else {
    "many"
  }

  rwr$RWRM_Results <- rwr$RWRM_Results %>% # nolint rwrmh output
    dplyr::mutate(
      num_in_network = num_in_network,
      num_seeds = num_seeds,
      num_leftout = num_leftout,
      networks = networks,
      fold = fold,
      modname = modname,
      geneset = geneset,
      seed = seed,
      leftout = leftout,
      method = method
    )

  rwr$RWRM_Results
}

# For 'loo' and 'kfold', each gene is leftout once, so gets ranked once.
RWR <- function(geneset,
                adjnorm,
                mpo,
                method,
                num_folds,
                chunks,
                restart = 0.7,
                tau = c(1, 1),
                name = "default",
                threads = 1,
                verbose = FALSE) {
  # This is the name of the combined networks.
  networks <- paste(names(mpo)[1:mpo$Number_of_Layers], collapse = "_")

  doParallel::registerDoParallel(cores = threads)
  res <- foreach::foreach(r = 1:num_folds) %dopar% {
    lo_seed_genelist <- extract_lo_and_seed_genes_cv(
      geneset,
      method,
      r,
      chunks
    )
    leftout <- lo_seed_genelist[[1]]
    seed_genes <- lo_seed_genelist[[2]]

    ### run RWR on a fold
    rwr <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(
      x = adjnorm,
      MultiplexObject = mpo,
      Seeds = seed_genes,
      r = restart,
      tau = tau
    )
    ranking_results <- create_rankings_cv(
      rwr,
      networks,
      r,
      name,
      geneset,
      method,
      seed_genes,
      leftout,
      mpo$Number_of_Nodes_Multiplex
    )
    ranking_results
  }
  doParallel::stopImplicitCluster()

  res
}

calc_metrics_cv <- function(res_combined, res_avg) {
  ### Singletons: each relevant gene (R) is ranked R-1 times.
  ###             - There are R folds, with R-1 left out per fold
  ### LOO:        each relevant gene (R) is ranked once.
  ###             - There are R folds, with 1 left out per fold
  ### KFOLD:      each relevant gene (R) is ranked once.
  ###             - There are K folds, with expectation of R/K
  ### .             left out per fold (but not necessarily!!)

  ### For KFOLD and Singletons::
  ### R-Precision is most appropriate metric for ranked retrieval where the
  ### number of relevant results (R) is known. It is simply RPREC = r/R,
  ### where r = number of hits, and R = number of leftout geneset genes.
  ### It is the same as Precision@R.
  ###
  ### For LOO this doesn't really work because R = 1, so we are looking
  ###  at Precision @ 1 (lol).

  # PR example: 5 genes left out need to be detected from 10 rankings
  # ranking:  1    2    3    4   5   6   7    8    9     10
  # truth:    1    0    1    0   0   1   0    0    1     1
  # recall:   0.2  0.2  0.4  0.4 0.4 0.6 0.6  0.6  0.8   1.0
  # Prec:     1    0.5  0.67 0.5 0.4 0.5 0.43 0.38 0.44  0.5
  # Avg Prec: 1 +  0  + 0.67 + 0 +0 +0.5 + 0 + 0 + 0.44 +0.5) / 5 = 0.62

  # n() no longer works if dplyr not attached #4062
  # https://github.com/tidyverse/dplyr/issues/4062
    # https://github.com/tidyverse/dplyr/blob/master/NEWS.md#breaking-changes
    if(res_combined$method[1] %in% c("kfold","singletons")) {  # cannot use opt$method here since that can be overridden in RWR

        ### calculate ROC/PRC per-fold
        res_combined <- res_combined %>%
            dplyr::mutate(fold = as.factor(fold)) %>%
            dplyr::group_by(fold) %>%
            dplyr::group_modify(~ calc_ROCPRC(df=.x, scorescol = "rank", labelscol = "InValset"))

        # 1. precision @ numleftout (aka R-PREC)
        output <- res_combined %>%
            dplyr::group_by(fold) %>%
            dplyr::filter(rank==num_leftout) %>%
            dplyr::summarise(value = PREC, measure="P@NumLeftOut")

        # 2. Avg PRC (uses Average Precision, not interpolated precision)
        output <- rbind(output, res_combined %>%
                            dplyr::group_by(fold) %>%
                            dplyr::filter(TP==1) %>%
                            dplyr::summarise(value = sum(PREC)/sum(InValset), measure="AvgPrec") )
        # 3. AUPRC (provide both the AUPRC per fold, and the AUPRC expected by an unskilled model)
        output <- rbind(output, res_combined %>%
                            dplyr::group_by(fold) %>%
                            dplyr::summarise(value = area_under_curve(REC, PREC, method="trapezoid", ties="max"), measure="AUPRC"))
        
        output <- rbind(output, res_combined %>%
                          dplyr::group_by(fold) %>% 
                          dplyr::summarise(value = dplyr::first(num_leftout)/dplyr::n(), measure="ExpectedAUPRC"))
        
        # 4. AUROC
        output <- rbind(output, res_combined %>%
            dplyr::group_by(fold) %>%
            dplyr::summarise(value = sum(REC)/dplyr::n(), measure="AUROC"))

        ### get metrics based on reranking of mean ranks
        res_avg <- calc_ROCPRC(res_avg, scorescol = "rerank", labelscol = "InValset")
        # 5. AvgPrec of mean ranks (uses Average Precision, not interpolated avg precision)
        output <- rbind(output, res_avg %>%
                            dplyr::filter(TP==1) %>%
                            dplyr::summarise(fold="meanrank", value = sum(PREC)/sum(InValset), measure="AvgPrec") )
        # 6. AUPRC of mean ranks
        output <- rbind(output, res_avg %>%
                            dplyr::summarise(fold="meanrank", value = area_under_curve(REC, PREC, method="trapezoid", ties="max") , measure="AUPRC") )
        # 7. AUROC of mean ranks
        output <- rbind(output, res_avg %>%
                            dplyr::summarise(fold="meanrank", value = sum(REC)/dplyr::n(), measure="AUROC") )
    }


    if(res_combined$method[1] == "loo") {
        ### DON'T get metrics per fold since only 1 gene is left out per fold!!

        # 1. rank of each leftout gene
        output <- res_combined %>%
            dplyr::group_by(fold) %>%
            dplyr::filter(InValset==1) %>%
            dplyr::summarise(value = rank, measure="leftout.rank")
        print("Rank of each left out gene: ")
        print(output)
        # 2. number of left-out genes with rank < NumGeneset
        output <- rbind(output, res_combined %>%
                            dplyr::filter(InValset==1 & rank <= sum(InValset)) %>%
                            dplyr::summarise(fold="all", value=dplyr::n(), measure="RankedBelowNumgeneset"))

        ### get metrics based on mean ranks
        res_avg <- calc_ROCPRC(res_avg, scorescol = "rerank", labelscol = "InValset")
        # 3. AvgPrec of mean ranks (uses Average Precision, not interpolated avg precision)
        output <- rbind(output, res_avg %>%
                            dplyr::filter(TP==1) %>%
                            dplyr::summarise(fold="meanrank", value = sum(PREC)/sum(InValset), measure="AvgPrec") )
        # 4. AUPRC of mean ranks (uses trapezoid method)
        output <- rbind(output, res_avg %>%
                            dplyr::summarise(fold="meanrank", value = area_under_curve(REC, PREC, method="trapezoid", ties="max"), measure="AUPRC"))
        # 5. AUROC of mean ranks
        output <- rbind(output, res_avg %>%
                            dplyr::summarise(fold="meanrank", value = sum(REC)/dplyr::n(), measure="AUROC") )
    }

    output$geneset <- res_combined$geneset[1]
    
    return(list(summary = output, res_combined = res_combined, res_avg = res_avg))

}

calculate_max_precision <- function(pr, metrics) {
  maxprec <- foreach::foreach(
    f = pr$fold@.Data,
    r = pr$recall,
    .combine = c
  ) %do% {
    prec <- metrics$res_combined %>%
      dplyr::filter(fold == f, REC >= r) %>%
      dplyr::pull(PREC)

    if (length(prec) > 0) {
      max(prec)
    } else {
      return(0)
    }
  }
  maxprec
}

save_plots_cv <- function(metrics, geneset, folds, dataPath, modname, outdirPath)
{
    message("Generating plots...\n")
    ggplot2::theme_set(ggplot2::theme_classic())
    
    thestats <- metrics$summary %>% 
        dplyr::filter(fold!="meanrank") %>% 
        dplyr::group_by(measure) %>% 
        dplyr::summarise(mean = round(mean(value),3)) %>% 
        tibble::column_to_rownames("measure")
    
    themeans <- metrics$summary %>% 
      dplyr::filter(fold=="meanrank") %>% 
      dplyr::mutate(value = round(value,3)) %>% 
        tibble::column_to_rownames("measure") %>% 
      dplyr::select(value)
      
    #updateFoldsByMethod_cv
    
    ####### KFOLD (few folds) ########

    if(metrics$res_combined$method[1] == "kfold")
    {
        ### ROC
        p1 <- ggplot2::ggplot() + 
            ggplot2::geom_line(data=metrics$res_combined %>% dplyr::group_by(FPR) %>% dplyr::summarise(mean = mean(REC)), 
                      ggplot2::aes(x=FPR, y=mean, col="mean of folds")) +
            ggplot2::geom_abline(ggplot2::aes(intercept=0, slope=1), linetype="dashed", col="darkgrey") +
            #ggplot2::annotate("text", label=paste0("mean AUROC = ",thestats["AUROC",]), x = 0.75, y = 0.01, size = 4, colour = "darkorange3") + 
            ggplot2::xlab("FPR") +
            ggplot2::ylab("TPR (Recall)") + 
            ggplot2::scale_color_manual(name = "ROC",
                               breaks = c("mean of folds"),
                               values = c("mean of folds" = "darkorange") ) +
            ggplot2::labs(title=paste0(metrics$res_combined$geneset[1]), 
                          caption = paste0("mean of folds AUROC = ", thestats["AUROC",])) +
            ggplot2::theme(legend.position="top", legend.box = "horizontal")
        
        ### PRC
        # p2 <- ggplot2::ggplot() + 
        #    # ggplot2::geom_line(data=metrics$res_combined, 
        #   #            ggplot2::aes(x=REC, y=PREC, group=fold, col="folds"), alpha=0.5) + 
        #     ggplot2::geom_line(data=metrics$res_combined %>% dplyr::group_by(REC) %>% dplyr::summarise(mean = mean(PREC)), 
        #                       ggplot2::aes(x=round(REC,2), y=mean, col="mean of folds")) +
        #    # ggplot2::geom_hline(yintercept=sum(metrics$res_avg$InValset)/nrow(metrics$res_avg), linetype="dotted", colour="darkred") +
        #     ggplot2::geom_hline(yintercept=thestats["ExpectedAUPRC",], linetype="dotted", colour="darkorange") +
        #     ggplot2::scale_color_manual(name = "PRC",
        #                        breaks = c("mean of folds"),
        #                        values = c("mean of folds" = "darkorange") ) +
        #     ggplot2::labs(title=paste0(metrics$res_combined$geneset[1]),caption = paste0("mean of folds AUPRC = ", thestats["AUPRC",])) +
        #     ggplot2::theme(legend.position="top", legend.box = "horizontal") +
        #     ggplot2::expand_limits(x=c(0,1), y=c(0, 1))
        
        ### Interpolated PRC
        # for each fold, interpolating the PR values at specific recall (x-axis) positions
        pr          <- expand.grid(fold=as.factor(seq(1,folds,by=1)), recall = seq(0,1,by=0.1))
        pr$maxprec  <- foreach::foreach(f=pr$fold, r = pr$recall, .combine = c) %do%
            {
                prec <- metrics$res_combined  %>% dplyr::filter(fold==f, REC >= r) %>% dplyr::pull(PREC)
                if(length(prec)>0)
                    max(prec)
                else
                    return(0)
            }
        # make it into a step curve for each fold (here we can show mean of folds line properly)
        lag1 <- pr %>% dplyr::group_by(fold) %>% dplyr::mutate(maxprec = lag(maxprec,1)) %>% dplyr::slice(-1)
        stepcurve <- rbind(pr, lag1) %>% dplyr::arrange(recall, dplyr::desc(maxprec))

        p3 <- ggplot2::ggplot(stepcurve) + 
           # ggplot2::geom_line(ggplot2::aes(x=recall, y=maxprec, group=fold, col="folds"), alpha=0.5) + 
            ggplot2::geom_line(data=stepcurve %>% dplyr::group_by(recall) %>% dplyr::summarise(meanprec=mean(maxprec)), 
                      ggplot2::aes(x=recall, y=meanprec, col="mean of folds")) + 
            ggplot2::geom_hline(yintercept=thestats["ExpectedAUPRC",], linetype="dotted", col="darkorange") +
            ggplot2::scale_color_manual(name = "Interpolated PRC",
                               breaks = c("mean of folds"),
                               values = c("mean of folds" = "darkorange")) +
            ggplot2::theme(legend.position="top", legend.box = "horizontal") + 
            ggplot2::expand_limits(x=c(0,1), y=c(0, 1)) + 
            ggplot2::labs(title=paste0(metrics$res_combined$geneset[1]),
                          caption = paste0("mean of folds AUPRC = ", thestats["AUPRC",]))


        # early recall (TPR)
        p6 <- ggplot2::ggplot() + 
            ggplot2::geom_line(data=metrics$res_combined %>% dplyr::filter(rank<=100) %>% 
                          dplyr::group_by(rank) %>% dplyr::summarise(REC = mean(REC)), 
                          ggplot2::aes(x=rank, y=REC, col="mean of folds")) + 
            ggplot2::xlab("rank position") + 
            ggplot2::ylab("TPR (Recall)") + 
            ggplot2::scale_color_manual(name = "Early Recall",
                               breaks = c("mean of folds"),
                               values = c("mean of folds" = "darkorange") ) +
            ggplot2::theme(legend.position="top", legend.box = "horizontal") + 
            ggplot2::expand_limits(y=c(0, 1))
        
        # rank distribution of true positive hits
        p7 <- ggplot2::ggplot(metrics$res_combined %>% dplyr::filter(InValset==1)) + 
          ggplot2::geom_density(ggplot2::aes(x=rank, group=fold), col="#AAAAAA1A") + 
          ggplot2::geom_hline( ggplot2::aes(col="uniform dist", yintercept=1/nrow(metrics$res_avg)), linetype="dashed") + 
          ggplot2::theme(legend.position="top", legend.box = "horizontal")
          
        plotcomp <- ggplot2::ggplot_build(p7)
        meandensity <- plotcomp$data[[1]] %>% dplyr::group_by(x) %>% dplyr::summarise(hit_density = mean(ymax))
        p7 <- p7 +  ggplot2::geom_line(data = meandensity, ggplot2::aes(x = x, y = hit_density, col="mean of folds")) + 
          ggplot2::scale_color_manual(name = "Recall Density",
                                      breaks = c("mean of folds", "uniform dist"),
                                      values = c("mean of folds" = "darkorange", "uniform dist"="black") )
          
        
        # p7 <- ggplot2::ggplot(metrics$res_combined %>% dplyr::filter(InValset==1)) + 
        #     ggplot2::geom_histogram(ggplot2::aes(x=rank), binwidth=100, fill="darkorange") + 
        #     ggplot2::xlab("CV rank (binwidth=100)") + 
        #     ggplot2::ylab("Geneset genes in bin")# ranking distribution of hits in bins of 100 for mean of folds

        
        fname = paste("RWR-CV", metrics$res_combined$geneset[1], get_base_name(dataPath), modname, sep="_")
        fname = paste(substr(fname, 1, 99), 'plots.png', sep='.')
        out_path = file.path(outdirPath, fname)
        
        png(out_path, width = 1200, height=1000)
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))

        print(p1, vp = vplayout(1, 1))
        print(p3, vp = vplayout(1, 2))
        print(p6, vp = vplayout(2, 1))
        print(p7, vp = vplayout(2, 2))
        
        dev.off()
        message('File saved:', out_path, "\n")
    }
    
    ######## SINGLETONS (many folds) #########

    if(metrics$res_combined$method[1] == "singletons") {
        ### ROC
        p1 <- ggplot2::ggplot() + 
          ggplot2::geom_line(data=metrics$res_combined %>% dplyr::group_by(FPR) %>% dplyr::summarise(mean = mean(REC)), 
                             ggplot2::aes(x=FPR, y=mean, col="mean of folds")) +
          ggplot2::geom_abline(ggplot2::aes(intercept=0, slope=1), linetype="dashed", col="darkgrey") +
          #ggplot2::annotate("text", label=paste0("mean AUROC = ",thestats["AUROC",]), x = 0.75, y = 0.01, size = 4, colour = "darkorange3") + 
          ggplot2::xlab("FPR") +
          ggplot2::ylab("TPR (Recall)") + 
          ggplot2::scale_color_manual(name = "ROC",
                                      breaks = c("mean of folds"),
                                      values = c("mean of folds" = "darkorange") ) +
          ggplot2::labs(title=paste0(metrics$res_combined$geneset[1]), 
                        caption = paste0("mean of folds AUROC = ", thestats["AUROC",])) +
          ggplot2::theme(legend.position="top", legend.box = "horizontal")
        
        ### PRC
        # p2 <- ggplot2::ggplot() + 
        #     ggplot2::geom_line(data=metrics$res_combined %>% dplyr::filter(InValset==1) %>% dplyr::group_by(REC) %>% dplyr::summarise(PREC = mean(PREC)), 
        #               ggplot2::aes(x=REC, y=PREC, col="mean of folds")) + 
        #     ggplot2::geom_line(data=metrics$res_avg, 
        #               ggplot2::aes(x=REC, y=PREC, col="mean ranks")) + 
        #     ggplot2::geom_hline(yintercept = sum(metrics$res_avg$InValset)/nrow(metrics$res_avg), linetype="dotted") + 
        #     ggplot2::scale_color_manual(name = "Precision / Recall",
        #                        breaks = c("mean of folds", "mean ranks"),
        #                        values = c("mean of folds" = "darkorange", "mean ranks" = "darkred") ) +
        #     ggplot2::theme(legend.position="top", legend.box = "horizontal") +
        #     ggplot2::labs(title=paste0(metrics$res_combined$geneset[1]),caption = paste0("mean of folds AUPRC = ", thestats["AUPRC",], " | expected = ", thestats["ExpectedAUPRC",])) +
        #     ggplot2::expand_limits(x=c(0,1), y=c(0, 1))
        
        ### Interpolated PRC
        # for each singleton, interpolating the PR values at specific recall (x-axis) positions
        pr          <- expand.grid(fold=seq(1, folds,by=10), recall = seq(0,1,by=0.1))
        pr$maxprec  <- foreach::foreach(f=pr$fold, r = pr$recall, .combine = c) %do%
            {
                prec <- metrics$res_combined  %>% dplyr::filter(fold==f, REC >= r) %>% dplyr::pull(PREC)
                if(length(prec)>0)
                    max(prec)
                else
                    return(0)
            }
        # make it into a step curve for each fold
        lag1 <- pr %>% dplyr::group_by(fold) %>% dplyr::mutate(maxprec = lag(maxprec,1)) %>% dplyr::slice(-1)
        stepcurve <- rbind(pr, lag1) %>% dplyr::arrange(recall, dplyr::desc(maxprec))
        p3 <- ggplot2::ggplot(stepcurve) + 
            ggplot2::geom_line(data=stepcurve %>% dplyr::group_by(recall) %>% dplyr::summarise(meanprec=mean(maxprec)), 
                      ggplot2::aes(x=recall, y=meanprec, col="mean of subset of folds")) + 
            ggplot2::scale_color_manual(name = "Interpolated PRC",
                               breaks = c("mean of subset of folds"),
                               values = c("mean of subset of folds" = "darkorange") ) +
            ggplot2::theme(legend.position="top", legend.box = "horizontal") + 
            ggplot2::expand_limits(x=c(0,1), y=c(0, 1)) + 
            ggplot2::labs(title=paste0(metrics$res_combined$geneset[1]),
                        caption = paste0("mean of folds AUPRC = ", thestats["AUPRC",]))
        
        # early recall (TPR)
        p5 <- ggplot2::ggplot() + 
          ggplot2::geom_line(data=metrics$res_combined %>% dplyr::filter(rank<=100) %>% 
                               dplyr::group_by(rank) %>% dplyr::summarise(REC = mean(REC)), 
                             ggplot2::aes(x=rank, y=REC, col="mean of folds")) + 
          ggplot2::xlab("rank position") + 
          ggplot2::ylab("TPR (Recall)") + 
          ggplot2::scale_color_manual(name = "Early Recall",
                                      breaks = c("mean of folds"),
                                      values = c("mean of folds" = "darkorange") ) +
          ggplot2::theme(legend.position="top", legend.box = "horizontal") + 
          ggplot2::expand_limits(y=c(0, 1))
        
        # rank distribution of true positive hits
        p7 <- ggplot2::ggplot(metrics$res_combined %>% dplyr::filter(InValset==1)) + 
          ggplot2::geom_density(ggplot2::aes(x=rank, group=fold),col="#AAAAAA1A") + 
          ggplot2::geom_hline( ggplot2::aes(col="uniform dist", yintercept=1/nrow(metrics$res_avg)), linetype="dashed") + 
          ggplot2::theme(legend.position="top", legend.box = "horizontal")
        
        plotcomp <- ggplot2::ggplot_build(p7)
        meandensity <- plotcomp$data[[1]] %>% dplyr::group_by(x) %>% dplyr::summarise(hit_density = mean(ymax))
        p7 <- p7 +  ggplot2::geom_line(data = meandensity, ggplot2::aes(x = x, y = hit_density, col="mean of folds")) + 
          ggplot2::scale_color_manual(name = "Recall Density",
                                      breaks = c("mean of folds", "uniform dist"),
                                      values = c("mean of folds" = "darkorange", "uniform dist"="black") )


        fname = paste("RWR-CV",metrics$res_combined$geneset[1],get_base_name(dataPath), modname, sep="_")
        fname = paste(substr(fname, 1, 99), 'plots.png', sep='.')
        out_path = file.path(outdirPath, fname)
        
        png(out_path, width = 1200, height=1000)
        # print( (p1+p2)/(p3+p4)/(p5) + plot_layout(heights=c(2,2,1)) )
          
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
        
        print(p1, vp = vplayout(1, 1))  # Top left
        print(p3, vp = vplayout(1, 2))  # top R
        print(p5, vp = vplayout(2, 1))  # Bottom L
        print(p7, vp = vplayout(2, 2))  # Bottom R
        
        dev.off()
        cat('File saved:', out_path,"\n")
    }
    

    ####### LOO ############

    if(metrics$res_combined$method[1] == "loo")
    {
        # early recall
        p1 <- ggplot2::ggplot() + 
          ggplot2::geom_line(data=metrics$res_avg %>% dplyr::filter(rerank <= 100) %>% 
                               dplyr::group_by(rerank) %>% dplyr::summarise(REC = mean(REC)), 
                             ggplot2::aes(x=rerank, y=REC, col="mean of folds")) + 
          ggplot2::xlab("rank position") + 
          ggplot2::ylab("TPR (Recall)") + 
          ggplot2::scale_color_manual(name = "Early Recall",
                                      breaks = c("mean of folds"),
                                      values = c("mean of folds" = "darkorange") ) +
          ggplot2::theme(legend.position="top", legend.box = "horizontal") + 
          ggplot2::expand_limits(y=c(0, 1))
        
       
        # PRC  step curve
        pr          <- expand.grid(fold="median", recall = seq(0,1,by=0.01))
        pr$maxprec  <- foreach::foreach(f=pr$fold, r = pr$recall, .combine = c) %do%
            {
                prec <- metrics$res_avg  %>% dplyr::filter(REC >= r) %>% dplyr::pull(PREC)
                if(length(prec)>0)
                    max(prec)
                else
                    return(0)
            }
        
        lag1 <- pr %>% dplyr::mutate(maxprec = lag(maxprec,1)) %>% dplyr::slice(-1)
        stepcurve <- rbind(pr, lag1) %>% dplyr::arrange(recall, dplyr::desc(maxprec))
        p2 <- ggplot2::ggplot(stepcurve) + 
            ggplot2::geom_line(data=stepcurve %>% dplyr::group_by(recall) %>% dplyr::summarise(meanprec=mean(maxprec)), 
                      ggplot2::aes(x=recall, y=meanprec, col="mean of folds")) + 
            ggplot2::scale_color_manual(name = "Interpolated PRC",
                               breaks = c("folds", "mean of folds"),
                               values = c("folds"="grey", "mean of folds" = "darkorange") ) +
            ggplot2::theme(legend.position="top", legend.box = "horizontal")
        
        # precision / recall
        # p2 <- ggplot(metrics$res_avg) + 
        #     geom_line(aes(x=REC, y=PREC , col="mean ranks")) +
        #     xlab("TPR (Recall)") + 
        #     ylab("Precision") +
        #     geom_hline(yintercept = sum(metrics$res_avg$InValset)/nrow(metrics$res_avg), linetype="dotted") + 
        #     scale_color_manual(name = "Precision / Recall",
        #                        breaks = c("mean ranks"),
        #                        values = c("mean ranks" = "darkred") ) +
        #     theme(legend.position="top", legend.box = "horizontal") + 
        #     expand_limits(x=c(0,1), y=c(0, 1))
        
            

        # ranking distribution of hits in bins of 100 
   # ranking distribution of hits in bins of 100 for meanof folds
        p4 <- ggplot2::ggplot(metrics$res_combined %>% dplyr::filter(InValset==1)) + 
            ggplot2::geom_histogram(ggplot2::aes(x=rank), binwidth=100, fill="darkred") + 
            ggplot2::xlab("CV Mean Rank (binwidth=100)") + 
            ggplot2::ylab("Geneset genes in bin")

        # p7 <- ggplot2::ggplot(in_valset_data) + 
        #   ggplot2::geom_density(ggplot2::aes(x=rank, group=fold),col="grey", alpha=0.5)
        # plotcomp <- ggplot2::ggplot_build(p7)
        # meandensity <- plotcomp$data[[1]] %>% dplyr::group_by(x) %>% dplyr::summarise(hit_density = mean(ymax))
        # p7 <- p7 +  ggplot2::geom_line(data = meandensity, ggplot2::aes(x = x, y = hit_density), col="darkorange")
        
        fname <- paste("RWR-CV", metrics$res_combined$geneset[1], get_base_name(dataPath), modname, sep = "_")
        fname = paste(substr(fname, 1, 99), 'plots.png', sep='.')
        out_path = file.path(outdirPath, fname)
        
        png(out_path, width = 1200, height=1000)
        # print( (p1 + p2)/(p3 + p4) + plot_layout(heights=c(1,1)) )
        # dev.off()
        # cat('File saved:', out_path,"\n")
        
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
        
        print(p1, vp = vplayout(1, 1))  # Top left
        print(p2, vp = vplayout(1, 2))  # Top right 
        print(p4, vp = vplayout(2, 1:2))  # Bottom 
        dev.off()
        message(paste('File saved:', out_path))
       
    }

}

post_process_rwr_output_cv <- function(res, extras, folds, nw.mpo) {
  ############# post process results  ########################################
  ## How to deal with missing geneset genes:
  # Way 1: Append all geneset genes that are missing from the multiplex (extras) to the end of folds, giving them worst rank.
  # and ensure they are spread evenly amongst the folds.
  # -- We do this so that metrics will be with respect to the full original geneset
  # -- For LOO, though, we need to append them to the end of the combined results, NOT to each fold.
  #
  # Problem with this approach: it makes a massive assumption that all the extras would rank at the bottom.
  # This ends up punishing the metrics very harshly (and gives whacky ROC plots!).
  # In reality, if the extras are functionally not relevant to the seeds (the null hyp) then they would have a
  # uniform distribution across the rankings. But this makes it's own assumptions. There is no best way to do this!

  res_method <- res[[1]]$method[1]
  extras_exist <- !is.null(extras)

  if (extras_exist & res_method %in% c("kfold", "singletons")) {
    for (i in 0:(nrow(extras) - 1))
    {
      res.tmp <- res[[(i %% folds) + 1]] # Get res[ i mod extras ]'th res
      extra_ith_row <- data.frame(
        NodeNames = extras$gene[i + 1],
        Score = 0,
        rank = max(res.tmp$rank) + 1,
        InValset = 1,
        num_in_network = dplyr::first(res.tmp$num_in_network),
        num_seeds = dplyr::first(res.tmp$num_seeds),
        num_leftout = dplyr::first(res.tmp$num_leftout),
        networks = dplyr::first(res.tmp$networks),
        fold = dplyr::first(res.tmp$fold),
        modname = dplyr::first(res.tmp$modname),
        geneset = dplyr::first(res.tmp$geneset),
        seed = "missing",
        leftout = "missing",
        method = dplyr::first(res.tmp$method)
      )
      res[[(i %% folds) + 1]] <- rbind(res.tmp, extra_ith_row)
    }
  }

  # combine results from the res list into one dataframe
  res_combined <- dplyr::bind_rows(res) %>% dplyr::arrange(rank)

  if (extras_exist & res_method == "loo") {
    res_combined <- rbind(res_combined, data.frame(
      NodeNames = extras$gene, Score = 0,
      rank = nw.mpo$Number_of_Nodes_Multiplex,
      InValset = 1,
      num_in_network = dplyr::first(res[[1]]$num_in_network),
      num_seeds = dplyr::first(res[[1]]$num_seeds),
      num_leftout = dplyr::first(res[[1]]$num_leftout),
      networks = dplyr::first(res[[1]]$networks),
      fold = seq(length(res) + 1, length(res) + nrow(extras)),
      modname = dplyr::first(res[[1]]$modname),
      geneset = dplyr::first(x = res[[1]]$geneset),
      seed = "", leftout = extras$gene, method = dplyr::first(res[[1]]$method)
    ))
  }

  # print(head(res_combined))

  res_combined
}

calculate_average_rank_across_folds_cv <- function(res_combined){
    # get the average ranks across all CV folds/runs
    res_avg <- res_combined %>%
        dplyr::group_by(NodeNames) %>%
        dplyr::summarise(meanrank = mean(rank), 
                         InValset = dplyr::first(InValset), 
                         geneset=dplyr::first(geneset),
                         num_in_network=dplyr::first(num_in_network)) %>%
                              dplyr::arrange(meanrank) %>%
                              dplyr::mutate(
                                rerank = rank(meanrank, ties.method = "min"),
                                .after=meanrank) 
                              # here we rerank the mean ranks of the folds to
                              #  go from 1:nrow again

  res_avg
}

########################################################################
# Main Function
########################################################################

#' RWR Cross Validation
#'
#' `RWR_CV` RWR Cross Validation performs K-fold cross validation on a single
#'  gene set, finding the RWR rank of the left-out genes.  Can choose: (1)
#'  leave-one-out (`loo`) to leave only one gene from the gene set out and find
#'  its rank, (2) cross-validation (`kfold`) to run k-fold cross-validation for
#'  a specified value of *k*, or (3) singletons (`singletons`) to use a single
#'  gene as a seed and find the rank of all remaining genes.
#'
#' @param data The path to the .Rdata file containing your multiplexed
#'                  functional networks. This file is produced by
#'                  RWR_make_multiplex. Default NULL
#' @param geneset_path The path to the gene set file. It must have the
#'                    following first two columns with no headers
#'                    tab-delimited: 
#'                    {<}setid{>} {<}gene{>} {<}weight{>}.   Default NULL
#' @param method Cross-validation method. Choice of: 'kfold', 'loo', or
#'               'singletons'. Default 'kfold'
#' @param folds Number (k) of folds to use in k-fold CV. Default 5
#' @param restart Set the restart parameter \[0,1). Higher value means the
#'                walker will jump back to seed node more often. Default 0.7
#' @param tau Comma-separated list of values between that MUST add up to the
#'            number of network layers in the .Rdata file. One value per 
#'            network layer that determines the probability that the random
#'            walker will restart in that layer. e.g. if there are three layers
#'            (A,B,C) in your multiplex network, then --tau '0.2,1.3,1.5' will
#'            mean that layer A is less likely to be walked on after a restart
#'            than layers B or C.  Default 1.0
#' @param numranked Proportion of ranked genes to return \[0,1\].  e.g. 0.1 will
#'                  return the top 10%. Default 1.0
#' @param outdir Path to the output directory. Both 'fullranks' and 
#'                   'medianranks' will be saved with auto-generated filenames.
#'                   Can be overridden by specifically setting 'out_full_ranks'
#'                   and 'out_mean_ranks' parameters. No defined path will 
#'                   output within the same directory from which the original
#'                   code was run.
#'                   Default NULL
#' @param modname String to include in output file name.  Default "default"
#' @param plot Output plots of ROC, PRC, etc. to file. Default FALSE
#' @param out_full_ranks Specify the full path for the full results. Ignores
#'                     outdir and modName, using this path instead.
#'                     Default NULL
#' @param out_mean_ranks Specify the full path for the mean results. Ignores
#'                       outdir and modName, using this path instead. 
#'                       Default NULL 
#' @param threads Specify the number of threads to use. Default for your system
#'                is all cores - 1.
#' @param verbose Verbose mode. Default FALSE
#' @param write_to_file Also write the result to a file. Default FALSE, however,
#'                      if output paths are included, the boolean is switched
#'                      to true.
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
RWR_CV <- function(data = NULL,
                   geneset_path = NULL,
                   method = "kfold",
                   folds = 5,
                   restart = 0.7,
                   tau = 1.0,
                   numranked = 1.0,
                   outdir = NULL,
                   modname = "default",
                   plot = FALSE,
                   out_full_ranks = NULL,
                   out_mean_ranks = NULL,
                   threads = 1,
                   verbose = FALSE,
                   write_to_file = FALSE) {

  ############# Initialize  ##################################################
  data_list <- load_multiplex_data(data)
  nw_mpo <- data_list$nw.mpo
  nw_adjnorm <- data_list$nw.adjnorm
  
  if ((!is.null(outdir) |
       !is.null(out_full_ranks) |
       !is.null(out_mean_ranks)
       ) & write_to_file == FALSE) {
    warning(sprintf("write_to_file was set to false, however, an output file path was set. write_to_file has been updated to TRUE.\n")) #nolint warning
    write_to_file <- TRUE
  }

  if (is.null(outdir)){
    outdir <- './'
  }

  tau <- get_or_set_tau(nw_mpo, tau)

  # Geneset is required.
  geneset_list <- load_geneset(geneset_path, nw_mpo, verbose = verbose)
  geneset <- geneset_list$geneset
  extras <- geneset_list$extras

  # check and update data if necessary 
  updated_data_list <- update_folds_by_method(geneset, method, folds)
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
             modname,
             threads,
             verbose)

  ############# post process results  ########################################
  res_combined <- post_process_rwr_output_cv(res, extras, folds, nw_mpo)
  # get the average ranks across all CV folds/runs
  res_avg <- calculate_average_rank_across_folds_cv(res_combined)
  # calculate ROC, PRC, P@k, AUROC, AUPRC, NDCG
  metrics <- calc_metrics_cv(res_combined, res_avg)
  ############# Save full results  ###########################################
  if (!is.null(out_full_ranks)) {
    out_path <- out_full_ranks
  } else {
    out_path <- get_file_path("RWR-CV_",
      res_combined$geneset[1],
      get_base_name(data),
      modname,
      outdir = outdir,
      ext = ".fullranks.tsv"
    )
  }

  if (write_to_file) {
  if (!file.exists(outdir)) {
    print("Creating directory")
    dir.create(outdir, recursive = TRUE)
  }

    combined <- res_combined %>%
      dplyr::group_by(fold) %>%
      dplyr::slice_head(prop = numranked)

    write_table(combined, out_path)
  }

  ############# Save averaged results  #######################################
  # Save the summary table of mean rank across CV folds/runs
  if (!is.null(out_mean_ranks)) {
    out_path <- out_mean_ranks
  } else {
    out_path <- get_file_path("RWR-CV_",
      res_avg$geneset[1],
      get_base_name(data),
      modname,
      outdir = outdir,
      ext = ".meanranks.tsv"
    )
  }

  if (write_to_file) {
    write_table(
      res_avg %>%
        dplyr::slice_head(prop = numranked),
      out_path
    )
  }

  ############# Save metrics  ################################################

  out_path <- get_file_path("RWR-CV_",
    metrics$res_combined$geneset[1],
    get_base_name(data),
    modname,
    outdir = outdir,
    ext = ".metrics.tsv"
  )
  if (write_to_file) {
    write_table(metrics$res_avg, out_path)
  }

  out_path <- get_file_path("RWR-CV_",
    metrics$res_combined$geneset[1],
    get_base_name(data),
    modname,
    outdir = outdir,
    ext = ".summary.tsv"
  )

  if (write_to_file) {
    write_table(metrics$summary, out_path)
  }

  ############# Save plots  ##################################################
  if (plot) {
    message("Saving plots ...")
    save_plots_cv(metrics, geneset, folds, data, modname, outdir)
  }


  return(
    list("fullranks" = res_combined,
         "meanranks" = res_avg,
         "metrics" = metrics$res_avg,
         "summary" = metrics$summary)
      )
}
