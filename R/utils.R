########################################################################
# Common functions for RWRtoolkit.
########################################################################
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



load_geneset <- function(path, nw.mpo = NULL, verbose = FALSE, select=NULL) {
  if (is.null(path)) {
    return(NULL)
  } else if (!file.exists(path)) {
    stop("ERROR: geneset file does not exist: ", path)
  } else {
    # geneset <- read.table(path, header = F, colClasses = c("character"))
    geneset <- data.table::fread(
                  path,
                  header =F,
                  select = select,
                  colClasses = c("character")
                )
    if (ncol(geneset) < 2) {
      stop("Your geneset file is incorrectly formatted. Please see documentation.") #nolint message
    }
    # Check if seed weights are included by user or not.
    numericV3 <- suppressWarnings(as.numeric(geneset$V3))
    if (!is.null(geneset$V3) && all(!is.null(numericV3)) && all(!is.na(numericV3))) {
      geneset <- dplyr::select(geneset, 1:3)
      geneset$V3 <- as.numeric(geneset$V3)
      colnames(geneset) <- c("setid", "gene", "weight")
    } else {
      # If there is a non-numeric third col, then just give all genes a weight of 1.
      geneset <- dplyr::select(geneset, 1:2) %>% dplyr::mutate(weight = 1)
      colnames(geneset) <- c("setid", "gene", "weight")
    }
  }

  # Filtering
  # Remove any duplicate genes
  geneset_orig <- geneset
  geneset <- geneset_orig %>% dplyr::distinct(gene, .keep_all = T)
  ngenes <- nrow(geneset)
  extras <- NULL

  # TODO: @kylesullivan in @GRINT.R has a duplicates written to file ~ll 416
  # TODO: Discussion on whether genes ought to be in multiplex (@kylesullivan - GRINT.R line 425)

  # Remove any genes not in multiplex
  if (!is.null(nw.mpo)) {
    geneset <- geneset %>% dplyr::filter(gene %in% nw.mpo$Pool_of_Nodes)

    # Warn user if some genes in the seed geneset are not in the multiplex
    if (nrow(geneset) < ngenes) {
      message(sprintf("%s genes from geneset (%s) are present in the multiplex \n", nrow(geneset), path))
      extras <- geneset_orig %>% dplyr::slice(which(!geneset_orig$gene %in% nw.mpo$Pool_of_Nodes))
      warning(sprintf("WARNING:: %s genes from geneset (%s) were not found in multiplex: %s\n ", nrow(extras), path, list(extras)))
      warning(sprintf("Please ensure your geneset files are formated: [ setids | genes | weights (if weights exist) ]\n"))
    } else {
      message(sprintf("All %s genes in the geneset (%s) are present in the multiplex\n", nrow(geneset), path))
    }
    # print(head(geneset))
  } else {
    message(sprintf("Geneset %s was not filtered (nw.mpo not passed to load_geneset utility function)\n", path))
  }

  if (verbose) {
    message(sprintf("Loaded gene set (%s genes):\n", nrow(geneset)))
    print(head(geneset))
  }

  return(list("geneset" = geneset, "extras" = extras))
}

get_or_set_tau <- function(nw.mpo, optTau) {
  # Ensure Tau param is appropriate (must be one value per layer, and must add up to NumLayers).

  tau <- if (typeof(optTau) != "double") as.numeric(unlist(strsplit(optTau, split = ","))) else optTau

  if (sum(tau) != nw.mpo$Number_of_Layers || length(tau) != nw.mpo$Number_of_Layers) {
    tau <- rep(1, nw.mpo$Number_of_Layers)
    warning(sprintf(
      "WARNING:: Your comma-delimited tau parameter values (%s) has the incorrect length or does not add up to the number of network layers (%s).  Automatically re-set tau to: %s \n",
      optTau, nw.mpo$Number_of_Layers, list(tau)
    ))
  }
  return(tau)
}

chunk <- function(x, n) {
  if (length(x) <= n) {
    warning(sprintf("Geneset Length:(%s) exceded by number of folds:(%s).\n", length(x), n))
  }
  split(x, cut(seq_along(x), n, labels = FALSE))
}

calc_ROCPRC <- function(df, scorescol=NULL, labelscol=NULL, totP = NULL) {
    ### assumes scorescol is already sorted in the order you want it (best to worst)
    df <- dplyr::mutate(df,
                 TP = dplyr::case_when(get(scorescol)>0  & get(labelscol)==1 ~ 1, TRUE ~ 0),
                 FP = dplyr::case_when(get(scorescol)>0  & get(labelscol)==0 ~ 1, TRUE ~ 0))
                 #TN = dplyr::case_when(get(scorescol)==0 & get(labelscol)==0 ~ 1, TRUE ~ 0),
                 #FN = dplyr::case_when(get(scorescol)==0 & get(labelscol)==1 ~ 1, TRUE ~ 0))
    
  if (is.null(totP)) {
    totP <- sum(get(labelscol, df) == 1)
  } # number of positives in gold set
    totN    <- sum(get(labelscol, df)==0) # number of negatives in gold set
    totTP   <- sum(df$TP==1, na.rm = T) # true positive predictions
    #totTN   <- sum(df$TN==1, na.rm = T) # true negative predictions
    totFP   <- sum(df$FP==1, na.rm = T) # false positive predictions (i.e. an edge/gene was predicted but should not have been )
    #totFN   <- sum(df$FN==1, na.rm = T) # false negative predictions (i.e. an edge/gene was not predicted but should have been)
    
    df <- dplyr::mutate(df, 
                 cum_TP = cumsum(TP), 
                 cum_FP = cumsum(FP))
                 #cum_TN = cumsum(TN), 
                 #cum_FN = cumsum(FN))
    
    df <- dplyr::mutate(df, 
                 FPR  = round(cum_FP / (totN),3),
                 PREC = round(cum_TP / (cum_TP + cum_FP), 3), # definitely correct
                 REC  = round(cum_TP / (totP),3) )            # recall = TPR.
    
    # df <- dplyr::mutate(df,
    #              dcg  = round(cumsum(TP/log2(dplyr::row_number()+1)),3),
    #              idcg = round(cumsum(1/log2(dplyr::row_number()+1)),3),
    #              ndcg = round(dcg/idcg,3))
    
    return(df)
}

area_under_curve <- function(x, y, from = min(x, na.rm = TRUE), to = max(x, na.rm = TRUE),
                             method = c("trapezoid", "step", "spline"), absolutearea = FALSE,
                             subdivisions = 100, na.rm = FALSE, ...) {
  ### Based on the DescTools `AUC` function.
  ### from here: https://github.com/cran/DescTools/blob/master/R/StatsAndCIs.r

  # calculates Area unter the curve
  # example:
  #   AUC( x=c(1,2,3,5), y=c(0,1,1,2))
  #   AUC( x=c(2,3,4,5), y=c(0,1,1,2))

  if (na.rm) {
    idx <- na.omit(cbind(x, y))
    x <- x[idx]
    y <- y[idx]
  }

  if (length(x) != length(y)) {
    stop("length x must equal length y")
  }

  if (length(x) < 2) {
    return(NA)
  }

  o <- order(x)
  x <- x[o]
  y <- y[o]

  ox <- x[o]
  oy <- y[o]

  method <- match.arg(method)

  if (method == "trapezoid") {
    # easy and short
    # , "trapezoid" = { a <- sum((apply( cbind(y[-length(y)], y[-1]), 1, mean))*(x[-1] - x[-length(x)])) }

    ## Default option
    if (!absolutearea) {
      values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
    } else { ## Absolute areas

      idx <- which(diff(oy >= 0) != 0)
      newx <- c(x, x[idx] - oy[idx] * (x[idx + 1] - x[idx]) / (y[idx + 1] - y[idx]))
      newy <- c(y, rep(0, length(idx)))
      values <- approx(newx, newy, xout = sort(unique(c(from, to, newx[newx > from & newx < to]))), ...)

      res <- 0.5 * sum(diff(values$x) * (abs(values$y[-1]) + abs(values$y[-length(values$y)])))
    }
  } else if (method == "step") {
    # easy and short
    # , "step" = { a <- sum( y[-length(y)] * (x[-1] - x[-length(x)])) }

    ## Default option
    if (!absolutearea) {
      values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)

      res <- sum(diff(values$x) * values$y[-length(values$y)])
      # res <- sum( y[-length(y)] * (x[-1] - x[-length(x)]))
    } else { ## Absolute areas

      idx <- which(diff(oy >= 0) != 0)
      newx <- c(x, x[idx] - oy[idx] * (x[idx + 1] - x[idx]) / (y[idx + 1] - y[idx]))
      newy <- c(y, rep(0, length(idx)))
      values <- approx(newx, newy, xout = sort(unique(c(from, to, newx[newx > from & newx < to]))), ...)

      res <- sum(diff(values$x) * abs(values$y[-length(values$y)]))
    }
  } else if (method == "spline") {
    if (absolutearea) {
      myfunction <- function(z) {
        abs(splinefun(x, y, method = "natural")(z))
      }
    } else {
      myfunction <- splinefun(x, y, method = "natural")
    }

    res <- integrate(myfunction, lower = from, upper = to, subdivisions = subdivisions)$value
  }

  return(res)
}

write_table <- function(table, path, row_names = F, col_names = T, verbose = FALSE) {
  if (length(table) == 0 || any(is.na(table))) {
    warning(sprintf("Table to be saved at %s is empty\n", path))
  }
  # Create out_dir if it doesn't exist (avoid warning message if out_dir exists).
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  # Save the table.
  write.table(table,
    path,
    sep = "\t",
    quote = F,
    col.names = col_names,
    row.names = row_names
  )
  if (verbose) {
    message(sprintf("Saved data to file: %s\n", path))
  }
}

get_base_name <- function(file_path){
  base_from_path <- basename(file_path)
  base_name <- sub("\\.[Rr][Dd][Aa][Tt][Aa].*", "", base_from_path)
  base_name
}

get_file_path <- function(..., outdir = NULL, ext = ".tsv") {
  filename <- paste(..., sep = "_")
  filename <- substr(filename, 1, 99)
  filename <- paste0(filename, ext)

  if (!is.null(outdir)) {
    filename <- file.path(outdir, filename)
  }
  return(filename)
}

dump_layers <- function(mpo, outdir = NULL) {
  # suppressPackageStartupMessages(library(igraph)) # as_data_frame
  # Helper for when you already have a multiplex object (eg, from a third-party),
  # and you want the underlying layers, eg to visualize in your favorite software.
  for (network in mpo) {
    # There are multiple types of objects in `nw.mpo`, but
    # we just want to dump the actual igraph networks.

    if (class(network) == "igraph") {
      df <- igraph::as_data_frame(network)
      name <- df$type[[1]]
      file_path <- get_file_path(name, outdir = outdir)
      write_table(df, file_path)
      message(sprintf("Saved network layer: %s\n", file_path))
    }
  }
}

dump_nodes <- function(mpo, outdir = NULL) {
  mpoNodes <- mpo$Pool_of_Nodes
  df <- as.data.frame(mpoNodes)

  if (length(df) == 0) stop("Pool of Nodes to be saved has 0 nodes")

  file_path <- get_file_path("pool-of-nodes", outdir = outdir)
  print("File path")
  print(file_path)
  print("actualdf")
  print(df)
  write_table(df, file_path)
  message(sprintf("Saved pool of nodes: %s\n", file_path))
}

generate_randomset <- function(mpo, n = 10, name = "RAND", outdir = NULL) {
  random_nodes <- data.frame(geneset = name, sample(mpo$Pool_of_Nodes, n))
  file_path <- get_file_path(paste0("geneset_", name), outdir = outdir)

  write.table(random_nodes,
              file_path,
              col.names = F,
              row.names = F, quote = F, sep = "\t")
  message(sprintf("Saved random ganeset: %s\n", file_path))
  random_nodes
}

## Layout function for plots:

vplayout <- function(x, y) {
  grid::viewport(layout.pos.row = x, layout.pos.col = y)
}


# End.

# Load the multiplex network and adjacency matrices
# (i.e. nw.mpo, nw.adj, nw.adjnorm)
load_multiplex_data <- function(filepath_or_url) {
  if (is.null(filepath_or_url)) {
    stop("ERROR: Mandatory arguement data is missing.")
  }

  is_url <- stringr::str_detect(filepath_or_url, pattern = "http")

  if (!is_url && !file.exists(filepath_or_url)) {
    stop("ERROR: Rdata input file does not exist: ", filepath_or_url)
  }

  updated_file_path <- if (is_url) url(filepath_or_url) else filepath_or_url
  # this contains the multiplex network layers and adj matrix
  load(updated_file_path)

  if (is.null(nw.mpo)) stop("ERROR: failed to load multiplex RData object") # nolint

  return(list(
    nw.mpo = nw.mpo,
    nw.adj = nw.adj,
    nw.adjnorm = nw.adjnorm
  ))
}
