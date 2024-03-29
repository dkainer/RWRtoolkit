% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RWR_network_aggregation.R
\name{RWR_network_aggregation}
\alias{RWR_network_aggregation}
\title{RWR_network_aggregation}
\usage{
RWR_network_aggregation(
  data = NULL,
  flist = NULL,
  outdir = NULL,
  merged_with_all_edges = F,
  merged_with_edgecounts = F,
  verbose = F
)
}
\arguments{
\item{merged_with_all_edges}{A boolean denoting a return of a merged
down multiplex network along with
network edge counts and vertex counts.
Default False}

\item{merged_with_edgecounts}{A boolean denoting a return of a merged
down multiplex, but simplified with
edge weights denoting the total number
of layers in which that edge existed.
Default False}
}
\description{
This function acts as an aggregator for RWR multiplex objects.
}
\examples{
# An example of aggregating multiplex networks 
# into a singular network 
extdata.dir <- system.file("example_data", package="RWRtoolkit")
mpo_path <- paste(extdata.dir, "string_interactions.Rdata", sep = "/")

output_aggregated_networks <- RWR_network_aggregation(
   data = mpo_path,
   merged_with_all_edges = T,
   merged_with_edgecounts = T,
   verbose = T
)

print(output_aggregated_networks)
}
