
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RWRtools

<!-- badges: start -->
<!-- badges: end -->

RWRtools enables easy use of RandomWalk with Restart on multiplex
networks. These functions are an extension to the
[RandomWalkRestartMH](https://github.com/alberto-valdeolivas/RandomWalkRestartMH)
R package. Also provided are scripts for use as command line tools.

## Installation

##### Dependencies

Installation of this R package requires R and r-devtools. If you use
prefer the use of conda you can create the base environment with
`conda create --name r-rwrtools -c conda-forge r-base=4.0.2 r-devtools`. You
can also install devtools from within a base R environment with
`install.packages("devtools")`.

##### Primary Method

You may clone this repo and install directly. This is particularly
useful to use the CLI scripts or for development purposes.

    git clone https://github.com/dkainer/RWRtools.git
    cd RWRtools
    git checkout packageAttempt1
    R
    devtools::install()

From a clean environment this may take a while (\~20 min).

##### Secondary Method (install as an R package directly)

*Note: this method is not yet available for general use as repository is
private.* You can install the released version of RWRtools from
[GitHub](https://github.com/dkainer/RWRtools/) with:

``` r
devtools::install_github("dkainer/RWRtools")
```

*Note: add in details about where to find the CLI scripts in
/bin/library ??*

## Running RWRtools

### General:

RWRtools enables RandomWalk with Restart (RWR) on both homogenous and
heterogeneous multiplex networks. A heterogeneous network is used when
integrating multiple network sources; for example in building a
multiplex network with a gene-to-gene network and a disease-to-disease
networks and combining them by defining a gene-to-disease network which
serves as the bi-partite edges. RWRtools provides functions for both
creating the muliplex networks and running RWR.

### Usage Options:

The tools provided by RWRtools can be used either directly in R or by
use of command line scripts. The R functions follow the convention of
`RWRtools::RWR_func` such as `RWRtools::RWR_make_multiplex`. View help
with `?RWRtools::RWR_make_multiplex`. The command line scripts are
available in `./inst/scripts` and can be used with `Rscript` such as
`Rscript run_make_multiplex.R`. Run `Rscript run_make_multiplex.R -h` to
view the help. You can use these scripts from any location, but remember
to either use complete paths or paths local to where you are running
when applicable.

### Inital Step:

The first step in RWRtools is to build the RData object that represents
the multiplex network using `RWR_make_multiplex`. This function requires
an `flist` (a **f**ile **list**) input file which represents the set of
networks to create the multiplex object. Each row in the flist is a
triple defining the network: {file\_path, name, group}. In a homogenous
network the group is all `1`. In a heterogeneous network, one set of
networks will use `1` (e.g. gene-to-gene), the other will use `2`
(e.g. disease-to-disease), and `3` for the connecting network
(e.g. gene-to-disease). An example flist for a homegenous networks looks
like (seperated by any of the following delimiters `,\t |;`):

|   **file\_path**   | **name**  | **group** |
|:------------------:|:---------:|:---------:|
| /path/to/file1.txt |    PPI    |     1     |
| /path/to/file2.txt | Co-Domain |     1     |

At this stage you also define values for delta and lambda. **Delta**
sets the probability to change between layers at the next step. If delta
= 0, the particle will always remain in the same layer after a
non-restart iteration. On the other hand, if delta = 1, the particle
will always change between layers, therefore not following the specific
edges of each layer. The default is 0.5. Note delta must be greater than
0 and less than or equal to 1.

**Lambda** is for heterogeneous networks only. When building a
heterogeneous network (i.e. multiple layer groups connected with
bipartite links), the walker can jump between layer groups with
probability = lambda when it is at a node with a bipartite link. If
lambda=1 then walker will oscillate between groups every time it is at a
node with a bipartite link. Default is 0.5.

Please note that for large networks or a large number of networks this
function may take a long time.

This function will not return anything, it will save the relevant
objects (the multiplex object *mpo*, adjacency matrix, and normalized
adjacency matrix) to file to be used in subsequent functions.

When using the CLI script, remember to use complete paths or paths local
to where you run `scripts/run_make_multiplex.R` in your `flist`.

Examples

-   **Running in R**

        RWRtools::RWR_make_multiplex(
          flist="./inst/example_data/flist.tsv",
          delta=0.25,
          lambda=0.75, 
          output="./outdir/myExampleNetwork.Rdata"
        )

-   **Running CLI**

        Rscript run_make_multiplex.R 
          --flist tests/testSTRINGDB/flist.tsv \
          --delta 0.25,
          --lambda 0.75,
          --out ./outdir/myExampleNetwork.Rdata

### Next Steps:

The choice of the next script depends on the type of analysis desired.
RWRtools provides several different workflows outlined below.

#### RWR\_CV.R

*RWR Cross Validation* performs K-fold cross validation on a single gene
set, finding the RWR rank of the left-out genes. Can choose between
three modes: (1) leave-one-out `loo` to leave only one gene from the
gene set out and find its rank, (2) cross-validation `kfold` to run
k-fold cross-validation for a specified value of *k*, or (3) singletons
`singletons` to use a single gene as a seed and find the rank of all
remaining genes.

-   **Input:** Pre-calculated interaction network (using
    `RWR_make_multiplex.R`), and a single geneset.
-   **Output:** Table/dataframe with the ranking of each gene in the
    gene set when left out, as well as AUPRC and AUROC curves.

Examples

-   **Running in R**

        RWRtools::RWR_CV(
          dataPath="./inst/example_data/string_interactions.Rdata",
          genesetPath="./inst/example_data/geneset1.txt",
          outdirPath="./outdir")

-   **Running CLI**

        Rscript run_cv.R
          --data ./inst/example_data/string_interactions.Rdata \
          --geneset ./inst/example_data/geneset1.txt \
          -o ./outdir

#### RWR\_LOE.R

*RWR Lines of Evidence* has two possible functions. Given one geneset of
seeds, rankings for all other genes in the network will be returned.
Given a second geneset of genes to be queried, rankings for just the
genes in that geneset will be returned. This can be used to build
multiple lines of evidence from the various input networks to relate the
two gene sets.

-   **Input:** Pre-calculated interaction network (using
    `RWR_make_multiplex`), and one or two genesets.
-   **Output:** Table/dataframe with a ranking of non-seed genes (either
    the rest of the genes in the network if only one input geneset is
    used, or just the genes in the second geneset if one is provided).

Examples

-   **Running in R**

        RWRtools::RWR_LOE(
          data="./inst/example_data/string_interactions.Rdata",
          seed_geneset="./inst/example_data/geneset1.txt",
          tau="1.0,1.0",
          outdir="./outdir")

-   **Running CLI**

        Rscript run_loe.R 
          --data            ./inst/example_data/string_interactions.Rdata \
          --seed_geneset    ./inst/example_data/geneset1.txt \
          --tau             "1.0,1.0" \
          -o                ./outdir

#### RWR\_netscore.R

*RWR Net Score* performs a network intersect between an input network
(`network`) and a gold truth network (`gold`), e.g. the GO network. It
will score the strength of the intersect with multiple metrics.

-   **Input:**
    -   A gold standard network as reference.
    -   Another network to compare to the gold standard.
-   **Output:** A table containing multiple metrics including the edge
    intersect between the input network and the gold standard network.

Examples

-   **Running in R**

        RWRtools::RWR_netscore(
          gold="./inst/example_data/netscore/combined_score-random-gold.tsv",
          network="./inst/example_data/netscore/combined_score-random-test.tsv",
          outdir="./outdir")

-   **Running CLI**

        Rscript run_loe.R \
          --gold    ./inst/example_data/netscore/combined_score-random-gold.tsv \
          --network ./inst/example_data/netscore/combined_score-random-test.tsv \
          --outdir  ./outdir

#### RWR\_shortestpaths.R

Find shortest paths between genes in gene sets. Given a single gene set,
find the shortest paths between the genes in that gene set. Given two
gene sets, find the shortest paths for pairs of genes between gene sets.

-   **Input:**
    -   Pre-calculated interaction network (`data`). The layers will be
        flattened into a single network to find the shortest paths.
    -   A file in TSV format containing genes of interest
        (`source-geneset`).
    -   Optional second file in TSV format containing genes of interest
        (`target-geneset`) to find pairs of paths to the
        `source-geneset`.
-   **Output:** Edge list table.

Examples

-   **Running in R**

        RWRtools::RWR_ShortestPaths(
            data="./inst/example_data/string_interactions.Rdata",
            source_geneset="./inst/example_data/geneset1.txt",
            target_geneset="./inst/example_data/geneset2.txt",
            outdir="./outdir"
        )

-   **Running CLI**

        Rscript run_shortestpaths.R \
            --data ./inst/example_data/string_interactions.Rdata \
            --source-geneset ./inst/example_data/geneset1.txt \
            --target-geneset ./inst/example_data/geneset2.txt \
            -o ./outdir
