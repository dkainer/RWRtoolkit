
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RWRtoolkit

<!-- badges: start -->
<!-- badges: end -->

RWRtoolkit enables easy use of RandomWalk with Restart on multiplex
networks. These functions are an extension to the
[RandomWalkRestartMH](https://github.com/alberto-valdeolivas/RandomWalkRestartMH)
R package. Also provided are scripts for use as command line tools.

## Installation

##### Dependencies

Installation of this R package requires R and r-devtools. If you use
prefer the use of conda you can create the base environment with
`conda create --name r-RWRtoolkit -c conda-forge r-base r-devtools r-irkernel`.
You can also install devtools from within a base R environment with
`install.packages("devtools")`.

**Note**: Depending on how your packages were installed, you may run
into an SSL issue when attempting to install devtools. This is due to
the installation of `gert`, which requires an installation of `libgit2`
(installable via the [binaries](https://libgit2.org/),
[conda](https://anaconda.org/conda-forge/libgit2),
[homebrew](https://formulae.brew.sh/formula/libgit2),
[yum](https://yum-info.contradodigital.com/view-package/epel/libgit2/),
or package manager of your choice).

##### Primary Method

You may clone this repo and install directly. This is particularly
useful to use the CLI scripts or for development purposes.

    git clone https://github.com/dkainer/RWRtoolkit.git
    cd RWRtoolkit
    R
    devtools::install()

From a clean environment this may take a while (\~20 min).

##### Secondary Method (install as an R package directly)

You can install the released version of RWRtoolkit from
[GitHub](https://github.com/dkainer/RWRtoolkit/) with:

``` r
devtools::install_github("dkainer/RWRtoolkit")
```

## Running RWRtoolkit

### Loading RWRtoolkit:

RWRtoolkit can be run as either an R package or a command line tool
depending on your preferences.

-   **R Package:** Simply loading the library with the `library`
    function in R loads RWRtoolkit:

    ``` r
    library(RWRtoolkit)
    ```

-   **Command Line Tool:**  
    If you have downloaded the code via GitHub, you can access the
    command line script code by navigating to the
    `RWRtoolkit/inst/scripts` directory.

    If you have downloaded the code via `devtools::install_github` open
    an R session and type:

    ``` r
    library(RWRtoolkit)
    .libPaths()
    ```

    Which ought to output a path similar to:

        /Library/Frameworks/R.framework/Versions/4.0/Resources/library/

    This is the directory in which your installed R libraries exist.

    From the above directory (hereby referred to as
    `<LIBPATHS_DIRECTORY>` ), the script files can be found on the path:

        <LIBPATHS_DIRECTORY>/RWRtoolkit/scripts

    Note: the paths are not the same as the GitHub repository due to the
    `devtools::install` function???s lifting of all directories within the
    `inst` directory during the build/installation phase.

    From the above path, all scripts can be accessed as:

    ``` bash
    Rscript <LIBPATHS_DIRECTORY>/RWRtoolkit/scripts/run_loe.R 
      --data            <LIBPATHS_DIRECTORY>/RWRtoolkit/example_data/string_interactions.Rdata \
      --seed_geneset    <LIBPATHS_DIRECTORY>/RWRtoolkit/example_data/geneset1.txt \
      --tau             "1.0,1.0" \
      -o                ./outdir
    ```

### General:

RWRtoolkit enables RandomWalk with Restart (RWR) on both homogenous and
heterogeneous multiplex networks. A heterogeneous network is used when
integrating multiple network sources; for example in building a
multiplex network with a gene-to-gene network and a disease-to-disease
networks and combining them by defining a gene-to-disease network which
serves as the bi-partite edges. RWRtoolkit provides functions for both
creating the muliplex networks and running RWR.

### Usage Options:

The tools provided by RWRtoolkit can be used either directly in R or by
use of command line scripts. The R functions follow the convention of
`RWRtoolkit::RWR_func` such as `RWRtoolkit::RWR_make_multiplex`. View
help with `?RWRtoolkit::RWR_make_multiplex`. The command line scripts
are available in `./inst/scripts` and can be used with `Rscript` such as
`Rscript run_make_multiplex.R`. Run `Rscript run_make_multiplex.R -h` to
view the help. You can use these scripts from any location, but remember
to either use complete paths or paths local to where you are running
when applicable.

### Initial Step:

The first step in RWRtoolkit is to build the RData object that
represents the multiplex network using `RWR_make_multiplex`. This
function requires an `flist` (a **f**ile **list**) input file which
represents the set of networks to create the multiplex object. Each row
in the flist is a triple defining the network: {file_path, name, group}.
In a homogeneous network the group is all `1`. In a heterogeneous
network, one set of networks will use `1` (e.g.??gene-to-gene), the other
will use `2` (e.g.??disease-to-disease), and `3` for the connecting
network (e.g.??gene-to-disease). An example flist for a homogeneous
networks looks like (separated by any of the following delimiters
`,\t |;`):

|   **file_path**    | **name**  | **group** |
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
heterogeneous network (i.e.??multiple layer groups connected with
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

-   **Running in R** The below code assumes an R session was initialized
    from within the `inst` directory of RWRtoolkit. Output will be
    within the `RWRtoolkit/inst` directory. (This is necessary due to
    the files within `flist.tsv` having relative paths)

    ``` r
    RWRtoolkit::RWR_make_multiplex(
      flist="./example_data/flist.tsv",
      delta=0.25,
      lambda=0.75, 
      output="./outdir/myExampleNetwork.Rdata"
    )
    ```

-   **Running CLI** If running the code from the cloned GitHub
    repository, the below code ought to be run from within the `inst`
    directory. If running from the `devtools::install_github` method,
    the below code ought to be run from with the RWRtoolkit directory
    located at `<LIBPATHS_DIRECTORY>/RWRtoolkit`. Output will be saved
    to your home directory.

    ``` bash
    Rscript scripts/run_make_multiplex.R \
      --flist example_data/flist.tsv \
      --delta 0.25 \
      --lambda 0.75 \
      --out ~/RWRtoolkitOutput/myExampleNetwork.Rdata
    ```

### Next Steps:

The choice of the next script depends on the type of analysis desired.
RWRtoolkit provides several different workflows outlined below.

#### RWR_CV.R

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

    ``` r
    # Can be run from anywhere so long as RWRtoolkit is installed. 
    extdata.dir <- system.file("example_data", package="RWRtoolkit")

    string.interactions.fp <- paste(extdata.dir, "string_interactions.Rdata", sep='/')
    geneset.path <- paste(extdata.dir, 'geneset1.tsv', sep='/')
    outdir.path <- '~/RWRtoolkitOutput/'

    RWRtoolkit::RWR_CV(
      dataPath = string.interactions.fp ,
      genesetPath = geneset.path,
      outdirPath = outdir.path)
    ```

-   **Running CLI** If running the code from the cloned GitHub
    repository, the below code ought to be run from within the `inst`
    directory. If running from the `devtools::install_github` method,
    the below code ought to be run from with the RWRtoolkit directory
    located at `<LIBPATHS_DIRECTORY>/RWRtoolkit`. Output will be saved
    to your home directory.

    ``` bash
    Rscript ./scripts/run_cv.R \
      --data ./example_data/string_interactions.Rdata \
      --geneset ./example_data/geneset1.tsv \
      -o ./outdircli
    ```

#### RWR_LOE.R

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

    ``` r
    # Can be run from anywhere so long as RWRtoolkit is installed. 
    extdata.dir <- system.file("example_data", package="RWRtoolkit")

    string.interactions.fp <- paste(extdata.dir, "string_interactions.Rdata", sep='/')
    geneset.path <- paste(extdata.dir, 'geneset1.tsv', sep='/')
    outdir.path <- '~/RWRtoolkitOutput/'

    RWRtoolkit::RWR_LOE(
      data= string.interactions.fp,
      seed_geneset= geneset.path,
      tau = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
      outdir= outdir.path )
    ```

-   **Running CLI**

    ``` bash
    Rscript scripts/run_loe.R \
      --data            ./example_data/string_interactions.Rdata \
      --seed_geneset    ./example_data/geneset1.tsv \
      --tau             "1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0" \
      -o                ./outdir
    ```

#### RWR_netscore.R

*RWR Net Score* performs a network intersect between an input network
(`network`) and a gold truth network (`gold`), e.g.??the GO network. It
will score the strength of the intersect with multiple metrics.

-   **Input:**
    -   A gold standard network as reference.
    -   Another network to compare to the gold standard.
-   **Output:** A table containing multiple metrics including the edge
    intersect between the input network and the gold standard network.

Examples

-   **Running in R**

    ``` r
    # Can be run from anywhere so long as RWRtoolkit is installed. 
    extdata.dir <- system.file("example_data", package="RWRtoolkit")

    gold.fp <- paste(extdata.dir, "netscore/combined_score-random-gold.tsv", sep='/')
    network.fp <- paste(extdata.dir, "netscore/combined_score-random-test.tsv", sep='/')
    outdir.path <- "~/RWRtoolkitOutput/"

    RWRtoolkit::RWR_netscore(
      gold = gold.fp,
      network = network.fp,
      outdir = outdir.path)
    ```

-   **Running CLI**

    ``` bash
    Rscript scripts/run_netscore.R \
      --gold    ./example_data/netscore/combined_score-random-gold.tsv \
      --network ./example_data/netscore/combined_score-random-test.tsv \
      --outdir  ./outdir
    ```

#### RWR_shortestpaths.R

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

    ``` r
    # Can be run from anywhere so long as RWRtoolkit is installed. 
    extdata.dir <- system.file("example_data", package="RWRtoolkit")

    string.interactions.fp <- paste(extdata.dir, "string_interactions.Rdata", sep='/')
    source.geneset.path <- paste(extdata.dir, 'geneset1.tsv', sep='/')
    target.geneset.path <- paste(extdata.dir, 'geneset1.tsv', sep='/')
    outdir.path <- '~/RWRtoolkitOutput/'



    RWRtoolkit::RWR_ShortestPaths(
        data = string.interactions.fp,
        source_geneset = source.geneset.path,
        target_geneset = target.geneset.path,
        outdir = outdir.path
    )
    ```

-   **Running CLI**

    ``` bash
    Rscript scripts/run_shortestpaths.R \
        --data ./example_data/string_interactions.Rdata \
        --source-geneset ./example_data/geneset1.tsv \
        --target-geneset ./example_data/geneset2.tsv \
        -o ./outdir
    ```
