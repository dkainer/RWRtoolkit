
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RWRtoolkit

<!-- badges: start -->

<!-- badges: end -->

RWRtoolkit enables easy use of RandomWalk with Restart on multiplex
networks. These functions are an extension to the
[RandomWalkRestartMH](https://github.com/alberto-valdeolivas/RandomWalkRestartMH)
R package. Also provided are scripts for use as command line tools.

## Installation

#### Dependencies

Installation of this R package requires R \>= 4.1.0 and devtools. If you
use prefer the use of conda you can create the base environment with
`conda create --name r-RWRtoolkit -c conda-forge "r>=4.1" "r-base>=4.1" r-devtools`
(`r-irkernel` is optional). You can also install devtools from within a
base R environment with `install.packages("devtools")`.

##### Installation Issues

###### Unable to access Bioconductor:

You may likely run into an issue with your R environment installing
packages via bioconductor:

``` r
devtools::install()
Error: Unknown remote type: bioc
  cannot open URL 'https://bioconductor.org/config.yaml'
```

To ensure the issue is a certificate issue, use another library to call
bioconductor:

``` r
httr::GET("https://bioconductor.org/config.yaml")
Error in curl::curl_fetch_memory(url, handle = handle) :
  SSL peer certificate or SSH remote key was not OK: [bioconductor.org] SSL certificate problem: self-signed certificate in certificate chain
```

This problem is an SSL error where you will need to update your SSL
certificate. To fix this, in your terminal, type the following:

``` bash
# 1. Get a certificate if you don't have one
curl -o ~/.ssh/cert.pem https://curl.se/ca/cacert.pem

# 2. Get a bioconductor specific certificate 
echo | openssl s_client -showcerts -servername bioconductor.org -connect bioconductor.org:443 2>/dev/null | openssl x509 -inform pem -outform pem > bioconductor_cert.pem

# 3. Append your bioconductor cert to your cert.pem file
cat bioconductor_cert.pem >> ~/.ssh/cert.pem

# 4. Add the cert path to your `.Renviron`
echo CURL_CA_BUNDLE=/Users/96v/.ssh/cert.pem > ~/.Renviron
```

In a newly restarted R environment, type:

``` r
httr::set_config(httr::config(cainfo = "/Users/96v/.ssh/cert.pem"))
response <- httr::GET("https://bioconductor.org/config.yaml")
print(response)
```

You ought to get an output similar to:

    Response [https://bioconductor.org/config.yaml]
      Date: 2024-08-14 10:57
      Status: 200
      Content-Type: <unknown>
      Size: 12.6 kB
    <BINARY BODY>

Now, `devtools::install()` ought to work.

###### devtools/r-devtools installation

It is possible you may run into issues installing `r-devtools` via conda
or `devtools` via R’s `install.packages()` function.

**textshaping** This might be due to a failure in the installation of
`textshaping`. `textshaping` requires the libraries `harfbuzz` and
`fribidi` libraries, yet uses the `pkg-config` command, which may be
external to your environment. There are multiple options for fixing
(linux/MacOS installation recommendations taken from R install.packages
ANTICONF): - Anaconda: conda install -c conda-forge pkg-config harfbuzz
fribidi - deb: libharfbuzz-dev libfribidi-dev (Debian, Ubuntu, etc) -
rpm: harfbuzz-devel fribidi-devel (Fedora, EPEL) - csw: libharfbuzz_dev
libfribidi_dev (Solaris) - brew: harfbuzz fribidi (OSX)

**libgit2** libgit2: Depending on how your packages were installed, you
may run into an SSL issue when attempting to install devtools. This is
due to the installation of `gert`, which requires an installation of
`libgit2` (installable via the [binaries](https://libgit2.org/),
[conda](https://anaconda.org/conda-forge/libgit2),
[homebrew](https://formulae.brew.sh/formula/libgit2),
[yum](https://yum-info.contradodigital.com/view-package/epel/libgit2/),
or package manager of your choice).

#### Package Installation

You may clone this repo and install directly. This is particularly
useful to use the CLI scripts or for development purposes.

    git clone https://github.com/dkainer/RWRtoolkit.git
    cd RWRtoolkit
    R
    devtools::install()

From a clean environment this may take a while (~20 min).

#### Secondary Method (install as an R package directly)

You can install the released version of RWRtoolkit from
[GitHub](https://github.com/dkainer/RWRtoolkit/) with:

``` r
devtools::install_github("dkainer/RWRtoolkit")
```

## Running RWRtoolkit

#### Loading RWRtoolkit:

RWRtoolkit can be run as either an R package or a command line tool
depending on your preferences.

- **R Package:** Simply loading the library with the `library` function
  in R loads RWRtoolkit:

  ``` r
  library(RWRtoolkit)
  ```

- **Command Line Tool:**  
  **If you have downloaded the code** via GitHub, you can access the
  command line script code by navigating to the
  `RWRtoolkit/inst/scripts` directory.

  **If you have downloaded the code** via `devtools::install_github`
  open an R session and type:

  ``` r
  library(RWRtoolkit)
  .libPaths()
  ```

  Which ought to output a path similar to:

      /Library/Frameworks/R.framework/Versions/4.0/Resources/library/

  This is the directory in which your installed R libraries exist.

  From the above directory (hereby referred to as `<LIBPATHS_DIRECTORY>`
  ), the script files can be found on the path:

      <LIBPATHS_DIRECTORY>/RWRtoolkit/scripts

  Note: the paths are not the same as the GitHub repository due to the
  `devtools::install` function’s lifting of all directories within the
  `inst` directory during the build/installation phase.

  From the above path, all scripts can be accessed as:

  ``` bash
  Rscript <LIBPATHS_DIRECTORY>/RWRtoolkit/scripts/run_loe.R 
    --data            <LIBPATHS_DIRECTORY>/RWRtoolkit/example_data/string_interactions.Rdata \
    --seed_geneset    <LIBPATHS_DIRECTORY>/RWRtoolkit/example_data/geneset1.txt \
    --tau             "1.0,1.0" \
    -o                ./outdir
  ```

#### Running

RWRtoolkit enables RandomWalk with Restart (RWR) on homogenous multiplex
networks. RWRtoolkit provides functions for both creating the muliplex
networks and running RWR.

##### Usage Options:

The tools provided by RWRtoolkit can be used either directly in R or by
use of command line scripts. The R functions follow the convention of
`RWRtoolkit::RWR_func` such as `RWRtoolkit::RWR_make_multiplex`. View
help with `?RWRtoolkit::RWR_make_multiplex`. The command line scripts
are available in `./inst/scripts` and can be used with `Rscript` such as
`Rscript run_make_multiplex.R`. Run `Rscript run_make_multiplex.R -h` to
view the help. You can use these scripts from any location, but remember
to either use complete paths or paths local to where you are running
when applicable.

##### Initial Step:

The first step in RWRtoolkit is to build the RData object that
represents the multiplex network using `RWR_make_multiplex`. This
function requires an `flist` (a **f**ile **list**) input file which
represents the set of networks to create the multiplex object. Each row
in the flist is a triple defining the network: {file_path, name, group}.
An example flist for a homogeneous networks looks like (separated by any
of the following delimiters `,\t |;`):

|   **file_path**    | **name**  |
|:------------------:|:---------:|
| /path/to/file1.txt |    PPI    |
| /path/to/file2.txt | Co-Domain |

At this stage you also define values for delta. **Delta** sets the
probability to change between layers at the next step. If delta = 0, the
particle will always remain in the same layer after a non-restart
iteration. On the other hand, if delta = 1, the particle will always
change between layers, therefore not following the specific edges of
each layer. The default is 0.5. Note delta must be greater than 0 and
less than or equal to 1. Please note that for large networks or a large
number of networks this function may take a long time.

This function will not return anything, it will save the relevant
objects (the multiplex object *mpo*, adjacency matrix, and normalized
adjacency matrix) to file to be used in subsequent functions.

When using the CLI script, remember to use complete paths or paths local
to where you run `scripts/run_make_multiplex.R` in your `flist`.

## RWRToolkit Examples

- **Running in R** The below code assumes an R session was initialized
  from within the `inst` directory of RWRtoolkit. Output will be within
  the `RWRtoolkit/inst` directory. (This is necessary due to the files
  within `flist.tsv` having relative paths)

  ``` r
  RWRtoolkit::RWR_make_multiplex(
    flist="./example_data/flist.tsv",
    delta=0.5,
    output="./RWRtoolkit_MPO_Output/myExampleNetwork.Rdata"
  )
  ```

- **Running CLI** If running the code from the cloned GitHub repository,
  the below code ought to be run from within the `inst` directory. If
  running from the `devtools::install_github` method, the below code
  ought to be run from with the RWRtoolkit directory located at
  `<LIBPATHS_DIRECTORY>/RWRtoolkit`. Output will be saved to your home
  directory.

  ``` bash
  Rscript scripts/run_make_multiplex.R \
    --flist example_data/flist.tsv \
    --delta 0.25 \
    --out ./RWRtoolkit_MPO_Output/myExampleNetwork.Rdata
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

- **Input:** Pre-calculated interaction network (using
  `RWR_make_multiplex.R`), and a single geneset.
- **Output:** Table/dataframe with the ranking of each gene in the gene
  set when left out, as well as AUPRC and AUROC curves.

Examples

- **Running in R**

  ``` r
  # Can be run from anywhere so long as RWRtoolkit is installed. 
  extdata.dir <- system.file("example_data", package="RWRtoolkit")

  string.interactions.fp <- paste(extdata.dir, "string_interactions.Rdata", sep='/')
  geneset.path <- paste(extdata.dir, 'geneset1.tsv', sep='/')
  outdir.path <- './RWRtoolkit_CV_Output/'

  RWRtoolkit::RWR_CV(
    data = string.interactions.fp ,
    genesetPath = geneset.path,
    outdirPath = outdir.path)
  ```

- **Running CLI** If running the code from the cloned GitHub repository,
  the below code ought to be run from within the `inst` directory. If
  running from the `devtools::install_github` method, the below code
  ought to be run from with the RWRtoolkit directory located at
  `<LIBPATHS_DIRECTORY>/RWRtoolkit`. Output will be saved to your home
  directory.

  ``` bash
  Rscript ./scripts/run_cv.R \
    --data ./example_data/string_interactions.Rdata \
    --geneset ./example_data/geneset1.tsv \
    -o ./RWRtoolkit_CV_Output/
  ```

#### RWR_LOE.R

*RWR Lines of Evidence* has two possible functions. Given one geneset of
seeds, rankings for all other genes in the network will be returned.
Given a second geneset of genes to be queried, rankings for just the
genes in that geneset will be returned. This can be used to build
multiple lines of evidence from the various input networks to relate the
two gene sets.

- **Input:** Pre-calculated interaction network (using
  `RWR_make_multiplex`), and one or two genesets.
- **Output:** Table/dataframe with a ranking of non-seed genes (either
  the rest of the genes in the network if only one input geneset is
  used, or just the genes in the second geneset if one is provided).

Examples

- **Running in R**

  ``` r
  # Can be run from anywhere so long as RWRtoolkit is installed. 
  extdata.dir <- system.file("example_data", package="RWRtoolkit")

  string.interactions.fp <- paste(extdata.dir, "string_interactions.Rdata", sep='/')
  geneset.path <- paste(extdata.dir, 'geneset1.tsv', sep='/')
  outdir.path <- './RWRtoolkitOutput_LOE/'

  RWRtoolkit::RWR_LOE(
    data= string.interactions.fp,
    seed_geneset= geneset.path,
    tau = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
    outdir= outdir.path )
  ```

- **Running CLI**

  ``` bash
  Rscript scripts/run_loe.R \
    --data            ./example_data/string_interactions.Rdata \
    --seed_geneset    ./example_data/geneset1.tsv \
    --tau             "1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0" \
    -o                ./RWRtoolkitOutput_LOE
  ```

#### RWR_netstats.R

*RWR Net Stats* performs offers a series of statistical methods for
extracting metrics for networks and multiplex layers. There are multiple
options within netstats:

- **Input:**
  - A multiplex object (from RWR_make_multiplex) or an flist.
  - A reference network: Optional (Depending on methods chosen)
  - A network of interest: Optional (Depending on methods chosen)
  - Network Scoring Metric: (“jaccard”, “overlap”, or “both”)
- **Output:** In R: a list containing tables of metrics flagged from
  input parameters. Files for each table can be saved by supplying an
  `output_dir`

Examples

- **Running in R**

  ``` r
  # Can be run from anywhere so long as RWRtoolkit is installed. 
  extdata.dir <- system.file("example_data", package="RWRtoolkit")

  mpo_path <- paste(extdata.dir, "string_interactions.Rdata", sep = "/") 
  gold.fp <- paste(extdata.dir, "netstat/combined_score-random-gold.tsv", sep='/')
  network.fp <- paste(extdata.dir, "netstat/combined_score-random-test.tsv", sep='/')
  outdir.path <- "~/RWRtoolkitOutput/"

  RWRtoolkit::RWR_netstats(
        data = mpo_path,
        network_1 = gold.fp,
        network_2 = network.fp,
        basic_statistics = T,
        scoring_metric = "both",
        pairwise_between_mpo_layer = T,
        multiplex_layers_to_refnet = T,
        net_to_net_similarity = T,
        calculate_tau_for_mpo = T,
        merged_with_all_edges = T,
        merged_with_edgecounts = T,
        calculate_exclusivity_for_mpo = T,
        outdir = "./",
        verbose = T
   )
  ```

- **Running CLI**

  ``` bash
  Rscript scripts/run_netstats.R \
    --data ./example_data/string_interactions.Rdata  \
    --network_1 ./example_data/netstat/combined_score-random-gold.tsv \
    --network_2 ./example_data/netstat/combined_score-random-test.tsv \
    --scoring_metric both \
    --outdir ./RWRtoolkitOutput_Netstats \
    --basic_statistics  \
    --pairwise_between_mpo_layer  \
    --multiplex_layers_to_refnet  \
    --net_to_net_similarity  \
    --calculate_tau_for_mpo  \
    --merged_with_all_edges  \
    --merged_with_edgecounts  \
    --calculate_exclusivity_for_mpo  \
    --verbose 
  ```

#### RWR_shortestpaths.R

Find shortest paths between genes in gene sets. Given a single gene set,
find the shortest paths between the genes in that gene set. Given two
gene sets, find the shortest paths for pairs of genes between gene sets.

- **Input:**
  - Pre-calculated interaction network (`data`). The layers will be
    flattened into a single network to find the shortest paths.
  - A file in TSV format containing genes of interest
    (`source-geneset`).
  - Optional second file in TSV format containing genes of interest
    (`target-geneset`) to find pairs of paths to the `source-geneset`.
- **Output:** Edge list table.

Examples

- **Running in R**

  ``` r
  # Can be run from anywhere so long as RWRtoolkit is installed. 
  extdata.dir <- system.file("example_data", package="RWRtoolkit")

  string.interactions.fp <- paste(extdata.dir, "string_interactions.Rdata", sep='/')
  source.geneset.path <- paste(extdata.dir, 'geneset1.tsv', sep='/')
  target.geneset.path <- paste(extdata.dir, 'geneset1.tsv', sep='/')
  outdir.path <- './RWRtoolkitOutput_SP/'



  RWRtoolkit::RWR_ShortestPaths(
      data = string.interactions.fp,
      source_geneset = source.geneset.path,
      target_geneset = target.geneset.path,
      outdir = outdir.path
  )
  ```

- **Running CLI**

  ``` bash
  Rscript scripts/run_shortestpaths.R \
      --data ./example_data/string_interactions.Rdata \
      --source_geneset ./example_data/geneset1.tsv \
      --target_geneset ./example_data/geneset2.tsv \
      -o ./RWRtoolkitOutput_SP/
  ```
