#!/usr/bin/env bash
# Usage:
#   $ cd RWRtoolkit
#   $ bash quickstart.sh
# Or:
#   $ bash path/to/quickstart.sh path/to/RWRtoolkit

set -e

export PATH_TO_RWRtoolkit="${1:-$PWD}"
echo "PATH_TO_RWRtoolkit=$PATH_TO_RWRtoolkit"

rwr_cv="$PATH_TO_RWRtoolkit/inst/scripts/run_cv.R"
rwr_loe="$PATH_TO_RWRtoolkit/inst/scripts/run_loe.R"
rwr_make_multiplex="$PATH_TO_RWRtoolkit/inst/scripts/run_make_multiplex.R"
rwr_netscore="$PATH_TO_RWRtoolkit/inst/scripts/run_netscore.R"
rwr_shortestpaths="$PATH_TO_RWRtoolkit/inst/scripts/run_shortestpaths.R"

out_dir="$PATH_TO_RWRtoolkit/tmp"

echo "Checking help ..."
Rscript $rwr_make_multiplex -h > /dev/null
Rscript $rwr_loe -h > /dev/null
Rscript $rwr_cv -h > /dev/null
Rscript $rwr_netscore -h > /dev/null
Rscript $rwr_shortestpaths -h > /dev/null
echo "Done."

echo "Creating multiplex ..."
Rscript $rwr_make_multiplex \
    --flist tests/testSTRINGDB/flist.tsv \
    --out $out_dir/multiplex.rdata \
    &> /dev/null
echo "Done."

echo "Running RWR_LOE ..."
Rscript $rwr_loe \
    --data tmp/multiplex.rdata \
    --seed_geneset tests/testSTRINGDB/geneset1.tsv \
    --outdir $out_dir \
    &> /dev/null
echo "Done."

echo "Running RWR_CV ..."
Rscript $rwr_cv \
    --data tmp/multiplex.rdata \
    --geneset tests/testSTRINGDB/geneset1.tsv \
    --outdir $out_dir \
    --plot \
    &> /dev/null
echo "Done."

echo "Running RWR_netscore ..."
Rscript $rwr_netscore \
    --gold tests/testSTRINGDB/layers/combined_score-random-gold.tsv \
    --network tests/testSTRINGDB/layers/combined_score-random-test.tsv \
    --outdir $out_dir \
    &> /dev/null
echo "Done."

echo "Running RWR_shortestpaths ..."
Rscript $rwr_shortestpaths \
    --data tmp/multiplex.rdata \
    --source-geneset tests/testSTRINGDB/geneset1.tsv \
    --target-geneset tests/testSTRINGDB/geneset2.tsv \
    --outdir $out_dir \
    &> /dev/null
echo "Done."

exit 0
