#!/bin/bash
# Run the superin pipeline on the small test data set (2 samples) located in this directory.
# Usage: bash run.sh [extra sequana_superin options, e.g. --apptainer-prefix PATH]
set -e

DATA=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

sequana_superin \
    --input-directory "$DATA/fastq" \
    --input-bam-pattern "$DATA/bam" \
    --samplesheet "$DATA/samplesheet.csv" \
    --plasmids-directory "$DATA/references" \
    --SI-reference-file "$DATA/references/SI.fasta" \
    --working-directory superin_test \
    --monitor \
    --force "$@"

cd superin_test
bash superin.sh
