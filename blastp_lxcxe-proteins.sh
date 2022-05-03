
# Script used to blast LxCxE-containing proteins against the
# Arabidopsis proteome using the Advanced Genomics Unit MAZORKA cluster


#PBS -q default
#PBS -l nodes=1:ppn=16,walltime=100:00:00,vmem=60gb
#PBS -N BLASTp-algae
#PBS -o blastp21jan.out
#PBS -e blastp21jan.err

cd $PBS_O_WORKDIR

module load ncbi-blast+/2.6.0


#!/bin/bash
set -euo pipefail

for filename in *.fasta
do
blastp \
-query $filename \
-db athdb \
-max_target_seqs 5 \
-outfmt "10 qseqid sseqid evalue bitscore" \
-evalue 1e-5 \
-out $filename.csv
done
