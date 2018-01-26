#!/bin/bash
#SBATCH --ntasks 8 --nodes 1 -p short 
module load hmmer/3

CPUS=8
DB=Rhizopus_peps.aa.fasta
for file in *.pep
do
 out=$(basename $file .pep)".phmmer"
 tbl=$(basename $file .pep)".tbl"
 phmmer -E 1e-5 --cpu $CPUS --domtblout $tbl -o $out $file $DB
done
