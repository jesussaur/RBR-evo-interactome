# Script used to find genome-wide orthologs of analyzed species
# 28 species

#PBS -q default
#PBS -l nodes=1:ppn=16,walltime=400:00:00,vmem=60gb
#PBS -N orthofinding-V4_BLAST
#PBS -o orthofinder30112020.out
#PBS -e orthofinder30112020.err

cd $PBS_O_WORKDIR

module load orthofinder/2.4.0
module load diamond/0.9.13
module load mafft/7.305b
module load ncbi-blast+/2.6.0
module load mcl/14-137


orthofinder -f . -t 16 -S blast -I 1.7 -p /scratch
