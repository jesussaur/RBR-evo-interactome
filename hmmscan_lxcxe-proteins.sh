# Script used to scan transcription associated domains
# in LxCxE-containing proteins against the TAPdb local database
# using the Advanced Genomics Unit MAZORKA cluster

#PBS -q default
#PBS -l nodes=1:ppn=16,walltime=10:00:00,vmem=60gb
#PBS -N HMMscan_TAPs
#PBS -o HMMERscan12-04-2021.out
#PBS -e HMMERscan12-04-2021.err
cd $PBS_O_WORKDIR

module load hmmer/3.1b1
# compresing TAPdb datbase
hmmpress TAPdb
# searching *.fasta sequences against TAPdb database
for filename in *.fasta
do
hmmscan \
--tblout $filename.csv \
TAPdb \
$filename
done
