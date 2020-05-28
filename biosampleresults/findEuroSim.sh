#!/bin/sh
#
#SBATCH -J findEuro # A single job name for the array
#SBATCH --ntasks=100
#SBATCH -N 3
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/cat7ep/slurmOut/05282020/findEuro.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/05282020/findEuro.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load samtools
#testing with UA_Cho_15_65

###This script will run through each of the dest_mapped European mel samples,
### then identify which ones have a higher contamination of sim data

cd /project/berglandlab/dest_mapped/

for dir in *_*; do
  #get only the 2L, 2R, 3L, 3R columns from mel.bam and sim.bam
  samtools idxstats ./$dir/mel.bam | \
    grep -E "2L|2R|3L|3R" | grep -vE "mapped|sim" \
    > /scratch/cat7ep/interDirEuro/"$dir".mel.txt

  samtools idxstats ./$dir/sim.bam | \
    grep -E "sim_2L|sim_2R|sim_3L|sim_3R" | grep -v "het" \
    > /scratch/cat7ep/interDirEuro/"$dir".sim.txt

  #sum up the mel and sim
  mel=$( cat /scratch/cat7ep/interDirEuro/"$dir".mel.txt | awk -F '\t' '{sum+=$3} END {print sum}' )
  sim=$( cat /scratch/cat7ep/interDirEuro/"$dir".sim.txt | awk -F '\t' '{sum+=$3} END {print sum}' )
  total="$(($mel+$sim))"
  #find proportion
  prop="$(echo $sim $total | awk '{ print $1/$2 }')"

  echo "$dir $prop" >> /scratch/cat7ep/simCline/biosampleresults/propsOut.txt
done
