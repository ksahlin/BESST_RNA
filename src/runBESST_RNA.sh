#!/bin/bash -l

#SBATCH -A b2010042
#SBATCH -J BESST_RNA
#SBATCH -p node
#SBATCH -n 8
#SBATCH -C mem72GB
#SBATCH -t 1:00:00
#SBATCH --mail-user kristoffer.sahlin@scilifelab.se
#SBATCH --mail-type=ALL


OUT="../output"

mkdir -p $OUT
module add python/2.7.1 
echo "saving results to directory $OUT"

CMD="time python2.7 Main.py 1 \
        -c /proj/b2010042/assembly/ASSEMBLYLOCK_2012_JULY/masterassembly/picea_abies.master.july2012.fa \
        -f /proj/b2010042/est/alignment/bwa/ASSEMBLYLOCK_2012_JULY_master/october2012/masterassembly/fq.normalized_K25_C80_pctSD100.bam \

        -e 3\
        -T 20000\
        -k 500\
        -d 1\
        -g 0\
        -z 1000\
        -o $OUT "

echo "running command: $CMD"

$CMD

