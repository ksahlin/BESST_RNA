#!/bin/bash -l

#SBATCH -A b2010042
#SBATCH -J BESSTtest0.7
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem72GB
#SBATCH -t 4-00:00:00
#SBATCH --mail-user kristoffer.sahlin@scilifelab.se
#SBATCH --mail-type=ALL




Project=b2010042
Contigs=CLCh003a

if [ ! "$Contigs" ] ; then
    echo "Usage: $0 contig-pool-name"   
    exit 1
else
    echo "$0: scaffolding WGS contigs $Contigs with BESST"
fi
OUT="../output"
#OUT="$POOL.180h+300h+650h+mp-BESSTv0.5-002"
mkdir -p $OUT
#module load python/2.7.1 
echo "saving results to directory $OUT"

CMD="time python2.7 Main.py 4 \
        -c /proj/$Project/assembly/CLC/CLCh003/CLCh003a/$Contigs.fa \
        -f /proj/$Project/assembly/mappings/Z4006s466.h.ipe.180-merged-ALL-to-$Contigs.bam \
           /proj/$Project/assembly/mappings/Z4006s466.h.ipe.300-merged-ALL-to-$Contigs.bam \
           /proj/$Project/assembly/mappings/Z4006s466.h.ipe.650-merged-ALL-to-$Contigs.bam \
           /proj/$Project/assembly/mappings/Z4006c01.g.imp.mp2k-merged-ALL-to-$Contigs.bam \

        -e 3 3 3 3 \
        -d 1\
        -g 0\
        -o $OUT "

echo "running command: $CMD"

$CMD

