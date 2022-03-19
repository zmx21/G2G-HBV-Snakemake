#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/src_gvcf/
#SBATCH --job-name 8_rsync
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/8_rsync.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/8_rsync.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
#SBATCH --account gr-fe
 

## copying files form scratch to main
rsync -avzh /scratch/rueger/gvcf_tmp/* /work/gr-fe/rueger/G2G-HBV/data/raw/gilead_20181126/wes_gvcf_to_plink_intermediate/

echo ended at `date`
exit
