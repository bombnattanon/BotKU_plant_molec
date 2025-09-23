#!/bin/bash

#SBATCH --job-name=GetOrganelle
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=0-6:00:00
#SBATCH -o stdout_getorganelle.log
#SBATCH -e stderr_getorganelle.log
#SBATCH --mail-user=nattanon.meep@live.ku.th
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

if [[ $# -eq 0 ]] ; then
    echo 'Error: Please provide a sample name as the first argument when running this script.'
    exit 0
fi

source activate /data/users/fscinnme/packages/miniforge3/envs/getorganelle

name=$1

echo $name

get_organelle_from_reads.py -1 /data/users/fscinnme/data/raw_reads/${name}/${name}*_1.fq.gz -2 /data/users/fscinnme/data/raw_reads/${name}/${name}*_2.fq.gz -o plastid_${name} --prefix ${name}_ -R 30 -k 35,85,127 -F embplant_pt
