#!/bin/bash
#SBATCH --account smontgom
#SBATCH --job-name=snakemake_slurm_submission
#SBATCH --output=snakemake_slurm_submission
#
#SBATCH --ntasks=1
#SBATCH --time=168:00:00
#SBATCH --mem=1G

snakemake --use-conda --cluster-config slurm_cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time} --mem {cluster.mem}" -j 200 -p
