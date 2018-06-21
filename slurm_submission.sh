
snakemake --use-conda --cluster-config slurm_cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time} --mem {cluster.mem}" -j 200 -p
