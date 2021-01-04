#! /bin/bash

conda activate bactopia
bactopia datasets --cpus 10 --species "Haemophilus influenzae" --include_genus
bactopia prepare ~/pojects/nthi/data/fastqs > fastqs.txt
mkdir bactopia
cd bactopia
bactopia --fastqs ../fastqs.txt --datasets ../datasets/ --species "Haemophilus influenzae" --genome_size median --cpus 4 -profile slurm
rm -rf work/
cd ..
bactopia search "Haemophilus influenzae" --prefix hflu --min_coverage 40 --genome_size 1860196 --min_read_length 75
mkdir bactopia-ena
cd bactopia-ena
bactopia --accessions ../hflu-accessions.txt --datasets ../datasets/ --species "Haemophilus influenzae" --genome_size median --cpus 4 -profile slurm --disable_auto_variants
bactopia tools summary --bactopia bactopia/ --prefix nthi
