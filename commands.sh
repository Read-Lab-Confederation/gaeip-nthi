#! /bin/bash

# Process Genomes With Bactopia

### Install
```
conda create -y -n bactopia -c conda-forge -c bioconda bactopia
conda activate bactopia
```

### Setup
```
bactopia datasets --cpus 10 --species "Haemophilus influenzae" --include_genus
bactopia prepare ~/pojects/nthi/data/fastqs > fastqs.txt
bactopia pull --default --include_tools
```

### Run Local Samples
```
mkdir bactopia
cd bactopia
bactopia --fastqs ../fastqs.txt \
         --datasets ../datasets/ \
         --species "Haemophilus influenzae" \
         --genome_size median \
         --cpus 4 \
         -profile slurm
rm -rf work/
```

### Run completed genome
```
mkdir GA81666-completed
cd GA81666-completed
bactopia --assembly ../../../../../data/completed-genomes/GA81666.fasta \
         --sample GA81666 \
         --datasets ../datasets/ \
         --species "Haemophilus influenzae" \
         --genome_size median \
         --cpus 4 \
         -profile slurm
rm -rf work/
```

### Process Public Samples
This was done in December 2020
```
cd ..
bactopia search "Haemophilus influenzae" --prefix hflu --min_coverage 40 --genome_size 1860196 --min_read_length 75
mkdir bactopia-ena
cd bactopia-ena
bactopia --accessions ../hflu-accessions.txt \
         --datasets ../datasets/ \
         --species "Haemophilus influenzae" \
         --genome_size median \
         --cpus 4 \
         -profile slurm \
         --disable_auto_variants
rm -rf work/
```

### Summary of Genomes
```
bactopia tools summary --bactopia bactopia/ --prefix nthi
mkdir ../results/bactopia-summary
cp -r bactopia-tools/summary/nthi/* ../results/bactopia-summary/
```
This generates a list samples to exclude from further analysis.

### FastANI
Used the GA81666 sample as a reference for FastANI. 
```
bactopia tools fastani --bactopia bactopia/ \
                       --species "Haemophilus influenzae" \
                       --reference ../../../../data/completed-genomes/GA81666.fasta \
                       --skip_pairwise \
                       --local_reference_only \
                       --exclude ../results/bactopia-summary/nthi-exclude.txt \
                       --cpus 20 \
                       -profile slurm
```

### Representative Set
```
python3 bin/generate-representative-set.py \
    results/bactopia-summary/nthi-report.txt \
    results/bactopia-summary/nthi-exclude.txt \
    results/fastani.tsv \
    data/st164-c1.txt \
    data/st1714-c2.txt
mv excludes.txt representatives.txt results/
```
This produces two files `excludes.txt` and `representatives.txt`

### Mashtree
```
bactopia tools mashtree --bactopia bactopia/ \
                        --species "Haemophilus influenzae" \
                        --exclude ../results/excludes.txt \
                        --max_time 720 \
                        -profile slurm
```

### PIRATE
#### C1 (ST164) only
```
bactopia tools pirate --bactopia bactopia/ \
                      --prefix st164-c1 \
                      --include st164-c1.txt \
                      --cpus 20 \
                      -profile slurm 
```

#### C2 (ST1714) only
```
bactopia tools pirate --bactopia bactopia/ \
                      --prefix st1714-c2 \
                      --include st1714-c2.txt \
                      --cpus 20 \
                      -profile slurm 
```

#### C1+C2
```
bactopia tools pirate --bactopia bactopia/ \
                      --prefix c1-c2 \
                      --include c1-c2.txt \
                      --cpus 20 \
                      -profile slurm 
```

#### Georgia Samples
```
bactopia tools pirate --bactopia bactopia/ \
                      --prefix ga-samples \
                      --include ga-samples.txt \
                      --cpus 20 \
                      -profile slurm 
```


#### Representative Set of Samples
```
bactopia tools pirate --bactopia bactopia/ \
                      --prefix representative-set \
                      --include ../results/representatives.txt \
                      --cpus 39 \
                      -profile slurm 
```

### Insertion Sequences
#### IS1016
```
bactopia tools ismapper --bactopia bactopia/ \
                        --insertions datasets/species-specific/haemophilus-influenzae/optional/insertion-sequences/IS1016.fasta \
                        --reference GA81666-completed/GA81666/annotation/GA81666.gbk \
                        --prefix IS1016 \
                        --exclude ../results/excludes.txt \
                        -profile slurm
```

### Capsule
```
bactopia tools hicap --bactopia bactopia/ \
                     --exclude ../results/excludes.txt \
                     -profile slurm
```
