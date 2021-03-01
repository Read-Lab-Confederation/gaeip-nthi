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
bactopia prepare ~/pojects/nthi/data/cdc/fastqs >> fastqs.txt
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
                       --exclude excludes.txt \
                       --cpus 20 \
                       -profile slurm
```

### Mashtree
```
bactopia tools mashtree --bactopia bactopia/ \
                        --species "Haemophilus influenzae" \
                        --exclude excludes.txt \
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

#### EIP Samples
```
bactopia tools pirate --bactopia bactopia/ \
                      --prefix ga-samples \
                      --include ga-samples.txt \
                      --cpus 20 \
                      -profile slurm 
```


### Insertion Sequences
#### IS1016
```
bactopia tools ismapper --bactopia bactopia/ \
                        --insertions datasets/species-specific/haemophilus-influenzae/optional/insertion-sequences/IS1016.fasta \
                        --reference GA81666-completed/GA81666/annotation/GA81666.gbk \
                        --prefix IS1016 \
                        --exclude excludes.txt \
                        -profile slurm
```
