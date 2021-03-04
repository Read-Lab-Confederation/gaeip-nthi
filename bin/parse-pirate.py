#! /usr/bin/env python3
"""

"""
from collections import OrderedDict
CLUSTER_COLS = [
    'allele_name', 'gene_family', 'consensus_gene_name', 'consensus_product', 'threshold', 
    'alleles_at_maximum_threshold', 'number_genomes', 'average_dose', 'min_dose', 'max_dose',
    'genomes_containing_fissions', 'genomes_containing_duplications', 'number_fission_loci',
    'number_duplicated_loci', 'no_loci', 'products', 'gene_names', 'min_length(bp)', 'max_length(bp)',
    'average_length(bp)', 'cluster', 'cluster_order'
]


def read_samples(c):
    samples = []
    with open(c, 'rt') as c_fh:
        for line in c_fh:
            samples.append(line.strip().split('\t')[0])
    return samples


def read_pirate(pirate):
    samples = OrderedDict()
    clusters = OrderedDict()
    pirate_cols = None
    total_samples = None
    with open(pirate, 'rt') as pirate_fh:
        for line in pirate_fh:
            vals = line.rstrip().split('\t')
            cluster_info = []
            cluster_dict = None
            if not pirate_cols:
                pirate_cols = vals
                total_samples = len(pirate_cols) - len(CLUSTER_COLS)
            else:
                for i, val in enumerate(vals):
                    if i < len(CLUSTER_COLS):
                        # first cols are cluster info
                        cluster_info.append(val)
                    else:
                        # remaining are sample presence/absence
                        if not cluster_dict:
                            # Add 
                            cluster_dict = dict(zip(CLUSTER_COLS, cluster_info))
                            clusters[cluster_dict['gene_family']] = cluster_dict

                        if cluster_dict['gene_family'] not in samples:
                            # Init sample
                            samples[cluster_dict['gene_family']] = {'samples': [], 'cluster': [], 'total': 0, 'percent': 0}

                        if val:
                            # add sample name (GA54827) and cluster name (GA54827_00366)
                            samples[cluster_dict['gene_family']]['samples'].append(pirate_cols[i])
                            samples[cluster_dict['gene_family']]['cluster'].append(val)
                samples[cluster_dict['gene_family']]['total'] = len(samples[cluster_dict['gene_family']]['samples'])
                samples[cluster_dict['gene_family']]['percent'] = samples[cluster_dict['gene_family']]['total'] / total_samples
    return(clusters, samples, total_samples)


if __name__ == '__main__':
    import argparse as ap
    
    import os
    import sys

    parser = ap.ArgumentParser(
        prog='parse-pirate',
        conflict_handler='resolve',
        description=(
            f'parse-pirate.py - Parse PIRATE output and give details about accessory genome.'
        )
    )

    parser.add_argument('pirate', metavar="STR", type=str, help='Presence/abscence file from PIRATE analyis')
    parser.add_argument('c1', metavar="STR", type=str, help='Text file with names of C1 samples')
    parser.add_argument('c2', metavar="STR", type=str, help='Text file with names of C2 samples')
    parser.add_argument('--full_detail', action='store_true', help='Only print for detailed output.')
    parser.add_argument('--exclusive', action='store_true', help='Gene is found in only one cluster (C1 or C2).')
    parser.add_argument('--accessory_only', action='store_true', help='Only print info about accessory genes.')
    parser.add_argument('--verbose', action='store_true', help='Print debug related text.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Make sure files exist
    error = False
    if not os.path.exists(args.pirate):
        print(f"cannot access '{args.pirate}': No such file or directory\n", file=sys.stderr)
        error = True
    if not os.path.exists(args.c1):
        print(f"cannot access '{args.c1}': No such file or directory\n", file=sys.stderr)
        error = True
    if not os.path.exists(args.c2):
        print(f"cannot access '{args.c2}': No such file or directory\n", file=sys.stderr)
        error = True

    if error:
        print("Please make sure the correct path is given.", file=sys.stderr)
        sys.exit(1)

    c1 = read_samples(args.c1)
    c2 = read_samples(args.c2)
    clusters, samples, total_samples = read_pirate(args.pirate)
    total_neither = total_samples - len(c1) - len(c2)

    colnames = [
        'gene_family', 'percent_all', 'percent_c1', 'percent_c2', 'percent_neither', 'info', 'product'
    ]

    if args.full_detail:
        colnames = [
            'gene_family', 'percent_all', 'percent_c1', 'percent_c2', 'percent_neither', 'info', 'product', 
            'c1_samples', 'c2_camples', 'neither_samples', 'clusters'
        ]
    print('\t'.join(colnames))

    for gene_family, vals in samples.items():
        if vals['percent'] == 1 and args.accessory_only:
            continue
        else:
            c1_samples = []
            c2_samples = []
            neither_samples = []
            for sample in vals['samples']:
                if sample in c1:
                    c1_samples.append(sample)
                elif sample in c2:
                    c2_samples.append(sample)
                else:
                    neither_samples.append(sample)
            
            # print details: gene_family, percent_all, percent_c1, percent_c2, percent_neither, c1_samples, c2_camples, neither_samples, clusters
            cols = [
                gene_family,
                f'{(vals["percent"] * 100):.2f}',
                f'{(len(c1_samples) / len(c1) * 100):.2f}', 
                f'{(len(c2_samples) / len(c2) * 100):.2f}', 
                f'{(len(neither_samples) / total_neither * 100):.2f}' if total_neither else"0",
            ]

            cluster = 'none'
            if len(c1_samples) and len(c2_samples):
                cluster = 'both'
            elif len(c1_samples): 
                cluster = 'c1 only'
            else:
                cluster = 'c2 only'
            cols.append(cluster)
            cols.append(clusters[gene_family]['products'],)

            if args.full_detail:
                cols.append(str(','.join(c1_samples)))
                cols.append(str(','.join(c2_samples)))
                cols.append(str(','.join(neither_samples)))
                cols.append(str(','.join(vals['cluster'])))

            if args.exclusive and (cluster == 'none' or cluster == 'both'):
                continue

            print('\t'.join(cols))
