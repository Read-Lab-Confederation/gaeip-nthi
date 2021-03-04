#! /usr/bin/env python3
"""

"""
import os
from collections import OrderedDict


def check_exists(files):
    has_error = False
    for f in files:
        if not os.path.exists(f):
            print(f"cannot access '{f}': No such file or directory\n", file=sys.stderr)
            has_error = True
    return has_error


def read_samples(c):
    samples = []
    with open(c, 'rt') as c_fh:
        for line in c_fh:
            samples.append(line.strip().split('\t')[0])
    return samples


def read_summary(summary_file):
    """
    So many columns!

    0: sample	                      32: original_total_bp	    64: final_qual_25th        	         96: mlst_blast_st
    1: is_paired	                  33: original_coverage	    65: final_qual_75th	                 97: mlst_blast_loci
    2: rank	                          34: original_read_total	66: CDS	                             98: mlst_blast_perfect_matches
    3: reason	                      35: original_read_min	    67: gene	                         99: mlst_ariba_st
    4: estimated_genome_size	      36: original_read_mean	68: rRNA	                        100: mlst_ariba_adk
    5: contig_non_acgtn	              37: original_read_std	    69: tRNA	                        101: mlst_ariba_atpg
    6: contig_percent_a	              38: original_read_median	70: tmRNA	                        102: mlst_ariba_frdb
    7: contig_percent_c	              39: original_read_max	    71: refseq_k21_id	                103: mlst_ariba_fuck
    8: contig_percent_g	              40: original_read_25th	72: refseq_k21_identity	            104: mlst_ariba_mdh
    9: contig_percent_n	              41: original_read_75th	73: refseq_k21_shared_hashes	    105: mlst_ariba_pgi
    10: contig_percent_t	          42: original_qual_min	    74: refseq_k21_median_multiplicity	106: mlst_ariba_reca
    11: contigs_greater_100k	      43: original_qual_mean	75: refseq_k21_p_value	            107: repeat_region
    12: contigs_greater_10k	          44: original_qual_std	    76: refseq_k21_comment	
    13: contigs_greater_1k	          45: original_qual_max	    77: refseq_k21_total	
    14: contigs_greater_1m	          46: original_qual_median	78: genbank_k21_match	
    15: l50_contig_count	          47: original_qual_25th	79: genbank_k21_overlap	
    16: max_contig_length	          48: original_qual_75th	80: genbank_k21_p_query	
    17: mean_contig_length	          49: final_total_bp	    81: genbank_k21_p_match	
    18: median_contig_length	      50: final_coverage	    82: genbank_k21_no_assignment	
    19: min_contig_length	          51: final_read_total	    83: genbank_k21_total	
    20: n50_contig_length	          52: final_read_min	    84: genbank_k31_match	
    21: num_contig_non_acgtn	      53: final_read_mean	    85: genbank_k31_overlap	
    22: percent_contigs_greater_100k  54: final_read_std	    86: genbank_k31_p_query	
    23: percent_contigs_greater_10k	  55: final_read_median	    87: genbank_k31_p_match	
    24: percent_contigs_greater_1k	  56: final_read_max	    88: genbank_k31_no_assignment	
    25: percent_contigs_greater_1m	  57: final_read_25th	    89: genbank_k31_total	
    26: total_contig	              58: final_read_75th	    90: genbank_k51_match	
    27: total_contig_length	          59: final_qual_min	    91: genbank_k51_overlap	
    28: checkm_lineage	              60: final_qual_mean	    92: genbank_k51_p_query	
    29: checkm_completeness	          61: final_qual_std	    93: genbank_k51_p_match	
    30: checkm_contamination	      62: final_qual_max	    94: genbank_k51_no_assignment	
    31: checkm_heterogeneity	      63: final_qual_median	    95: genbank_k51_total	

    """
    ga_samples = []
    excluded = {}
    summary = {}
    st_included = {}
    sequence_types = {}
    summary_cols = None
    with open(summary_file, 'rt') as summary_fh:
        for line in summary_fh:
            vals = line.rstrip().split('\t')
            if not summary_cols:
                summary_cols = vals
            else:
                if vals[1] == "True" and vals[2] != "exclude":
                    # Only include paired-end samples that passed
                    summary_info = dict(zip(summary_cols, vals))
                    st = vals[99].replace("*", "")
                    if summary_info['genbank_k31_match'] != "Haemophilus influenzae":
                        excluded[vals[0]] = f'sourmash (k31) matches {summary_info["genbank_k31_match"]}'
                    elif st == "ND" or st == "Novel":
                        excluded[vals[0]] = f'sequence type ({st}) was not an integer (e.g. known ST)'
                    else:
                        summary[vals[0]] = summary_info
                        summary[vals[0]]['st'] = st
                        if st not in sequence_types:
                            sequence_types[st] = []
                        sequence_types[st].append(vals[0])
                        st_included[st] = False
                        if vals[0].startswith("GA"):
                            ga_samples.append(vals[0])
                else:
                    excluded[vals[0]] = vals[3]
    return [summary, st_included, sequence_types, ga_samples, excluded]


def read_fastani(fastani_file):
    """
    0: query 
    1: reference
    2: ani
    3: mapped_fragments
    4: total_fragments
    """
    fastani = {}
    fastani_cols = None
    with open(fastani_file, 'rt') as fastani_fh:
        for line in fastani_fh:
            vals = line.rstrip().split('\t')
            if not fastani_cols:
                fastani_cols = vals
            else:
                if not vals[0].endswith('reference'):
                    fastani[vals[0]] = dict(zip(fastani_cols, vals))
    return fastani


if __name__ == '__main__':
    import argparse as ap
    import math
    import os
    import sys

    parser = ap.ArgumentParser(
        prog='parse-pirate',
        conflict_handler='resolve',
        description=(
            f'parse-pirate.py - Parse PIRATE output and give details about accessory genome.'
        )
    )

    parser.add_argument('summary', metavar="STR", type=str, help='Summary file from bactopia summary')
    parser.add_argument('exclude', metavar="STR", type=str, help='Exclude file from bactopia summary')
    parser.add_argument('fastani', metavar="STR", type=str, help='Output from bactopia fastani')
    parser.add_argument('c1', metavar="STR", type=str, help='Text file with names of C1 samples')
    parser.add_argument('c2', metavar="STR", type=str, help='Text file with names of C2 samples')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    # Make sure files exist
    if check_exists([args.summary, args.exclude, args.fastani, args.c1, args.c2]):
        print("Please make sure the correct path is given.", file=sys.stderr)
        sys.exit(1)

    summary, st_included, sequence_types, ga_samples, excluded = read_summary(args.summary)
    fastani = read_fastani(args.fastani)
    c1 = read_samples(args.c1)
    c2 = read_samples(args.c2)
    cutoff = math.floor(min([float(fastani[c]['ani']) for c in c2]))

    # Append to exluded
    with open(args.exclude, 'rt') as exclude_fh:
        for line in exclude_fh:
            vals = line.rstrip().split('\t')
            if vals[0] not in excluded:
                excluded[vals[0]] = vals[2]
 
    # Include all Georgia samples that are Haemophilus influenzae
    samples = OrderedDict()
    for sample in sorted(ga_samples):
        st_included[summary[sample]['st']] = True
        samples[sample] = "Georgia sample"

    # Include samples with ANI > cutoff
    for sample, vals in fastani.items():
        if sample in summary:
            if float(vals['ani']) >= cutoff:
                if sample not in samples:
                    samples[sample] = f"ANI ({vals['ani']})to {vals['query']} greater than cutoff ({cutoff})"
                    st_included[summary[sample]['st']] = True

    # Get a representative of the remaining STs
    for st in sorted(st_included.keys(), key=int):
        is_included = st_included[st]
        if not is_included:
            representative = None
            for sample in sequence_types[st]:
                if sample not in samples:
                    # Prioritorize gold samples first
                    if representative:
                        if summary[sample]['rank'] == "gold":
                            # check if other is gold as well
                            if summary[representative]['rank'] == "gold":
                                # Both gold, take the one with highest checkm completeness                       
                                if float(summary[representative]['checkm_completeness']) == float(summary[sample]['checkm_completeness']):
                                    # Take the one with fewest contigs
                                     if int(summary[representative]['total_contig']) < int(summary[sample]['total_contig']):
                                         # New representative
                                        representative = sample
                                elif float(summary[representative]['checkm_completeness']) < float(summary[sample]['checkm_completeness']):
                                    # New representative
                                    representative = sample
                            else:
                                # New representative
                                representative = sample
                        elif summary[sample]['rank'] == "silver":
                            # check if other is not gold
                            if summary[representative]['rank'] != "gold":
                                # check if other is silver as well
                                if summary[representative]['rank'] == "silver":
                                    # Both silver, take the one with highest checkm completeness                       
                                    if float(summary[representative]['checkm_completeness']) == float(summary[sample]['checkm_completeness']):
                                        # Take the one with fewest contigs
                                        if int(summary[representative]['total_contig']) < int(summary[sample]['total_contig']):
                                            # New representative
                                            representative = sample
                                    elif float(summary[representative]['checkm_completeness']) < float(summary[sample]['checkm_completeness']):
                                        # New representative
                                        representative = sample
                                else:
                                    # New representative
                                    representative = sample
                        elif summary[sample]['rank'] == "silver":
                            # check if other is not gold or silver
                            if summary[representative]['rank'] != "gold" and summary[representative]['rank'] != "silver":
                                # check if other is bronze as well
                                if summary[representative]['rank'] == "silver":
                                    # Both silver, take the one with highest checkm completeness                       
                                    if float(summary[representative]['checkm_completeness']) == float(summary[sample]['checkm_completeness']):
                                        # Take the one with fewest contigs
                                        if int(summary[representative]['total_contig']) < int(summary[sample]['total_contig']):
                                            # New representative
                                            representative = sample
                                    elif float(summary[representative]['checkm_completeness']) < float(summary[sample]['checkm_completeness']):
                                        # New representative
                                        representative = sample
                                else:
                                    # New representative
                                    representative = sample
                    else:
                        representative = sample
            st_included[st] = True
            samples[sample] = f"Representative for ST {st}"

    # Write updated excludes files
    with open('excludes.txt', 'wt') as excludes_fh:
        excludes_fh.write(f'sample\treason\n')
        for sample, reason in sorted(excluded.items()):
            if sample in fastani:
                reason = f'{reason};ANI against C1 ({fastani[sample]["reference"]}) is {fastani[sample]["ani"]}'
            excludes_fh.write(f'{sample}\t{reason}\n')
    
    # Write updated representative set
    with open('representatives.txt', 'wt') as representatives_fh:
        cols = [
            'sample',
            'st',
            'is_georgia',
            'is_c1',
            'is_c2',
            'is_ani_cutoff'
            'inclusion_reason'
        ]
        col_string = '\t'.join(cols)
        representatives_fh.write(f'{col_string}\n')
        for sample, reason in samples.items():
            row = [sample, summary[sample]['st']]
            row.append(True if sample in ga_samples else False)
            row.append(True if sample in c1 else False)
            row.append(True if sample in c2 else False)
            row.append(True if float(fastani[sample]['ani']) > cutoff else False)
            row.append(reason)
            row_string = '\t'.join([str(r) for r in row])
            representatives_fh.write(f'{row_string}\n')
