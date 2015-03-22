#!/usr/bin/env python
"""
This script creates a FASTA file of representative sequences
with staggered abundances of each OTU representative sequence
in the even community. In addition, the script outputs the
abundance file to reproduce the dataset.

usage: generate_abundance_file.py otu_map.txt min_num_reads_per_seq \
                                  max_num_reads_per_seq num_reads_to_generate \
                                  original_fasta_repset_fp out_1_abundance_file.txt \
                                  out_2_repset_staggered_art.fasta
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


import random
import sys
from collections import defaultdict
from skbio.parse.sequences import parse_fasta

def output_fasta_abundance_files(abund,
                                 num_reads_to_generate,
                                 rep_set_fp,
                                 output_fasta_fp,
                                 output_abundance_fp):
    """
    Extract all FASTA sequences matching to reference IDs

    Parameters
    ----------
    abund : dictionary
        keys are reference IDs and values are the fractional
        abundance of number of reads to generate from the
        sequence with reference ID
    num_reads_to_generate : integer
        number of total reads to generate
    rep_set_fp : string
        filepath to FASTA file storing all sequences with
        reference IDs in abund keys
    output_fasta_fp : string
        filepath to resulting FASTA file of staggered sequences,
        representing all reference IDs in abund with the
        fractional number of occurrences specified by abund
    output_abundance_fp : string
        filepath to resulting abundance file
    """

    out_1 = open(output_abundance_fp, 'w')
    out_2 = open(output_fasta_fp, 'w')
    with open(rep_set_fp, 'U') as rep_set_fp:
        for label, seq in parse_fasta(rep_set_fp):
            # sequence to be amplified
            if label in abund:
                generate = int(abund[label] * num_reads_to_generate)
                out_1.write("%s\t%s\n" % (label, abund[label]))
                if generate == 0:
                    generate = 1
                for i in range(0,generate):
                    out_2.write(">%s\n%s\n" % (label, seq))
    out_1.close()
    out_2.close()

def create_abundance(otu_map_fp, minimum, maximum):
    """
    Parameters
    ----------
    otu_map_fp : string
        filepath to OTU map
    minimum : integer
        minimum number of reads to generate per reference sequence
    maximum : integer
        maximum number of reads to generate per reference sequence

    Return
    ------
    abundance_file : dictionary
        a dict containing the fractional abundance for each seq
    """
    # OTU map IDs
    ids = []
    with open(otu_map_fp, 'U') as otu_map_f:
        for line in otu_map_f:
            _id = line.strip().split("\t")[0]
            ids.append(_id)

    ids_dict = defaultdict(list)
    total = 0.0
    for i in ids:
        r = random.randint(minimum, maximum)
        ids_dict[i] = r
        total += r
    return {key: ids_dict[key] / total for key in ids_dict}

def main(argv):

    # return dict of OTU ids and fractional abundance for each seq 
    abund = create_abundance(otu_map_fp=sys.argv[1],
                             minimum=int(sys.argv[2]),
                             maximum=int(sys.argv[3]))

    output_fasta_abundance_files(
        abund=abund,
        num_reads_to_generate=int(sys.argv[4]),
        rep_set_fp=sys.argv[5],
        output_fasta_fp=sys.argv[6],
        output_abundance_fp=sys.argv[7])


if __name__ == "__main__":
    main(sys.argv[1:])

