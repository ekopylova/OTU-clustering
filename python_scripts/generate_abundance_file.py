#!/opt/python-2.7.3/bin/python 

# usage   : generate_abundance_file.py otu_map.txt min_num_reads_per_seq max_num_reads_per_seq num_reads_to_generate original_fasta_repset_fp out_1_abundance_file.txt out_2_repset_staggered_art.fasta
# purpose : create a FASTA file of staggered abundances of each OTU representative sequence in the even community
#           (also output the abundance file to reproduce the dataset)
# date    : 24 Nov 2014
# author  : Evguenia Kopylova (jenya.kopylov@gmail.com)

import random
import sys
from collections import defaultdict
from skbio.parse.sequences import parse_fasta

def create_abundance (ids, minimum, maximum):
    '''
    ids: a list of seq IDs
    min: the min number of reads from any seq
    max: the max number of reads from any seq

    return a dict contain the fractional abundance for each seq
    '''
    ids_dict = defaultdict(list)
    total = 0.0
    for i in ids:
        r = random.randint(minimum, maximum)
        ids_dict[i] = r
        total += r
    return {key: ids_dict[key] / total for key in ids_dict}

def main(argv):
    # OTU map IDs
    ids = []
    with open(sys.argv[1], 'U') as otu_map_fp:
        for line in otu_map_fp:
            _id = line.strip().split("\t")[0]
            ids.append(_id)

    # return dict of OTU ids and fractional abundance for each seq 
    abund = create_abundance(ids, int(sys.argv[2]), int(sys.argv[3]))

    num_reads_to_generate = int(sys.argv[4])

    # extract all FASTA sequences matching to reference IDs
    out_1 = open(sys.argv[6], 'w')
    out_2 = open(sys.argv[7], 'w')
    with open(sys.argv[5], 'U') as rep_set_fp:
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

if __name__ == "__main__":
    main(sys.argv[1:])

