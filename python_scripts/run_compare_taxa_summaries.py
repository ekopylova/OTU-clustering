#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Compute Pearson correlation for summarized taxonomies for all tools
vs. all tools.

Dependencies: QIIME 1.9.0
usage:        python run_compare_taxa_summaries.py
"""

import sys
import os
from subprocess import Popen, PIPE
from glob import glob

if __name__ == '__main__':

    # output results directory
    outdir_root = "/scratch/Users/evko1434/working_dir/compare_taxa_summaries_mc2"

    # summarized taxonomies for all tools (output of run_summarize_taxa.py)
    rootdir = "/scratch/Users/evko1434/working_dir/summarize_taxa_mc2"

    # mapping files directory
    mapping_dir = "/scratch/Users/evko1434/supplemental_otu_clustering_datasets/mapping_files"

    # genes
    datatypes = ['16S', '18S']

    # OTU picking methods
    picking = ['de_novo', 'closed_ref', 'open_ref']

    # studies per genes
    studies = {'16S': ['1685', '1686', '1688', '449', '632'],
               '18S': ['nematodes', '2107']}

    # tools per OTU picking method
    tools = {'de_novo': ['swarm', 'uclust', 'sumaclust', 'usearch', 'usearch61', 'uparse_q3', 'uparse_q16'],
             'closed_ref': ['sortmerna', 'uclust', 'usearch', 'usearch61'],
             'open_ref': ['sortmerna_sumaclust', 'uclust', 'usearch61']}

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in picking:
            # ex. 1685
            for study in studies[datatype]:
                i = 1
                for tool1 in tools[method]:
                    search_dir_1 = os.path.join(rootdir, datatype, method, "%s_%s" % (tool1, study))
                    if not os.path.exists(search_dir_1):
                        print "%s does not exist" % search_dir_1
                        continue
                    otu_table = "otu_table_mc2_L6.txt"
                    for tool2 in tools[method][i:]:
                        search_dir_2 = os.path.join(rootdir, datatype, method, "%s_%s" % (tool2, study))
                        if not os.path.exists(search_dir_2):
                            print "%s does not exist" % search_dir_2
                            continue
                        outdir = os.path.join(outdir_root, datatype, method, "%s_%s_%s" % (tool1, tool2, study))
                        # run QIIME's compare_taxa_summaries.py script
                        if not os.path.exists(outdir):
                            taxonomies_to_compare = "%s,%s" % (
                                os.path.join(search_dir_1, otu_table), os.path.join(search_dir_2, otu_table))
                            compare_taxa_summaries_command = ["compare_taxa_summaries.py",
                                                               "-i",
                                                               taxonomies_to_compare,
                                                               "-m",
                                                               "paired",
                                                               "-o",
                                                               outdir]
                            print "command = ", compare_taxa_summaries_command
                            proc = Popen(compare_taxa_summaries_command,
                                         stdout=PIPE,
                                         stderr=PIPE,
                                         close_fds=True)
                            proc.wait()
                            stdout, stderr = proc.communicate()
                            if stderr:
                                print stderr
                    i = i+1
            
                print "STUDY = ", study
                # output the results of compare_taxa_summaries.py
                i = 1
                for tool1 in tools[method]:
                    for tool2 in tools[method][i:]:
                        search_dir = os.path.join(
                            outdir_root, datatype, method, "%s_%s_%s" % (tool1, tool2, study))
                        if os.path.isfile(os.path.join(search_dir, "overall_comparison.txt")):
                            with open(os.path.join(search_dir, "overall_comparison.txt")) as in_file:
                                for line in in_file:
                                    if (line.startswith("#") or line.startswith("Correlation")):
                                        continue
                                    sys.stdout.write("%s\t%s\t%s" % (tool1, tool2, line))
                    i = i+1               

        
            
                                    
                                
