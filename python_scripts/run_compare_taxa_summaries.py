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
    outdir_root = sys.argv[1]

    # summarized taxonomies for all tools (output of run_summarize_taxa.py)
    rootdir = sys.argv[2]

    # mapping files directory
    mapping_dir = sys.argv[3]

    # studies 16S
    studies_bac = sys.argv[4].split()

    # studies 18S
    studies_euk = sys.argv[5].split()

    # genes
    datatypes = ['16S', '18S']

    # OTU picking methods
    methods = ['de_novo', 'closed_ref', 'open_ref']

    # list of studies for each gene type
    studies = {'16S': [], '18S': []}

    for study in studies_bac:
      if study not in studies['16S']:
        studies['16S'].append(study)

    for study in studies_euk:
      if study not in studies['18S']:
        studies['18S'].append(study)

    # tools
    tools_denovo = sys.argv[6].split()
    tools_closed_ref = sys.argv[7].split()
    tools_open_ref = sys.argv[8].split()

    # list of tools for each OTU picking method
    tools = {'de_novo': [], 'closed_ref': [], 'open_ref': []}

    for tool in tools_denovo:
      if tool not in tools['de_novo']:
        tools['de_novo'].append(tool)
    for tool in tools_closed_ref:
      if tool not in tools['closed_ref']:
        tools['closed_ref'].append(tool)
    for tool in tools_open_ref:
      if tool not in tools['open_ref']:
        tools['open_ref'].append(tool)

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in methods:
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
                            proc = Popen(compare_taxa_summaries_command,
                                         stdout=PIPE,
                                         stderr=PIPE,
                                         close_fds=True)
                            proc.wait()
                            stdout, stderr = proc.communicate()
                            if stderr:
                                print stderr
                        else:
                            print "skipping %s already exists" % (outdir)
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

        
            
                                    
                                
