#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Generate table summaries for all results.
Table summaries are used to choose sampling depths
during rarefaction (for run_single_rarefaction_and_plot.py)
and to obtain the OTU count.

Dependencies: QIIME 1.9.0, BIOM format 2.1
usage:        python run_summarize_tables.py
"""


import sys
import os
from subprocess import Popen, PIPE
from glob import glob

if __name__ == '__main__':

    # original BIOM tables excl. singletons directory
    # (same outdir_root path as in run_filter_singleton_otus.py)
    rootdir = "/scratch/Users/evko1434/working_dir/otu_tables_mc2"

    # output results directory 
    outdir_root = "/scratch/Users/evko1434/working_dir/summarize_tables_mc2"

    # genes
    datatypes = ['16S', '18S']

    # OTU picking methods
    picking = ['closed_ref', 'de_novo', 'open_ref']

    # list of studies per gene
    studies = {'16S': ['even', 'staggered', '1685', '1686', '1688', '449'],
               '18S': ['nematodes', '2107']}

    # list of tools
    tools = {'de_novo': ['uclust', 'usearch', 'usearch61', 'swarm', 'sumaclust', 'uparse_q3', 'uparse_q16'],
             'closed_ref': ['sortmerna', 'uclust', 'usearch', 'usearch61'],
             'open_ref': ['sortmerna_sumaclust', 'uclust', 'usearch61']}

    # ex. 16S
    for datatype in datatypes:
        # ex. de_novo
        for method in picking:
            # ex. swarm
            for tool in tools[method]:
                # ex. 1685
                for study in studies[datatype]:
                    search_dir = os.path.join(rootdir, datatype, method, "%s_%s" % (tool, study))
                    #print "search_dir = ",search_dir
                    outdir = os.path.join(outdir_root, datatype, method, "%s_%s" % (tool, study))
                    #print "\toutdir = ", outdir
                    if not os.path.exists(outdir):
                        os.makedirs(outdir)
                        otu_table = "otu_table_mc2.biom"

                        # run BIOM's summarize table command
                        summarize_tbl_command = ["biom",
                                                 "summarize-table",
                                                 "-i",
                                                 os.path.join(search_dir, otu_table),
                                                 "-o",
                                                 os.path.join(outdir, "table_summary.txt")]
                        #print "command = ", summarize_tbl_command
                        proc = Popen(summarize_tbl_command,
                                     stdout=PIPE,
                                     stderr=PIPE,
                                     close_fds=True)
                        proc.wait()
                        stdout, stderr = proc.communicate()
                        if stderr:
                            print stderr
               


        
            
                                    
                                
