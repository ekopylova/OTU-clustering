#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Generate taxa summaries for all BIOM tables (excl. singletons).

Dependencies: QIIME 1.9.0
usage:        python run_summarize_taxa.py
note:         This script will need to be modified if using original
              BIOM tables incl. singleton OTUs. The major change is the
              assignation "otu_table = 'otu_table_mc2.biom'" which have
              slightly different names for closed_ref/de_novo vs. open_ref
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
    outdir_root = "/scratch/Users/evko1434/working_dir/summarize_taxa_mc2"

    # genes
    datatypes = ['16S', '18S']
    
    # OTU picking methods
    picking = ['de_novo', 'closed_ref', 'open_ref']

    # studies per gene
    studies = {'16S': ['even', 'staggered', '1685', '1686', '1688', '449'],
               '18S': ['nematodes', '2107']}

    # tools per OTU picking method
    tools = {'de_novo': ['swarm', 'sumaclust', 'usearch', 'usearch61', 'uclust', 'uparse_q16', 'uparse_q3'],
             'closed_ref': ['sortmerna', 'uclust', 'usearch', 'usearch61'],
             'open_ref': ['sortmerna_sumaclust', 'uclust', 'usearch61']}

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in picking:
            # ex. swarm
            for tool in tools[method]:
                # ex. 1685
                for study in studies[datatype]:
                    search_dir = os.path.join(rootdir, datatype, method, "%s_%s" % (tool, study))
                    outdir = os.path.join(outdir_root, datatype, method, "%s_%s" % (tool, study))
                    if not os.path.exists(outdir):
                        otu_table = "otu_table_mc2.biom"

                        if os.path.isfile(os.path.join(search_dir, otu_table)):
                            # run QIIME's summarize_taxa.py script
                            summarize_taxa_command = ["summarize_taxa.py",
                                                      "-i",
                                                      os.path.join(search_dir, otu_table),
                                                      "-o",
                                                      outdir]
                            #print "command = ", summarize_taxa_command
                            proc = Popen(summarize_taxa_command,
                                         stdout=PIPE,
                                         stderr=PIPE,
                                         close_fds=True)
                            proc.wait()
                            stdout, stderr = proc.communicate()
                            if stderr:
                                print stderr
                        else:
                            print "%s does not exist" % os.path.join(search_dir, otu_table)

                    else:
                        print "%s already exists" % outdir


        
            
                                    
                                
