#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Filter singleton OTUs from original BIOM table for multiple tools.

Dependencies: QIIME 1.9.0
usage:        python run_filter_singleton_otus.py

"""

import sys
import os
from subprocess import Popen, PIPE
from glob import glob

if __name__ == '__main__':
    
    # program results directory
    # must follow the structure (default output of commands_16S.sh and commands_18S.sh):
    # $rootdir/16S/
    #             /closed_ref/
    #                        ..
    #             /open_ref/
    #                       ..
    #             /de_novo/
    #                     ..
    #                     /swarm_1685/
    #                                /(output from OTU-picking)
    #                     /sumaclust_1685/
    #                     ..
    rootdir = "/scratch/Users/evko1434/working_dir/program_results"

    # output results directory
    outdir_root = "/scratch/Users/evko1434/working_dir/otu_tables_mc2"

    # genes 
    datatypes = ['16S', '18S']

    # OTU picking methods
    picking = ['de_novo', 'closed_ref', 'open_ref']

    # list of studies for each gene type
    studies = {'16S': ['1685', '1686', '1688', '449', '632', 'even', 'staggered'],
               '18S': ['2107', 'nematodes']}

    # list of tools for each OTU picking method
    tools = {'de_novo': ['uclust', 'usearch', 'usearch61', 'swarm', 'sumaclust', 'uparse_q3', 'uparse_q16'],
             'closed_ref': ['sortmerna', 'uclust', 'usearch61', 'usearch'],
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
                    if os.path.isfile(os.path.join(search_dir, "otu_table.biom")):
                        outdir = os.path.join(outdir_root, datatype, method, "%s_%s" % (tool, study))
                        if not os.path.exists(outdir):
                            os.makedirs(outdir)

                            # run QIIME's filter_otus_from_otu_table.py script
                            filter_otus_from_otu_table_command = ["filter_otus_from_otu_table.py",
                                                                  "-i",
                                                                  os.path.join(search_dir, "otu_table.biom"),
                                                                  "-o",
                                                                  os.path.join(outdir, "otu_table_mc2.biom"),
                                                                  "--min_count",
                                                                  "2"]
                            #print "command = ", filter_otus_from_otu_table_command
                            proc = Popen(filter_otus_from_otu_table_command,
                                         stdout=PIPE,
                                         stderr=PIPE,
                                         close_fds=True)
                            proc.wait()
                            stdout, stderr = proc.communicate()
                            if stderr:
                                print stderr
                    else:
                        print "%s does not exist" % os.path.join(search_dir, "otu_table.biom")


        
            
                                    
                                
