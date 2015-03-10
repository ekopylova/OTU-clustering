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
    rootdir = sys.argv[1]

    # output results directory
    outdir_root = sys.argv[2]

    # studies 16S
    studies_bac = sys.argv[3].split()

    # studies 18S
    studies_euk = sys.argv[4].split()

    # list of studies for each gene type
    studies = {'16S': [], '18S': []}

    for study in studies_bac:
      if study not in studies['16S']:
        studies['16S'].append(study)

    for study in studies_euk:
      if study not in studies['18S']:
        studies['18S'].append(study)

    # tools
    tools_denovo = sys.argv[5].split()
    tools_closed_ref = sys.argv[6].split()
    tools_open_ref = sys.argv[7].split()

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

    # genes 
    datatypes = ['16S', '18S']

    # OTU picking methods
    methods = ['de_novo', 'closed_ref', 'open_ref']

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in methods:
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
                            print "skipping %s does not exist" % os.path.join(search_dir, otu_table)

                    else:
                        print "%s already exists" % outdir


        
            
                                    
                                
