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

    print "tools = ", tools

    exit()

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in methods:
            # ex. swarm
            for tool in tools[method]:
                # ex. 1685
                for study in studies[datatype]:
                  if method == "open_ref":
                      biom_table = "otu_table_mc2_w_tax_no_pynast_failures.biom"
                  else:
                      biom_table = "otu_table.biom"
                  search_dir = os.path.join(rootdir, datatype, method, "%s_%s" % (tool, study))
                  if os.path.isfile(os.path.join(search_dir, biom_table)):
                      outdir = os.path.join(outdir_root, datatype, method, "%s_%s" % (tool, study))
                      if not os.path.exists(outdir):
                          os.makedirs(outdir)

                          # run QIIME's filter_otus_from_otu_table.py script
                          filter_otus_from_otu_table_command = ["filter_otus_from_otu_table.py",
                                                                "-i",
                                                                os.path.join(search_dir, biom_table),
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
                      print "skipping %s: does not exist" % os.path.join(search_dir, biom_table)


        
            
                                    
                                
