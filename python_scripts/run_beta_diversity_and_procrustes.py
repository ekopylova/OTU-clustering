#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Run QIIME's beta-diversity and Procrustes analysis for all tools vs. UCLUST

Dependencies: QIIME 1.9.0, Emperor
usage:        python run_beta_diversity_and_procrustes.py
"""

import sys
import os
from subprocess import Popen, PIPE
from glob import glob

if __name__ == '__main__':

    # original BIOM tables excl. singletons directory
    # (same outdir_root path as in run_filter_singleton_otus.py)
    rootdir = sys.argv[1]

    # output results for beta diversity analysis
    outdir_beta = sys.argv[2]

    # output results for Procrustes analysis
    outdir_procrustes = sys.argv[3]

    # trees (generated during de novo/open ref OTU picking, will be in
    # same directory as all results)
    treedir = sys.argv[4]

    # trees for closed-reference analysis (16S distributed with Greengenes, 18S built)
    tree_fp = {'16S': "", '18S': ""}
    tree_fp['16S'] = sys.argv[5]
    tree_fp['18S'] = sys.argv[6]

    # mapping files directory
    mapping_dir = sys.argv[7]

    # studies 16S
    studies_bac = sys.argv[8].split()

    # studies 18S
    studies_euk = sys.argv[9].split()

    # list of studies for each gene type
    studies = {'16S': [], '18S': []}

    for study in studies_bac:
      if study not in studies['16S']:
        studies['16S'].append(study)

    for study in studies_euk:
      if study not in studies['18S']:
        studies['18S'].append(study)

    # tools
    tools_denovo = sys.argv[10].split()
    tools_closed_ref = sys.argv[11].split()
    tools_open_ref = sys.argv[12].split()

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

    # coordinate matrices
    coordinate_matrices = sys.argv[13].split()

    # genes
    datatypes = ['16S', '18S']

    # OTU picking methods
    methods = ['closed_ref', 'de_novo', 'open_ref']

    # depth of coverage for even sampling (based on results of run_summarize_tables.py)
    sampling_depth_cr = {'449': '380',
                         '632': '1000',
                         '2105': '1000',
                         '2107': '1000',
                         '1919': '100',
                         '2020': '1000'}
    sampling_depth_dn = {'449': '100',
                         '632': '1000',
                         '2105': '1000',
                         '2107': '1000',
                         '1919': '100',
                         '2020': '1000'}
    sampling_depth_or = {'449': '413',
                         '632': '1000',
                         '2105': '1000',
                         '2107': '1000',
                         '1919': '100',
                         '2020': '1000'}

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in methods:
            # ex. 449
            for study in studies[datatype]:
                # ex. swarm
                for tool in tools[method]:
                    search_dir = os.path.join(rootdir, datatype, method, "%s_%s" % (tool, study))
                    outdir = os.path.join(outdir_beta, datatype, method, "%s_%s" % (tool, study))
                    if not os.path.exists(outdir):
                        sample_depth = ''
                        otu_table = "otu_table_mc2.biom"

                        if (method == "closed_ref" and sampling_depth_cr[study] is not ''):
                            sample_depth = sampling_depth_cr[study]
                            tree = tree_dir[datatype]
                                
                        elif (method == "open_ref" and sampling_depth_or[study] is not ''):
                            sample_depth = sampling_depth_or[study]
                            tree = os.path.join(
                                outdir_root, "program_results", datatype, method, "%s_%s" % (tool, study), "rep_set.tre")

                        elif (method == "de_novo" and sampling_depth_dn[study] is not ''):
                            sample_depth = sampling_depth_dn[study]
                            tree = os.path.join(
                                outdir_root, "program_results", datatype, method, "%s_%s" % (tool, study), "rep_set.tre")

                        # BIOM table exists
                        if os.path.isfile(os.path.join(search_dir, otu_table)):
                                beta_div_command = ["beta_diversity_through_plots.py",
                                                    "-i",
                                                    os.path.join(search_dir, otu_table),
                                                    "-o",
                                                    outdir,
                                                    "-t",
                                                    tree,
                                                    "-m",
                                                    os.path.join(mapping_dir, "study_%s_mapping_file.txt" % study),
                                                    "-e",
                                                    sample_depth]
                            print "command = ", beta_div_command
                            proc = Popen(beta_div_command,
                                         stdout=PIPE,
                                         stderr=PIPE,
                                         close_fds=True)
                            proc.wait()
                            stdout, stderr = proc.communicate()
                            if stderr:
                                print stderr
                        else:
                            print "skpping %s not found" % os.path.join(search_dir, otu_table)
            # run procrustes analysis
            # loop through all the study output folders
            for study in studies[datatype]:
                # loop through all the tools
                for tool in tools[method]:
                    search_dir = os.path.join(outdir_beta, datatype, method, "%s_%s" % (tool, study))
                    # beta diversity was not computed for this tool_study
                    if not os.path.exists(search_dir):
                        continue
                    outdir = os.path.join(outdir_procrustes, datatype, method, "uclust_vs_%s_%s" % (tool, study))
                    if not os.path.exists(outdir):
                            for matrix in coordinate_matrices:
                                matrices = "%s,%s" % (
                                    os.path.join(
                                        outdir_beta, datatype, method, "uclust_%s" % study, "%s_unifrac_pc.txt" % matrix), os.path.join(search_dir, "%s_unifrac_pc.txt" % matrix))

                                # transform coordiate matrices
                                transform_coordinate_matrices_command = ["transform_coordinate_matrices.py",
                                                                         "-i",
                                                                         matrices,
                                                                         "-o",
                                                                         os.path.join(outdir, matrix)]
                                proc = Popen(transform_coordinate_matrices_command,
                                             stdout=PIPE,
                                             stderr=PIPE,
                                             close_fds=True)
                                proc.wait()
                                stdout, stderr = proc.communicate()
                                if stderr:
                                    print stderr
                    
                                # plot emperor plots
                                make_emperor_command = ["make_emperor.py",
                                                        "-c",
                                                        "-i",
                                                        os.path.join(outdir, matrix),
                                                        "-o",
                                                        os.path.join(outdir, matrix, "plots"),
                                                        "-m",
                                                        os.path.join(mapping_dir, "study_%s_mapping_file.txt" % study)]
                                proc = Popen(make_emperor_command,
                                             stdout=PIPE,
                                             stderr=PIPE,
                                             close_fds=True)
                                proc.wait()
                                stdout, stderr = proc.communicate()
                                if stderr:
                                    print stderr
