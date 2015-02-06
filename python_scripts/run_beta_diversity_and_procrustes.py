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

    # output results for beta diversity analysis
    outdir_beta = "/scratch/Users/evko1434/working_dir/beta_diversity_through_plots"

    # output results for Procrustes analysis
    outdir_procrustes = "/scratch/Users/evko1434/working_dir/procrustes_analyses"

    # original BIOM tables excl. singletons directory
    # (same outdir_root path as in run_filter_singleton_otus.py)
    rootdir = "/scratch/Users/evko1434/working_dir/otu_tables_mc2"

    # trees for closed-reference analysis (16S distributed with Greengenes, 18S built)
    tree_fp = {'16S': "/scratch/Users/evko1434/reference/gg_13_8_otus/trees/97_otus.tree",
               '18S': "/scratch/Users/evko1434/reference/Silva_111_post/trees/97_Silva_111_rep_set_pfiltered.tre"}

    # mapping files directory
    mapping_dir = "/scratch/Users/evko1434/supplemental_otu_clustering_datasets/mapping_files"

    # genes
    datatypes = ['16S', '18S']

    # OTU picking methods
    picking = ['closed_ref', 'de_novo', 'open_ref']

    # list of studies per gene type
    studies = {'16S': ['632', '449'],
               '18S': ['2107']}

    # list of tools per OTU picking method
    tools = {'de_novo': ['sumaclust', 'swarm', 'uclust', 'usearch61', 'usearch', 'uparse_q3', 'uparse_q16'],
             'closed_ref': ['sortmerna', 'uclust', 'usearch', 'usearch61'],
             'open_ref': ['sortmerna_sumaclust', 'uclust', 'usearch61']}

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

    # coordinate matrices for UniFrac
    coordinate_matrices = ['weighted', 'unweighted']

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in picking:
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
                            print "OTU table %s not found" % os.path.join(search_dir, otu_table)
            # run procrustes analysis
            # loop through all the study output folders
            print "Run procrustes analysis .."
            for study in studies[datatype]:
                # loop through all the tools
                for tool in tools[method]:
                    search_dir = os.path.join(outdir_beta, datatype, method, "%s_%s" % (tool, study))
                    # beta diversity was not computed for this tool_study
                    if not os.path.exists(search_dir):
                        continue
                    outdir = os.path.join(outdir_procrustes, datatype, method, "uclust_vs_%s_%s" % (tool, study))
                    print "outdir = ", outdir
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
                                print "command = ",transform_coordinate_matrices_command
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
                                print "command = ",make_emperor_command
                                proc = Popen(make_emperor_command,
                                             stdout=PIPE,
                                             stderr=PIPE,
                                             close_fds=True)
                                proc.wait()
                                stdout, stderr = proc.communicate()
                                if stderr:
                                    print stderr
