#!/usr/bin/env python

"""
Use the data (taxonomy_mean_f and taxonomy_stdev_f) output by
run_compute_precision_recall.py to generate graphs.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import brewer2mpl
import sys
from os.path import join, isfile

def plot_tp_fp_distribution(
    taxonomy_summary_dir,
    datatypes,
    methods,
    tools,
    studies,
    output_dir):
    """
    Using the taxonomy mean and taxonomy stdev output of run_compute_precision_recall.py,
    plot the mean number of reads per assigned taxonomy vs. a given measure (TP, FP-known,
    FP-chimeric, FP-other).

    Parameters
    ----------
    taxonomy_summary_dir : string
        path to directory storing taxonomy mean / stdev files for each
        datatype, study and method benchmarked
    datatypes : list
        list of genes (16S, 18S) to examine
    methods : list
        OTU clustering methods (closed_ref, de_novo, open_ref)
    tools : dictionary
        for each method (key), the tools benchmarked as a list (value)
    studies : dictionary
        studies benchmarked for each datatype
    output_dir : string
        dirpath for storing the resulting bar graph

    Return
    ------
    None

    """
    # ex. 16S
    for datatype in datatypes:
        # ex. 1685
        for study in studies[datatype]:
            # ex. de novo
            for method in methods:
                taxonomy_mean_fp = join(taxonomy_summary_dir, datatype, method, "%s_taxonomy_mean.txt" % study)
                taxonomy_stdev_fp = join(taxonomy_summary_dir, datatype, method, "%s_taxonomy_stdev.txt" % study)

                taxonomy_mean = {}
                taxonomy_stdev = {}

                # load taxonomy_mean array
                if isfile(taxonomy_mean_fp):
                    with open(taxonomy_mean_fp, 'U') as taxonomy_mean_f:
                        for line in taxonomy_mean_f:
                            line = line.strip().split("\t")
                            if line[0] not in taxonomy_mean:
                                taxonomy_mean[line[0]] = []
                                for value in line[1:]:
                                    taxonomy_mean[line[0]].append(float(value))
                            else:
                                print "ERROR: %s is already in taxonomy_mean" % line[0]
                else:
                    print "skipping %s does not exist" % taxonomy_mean_fp
                    continue

                # load taxonomy_stdev array
                skip_stddev = False
                if isfile(taxonomy_stdev_fp):
                    with open(taxonomy_stdev_fp, 'U') as taxonomy_stdev_f:
                        for line in taxonomy_stdev_f:
                            line = line.strip().split("\t")
                            if line[0] not in taxonomy_stdev:
                                taxonomy_stdev[line[0]] = []
                                for value in line[1:]:
                                    taxonomy_stdev[line[0]].append(float(value))
                            else:
                                print "ERROR: %s is already in taxonomy_stdev" % line[0]
                else:
                    print "skipping %s does not exist" % taxonomy_stdev_fp
                    skip_stddev = True

                N = 4
                ind = np.arange(N)
                width = 0.1
                fig, ax = plt.subplots()
                colors = brewer2mpl.get_map('Set2', 'qualitative', len(tools[method])).mpl_colors

                i = 0
                plots = []
                tools_lgd = []
                for tool in tools[method]:
                    if tool not in taxonomy_mean:
                        print "WARNING: %s is not in taxonomy_mean tools" % tool
                        continue
                    #if not skip_stddev:
                    #   p = ax.bar(ind+(i*width), taxonomy_mean[tool], width, color=colors[i], yerr=taxonomy_stdev[tool])
                    #else:
                    #   p = ax.bar(ind+(i*width), taxonomy_mean[tool], width, color=colors[i])
                    p = ax.bar(ind+(i*width), taxonomy_mean[tool], width, color=colors[i])
                    plots.append(p)
                    tools_lgd.append(tool)
                    i += 1

                ax.set_ylabel("Mean number of reads per assigned taxonomy")
                ax.set_xlabel("Measure")
                ax.set_xticks(ind+4*width)
                ax.set_xticklabels(('TP', 'FP-known', 'FP-other', 'FP-chimeric'))
                legend=[]
                for i in range(0, len(plots)):
                    legend.append(plots[i][0])
                ax.legend(legend, tools_lgd)

                plt.savefig(join(output_dir, "tp_plot_%s_%s.png" % (method, study)), bbox_inches="tight", bbox_extra_artist=[ax.legend])
                plt.clf()

def main(argv):

    # output of run_compute_precision_recall.py
    taxonomy_summary_dir = sys.argv[1]

    # output directory
    output_dir = sys.argv[2]

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

    plot_tp_fp_distribution(
        taxonomy_summary_dir=taxonomy_summary_dir,
        datatypes=datatypes,
        methods=methods,
        tools=tools,
        studies=studies,
        output_dir=output_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
