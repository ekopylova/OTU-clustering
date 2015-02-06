#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Generate layered taxonomy barcharts for each tool.

usage:        python run_generate_taxa_barcharts [16S, 18S, ITS] study \
              tool1_summarize_taxa_L6.txt tool2_summarize_taxa_L6.txt .. tooln_summarize_taxa_L6.txt
"""

import sys
import os
import random
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import brewer2mpl

if __name__ == '__main__':
    rootdir = "/scratch/Users/evko1434/working_dir/taxonomy_barcharts"
    num_tools = len(sys.argv)-3
    gene = sys.argv[1]
    study = sys.argv[2]

    top_N_taxa = 20
    abundances = {}
    tools = {"de_novo": ['sumaclust', 'swarm', 'uclust', 'usearch', 'usearch61', 'uparse_q3', 'uparse_q16'],
             "closed_ref": ['sortmerna', 'uclust', 'usearch', 'usearch61'],
             "open_ref": ['sortmerna_sumaclust', 'uclust', 'usearch61']}
    labels = []

    for i in range(3, num_tools+3):
        for method in tools:
            if method in sys.argv[i]:
                for tool in tools[method]:
                    if "%s_" % tool in sys.argv[i]:
                        if (method == "closed_ref" and (tool in tools["de_novo"] or tool in tools["open_ref"])):
                            tool1 = "%s_cr" % tool
                        elif (method == "open_ref" and (tool in tools["closed_ref"] or tool in tools["de_novo"])):
                            tool1 = "%s_or" % tool
                        elif (method == "de_novo" and (tool in tools["closed_ref"] or tool in tools["open_ref"])):
                            tool1 = "%s_dn" % tool
                        else:
                            tool1 = tool
                        labels.append(tool1)

        with open(sys.argv[i], 'U') as taxa_summary_f:
            for line in taxa_summary_f:
                if line.startswith("#"):
                    continue
                line = line.strip().split('\t')
                # sum abundance across all samples
                all_abundances = map(float, line[1:])
                num_samples = len(line[1:])
                total_abundance = sum(all_abundances)
                if not line[0] in abundances:
                    abundances[line[0]] = {}
                    # new taxa not found previously by any tool,
                    # add to dictionary and set all abundances for
                    # other tools to 0.0
                    for j in range(0, i-3):
                        abundances[line[0]][j] = 0.0
                abundances[line[0]][i-3] = total_abundance/float(num_samples)
        # current tool is missing some taxa previously found
        # by other tools, set abundance for this taxa to 0.0
        for taxa in abundances:
            if i-3 not in abundances[taxa]:
                abundances[taxa][i-3] = 0.0

#    for taxa in abundances:
#        sys.stdout.write("%s\t" % taxa)
#        for tool_abu in abundances[taxa]:
#            sys.stdout.write("%s\t" % abundances[taxa][tool_abu])
#        sys.stdout.write("\n")

    # find the min abundance for storing top N taxa
    min_abundance = {}
    for i in range(0, num_tools):
        tmp_abundance = []
        for taxa in abundances:
            tmp_abundance.append(abundances[taxa][i])
        tmp_abundance.sort(reverse=True)
        if (top_N_taxa <= len(tmp_abundance)):
            min_abundance[i] = tmp_abundance[top_N_taxa-1]
        else:
            min_abundance[i] = 0.0

    for i in range(0, num_tools):
        for taxa in abundances:
            if abundances[taxa][i] < min_abundance[i]:
                abundances[taxa][i] = 0.0

    total_colors = len(abundances)
    for taxa in abundances:
        tmp_taxa_abu_per_tool = []
        for tool_abu in abundances[taxa]:
            tmp_taxa_abu_per_tool.append(abundances[taxa][tool_abu])

        if all (v == 0.0 for v in tmp_taxa_abu_per_tool):
            total_colors -= 1

    if len(labels) % 2 != 0:
        end = len(labels)/2 + 0.5
    else:
        end = len(labels)/2

    ind = np.arange(1,end+1,0.5)
    width = 0.25
    plots = []
    all_taxa = []
    colors = cm.rainbow(np.linspace(0,1,total_colors))
    colors = np.random.permutation(colors)

    k = 0
    yoff = np.array([0.0] * len(labels))
    for taxa in abundances:
        taxa_abu_per_tool = []
        for tool_abu in abundances[taxa]:
            taxa_abu_per_tool.append(abundances[taxa][tool_abu]) 

        if not all(v == 0.0 for v in taxa_abu_per_tool):
            all_taxa.append(taxa)
            p = plt.bar(ind, taxa_abu_per_tool, width, bottom=yoff, color=colors[k])
            yoff+=taxa_abu_per_tool
            plots.append(p)
            k+=1

    plt.ylabel('Relative abundance [0,1]')
    plt.xlabel('Software')
    plt.xticks(ind+width/2., labels, size="small", rotation=60)
    plt.yticks(size="small")
    plt.grid(True, zorder=0, axis="y")
    legend = []
    for i in range(0,len(plots)): 
        legend.append(plots[i][0])

    lgd = plt.legend(legend, all_taxa, fontsize=8, bbox_to_anchor=(0.5, -1.3), loc='lower center')

    outdir = os.path.join(rootdir, gene)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    plt.savefig(os.path.join(
        outdir, "barchart_%s_%s_top_%s_.png" % (tool, study, top_N_taxa)), bbox_inches="tight", bbox_extra_artist=[lgd])
