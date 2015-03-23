#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Generate layered taxonomy barcharts for each OTU clustering tool
"""

import sys
import os
from os.path import isfile
import random
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
# requires matplot version 1.2.1 (or one that supports fontsize in the legend)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import brewer2mpl

def generate_taxa_barcharts(output_dir,
                            taxa_summary_dir,
                            datatypes,
                            studies,
                            studies_bac_mock,
                            studies_euk_mock,
                            methods,
                            tools,
                            top_N_taxa_mock,
                            top_N_taxa_env):
    """
    This method draws a layered bargraph for summarizing the top N abundant
    taxa recovered by each OTU clustering tool

    Parameters
    ----------
    output_dir : string
        dirpath to store output
    taxa_summary_dir : string
        dirpath to summarized taxonomies (not including singleton OTUs)
    datatypes : list
        list of genes to analyze (16S, 18S)
    studies : dictionary
        dictionary with keys being a datatype (16S or 18S) and the values
        being lists of studies to analyze
    studies_bac_mock : list
        list of mock bacterial study IDs
    studies_euk_mock : list
        list of mock eukaryotic study IDs
    methods : list
        list of OTU clustering methods (closed_ref, de_novo, open_ref)
    tools : dictionary
        dictionary with keys being a method and the values being a list
        of tools designed to perform that method
    top_N_taxa_mock : integer
        number of top abundant taxonomies to illustrate for mock communities
    top_N_taxa_env : integer
        number of top abundant taxonomies to illustrate for environmental
        communities

    Returns
    -------
    None
    """

    # ex. 16S
    for datatype in datatypes:
        # ex. 1685
        for study in studies[datatype]:
            if study == "nematodes":
                summary_level = "L5"
            else:
                summary_level = "L6"
            abundances = {}
            num_tools = 0
            i=0
            labels = []
            # ex. closed_ref
            for method in methods:
                # ex. swarm
                for tool in tools[method]:
                    tax_sum_fp = os.path.join(taxa_summary_dir, datatype, method, "%s_%s" % (tool, study), "otu_table_mc2_%s.txt" % summary_level)
                    if isfile(tax_sum_fp):
                        num_tools+=1
                        with open(tax_sum_fp, 'U') as taxa_summary_f:
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
                                    for j in range(0, i):
                                        abundances[line[0]][j] = 0.0
                                abundances[line[0]][i] = total_abundance/float(num_samples)
                        # current tool is missing some taxa previously found
                        # by other tools, set abundance for this taxa to 0.0
                        for taxa in abundances:
                            if i not in abundances[taxa]:
                                abundances[taxa][i] = 0.0

                        i+=1
                    else:
                        print "skipping %s does not exist" % tax_sum_fp
                        continue

                    if (method == "closed_ref" and (tool in tools["de_novo"] or tool in tools["open_ref"])):
                        tool1 = "%s_cr" % tool
                    elif (method == "open_ref" and (tool in tools["closed_ref"] or tool in tools["de_novo"])):
                        tool1 = "%s_or" % tool
                    elif (method == "de_novo" and (tool in tools["closed_ref"] or tool in tools["open_ref"])):
                        tool1 = "%s_dn" % tool
                    else:
                        tool1 = tool
                    labels.append(tool1)

            # go to next study
            if num_tools == 0:
                continue

#            for taxa in abundances:
#                sys.stdout.write("%s\t" % taxa)
#                for tool_abu in abundances[taxa]:
#                    sys.stdout.write("%s\t" % abundances[taxa][tool_abu])
#                sys.stdout.write("\n")

            if (study in studies_bac_mock or study in studies_euk_mock):
                top_N_taxa = top_N_taxa_mock
            else:
                top_N_taxa = top_N_taxa_env

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

            outdir = os.path.join(output_dir, datatype)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            plt.savefig(os.path.join(
                outdir, "barchart_%s_top_%s.pdf" % (study, top_N_taxa)), bbox_inches="tight", bbox_extra_artist=[lgd])
            plt.clf()


def main(argv):

    # output directory
    output_dir = sys.argv[1]

    # summarized taxa (BIOM tables) directory
    taxa_summary_dir = sys.argv[2]

    # tools
    tools_denovo = sys.argv[3].split()
    tools_closed_ref = sys.argv[4].split()
    tools_open_ref = sys.argv[5].split()

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

    # studies 16S
    studies_bac_mock = sys.argv[6].split()
    studies_bac_env = sys.argv[7].split()
    studies_bac = studies_bac_mock + studies_bac_env

    # studies 18S
    studies_euk_mock = sys.argv[8].split()
    studies_euk_env = sys.argv[9].split()
    studies_euk = studies_euk_mock + studies_euk_env

    # list of studies for each gene type
    studies = {'16S': [], '18S': []}

    for study in studies_bac:
      if study not in studies['16S']:
        studies['16S'].append(study)

    for study in studies_euk:
      if study not in studies['18S']:
        studies['18S'].append(study)

    # genes
    datatypes = ['16S', '18S']

    # OTU picking methods
    methods = ['de_novo', 'closed_ref', 'open_ref']

    top_N_taxa_mock = int(sys.argv[10])
    top_N_taxa_env = int(sys.argv[11])

    generate_taxa_barcharts(
        output_dir=output_dir,
        taxa_summary_dir=taxa_summary_dir,
        datatypes=datatypes,
        studies=studies,
        studies_bac_mock=studies_bac_mock,
        studies_euk_mock=studies_euk_mock,
        methods=methods,
        tools=tools,
        top_N_taxa_mock=top_N_taxa_mock,
        top_N_taxa_env=top_N_taxa_env)


if __name__ == "__main__":
    main(sys.argv[1:])
