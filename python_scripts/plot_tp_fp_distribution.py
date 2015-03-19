#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Use the data (taxonomy_mean_f and taxonomy_stdev_f) output by
run_compute_precision_recall.py to generate graphs.

usage:        python plot_tp_fp_distribution.py taxonomy_mean_f \
			  taxonomy_stdev_f [de_novo, closed_ref, open_ref] study
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import brewer2mpl
import sys

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

    # ex. 16S
    for datatype in datatypes:
    	# ex. 1685
    	for study in studies[datatype]:
    		# ex. de novo
    		for method in methods:
    			taxonomy_mean_fp = os.path.join(taxonomy_summary_dir, datatype, method, "%s_taxonomy_mean.txt" % study)
    			taxonomy_stdev_fp = os.path.join(taxonomy_summary_dir, datatype, method, "%s_taxonomy_stdev.txt" % study)

				taxonomy_mean = {}
				taxonomy_stdev = {}

				# load taxonomy_mean array
				if os.path.exists(taxonomy_mean_fp):
					with open (taxonomy_mean_f, 'U') as taxonomy_mean_fp:
						for line in taxonomy_mean_fp:
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
				if os.path.exists(taxonomy_stdev_fp):
					with open (taxonomy_stdev_f, 'U') as taxonomy_stdev_fp:
						for line in taxonomy_stdev_fp:
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
				    #	p = ax.bar(ind+(i*width), taxonomy_mean[tool], width, color=colors[i], yerr=taxonomy_stdev[tool])
				    #else:
					#	p = ax.bar(ind+(i*width), taxonomy_mean[tool], width, color=colors[i])
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

				plt.savefig(os.path.join(output_dir, "tp_plot_%s_%s.png" % (method, study)), bbox_inches="tight", bbox_extra_artist=[ax.legend])
				plt.clf()


if __name__ == "__main__":
    main(sys.argv[1:])
