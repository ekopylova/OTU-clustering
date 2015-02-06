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

	taxonomy_mean_f = sys.argv[1]
	taxonomy_stdev_f = sys.argv[2]
	method = sys.argv[3]
	study = sys.argv[4]

	tools = {"de_novo": ['sumaclust', 'uclust', 'usearch61', 'usearch', 'swarm', 'uparse_q16', 'uparse_q3'],
	         "closed_ref": ['sortmerna', 'uclust', 'usearch61', 'usearch'],
	         "open_ref": ['sortmerna_sumaclust', 'uclust', 'usearch61']}

	taxonomy_mean = {}
	taxonomy_stdev = {}

	# load taxonomy_mean array
	with open (taxonomy_mean_f, 'U') as taxonomy_mean_fp:
		for line in taxonomy_mean_fp:
			line = line.strip().split("\t")
			if line[0] not in taxonomy_mean:
				taxonomy_mean[line[0]] = []
				for value in line[1:]:
					taxonomy_mean[line[0]].append(float(value))
			else:
				print "ERROR: %s is already in taxonomy_mean" % line[0]

	# load taxonomy_stdev array
	with open (taxonomy_stdev_f, 'U') as taxonomy_stdev_fp:
		for line in taxonomy_stdev_fp:
			line = line.strip().split("\t")
			if line[0] not in taxonomy_stdev:
				taxonomy_stdev[line[0]] = []
				for value in line[1:]:
					taxonomy_stdev[line[0]].append(float(value))
			else:
				print "ERROR: %s is already in taxonomy_stdev" % line[0]

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
	        #p = ax.bar(ind+(i*width), taxonomy_mean[tool], width, color=colors[i], yerr=taxonomy_stdev[tool])
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

	plt.savefig("tp_plot_%s_%s.png" % (method, study) , bbox_inches="tight", bbox_extra_artist=[ax.legend])

if __name__ == "__main__":
    main(sys.argv[1:])
