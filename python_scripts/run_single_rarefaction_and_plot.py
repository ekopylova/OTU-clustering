#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
For studies with multiple samples, rarefy the OTU tables at different depths
(using information from run_summarize_tables.py), run alpha diversity for
each new OTU table/depth and plot the results.

Dependencies: QIIME 1.9.0
usage:        python run_single_rarefaction_and_plot.py
"""

import sys
import os
from subprocess import Popen, PIPE
from glob import glob

from numpy import arange, array, mean, zeros 
import matplotlib as mpl
import sys

mpl.use('agg')

import matplotlib.pyplot as plt
import brewer2mpl
 

def main():

    # original BIOM tables excl. singletons directory
    rootdir = sys.argv[1]

    # output results directory 
    outdir_root = sys.argv[2]
    
    # trees (generated during de novo/open ref OTU picking, will be in
    # same directory as all results)
    treedir = sys.argv[3]

    # trees for closed-reference analysis (16S distributed with Greengenes, 18S built)
    tree_fp = {'16S': "", '18S': ""}
    tree_fp['16S'] = sys.argv[4]
    tree_fp['18S'] = sys.argv[5]

    # mapping files directory
    mapping_dir = sys.argv[6]

    # studies 16S
    studies_bac = sys.argv[7].split()

    # studies 18S
    studies_euk = sys.argv[8].split()

    # list of studies for each gene type
    studies = {'16S': [], '18S': []}

    for study in studies_bac:
      if study not in studies['16S']:
        studies['16S'].append(study)

    for study in studies_euk:
      if study not in studies['18S']:
        studies['18S'].append(study)

    # tools
    tools_denovo = sys.argv[9].split()
    tools_closed_ref = sys.argv[10].split()
    tools_open_ref = sys.argv[11].split()

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

    # OTU picking method
    methods = ['de_novo', 'closed_ref', 'open_ref']

    # sampling depths for each study (determined using output of run_summarize_tables.py)
    sample_depth_dn = {'449': ['100', '200', '380'],
                       '632': ['1000', '80000', '170000'],
                       '1685': ['1000', '539800', '1000000'],
                       '1686': ['1000', '471000', '800000'],
                       '2107': ['1000', '100000', '200000']}
    sample_depth_cr = {'449': ['300', '1000', '2000'],
                       '632': ['1000', '140000', '590000'],
                       '1685': ['1000','500000', '1000000'],
                       '1686': ['1000', '471000', '1000000'],
                       '2107': ['100', '1000', '170000']}
    sample_depth_or = {'449': ['300', '1000', '2000'],
                       '632': ['1000', '80000', '170000'],
                       '1685': ['1000', '539800', '1000000'],
                       '1686': ['1000', '471000', '800000'],
                       '2107': ['1000', '200000', '400000']}

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in methods:
            otu_table = "otu_table_mc2.biom"
            # ex. 1685
            for study in studies[datatype]:
                if (method == "closed_ref" and study not in sample_depth_cr[study]):
                    continue
                elif (method == "de_novo" and study not in sample_depth_dn[study]):
                    continue
                elif (method == "open_ref" and study not in sample_depth_or[study]):
                    continue
                # ex. swarm
                for tool in tools[method]:
                    search_dir = os.path.join(rootdir, datatype, method, "%s_%s" % (tool, study)) 
                    # OTU table doesn't exist (tool probably failed on this study)
                    if not os.path.isfile(os.path.join(search_dir, otu_table)):
                        print "skipping %s does not exist" % os.path.join(search_dir, otu_table)
                        continue
                    outdir = os.path.join(outdir_root, datatype, method, "%s_%s" % (tool, study))
                    if not os.path.exists(outdir):
                        os.makedirs(outdir)

                        if (method == "closed_ref" and sample_depth_cr[study] is not ''):
                            sample_depths = sample_depth_cr[study]
                        elif (method == "open_ref" and sample_depth_or[study] is not ''):
                            sample_depths = sample_depth_or[study]
                        elif (method == "de_novo" and sample_depth_dn[study] is not ''):
                            sample_depths = sample_depth_dn[study]

                        for depth in sample_depths:
                            single_rare_command = ["single_rarefaction.py",
                                                "-i",
                                                os.path.join(search_dir, otu_table),
                                                "-o",
                                                os.path.join(outdir, "otu_table_even%s.biom" % depth),
                                                "-d",
                                                depth]
                            #print "command = ", single_rare_command
                            proc = Popen(single_rare_command,
                                         stdout=PIPE,
                                         stderr=PIPE,
                                         close_fds=True)
                            proc.wait()
                            stdout, stderr = proc.communicate()
                            if stderr:
                                print stderr
                
                # run alpha diversity on rarefied tables
                for tool in tools[method]:
                    search_dir = os.path.join(outdir_root, datatype, method, "%s_%s" % (tool, study))
                    if not os.path.exists(search_dir):
                        continue

                    #print "outdir = ", search_dir

                    if (method == "closed_ref" and sample_depth_cr[study] is not ''):
                        sample_depths = sample_depth_cr[study]
                        tree = tree_fp[datatype]
                    elif (method == "open_ref" and sample_depth_or[study] is not ''):
                        sample_depths = sample_depth_or[study]
                        tree = os.path.join(treedir, datatype, method, "%s_%s" % (tool, study), "rep_set.tre")
                    elif (method == "de_novo" and sample_depth_dn[study] is not ''):
                        sample_depths = sample_depth_dn[study]
                        tree = os.path.join(treedir, datatype, method, "%s_%s" % (tool, study), "rep_set.tre")

                    for depth in sample_depths:
                        # output file doesn't already exist
                        if not os.path.isfile(os.path.join(search_dir, "alpha_div_even%s.txt" % depth)):
                            alpha_div_command = ["alpha_diversity.py",
                                                 "-i",
                                                 os.path.join(search_dir, "otu_table_even%s.biom" % depth),
                                                 "-o",
                                                 os.path.join(search_dir, "alpha_div_even%s.txt" % depth),
                                                 "-t",
                                                 tree]
                            proc = Popen(alpha_div_command,
                                         stdout=PIPE,
                                         stderr=PIPE,
                                         close_fds=True)
                            proc.wait()
                            stdout, stderr = proc.communicate()
                            if stderr:
                                print stderr

                # For each depth, plot graphs
                if (method == "closed_ref" and sample_depth_cr[study] is not ''):
                    sample_depths = sample_depth_cr[study]
                elif (method == "open_ref" and sample_depth_or[study] is not ''):
                    sample_depths = sample_depth_or[study]
                elif (method == "de_novo" and sample_depth_dn[study] is not ''):
                    sample_depths = sample_depth_dn[study]

                # info to generate separated box plot
                tools_not_to_plot = []
                for tool in tools[method]:
                    for depth in sample_depths:
                        file_s = os.path.join(
                            outdir_root, datatype, method, "%s_%s" % (tool, study), "alpha_div_even%s.txt" % depth)
                        if not os.path.isfile(file_s):
                            print "skipping %s does not exist" % file_s
                            tools_not_to_plot.append(tool)
                            break
                    
                for tool in tools_not_to_plot:
                    tools[method].remove(tool)

                num_groups = len(tools[method])
                num_time_points = len(sample_depths)
                labels = tools[method]
                #print "labels = ",labels
                labels_to_remove = []
                xticklabels = sample_depths
                data_oo = zeros((num_groups, num_time_points), dtype=list)
                data_pd = zeros((num_groups, num_time_points), dtype=list)

                for i in xrange(len(tools[method])):
                    for j in xrange(len(sample_depths)):
                        tool = tools[method][i]
                        depth = sample_depths[j]
                        search_dir = os.path.join(outdir_root, datatype, method, "%s_%s" % (tool, study))
                        file_s = os.path.join(search_dir, "alpha_div_even%s.txt" % depth)
                        #print "file = %s" % file_s
                        pd_list = []
                        oo_list = []
                        with open(file_s, 'U') as collection:
                            next(collection)
                            for line in collection:
                                pd_whole_tree = float(line.strip().split()[1])
                                observed_otus = float(line.strip().split()[3])
                                pd_list.append(pd_whole_tree)
                                oo_list.append(observed_otus)
                            data_oo[i][j] = oo_list
                            data_pd[i][j] = pd_list

                # observed OTUs
                fig = plt.figure()
                ax = fig.add_subplot(111)
                #ax.set_title('Study %s: Alpha rarefaction at various depths (sequences/sample)' % study)
                ax.set_xticklabels(xticklabels, rotation=40, ha='center')
                ax.set_ylabel('observed_otus')
                ax.set_xlabel('sampling depths (# of samples)')
                lgd = make_separated_box(ax, data_oo, labels, xticklabels=xticklabels, legend_pos='lower right')
                lgd = None

                plots_dir = os.path.join(outdir_root, datatype, method, "plots_%s" % study)
                if not os.path.exists(plots_dir):
                    os.makedirs(plots_dir)
                fig.savefig(os.path.join(plots_dir, 'observed_otus.png'), bbox_inches='tight')
                # pd_whole_tree
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_xticklabels(tools[method], rotation=40, ha='center')
                ax.set_ylabel("Faith's phylogenetic diversity (PD_whole_tree)")
                ax.set_xlabel('sampling depth (# of samples)')
                lgd = make_separated_box(ax, data_pd, labels, xticklabels=xticklabels, legend_pos='lower right')
                
                fig.savefig(os.path.join(plots_dir, 'pd_whole_tree.png'), bbox_extra_artists=(lgd,), bbox_inches='tight')

def color_bp(bp, color):
    c = array(color) * 0.5
    c = tuple(c)

    for x in bp['boxes']: 
        plt.setp(x, color=c)
        x.set_facecolor(color)
    for x in bp['medians']:
        plt.setp(x, color=c)
    for x in bp['whiskers']: 
        plt.setp(x, color=c)
    for x in bp['fliers']: 
        plt.setp(x, color=c)
    for x in bp['caps']:
        plt.setp(x, color=c)

# code (modified) from https://github.com/samfway/biotm/blob/master/plotting/grouped_box.py
def make_separated_box(ax, data, labels=None, colors=None,
                       xticklabels=[], width=0.8, legend_pos=0,
                       dot_mean=False, mean_color='w'):
    if labels and len(data) != len(labels):
        raise ValueError('Number of labels must match ',
                         'size of data matrix.')

    if colors and len(colors) != len(labels):
        raise ValueError('Number of colors must match ',
                         'size of data matrix.')

    num_groups = len(labels)
    num_points = data.shape[1]

    if not colors:
        num_colors = max(3, num_groups)
        colors = brewer2mpl.get_map('Set2', 
                                    'qualitative', 
                                    num_colors).mpl_colors
    current_pos = 0
    xticks = []
    xlabels = []

    for i in xrange(num_groups):
        color = colors[i]
        for j in xrange(num_points):
            bp = ax.boxplot(data[i][j], positions=[current_pos],
                            widths=[width], patch_artist=True)
            xticks.append(current_pos)
            xlabels.append(xticklabels[j] + " (%s)" % len(data[i][j])) 
            color_bp(bp, color)
            if dot_mean:
                means = [mean(data[i][j])]
                ax.plot([current_pos], means, linestyle='None', 
                    marker='o', markerfacecolor=mean_color,
                    markeredgecolor='k')
            current_pos += 1.6 
        current_pos += 2

    if labels:
        lgd = legend_hack(ax, labels, colors, legend_pos)
    else:
        lgd = None

    ax.set_xlim(-1,current_pos-2)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    return lgd

def legend_hack(ax, labels, colors, legend_pos):
    """ Hack a legend onto a plot. 
    """ 
    handles = []
    for i, l in enumerate(labels):
        temp = plt.Line2D(range(1), range(1), 
                          linewidth=2,
                          color=colors[i])
        handles.append(temp)
    #plt.legend(handles, labels, numpoints=1, loc=legend_pos)
    lgd = plt.legend(handles, labels, numpoints=1, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #lgd = plt.legend(handles, labels, numpoints=1, bbox_to_anchor=(0, -0.7), loc='lower left', borderaxespad=0.)
    for handle in handles:
        handle.set_visible(False)

    return lgd 
                                
if __name__ == '__main__':
    main()
