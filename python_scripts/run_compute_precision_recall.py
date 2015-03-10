#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
This script computes the true-positive (TP), false-positive (FP),
false-negative (FN), precision, recall, F-measure, FP-chimeric
(taxa fully comprising of chimeric OTUs), FP-known (taxa fully
comprising of OTUs mapping to BLAST's NT database with >=97% id
and coverage and FP-other (taxa fully comprising of OTUs mapping
to BLAST's NT database with <97% id and coverage.

Dependencies: QIIME 1.9.0, BIOM-format >=2.1.3, <2.2.0, BLAST 2.2.29+,
              Blast NT database indexed, USEARCH Uchime (7.0.1090)
usage:        python run_compute_precision_recall.py [16S, 18S] study \
              taxonomy_level expected_taxa.txt
"""

import sys
import os
from biom import load_table
from subprocess import Popen, PIPE
from skbio.parse.sequences import parse_fasta
import copy
import numpy as np 
import matplotlib.pyplot as plt
import brewer2mpl


def graph_abundance_func(true_positive_otus,
                         false_positive_known_otus,
                         false_positive_other_otus,
                         false_positive_chimeric_otus,
                         datatype,
                         tool,
                         study,
                         method,
                         results_dir,
                         taxonomy_mean,
                         taxonomy_stdev):
    '''Function to build taxonomy_mean and taxonomy_stdev dictionaries for
       storing TP, FP-chimeric, FP-known and FP-other mean number of reads
       and stdev. Accessed once per gene/method/tool/study. At the end of main(),
       these dictionaries are written to files which can be passed to
       OTU-picking/python_scripts/plot_tp_fp_distribution.py to generate
       graphs.
    '''

    total_taxa = len(false_positive_known_otus) + len(false_positive_other_otus) \
        + len(false_positive_chimeric_otus)

    if (method == "de_novo" and (tool == "uparse_q3" or tool == "uparse_q16")):
            otu_map_f = os.path.join(
                results_dir, datatype, method, "%s_%s" % (tool, study), "seqs_otus.txt")
    elif (method == "closed_ref" and (tool == "uclust" or tool == "usearch" or tool == "usearch61")):
            otu_map_f = os.path.join(
                results_dir, datatype, method, "%s_%s" % (tool, study), "%s_ref_picked_otus" % tool, "seqs_otus.txt")
    elif method == "open_ref":
            otu_map_f = os.path.join(
                results_dir, datatype, method, "%s_%s" % (tool, study), "final_otu_map_mc2.txt")
    else:
        otu_map_f = os.path.join(
            results_dir, datatype, method, "%s_%s" % (tool, study), "%s_picked_otus" % tool, "seqs_otus.txt")

    # load OTU map into dict
    otu_map_dict = {}
    with open (otu_map_f, 'U') as otu_map_fp:
        for line in otu_map_fp:
            line = line.strip().split("\t")
            if line[0] not in otu_map_dict:
                otu_map_dict[line[0]] = line[1:]
            else:
                print "ERROR: %s is already in dict" % line[0]
                exit(1)

    if tool not in taxonomy_mean:
        taxonomy_mean[tool] = []
    if tool not in taxonomy_stdev:
        taxonomy_stdev[tool] = []

    # TP
    total_reads = []
    for taxa in true_positive_otus:
        t = 0
        for otu_id in true_positive_otus[taxa]:
            if otu_id in otu_map_dict:
                t += len(otu_map_dict[otu_id])
        total_reads.append(t)
    arr = np.array(total_reads)
    taxonomy_mean[tool].append(np.rint(np.nan_to_num(np.mean(arr, axis=0))))
    taxonomy_stdev[tool].append(np.rint(np.nan_to_num(np.std(arr, axis=0))))

    # FP-known
    total_reads = []
    for taxa in false_positive_known_otus:
        t = 0
        for otu_id in false_positive_known_otus[taxa]:
            if otu_id in otu_map_dict:
                t += len(otu_map_dict[otu_id])            
        total_reads.append(t)
    arr = np.array(total_reads)
    taxonomy_mean[tool].append(np.rint(np.nan_to_num(np.mean(arr, axis=0))))
    taxonomy_stdev[tool].append(np.rint(np.nan_to_num(np.std(arr, axis=0))))

    # FP-other
    total_reads = []
    for taxa in false_positive_other_otus:
        t = 0
        for otu_id in false_positive_other_otus[taxa]:
            if otu_id in otu_map_dict:
                t += len(otu_map_dict[otu_id])            
        total_reads.append(t)
    arr = np.array(total_reads)
    taxonomy_mean[tool].append(np.rint(np.nan_to_num(np.mean(arr, axis=0))))
    taxonomy_stdev[tool].append(np.rint(np.nan_to_num(np.std(arr, axis=0))))

    # FP-chimeric
    total_reads = []
    for taxa in false_positive_chimeric_otus:
        t = 0
        for otu_id in false_positive_chimeric_otus[taxa]:
            if otu_id in otu_map_dict:
                t += len(otu_map_dict[otu_id])
        total_reads.append(t)
    arr = np.array(total_reads)
    taxonomy_mean[tool].append(np.rint(np.nan_to_num(np.mean(arr, axis=0))))
    taxonomy_stdev[tool].append(np.rint(np.nan_to_num(np.std(arr, axis=0))))

    print "taxonomy_mean = ", taxonomy_mean
    print "taxonomy_stdev = ", taxonomy_stdev


def compute_fp_other(results_dir,
                     out_dir,
                     filter_otus_dir,
                     chimera_db_18S,
                     chimera_db_16S,
                     taxonomy_mean,
                     taxonomy_stdev,
                     blast_nt_index,
                     actual_tax,
                     expected_tax,
                     tool,
                     study,
                     datatype,
                     method,
                     tax_level,
                     graph_abundance=True):
    '''This function completes the following steps,
         1. Load original OTU table (excl. singleton OTUs) into dict with taxonomies as
            keys and OTUs representing them as values in a list (only L5 and L6 supported)
         2. Load all true-positive taxonomies as keys in a dict and all OTUs representing
            them (known from dict in 1.) as values in a list
         3. Load all false-positive taxonomies as keys in a dict and all OTUs representing
            them (known from dict in 1.) as values in a list
         4. Using results from 3. output all false-positive OTUs to a FASTA file and
            run UCHIME chimera filter on them
         5. Remove all chimeric OTUs from 3., if all OTUs representing a taxa have
            been removed, count this taxa chimeric (FP-chimeric)
         6. Write all non-chimeric OTUs to a file and run MEGABLAST on them against
            BLAST's NT database
         7. Remove all OTUs from 3. if mapping with >=97% identity and coverage to BLAST's
            NT database. If all OTUs representing a taxa have been removed, count
            this taxa as a known species (FP-known)
         8. The remaining taxa in 3. comprise of OTUs mapping with <97% identity and coverage
            to BLAST's NT database (FP-other)
         9. Pass the true-positive dict, FP-chimeric dict, FP-known dict and FP-other dict
            to graph_abundance_func() to compute the mean number of reads representing each
            taxa in each of those dicts.
    '''

    fp_chimera = 0
    fp_known = 0
    fp_other = 0

    # load taxonomies from OTU table as keys into dictionary,
    # and all OTU ids that share that taxonomy as values in a list
    otu_table_dict = {}
    biom_table_f = os.path.join(
        filter_otus_dir, datatype, method, "%s_%s" % (tool, study), "otu_table_mc2.biom")
    if not os.path.exists(biom_table_f):
        print "%s does not exist, cannot search for contaminants" % biom_table_f
    else:
        print "loading %s" % biom_table_f
        biom_table = load_table(biom_table_f)
        obs_ids_list = biom_table._observation_ids
        for obs_id in obs_ids_list:
            obs_data = biom_table.metadata(obs_id,'observation')
            # taxonomy must be up to genus level (L6)
            if tax_level == "L6":
                if len(obs_data['taxonomy']) == 7:
                    del obs_data['taxonomy'][-1]
            # taxonomy must be up to family level (L5)
            elif tax_level == "L5":
                if len(obs_data['taxonomy']) > 5:
                    obs_data['taxonomy'] = obs_data['taxonomy'][0:5]
            else:
                print "ERROR: taxonomy level %s not supported" % tax_level
                exit(1)
            assignment = ";".join(obs_data['taxonomy'])
            if assignment not in otu_table_dict:
                otu_table_dict[assignment] = [obs_id]
            else:
                otu_table_dict[assignment].append(obs_id)

    # collect all OTUs representing summarized taxa true positive matches
    # (this info is for the graph)
    true_positive_otus = {}
    lis_tp = actual_tax & expected_tax
    for l in lis_tp:
        if l not in otu_table_dict:
            print "Error: TP %s not in OTU table" % l
            exit(1)
        else:
            if l not in true_positive_otus:
                true_positive_otus[l] = otu_table_dict[l]
            else:
                true_positive_otus[l].extend(otu_table_dict[l])

    # collect all OTUs representing summarized taxa false positive matches
    false_positive_otus = {}
    lis = actual_tax - expected_tax
    for l in lis:
        if ";Other" in l:
            l = l.replace(";Other", "")
        if l not in otu_table_dict:
            print "Error: FP %s not in OTU table" % l
            exit(1)
        else:
            if l not in false_positive_otus:
                false_positive_otus[l] = otu_table_dict[l]
            else:
                false_positive_otus[l].extend(otu_table_dict[l])

    print "Total true positive taxa = ", len(true_positive_otus)
    print "Total false positive taxa = ", len(false_positive_otus)
    false_positive_otus_count = 0
    for s in false_positive_otus:
        false_positive_otus_count += len(false_positive_otus[s])
    print "Total false positive otus = ", false_positive_otus_count
    true_positive_otus_count = 0
    for s in true_positive_otus:
        true_positive_otus_count += len(true_positive_otus[s])
    print "Total true positive otus = ", true_positive_otus_count

    # create list file of false positive OTUs
    fp_otus_ids_f = os.path.join(out_dir, datatype, method, "%s_%s_fp_ids.txt" % (tool, study))
    with open(fp_otus_ids_f, 'w') as out_file:
        for tax in false_positive_otus:
            for otu in false_positive_otus[tax]:
                out_file.write("%s\n" % otu)

    # create FASTA file of false positive OTUs                                                                                                                                                                                                                   
    if method == "closed_ref":
        rep_set_fasta = os.path.join(
            results_dir, datatype, method, "%s_%s" % (tool, study), "rep_set", "seqs_rep_set.fasta")
    elif method == "de_novo":
        if (tool == "uparse_q3" or tool == "uparse_q16"):
            rep_set_fasta = os.path.join(
                results_dir, datatype, method, "%s_%s" % (tool, study), "otus.fa")
        else:
            rep_set_fasta = os.path.join(
                results_dir, datatype, method, "%s_%s" % (tool, study), "rep_set", "seqs_rep_set.fasta")
    else:
        rep_set_fasta = os.path.join(
            results_dir, datatype, method, "%s_%s" % (tool, study), "rep_set.fna")

    otus_fasta = os.path.join(
        out_dir, datatype, method, "%s_%s_fp.fasta" % (tool, study))
    filter_fasta_command = ["filter_fasta.py",
                            "-f", rep_set_fasta, "-o", otus_fasta,
                            "-s", fp_otus_ids_f]
    print "command = ", filter_fasta_command
    proc = Popen(filter_fasta_command, stdout=PIPE, stderr=PIPE, close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if stderr:
        print stderr

    # search for chimeras in all false positive OTUs using UCHIME
    if datatype == "16S":
        chimera_db = chimera_db_16S
    elif datatype == "18S":
        chimera_db = chimera_db_18S
    else:
        raise ValueError("%s not supported" % datatype)

    chimeric_otus_fasta = os.path.join(
        out_dir, datatype, method, "%s_%s_fp_chimeras.fasta" % (tool, study))
    uchime_command = ["usearch70", "-uchime_ref", otus_fasta,
                      "-db", chimera_db,
                      "-strand", "plus", "-chimeras", chimeric_otus_fasta]
    print "command = ", uchime_command
    proc = Popen(uchime_command, stdout=PIPE, stderr=PIPE, close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    #if stderr:
    #    print stderr

    # get list of chimeric OTU ids
    chimeric_ids = []
    with open(chimeric_otus_fasta, "U") as identified_chimeras_fp:
        for label, seq in parse_fasta(identified_chimeras_fp):
            chimeric_ids.append(label)

    print "Total chimeric OTUs identified in false-positive set: %s" % len(chimeric_ids)

    false_positive_chimeric_otus = copy.deepcopy(false_positive_otus)
    chimeric_taxa = []

    # remove chimeric OTU ids from false_positive_otus and if list becomes
    # empty for a taxa, count taxa as formed by chimeric OTUs
    for chimera in chimeric_ids:
        for s in false_positive_otus:
            if chimera in false_positive_otus[s]:
                false_positive_otus[s].remove(chimera)
                if len(false_positive_otus[s]) == 0:
                    fp_chimera += 1
                    chimeric_taxa.append(s)
    print "Total chimeric taxa: ", fp_chimera

    for taxa in false_positive_otus:
        if len(false_positive_otus[taxa]) != 0:
            del false_positive_chimeric_otus[taxa]

    # write remaining non-chimeric false positive OTU ids into file
    fp_otus_non_chimeric_ids = os.path.join(
        out_dir, datatype, method, "%s_%s_fp_otus_non_chimeric.txt" % (tool, study))
    with open(fp_otus_non_chimeric_ids, "w") as fp_otus_non_chimeric_ids_fp:
        for s in false_positive_otus:
            for otu_id in false_positive_otus[s]:
                fp_otus_non_chimeric_ids_fp.write("%s\n" % otu_id)

    # get FASTA file of non-chimeric false positive OTU ids
    fp_otus_non_chimeric_fasta = os.path.join(
        out_dir, datatype, method, "%s_%s_fp_otus_non_chimeric.fasta" % (tool, study))
    filter_fasta_command = ["filter_fasta.py", "-f", rep_set_fasta,
                            "-o", fp_otus_non_chimeric_fasta, "-s",
                            fp_otus_non_chimeric_ids]
    print "command = ", filter_fasta_command
    proc = Popen(filter_fasta_command, stdout=PIPE, stderr=PIPE, close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if stderr:
        print stderr

    # megablast all OTUs against nt database
    otus_blast = os.path.join(
        out_dir, datatype, method, "%s_%s_fp_non_chimeric.blast" % (tool, study))
    blast_command = ["blastn",
                     "-task", "megablast",
                     "-db", blast_nt_index,
                     "-query", fp_otus_non_chimeric_fasta,
                     "-out", otus_blast,
                     "-evalue", "1e-5",
                     "-outfmt", "6 std qcovs",
                     "-max_target_seqs", "1",
                     "-num_threads", "50"]
    print "command = ", blast_command
    proc = Popen(blast_command, stdout=PIPE, stderr=PIPE, close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if stderr:
        print stderr                

    # parse Blast output to collect OTU ids that mapped with >=97% id and 97% coverage
    # to some sequence in the nt database
    contaminants = set()
    with open (otus_blast, "U") as blast_in:
        for line in blast_in:
            line = line.strip().split("\t")
            if (float(line[2]) >= 97.0 and float(line[12]) >= 97.0):
                if (tool == "uparse_q3"or tool== "uparse_q16"):
                    line[0] = "%s;" % line[0]
                contaminants.add(line[0])

    print "Total known OTUs (map with >=97%% id and 100%% coverage) = ", len(contaminants)
    print "Total other OTUs (map with <97%% id and <100%% coverage) = ", false_positive_otus_count-len(chimeric_ids)-len(contaminants)

    # Create a copy of FP - chimera. In the next step,
    # FP - chimera also removes known taxa from the dictionary,
    # we need the original FP - chimera to also keep *only*
    # the known taxa (FP-known)
    false_positive_known_otus = copy.deepcopy(false_positive_otus)

    # remove mapped OTU ids from false positive OTUs
    for id_otu in contaminants:
        for tax in false_positive_otus:
            if id_otu in false_positive_otus[tax]:
                false_positive_otus[tax].remove(id_otu)
                # remove taxa with empty OTU id lists from false positive taxa list
                if len(false_positive_otus[tax]) == 0:
                    fp_known+=1

    for taxa in false_positive_otus:
        if (len(false_positive_otus[taxa]) != 0 or taxa in chimeric_taxa):
            del false_positive_known_otus[taxa]

    false_positive_other_otus = copy.deepcopy(false_positive_otus)
    for taxa in false_positive_otus:
        if len(false_positive_otus[taxa]) == 0:
            del false_positive_other_otus[taxa]

    print "Total known taxa = ", fp_known

    for s in false_positive_otus:
        if len(false_positive_otus[s]) != 0:
            fp_other += 1
    print "Total other taxa = ", fp_other

    # make bar charts for TP/FP-known/FP-other 
    if graph_abundance:
        # TP, FP-known, FP-other
        graph_abundance_func(
            true_positive_otus, false_positive_known_otus, false_positive_other_otus, \
            false_positive_chimeric_otus, datatype, tool, study, method, \
            results_dir, taxonomy_mean, taxonomy_stdev)

    return fp_chimera, fp_known, fp_other

def main(argv):
    '''This function loads the expected taxonomy composition and the
       summarized taxonomies for the same level and computes
       the number of true-positive (TP), false-positive (FP),
       false-negative (FN), precision, recall, F-measure and
       number of OTUs and output the statistics to stdout for
       each method and tool.
    '''

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
    #
    results_dir = sys.argv[1]

    # summarized taxonomies directory
    # (same outdir_root path as in run_summarize_taxa.py)
    summarize_taxa_dir = sys.argv[2]

    # output directory
    out_dir = sys.argv[3]

    # chimera database for 18S
    chimera_db_18S = sys.argv[4]

    # chimera database for 16S
    chimera_db_16S = sys.argv[5]

    # path to blast NT indexed database
    blast_nt_index = sys.argv[6]

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

    # filepath to expected summarized taxonomies
    expected_fp = sys.argv[12]

    # OTU-picking methods
    methods = ['de_novo', 'closed_ref', 'open_ref']

    # genes 
    datatypes = ['16S', '18S']

    singletons_removed = True

    # ex. 16S
    for datatype in datatypes:
        # ex. closed_ref
        for method in methods:
            # ex. swarm
            for tool in tools[method]:
                sys.stdout.write("%s\n" % method)
                # ex. 1685
                for study in studies[datatype]:
                    sys.stdout.write("%s\n" % study)
                    if study == "nematodes":
                        tax_level = "L5"
                    else:
                        tax_level = "L6"
                    expected_f = os.path.join(expected_fp, study, "%s_%s.txt" % (study, tax_level))

                    # these lists will contain mean number of reads + stdev for TP, FP-known,
                    # FP-other, FP-chimeric groups ex. taxonomy_mean = {"sumaclust": [4,76,21,53],
                    # "swarm": [489,2,32,4,3], ..}
                    taxonomy_mean = {}
                    taxonomy_stdev = {}

                    expected_tax = set()

                    # do not need to find contaminants for simulated communities
                    if (study == "even" or study == "staggered"):
                        find_contaminants = False
                        graph_abundance = False

                    with open(expected_f, 'U') as expected:
                        for line in expected:
                            if line.startswith('#'):
                                continue
                            else:
                                tax = line.split()[0]
                                expected_tax.add(tax)

                    actual_tax = set()
                    if singletons_removed:
                        otu_table = "otu_table_mc2_%s.txt" % tax_level
                    else:
                        otu_table = "otu_table_%s.txt" % tax_level

                    # create output directory for 
                    if not os.path.exists(os.path.join(out_dir, datatype, method)):
                        os.makedirs(os.path.join(out_dir, datatype, method))

                    # skip analysis if summarized taxa file does exist
                    if not os.path.exists(os.path.join(
                        summarize_taxa_dir, datatype, method, "%s_%s" % (tool, study), otu_table)):
                        print "skipping %s does not exist" % os.path.join(
                            summarize_taxa_dir, datatype, method, "%s_%s" % (tool, study), otu_table)
                        continue

                    with open(os.path.join(
                        summarize_taxa_dir, datatype, method, "%s_%s" % (tool, study), otu_table), 'U') as actual:
                        for line in actual:
                            if line.startswith('#'):
                                continue
                            else:
                                # remove ";Other" strings appended by summarize_taxa.py to extend
                                # taxonomy up to specified taxonomy level
                                if ";Other" in line:
                                    line = line.replace(";Other","")
                                tax = line.strip().split("\t")[0]
                                actual_tax.add(tax)

                    if out_dir != "None":
                        # output list of false-positive taxa
                        with open(os.path.join(out_dir, datatype, method, "%s_%s_fp.txt" % (tool, study)), 'w') as of:
                            lis = actual_tax - expected_tax
                            for l in lis:
                                of.write("%s\n" % l)
                            
                        # output list of false-negative taxa
                        with open(os.path.join(out_dir, datatype, method, "%s_%s_fn.txt" % (tool, study)), 'w') as fn_o:
                            lis = expected_tax - actual_tax
                            for l in lis:
                                fn_o.write("%s\n" % l)
                    
                        # output list of true-positive taxa
                        with open(os.path.join(out_dir, datatype, method, "%s_%s_tp.txt" % (tool, study)), 'w') as tp_o:
                            lis = actual_tax & expected_tax
                            for l in lis:
                                tp_o.write("%s\n" % l)

                    fp_chimera = 0
                    fp_known = 0
                    fp_other = 0

                    if (study is not "even" and study is not "staggered"):
                        fp_chimera, fp_known, fp_other = compute_fp_other(results_dir, out_dir,
                            filter_otus_dir, chimera_db_18S, chimera_db_16S, taxonomy_mean, taxonomy_stdev,
                            blast_nt_index, actual_tax, expected_tax, tool, study, datatype,
                            method, tax_level, graph_abundance)

                    tp = len(actual_tax & expected_tax)
                    fp = len(actual_tax - expected_tax)
                    fn = len(expected_tax - actual_tax)

                    p = tp / float(tp + fp)
                    r = tp / float(tp + fn)
                    f = float(2 * p * r) / float(p + r)

                    sys.stdout.write("%s\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tool, p, r, f, tp, fn, fp, fp_chimera, fp_known, fp_other))

                    if (study is not "even" and study is not "staggered"):
                        # output taxonomy_mean and taxonomy_stdev values to file
                        with open (os.path.join(out_dir, datatype, method, "%s_taxonomy_mean.txt" % study), 'w') as out_fp:
                            for tool in taxonomy_mean:
                                out_fp.write("%s\t" % tool)
                                for value in taxonomy_mean[tool]:
                                    out_fp.write("%s\t" % value)
                                out_fp.write("\n")

                        with open (os.path.join(out_dir, datatype, method, "%s_taxonomy_stdev.txt" % study), 'w') as out_fp:
                            for tool in taxonomy_stdev:
                                out_fp.write("%s\t" % tool)
                                for value in taxonomy_stdev[tool]:
                                    out_fp.write("%s\t" % value)
                                out_fp.write("\n")


if __name__ == "__main__":
    main(sys.argv[1:])

