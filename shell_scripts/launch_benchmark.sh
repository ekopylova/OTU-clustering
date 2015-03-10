#!/bin/bash

######### USER DEFINED VARIABLES #########
# root output directory
output_dir=/scratch/Users/evko1434/working

# filepath to Greengenes 97% representative set
gg_reference=/scratch/Users/evko1434/reference/gg_13_8_otus/rep_set/97_otus.fasta

# filepath to SILVA reference 97% representative set
silva_reference=/scratch/Users/evko1434/reference/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta

# filepath to Greengenes 97% taxonomy
gg_taxonomy=/scratch/Users/evko1434/reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt

# filepath to SILVA 97% taxonomy (RDP compatible)
silva_taxonomy=/scratch/Users/evko1434/reference/Silva_111_post/taxonomy/97_Silva_111_taxa_map_RDP_6_levels.txt

# filepath to 16S chimera database
gold_fp=/scratch/Users/evko1434/reference/gold.fa

# filepath to PyNast template 16S core set
template_fp_bac=/scratch/Users/evko1434/reference/core_set_aligned.fasta.imputed

# filepath to PyNast template 18S core set
template_fp_euk=/scratch/Users/evko1434/reference/core_Silva119_alignment.fna

# directory path to OTU-clustering directory
otu_clustering=/scratch/Users/evko1434/OTU-clustering

# directory path to otu_clustering_datasets directory
datasets=/scratch/Users/evko1434/otu_clustering_datasets

# number of threads to use (all of the $otu_clustering/param_files/* files should
# also have threads set to use this number)
num_threads=10

# number of jobs to launch
num_jobs=1

# list of 16S studies to analyze
#studies="even staggered 1685 1686 1688 449 632"
studies_bac="staggered"

# list of 18S studies to analyze
#studies_euk="nematodes 2107"
studies_euk="nematodes"

# qsub params
qsub_params="-k oe -q long8gb -l nodes=1:ppn=$num_threads -l walltime=120:00:00"

###########################################


# filepath to studies
# each study must follow the directory structure:
#
# QIIME filtered reads:
# otu_clustering_datasets/QIIME_filtered/16S/
#                                           /1685/
#                                                /seqs.fna
#
# UPARSE non-filtered reads in QIIME label format (used by UPARSE only):
# otu_clustering_datasets/UPARSE_not_filtered_QIIME_label_format/16S/
#                                                                   /1685/
#                                                                        /seqs.fastq
#
studies_path_qiime=$datasets/QIIME_filtered
studies_path_uparse=$datasets/UPARSE_not_filtered_QIIME_label_format

mkdir $output_dir
mkdir $output_dir/program_results

# simulate even / staggered reads 
#mkdir $output_dir/simulated_reads
#bash $otu_clustering/shell_scripts/simulate_reads.sh $output_dir/simulated_reads $gg_reference $gg_taxonomy $otu_clustering $datasets/mapping_files

# launch software on 16S data
#bash $otu_clustering/shell_scripts/commands_16S.sh $gg_reference \
#    $gg_taxonomy \
#    $gold_fp \
#    $template_fp_bac \
#    $studies_path_qiime/16S \
#    $studies_path_uparse/16S \
#    $output_dir/program_results \
#    $otu_clustering/param_files/16S \
#    $otu_clustering/shell_scripts \
#    $otu_clustering/python_scripts \
#    $num_threads \
#    $num_jobs \
#    "${studies_bac}" \
#    "${qsub_params}"

# launch software on 18S data
bash $otu_clustering/shell_scripts/commands_18S.sh $silva_reference \
    $silva_taxonomy \
    $silva_reference \
    $template_fp_euk \
    $studies_path_qiime/18S \
    $studies_path_uparse/18S \
    $output_dir/program_results \
    $otu_clustering/param_files/18S \
    $otu_clustering/shell_scripts \
    $otu_clustering/python_scripts \
    $num_threads \
    $num_jobs \
    "${studies_euk}" \
    "${qsub_params}"

# remove singleton OTUs (OTUs comprising of only 1 read) from the final OTU tables
#mkdir $output_dir/run_filter_singleton_otus
#python $otu_clustering/python_scripts/run_filter_singleton_otus.py $output_dir/program_results $output_dir/run_filter_singleton_otus

