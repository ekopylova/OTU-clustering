#!/bin/bash

#
#
#
# Note: swarm v 1.2.19 will not execute if there are ambiguous/non-iupac nucleotides in any of the
#       sequences. The tool prinseq-lite was used to filter out these foreign nucleotides from
#       the even & staggered datasets using the command:
#          prinseq-lite.pl -fasta seqs.fna -ns_max_n 0 -noniupac -out_good seqs_swarm -out_bad null
#       This explains the additional seqs_swarm.fna files in both even & staggered distributed
#       datasets 
#

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

# filepath to Greengees 97% tree
gg_tree=/scratch/Users/evko1434/reference/gg_13_8_otus/trees/97_otus.tree

# filepath to SILVA 97% tree
silva_tree=/scratch/Users/evko1434/reference/Silva_111_post/trees/97_Silva_111_rep_set_pfiltered.tre

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

# filepath to Blast's NT database
blast_nt=/scratch/Users/evko1434/reference/blast_databases/nt

# number of threads to use (all of the $otu_clustering/param_files/* files should
# also have threads set to use this number)
num_threads=10

# number of jobs to launch
num_jobs=1

# list of 16S studies to analyze (each study must be separated by a space)
#studies="even staggered 1685 1686 1688 449 632"
studies_bac="staggered 1688"

# subset of $studies_bac that are simulated or mock (to be passed to
# run_compute_precision_recall.py)
#simulated_mock_studies_bac="even staggered 1685 1688"
simulated_mock_studies_bac="staggered 1688"

# subset of $studies_bac that are environmental (to be passed to
# run_beta_diversity_and_procrustes.py)
env_studies_bac="632 449"

# list of 18S studies to analyze (each study must be separated by a space)
#studies_euk="nematodes 2107"
studies_euk="nematodes"

# subset of $studies_euk that are simulated or mock (to be passed to
# run_compute_precision_recall.py)
simulated_mock_studies_euk="nematodes"

# subset of $studies_euk that are environmental (to be passed to
# run_beta_diversity_and_procrustes.py)
env_studies_euk="2107"

# lists of tools for each method
# these lists will only be used in subsequent scripts after commands_16.sh
# and commands_18.sh. If the user would like to add other tools to any of
# the methods, they must also update the commands in the two previously
# mentioned files
tools_denovo="uclust usearch usearch61 swarm sumaclust uparse_q3 uparse_q16"
tools_closed_ref="sortmerna uclust usearch61 usearch"
tools_open_ref="sortmerna_sumaclust uclust usearch61"

# qsub params
qsub_params="-k oe -q long8gb -l nodes=1:ppn=$num_threads -l walltime=120:00:00"

# coordinate matrices to be used in run_beta_diversity_and_procrustes.py
coordinate_matrices="weighted unweighted"

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

#mkdir $output_dir
#mkdir $output_dir/program_results

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
#bash $otu_clustering/shell_scripts/commands_18S.sh $silva_reference \
#    $silva_taxonomy \
#    $silva_reference \
#    $template_fp_euk \
#    $studies_path_qiime/18S \
#    $studies_path_uparse/18S \
#    $output_dir/program_results \
#    $otu_clustering/param_files/18S \
#    $otu_clustering/shell_scripts \
#    $otu_clustering/python_scripts \
#    $num_threads \
#    $num_jobs \
#    "${studies_euk}" \
#    "${qsub_params}"

# remove singleton OTUs (OTUs comprising of only 1 read) from the final OTU tables
echo "Step 1: remove singleton OTUs from the final OTU tables"
mkdir $output_dir/run_filter_singleton_otus
python $otu_clustering/python_scripts/run_filter_singleton_otus.py \
    $output_dir/program_results $output_dir/run_filter_singleton_otus \
    "${studies_bac}" "${studies_euk}" "${tools_denovo}" "${tools_closed_ref}" "${tools_open_ref}"

# summarize taxa for all BIOM tables (not including singletons)
echo "Step 2: summarize taxonomies for all BIOM tables (not including singletons)"
mkdir $output_dir/run_summarize_taxa
python $otu_clustering/python_scripts/run_summarize_taxa.py \
    $output_dir/run_filter_singleton_otus $output_dir/run_summarize_taxa \
    "${studies_bac}" "${studies_euk}" "${tools_denovo}" "${tools_closed_ref}" "${tools_open_ref}"

# summarize all BIOM tables (not including singletons)
echo "Step 3: summarize all BIOM tables (not including singletons)"
mkdir $output_dir/run_summarize_tables
python $otu_clustering/python_scripts/run_summarize_tables.py \
    $output_dir/run_filter_singleton_otus $output_dir/run_summarize_tables \
    "${studies_bac}" "${studies_euk}" "${tools_denovo}" "${tools_closed_ref}" "${tools_open_ref}"

# compute true positive, false positive, false negative, precision, recall,
# F-measure and FP-chimera, FP-known, FP-other metrics using the summarized
# taxonomy results
echo "Step 4: compute TP, FP, FN, precision, recall & other metrics using the summarized taxonomies"
mkdir $output_dir/run_compute_precision_recall
python $otu_clustering/python_scripts/run_compute_precision_recall.py \
    $output_dir/program_results $output_dir/run_filter_singleton_otus $output_dir/run_summarize_taxa \
    $output_dir/run_summarize_tables $output_dir/run_compute_precision_recall $silva_reference $gold_fp \
    $blast_nt "${simulated_mock_studies_bac}" "${simulated_mock_studies_euk}" \
    "${tools_denovo}" "${tools_closed_ref}" "${tools_open_ref}" \
    $datasets/expected_taxonomies

# generate alpha diversity plots (if additional studies are added to the
# benchmark, their sampling depths should be also added to the
# run_single_rarefaction_and_plot.py)
echo "Step 5: Generate alpha diversity plots"
mkdir $output_dir/run_single_rarefaction_and_plot
python $otu_clustering/python_scripts/run_single_rarefaction_and_plot.py \
    $output_dir/run_filter_singleton_otus $output_dir/run_single_rarefaction_and_plot \
    $output_dir/program_results $gg_tree $silva_tree $datasets/mapping_files \
    "${studies_bac}" "${studies_euk}" "${tools_denovo}" "${tools_closed_ref}" "${tools_open_ref}"

# Run QIIME's beta-diversity and Procrustes analysis for all tools vs. UCLUST (default)
echo "Step 6: Run QIIME's beta-diversity and Procrustes analysis"
mkdir $output_dir/run_beta_diversity_and_procrustes
mkdir $output_dir/run_beta_diversity_and_procrustes/beta_diversity
mkdir $output_dir/run_beta_diversity_and_procrustes/procrustes
python $otu_clustering/python_scripts/run_beta_diversity_and_procrustes.py \
    $output_dir/run_filter_singleton_otus \
    $output_dir/run_beta_diversity_and_procrustes/beta_diversity \
    $output_dir/run_beta_diversity_and_procrustes/procrustes \
    $output_dir/program_results $gg_tree $silva_tree \
    $datasets/mapping_files "${env_studies_bac}" "${env_studies_euk}" \
    "${tools_denovo}" "${tools_closed_ref}" "${tools_open_ref}" \
    $coordinate_matrices

