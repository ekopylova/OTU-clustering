#!/bin/bash

# usage   : bash commands_16.sh
# purpose : launch all OTU picking methods (via QIIME, or run_uparse.sh for UPARSE) for 16S studies
#           in a qsub environment
# author  : Jose A. Navas-Molina (josenavasmolina@gmail.com)
#           script generated using https://github.com/josenavas/QIIME-Scaling

########### USER EDIT PATHS #############

# root dir
root=/scratch/Users/evko1434

# Greengenes 13.8 directory
gg_path=$root/reference/scratch/Users/evko1434/reference/gg_13_8_otus

# Greengenes 97% OTUs
gg_rep_set=$gg_path/rep_set/97_otus.fasta

# Greengenes 97% OTUs taxonomy
gg_tax=$gg_path/taxonomy/97_otu_taxonomy.txt

# 16S gold database for chimera checking 
gold_fp=$root/reference/gold.fa

# 16S PyNAST template
template_fp=$root/reference/core_set_aligned.fasta.imputed

# Define the studies path
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
studies_path_qiime=$root/supplemental_otu_clustering_datasets/QIIME_filtered/16S
studies_path_uparse=$root/supplemental_otu_clustering_datasets/UPARSE_not_filtered_QIIME_label_format/16S

# Define the output directory
output_dir=$root/program_results
mkdir $output_dir

# Define the parameter folder
param_dir=$root/OTU-clustering/param_files/16S

# UPARSE script (run_uparse.sh)
u=$root/OTU-clustering/shell_scripts

# Number of threads per job
procs=40

# Number of jobs 
jobs=10

########### USER END EDIT PATHS ##########

mkdir $output_dir/16S

# list of studies to analyze
studies=(even staggered 1685 1686 1688 449 632)

# Run de-novo OTU picking on all the studies
out_denovo_dir=$output_dir/de_novo
qsub_params="-k oe -M jenya.kopylov@gmail.com -q mem8gbq -l nodes=1:ppn=$procs -l walltime=230:00:00"

mkdir $out_denovo_dir
for i in ${studies[@]}
do
    # Run with SUMACLUST
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/sumaclust_$i -a -O $jobs -p $param_dir/DN_sumaclust_params.txt" | qsub -N 16DN_SC_$i $qsub_params; sleep 2
    # Run with UCLUST                                                                                                                                
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/uclust_$i -a -O $jobs -p $param_dir/DN_uclust_params.txt" | qsub -N 16DN_UC_$i $qsub_params; sleep 2
    # Run with USEARCH61                                                                                                                                                 
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/usearch61_$i -a -O $jobs -p $param_dir/DN_usearch61_params.txt" | qsub -N 16DN_US61_$i $qsub_params; sleep 2
    # Run with Swarm                                                                                                                                                                                                                      
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/swarm_$i -a -O $jobs -p $param_dir/DN_swarm_params.txt" | qsub -N 16DN_SW_$i $qsub_params; sleep 2
    # Run with USEARCH52 - no chimera detection
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/usearch_$i -a -O $jobs -p $param_dir/DN_usearch_params.txt" | qsub -N 16DN_US_$i $qsub_params; sleep 2
    # Run UPARSE q16  
    echo "bash $u/run_uparse.sh $studies_path_uparse/$i/seqs.fastq $out_denovo_dir/uparse_q16_${i} $gold_fp $procs $jobs $gg_rep_set $gg_tax $i $template_fp 16" | qsub -N 16DN_UP16_$i $qsub_params; sleep 2
    # Run UPARSE q3
    echo "bash $u/run_uparse.sh $studies_path_uparse/$i/seqs.fastq $out_denovo_dir/uparse_q3_${i} $gold_fp $procs $jobs $gg_rep_set $gg_tax $i $template_fp 3" | qsub -N 16DN_UP3_$i $qsub_params; sleep 2
done

# Run closed-reference OTU picking on all the studies
out_closed_dir=$output_dir/closed_ref
qsub_params="-k oe -M jenya.kopylov@gmail.com -q mem8gbq -l nodes=1:ppn=$procs -l walltime=230:00:00"

mkdir $out_closed_dir
for i in ${studies[@]}
do
    # Run with SortMeRNA
    echo "pick_closed_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_closed_dir/sortmerna_$i -t $gg_tax -s -p $param_dir/CR_sortmerna_params.txt" | qsub -N 16CR_SMR_$i $qsub_params; sleep 2
    # Run with UCLUST
    echo "pick_closed_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_closed_dir/uclust_$i -t $gg_tax -s -a -O $jobs -p $param_dir/CR_uclust_params.txt" | qsub -N 16CR_UC_$i $qsub_params; sleep 2
    # Run with USEARCH52
    echo "pick_closed_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_closed_dir/usearch_$i -t $gg_tax -p $param_dir/CR_usearch_params.txt -s -a -O $jobs" | qsub -N 16CR_US_$i $qsub_params; sleep 2
    # Run with USEARCH61
    echo "pick_closed_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_closed_dir/usearch61_$i -t $gg_tax -p $param_dir/CR_usearch61_params_449.txt -s" | qsub -N 16CR_US61_$i $qsub_params; sleep 2
done

# Run open-reference OTU picking on all the studies
out_open_dir=$output_dir/open_ref
qsub_params="-k oe -M jenya.kopylov@gmail.com -q mem16gbq -l nodes=1:ppn=$procs -l walltime=230:00:00"

mkdir $out_open_dir
for i in ${studies[@]}
do
    # Run with SortMeRNA and SUMACLUST
    echo "pick_open_reference_otus.py -m sortmerna_sumaclust -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_open_dir/sortmerna_sumaclust_$i -p $param_dir/OR_sortmerna_sumaclust_params.txt -a -O $jobs" | qsub -N 16OR_SMR_SC_$i $qsub_params; sleep 2
    # Run with UCLUST
    echo "pick_open_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_open_dir/uclust_$i -p $param_dir/OR_params.txt -a -O $jobs" | qsub -N 16OR_UC_$i $qsub_params; sleep 2
    # Run with USEARCH61 (even and staggered datasets)
    echo "pick_open_reference_otus.py -m usearch61 -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_open_dir/usearch61_$i -p $param_dir/OR_params.txt -a -O $jobs --percent_subsample=0.1" | qsub -N 16OR_US61_$i $qsub_params; sleep 2
    # Run with USEARCH61
    echo "pick_open_reference_otus.py -m usearch61 -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_open_dir/usearch61_$i -p $param_dir/OR_params.txt -a -O $jobs" | qsub -N 16OR_US61_$i $qsub_params; sleep 2
done



