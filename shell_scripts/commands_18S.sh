#!/bin/bash

# usage   : bash commands_18.sh
# purpose : launch all OTU picking methods (via QIIME, or run_uparse.sh for UPARSE) for 18S studies
#           in a qsub environment
# author  : Jose A. Navas-Molina (josenavasmolina@gmail.com)
#           script generated using https://github.com/josenavas/QIIME-Scaling

########### USER EDIT PATHS #############

# root dir
root=/scratch/Users/evko1434

# Silva directory 
si_path=$root/reference/Silva_111_post

# Silva 97% OTUs
si_rep_set=$si_path/rep_set/97_Silva_111_rep_set.fasta

# Silva 97% OTUs taxonomy
si_tax=$si_path/taxonomy/97_Silva_111_taxa_map_RDP_6_levels.txt

# 18S database for chimera checking
chimera_fp=$si_rep_set

# 18S PyNAST template
template_fp=$root/reference/core_Silva119_alignment.fna

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
studies_path_qiime=$root/otu_clustering_datasets/QIIME_filtered/18S
studies_path_uparse=$root/otu_clustering_datasets/UPARSE_not_filtered_QIIME_label_format/18S

# Define the output directory
output_dir=$root/program_results
mkdir $output_dir

# Define the parameter folder
param_dir=$root/OTU-clustering/param_files/18S

# UPARSE script (run_uparse.sh)
u=$root/OTU-clustering/shell_scripts

# Number of threads per job
procs=40

# Number of jobs
jobs=10

########### USER END EDIT PATHS ########## 

mkdir $output_dir/18S

# list of studies to analyze
studies=(nematodes 2107)

# Run de-novo OTU picking on all the studies                      
out_denovo_dir=$output_dir/de_novo
qsub_params="-k oe -M jenya.kopylov@gmail.com -q mem8gbq -l nodes=1:ppn=$procs -l walltime=230:00:00"

mkdir $out_denovo_dir
for i in ${studies[@]}
do
    # Run with SUMACLUST
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/sumaclust_$i -a -O $jobs -p $param_dir/DN_sumaclust_params.txt" | qsub -N 18DN_SC_$i $qsub_params; sleep 2
    # Run with UCLUST                                                                                                                                
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/uclust_$i -a -O $jobs -p $param_dir/DN_uclust_params.txt" | qsub -N 18DN_UC_$i $qsub_params; sleep 2
    # Run with USEARCH61                                                                                                                                                 
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/usearch61_$i -a -O $jobs -p $param_dir/DN_usearch61_params.txt" | qsub -N 18DN_US61_$i $qsub_params; sleep 2
    # Run with Swarm
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/swarm_$i -a -O $jobs -p $param_dir/DN_swarm_params.txt" | qsub -N 18DN_SW_$i $qsub_params; sleep 2
    # Run with USEARCH52 - no chimera detection            
    echo "pick_de_novo_otus.py -i $studies_path_qiime/$i/seqs.fna -o $out_denovo_dir/usearch_$i -a -O $jobs -p $param_dir/DN_usearch_params.txt" | qsub -N 18DN_US_$i $qsub_params; sleep 2
    # Run UPARSE q16
    echo "bash $u/run_uparse.sh $studies_path_uparse/$i/seqs.fastq $out_denovo_dir/uparse_q16_$i $chimera_fp $procs $jobs $si_rep_set $si_tax $i $template_fp 16" | qsub -N 18DN_UP_$i $qsub_params; sleep 2
    # Run UPARSE_q3
    echo "bash $u/run_uparse.sh $studies_path_uparse/$i/seqs.fastq $out_denovo_dir/uparse_q3_$i $chimera_fp $procs $jobs $si_rep_set $si_tax $i $template_fp 3" | qsub -N 18DN_UP_$i $qsub_params; sleep 2
done

# Run closed-reference OTU picking on all the studies
out_closed_dir=$output_dir/closed_ref
qsub_params="-k oe -M jenya.kopylov@gmail.com -q mem8gbq -l nodes=1:ppn=$procs -l walltime=230:00:00"

mkdir $out_closed_dir
for i in ${studies[@]}
do
    # Run with SortMeRNA
    echo "pick_closed_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $si_rep_set -o $out_closed_dir/sortmerna_$i -t $si_tax -s -p $param_dir/CR_sortmerna_params.txt" | qsub -N 18CR_SMR_$i $qsub_params; sleep 2
    # Run with UCLUST
    echo "pick_closed_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $si_rep_set -o $out_closed_dir/uclust_$i -t $si_tax -s -p $param_dir/CR_uclust_params.txt -a -O $jobs" | qsub -N 18CR_UC_$i $qsub_params; sleep 2
    # Run with USEARCH52
    echo "pick_closed_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $si_rep_set -o $out_closed_dir/usearch_$i -t $si_tax -s -p $param_dir/CR_usearch_params.txt -a -O $jobs" | qsub -N 18CR_US_$i $qsub_params; sleep 2
    # Run with USEARCH61
    echo "pick_closed_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $si_rep_set -o $out_closed_dir/usearch61_$i -t $si_tax -s -p $param_dir/CR_usearch61_params.txt -a -O $jobs" | qsub -N 18CR_US61_$i $qsub_params; sleep 2
done

# Run open-reference OTU picking on all the studies
out_open_dir=$output_dir/open_ref
qsub_params="-k oe -M jenya.kopylov@gmail.com -q mem16gbq -l nodes=1:ppn=$procs -l walltime=230:00:00"

mkdir $out_open_dir
for i in ${studies[@]}
do
    # Run with SortMeRNA and SUMACLUST
    echo "pick_open_reference_otus.py -m sortmerna_sumaclust -i $studies_path_qiime/$i/seqs.fna -r $gg_rep_set -o $out_open_dir/sortmerna_sumaclust_$i -p $param_dir/OR_sortmerna_sumaclust_params.txt -a -O $jobs" | qsub -N 18OR_SMR_SC_$i $qsub_params; sleep 2
    # Run with UCLUST
    echo "pick_open_reference_otus.py -i $studies_path_qiime/$i/seqs.fna -r $si_rep_set -o $out_open_dir/uclust_$i -p $param_dir/OR_params.txt -a -O $jobs" | qsub -N 18OR_UC_$i $qsub_params; sleep 2
    # Run with USEARCH61
    echo "pick_open_reference_otus.py -m usearch61 -i $studies_path_qiime/$i/seqs.fna -r $si_rep_set -o $out_open_dir/usearch61_$i -p $param_dir/OR_params.txt -a -O $jobs" | qsub -N 18OR_US61_$i $qsub_params; sleep 2
done



