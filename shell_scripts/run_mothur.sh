#!/bin/bash
# This script closely follows the MiSeq SOP outlined here: http://www.mothur.org/wiki/MiSeq_SOP
# Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013): Development of a
# dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence
# data on the MiSeq Illumina sequencing platform. Applied and Environmental Microbiology.
# 79(17):5112-20.

reads=$1                      # from split_libraries_fastq.py (demultiplexed and barcode-primer removed FASTQ file)
output_dir=$2
threads=$3
jobs_to_start=$4              # for parallel_assign_taxonomy_rdp.py
reference_fasta=$5            # for parallel_assign_taxonomy_rdp.py
reference_fasta_aligned=$6    # for mothur alignment
reference_taxonomy=$7         # for parallel_assign_taxonomy_rdp.py
study=$8
trainset_tax=$9               # for mothur classify.seqs
trainset_fasta=${10}          # for mothur classify.seqs
taxon_to_remove=${11}
algorithm=${12}               # nearest, furthest, average
template_str=${13}            # PyNAST

#declare -A mothur_trimlen=( ["1688"]="150" ["1686"]="151" ["1685"]="251" ["even"]="150" ["staggered"]="150" ["632"]="100" ["2107"]="151" )
declare -A mothur_trimlen=( ["10000"]="101" ["100000"]="101" ["1000000"]="101" ["10000000"]="101" ["100000000"]="101" )

algorithm_abr=""
if [ "$algorithm" == "average" ]; then
	algorithm_abr="an"
elif [ "$algorithm" == "nearest" ]; then
	algorithm_abr="nn"
elif [ "$algorithm" == "furthest" ]; then
	algorithm_abr="fn"
else
	echo "$algorithm not supported"
	exit
fi

mkdir $output_dir

reads_file=$(basename $reads)
reads_name="${reads_file%.*}"
#
# The first function of mothur's MiSeq SOP is make.contigs() which requires
# a 3-column file with the first column representing sample IDs, second
# being the forward reads and third the reverse reads. Since we're only
# working with forward reads, here we make use of seqtk, QIIME's
# split_sequence_file_on_sample_ids.py and bash commands to create the
# reverse reads (reverse-complement of the forward reads) and the input
# file for make.contigs().
#
mothur_output=$output_dir/mothur_output
mkdir $mothur_output
touch $mothur_output/"mothur.file"
## Reverse-complement FASTQ reads
seqtk seq -r $reads > $output_dir/${reads_name}_rc.fastq
## Split FASTA files on sample IDs
## Forward (rename output files to include 'R1')
split_sequence_file_on_sample_ids.py -i $reads -o $output_dir/split --file_type fastq
for f in $output_dir/split/*.fastq
do
	mv $f ${f%.fastq}_R1.fastq
	# Create entry for sample in mothur mapping file for make.contigs()
	file=$(basename $f .fastq)
	printf "$file ${file}_R1.fastq ${file}_R2.fastq\n" >> $mothur_output/"mothur.file"
done
## Reverse (rename output files to include 'R2')
split_sequence_file_on_sample_ids.py -i $output_dir/${reads_name}_rc.fastq -o $output_dir/split_rc --file_type fastq
for f in $output_dir/split_rc/*.fastq
do
	mv $f ${f%.fastq}_R2.fastq
done

## Move preprocessed sequences to output directory
mv $output_dir/split/*.fastq $mothur_output
mv $output_dir/split_rc/*.fastq $mothur_output
cd $mothur_output

## Overlap forward and reverse reads
## Output File Names: 
##  + mothur.trim.contigs.fasta
##  + mothur.contigs.qual
##  + mothur.contigs.report
##  + mothur.scrap.contigs.fasta
##  + mothur.scrap.contigs.qual
##  + mothur.contigs.groups
mothur "#make.contigs(file=mothur.file)" > step_1_make_contigs.log

## Filter ambiguous sequences and sequences longer than ${mothur_trimlen["$study"]}
## Output File Names: 
##  + mothur.trim.contigs.good.fasta
##  + mothur.trim.contigs.bad.accnos
##  + mothur.contigs.good.groups
mothur "#screen.seqs(fasta=mothur.trim.contigs.fasta, group=mothur.contigs.groups, maxambig=0, maxlength=${mothur_trimlen["$study"]}, processors=$threads)" > step_2_screen_seqs.log

## Dereplication
## Output File Names: 
##  + mothur.trim.contigs.good.names
##  + mothur.trim.contigs.good.unique.fasta
mothur "#unique.seqs(fasta=mothur.trim.contigs.good.fasta)" > step_3_unique_seqs.log

## Calculate the frequencies of each sequence in each sample
## Output File Names: 
##  + mothur.trim.contigs.good.count_table
mothur "#count.seqs(name=mothur.trim.contigs.good.names, group=mothur.contigs.good.groups)" > step_4_count_seqs.log

## Align sequences
## Output File Names:
##  + mothur.trim.contigs.good.unique.align
##  + mothur.trim.contigs.good.unique.align.report
mothur "#align.seqs(fasta=mothur.trim.contigs.good.unique.fasta, reference=${reference_fasta_aligned}, processors=$threads)" > step_5_align_seqs.log

## Run summary.seqs and use output file mothur.trim.contigs.good.unique.summary to
## determine the start and end positions of alignments
mothur "#summary.seqs(fasta=mothur.trim.contigs.good.unique.align, count=${mothur_output}/mothur.trim.contigs.good.count_table, processors=$threads)" > step_6_summary_seqs.log

start_pos=$(awk '{if ($1=="Median:") printf $2}' step_6_summary_seqs.log)
end_pos=$(awk '{if ($1=="Median:") printf $3}' step_6_summary_seqs.log)

if [ -z "$start_pos" ]; then
	echo "Start position of alignment could not be determined, exiting."
	exit
elif [ -z "$end_pos" ]; then
	echo "End position of alignment could not be determined, exiting."
	exit
fi

## Filter sequences that fail to overlap from start_pos to end_pos
## Output File Names: 
##  + mothur.trim.contigs.good.unique.good.summary
##  + mothur.trim.contigs.good.unique.good.align
##  + mothur.trim.contigs.good.unique.bad.accnos
##  + mothur.trim.contigs.good.good.count_table
mothur "#screen.seqs(fasta=mothur.trim.contigs.good.unique.align, count=mothur.trim.contigs.good.count_table, summary=mothur.trim.contigs.good.unique.summary, start=${start_pos}, end=${end_pos}, maxhomop=8, processors=$threads)" > step_7_screen_seqs.log
## Output File Names: 
##  + mothur.filter
##  + mothur.trim.contigs.good.unique.good.filter.fasta
mothur "#filter.seqs(fasta=mothur.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=$threads)" > step_8_filter_seqs.log

## Dereplicate aligned sequences
## Output File Names: 
## + mothur.trim.contigs.good.unique.good.filter.count_table
## + mothur.trim.contigs.good.unique.good.filter.unique.fasta
mothur "#unique.seqs(fasta=mothur.trim.contigs.good.unique.good.filter.fasta, count=mothur.trim.contigs.good.good.count_table)" > step_9_unique_seqs.log
mothur "#unique.seqs(fasta=mothur.trim.contigs.good.unique.good.filter.fasta, name=mothur.trim.contigs.good.names)" > step_9b_unique_seqs_names.log

## Pre-cluster the aligned sequences
## Output File Names: 
## + mothur.trim.contigs.good.unique.good.filter.unique.precluster.fasta
## + mothur.trim.contigs.good.unique.good.filter.unique.precluster.count_table
## + mothur.trim.contigs.good.unique.good.filter.unique.precluster.s1.map
mothur "#pre.cluster(fasta=mothur.trim.contigs.good.unique.good.filter.unique.fasta, count=mothur.trim.contigs.good.unique.good.filter.count_table, diffs=2, processors=$threads)" > step_10_precluster.log
cat *.map > single.map

## Identify chimeric sequences and remove them
## Output File Names: 
## + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
mothur "#chimera.uchime(fasta=mothur.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=mothur.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t, processors=$threads); remove.seqs(fasta=current, accnos=current)" > step_11_uchime.log

## Assign taxonomy to clean sequences
## Output File Names: 
## + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy
## + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
mothur "#classify.seqs(fasta=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, reference=${trainset_fasta}, taxonomy=${trainset_tax}, cutoff=80)" > step_12_classify_seqs.log

## Remove unexpected lineages
## Output File Names: 
##  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
##  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta
##  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table
mothur "#remove.lineage(fasta=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, taxonomy=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=${taxon_to_remove})" > step_13_remove_lineage.log

# Cluster
# Output File Names: 
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list
mothur "#cluster.split(fasta=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.15, method=$algorithm, processors=$threads)" > step_14_cluster_split.log

# Compute number of sequences in each OTU
# Output File Names: 
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.shared
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.s1.rabund
mothur "#make.shared(list=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, count=mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, label=0.03)" > step_15_make_shared.log

# Create BIOM table
# Output File Names: 
#   + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.0.03.biom
mothur "#make.biom(shared=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.shared, label=0.03)" > step_16_make_biom.log
mv mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.0.03.biom ${output_dir}/${reads_name}.biom

## Create representative FASTA file for OTUs
## Output File Names:
##   + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.0.03.rep.names
##   + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.0.03.rep.fasta
mothur "#get.oturep(fasta=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, list=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, group=mothur.contigs.good.groups, name=mothur.trim.contigs.good.names, method=abundance, label=0.03)" > step_17_get_oturep.log

## Get OTU map
## Output File Names:
##  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.otu
mothur "#get.otulist(list=mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, label=0.03)" > step_18_get_otulist.log
#python ~/OTU-clustering/python_scripts/get_otu_map_mothur_list.py --get-otulist-fp mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.otu --names-unique-seqs-fp mothur.trim.contigs.good.unique.good.filter.names --pre-cluster-map-fp mothur.trim.contigs.good.unique.good.filter.unique.precluster.s1.map --output-otu-map-fp seqs_otus.txt

## Convert FASTA alignment to FASTA file
#sed 's|[-.]||g' mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.0.03.rep.fasta | awk -F"[>|\t]" '{if (NF>1) print ">"$3; else print $1;}' > ${output_dir}/${reads_name}_rep.fa

## Assign taxonomy
#echo "parallel_assign_taxonomy_rdp.py"
#parallel_assign_taxonomy_rdp.py -i $output_dir/${reads_name}_rep.fa -o $output_dir/rdp_assigned_taxonomy \
#-T --jobs_to_start $jobs_to_start --reference_seqs_fp $reference_fasta --id_to_taxonomy_fp $reference_taxonomy --confidence 0.8
#awk -F "[\t|]" '{print $2"\t"$5"\t"$6}' $output_dir/rdp_assigned_taxonomy/${reads_name}_rep_tax_assignments.txt > $output_dir/rdp_assigned_taxonomy/${reads_name}_rep_tax_assignments_2.txt
#mv $output_dir/rdp_assigned_taxonomy/${reads_name}_rep_tax_assignments_2.txt $output_dir/rdp_assigned_taxonomy/${reads_name}_rep_tax_assignments.txt

## Add taxonomy to BIOM table
#echo "biom add-metadata"
#biom add-metadata -i ${output_dir}/${reads_name}.biom -o ${output_dir}/${reads_name}_wtax.biom --observation-metadata-fp $output_dir/rdp_assigned_taxonomy/${reads_name}_rep_tax_assignments.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy
#v ${output_dir}/${reads_name}_wtax.biom ${output_dir}/otu_table.biom
#rm ${output_dir}/${reads_name}.biom

## Align sequences command 
#echo "parallel_align_seqs_pynast.py"
#parallel_align_seqs_pynast.py -i $output_dir/${reads_name}_rep.fa -o $output_dir/pynast_aligned_seqs \
#-T --jobs_to_start $jobs_to_start --template_fp $template_str

## Filter alignment command 
#echo "filter_alignment.py"
#if [ "$study" == "nematodes" ] || [ "$study" == "2107" ]; then
#    filter_alignment.py -o $output_dir/pynast_aligned_seqs -i $output_dir/pynast_aligned_seqs/${reads_name}_rep_aligned.fasta -e 0.10 -g 0.80
#else
#    filter_alignment.py -o $output_dir/pynast_aligned_seqs -i $output_dir/pynast_aligned_seqs/${reads_name}_rep_aligned.fasta
#fi

## Build phylogenetic tree command 
#echo "make_phylogeny.py"
#make_phylogeny.py -i $output_dir/pynast_aligned_seqs/${reads_name}_rep_aligned_pfiltered.fasta -o $output_dir/rep_set.tre
#awk '{if ($1 ~ /^>/) print $0; }' $output_dir/${reads_name}_rep.fa | awk -F"[>|\t]" '{print $2"\t"$3}' > $output_dir/tree_ids.txt

