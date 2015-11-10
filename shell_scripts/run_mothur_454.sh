# This script closely follows the 454 SOP outlined here: http://www.mothur.org/wiki/454_SOP
# Schloss PD, Gevers D, Westcott SL. (2011). Reducing the effects of PCR
# amplification and sequencing artifacts on 16S rRNA-based studies. PloS ONE. 6:e27310.

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
declare -A mothur_trimlen=( ["449"]="200" ["nematodes"]="200" )

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

mothur_output=$output_dir/mothur_output
mkdir $mothur_output
# Here we split a FASTQ file based on sample IDs (this file has already
# had barcode and primer sequences removed)
convert_fastaqual_fastq.py -f $reads -o $mothur_output/fastaqual -c fastq_to_fastaqual
mv $mothur_output/fastaqual/* $mothur_output/
cd $mothur_output

# Clean sequences
# Output File Names: 
#  + seqs.trim.fasta
#  + seqs.scrap.fasta
#  + seqs.trim.qual
#  + seqs.scrap.qual
mothur "#trim.seqs(fasta=${reads_name}.fna, qfile=${reads_name}.qual, flip=T, maxambig=0, maxhomop=8, qwindowaverage=35, qwindowsize=50, processors=$threads)" > step_1_trim_seqs.log

# Create a group file (using QIIME's formatted labels)
grep '>' ${reads_name}.trim.fasta | awk -F"[>,_]" '{printf $2"_"$3"\t"$2"\n";}' | sort -k 2,2 > ${reads_name}.trim.groups

# Dereplication
# Output File Names: 
# + seqs.trim.names
# + seqs.trim.unique.fasta
mothur "#unique.seqs(fasta=${reads_name}.trim.fasta)" > step_2_unique_seqs.log

# Align sequences
# Output File Names: 
# + seqs.trim.unique.align
# + seqs.trim.unique.align.report
# + seqs.trim.unique.flip.accnos
mothur "#align.seqs(fasta=${reads_name}.trim.unique.fasta, reference=${reference_fasta_aligned}, processors=$threads)" > step_3_align_seqs.log

# Run summary.seqs and use the summary to determine the
# start and end positions of alignments
mothur "#summary.seqs(fasta=${reads_name}.trim.unique.align, name=${reads_name}.trim.names, processors=$threads)" > step_4_summary_seqs.log

end_pos=$(awk '{if ($1=="Median:") printf $3}' step_4_summary_seqs.log)

elif [ -z "$end_pos" ]; then
	echo "End position of alignment could not be determined, exiting."
	exit
fi

# Filter sequences that fail to overlap from start_pos to end_pos
# Output File Names: 
# + seqs.trim.unique.good.align
# + seqs.trim.unique.bad.accnos
# + seqs.trim.good.names
# + seqs.trim.good.groups
mothur "#screen.seqs(fasta=${reads_name}.trim.unique.align, name=${reads_name}.trim.names, group=${reads_name}.trim.groups, end=${end_pos}, optimize=start, criteria=95, processors=$threads)" > step_5_screen_seqs.log

# Output File Names: 
# + seqs.filter
# + seqs.trim.unique.good.filter.fasta
mothur "#filter.seqs(fasta=${reads_name}.trim.unique.good.align, vertical=T, trump=., processors=$threads)" > step_6_filter_seqs.log

# Dereplicate aligned sequences
# Output File Names: 
#  + seqs.trim.unique.good.filter.names
#  + seqs.trim.unique.good.filter.unique.fasta
mothur "#unique.seqs(fasta=${reads_name}.trim.unique.good.filter.fasta, name=${reads_name}.trim.good.names)" > step_7_unique_seqs.log

# Pre-cluster the aligned sequences
# Output File Names: 
# + seqs.trim.unique.good.filter.unique.precluster.fasta
# + seqs.trim.unique.good.filter.unique.precluster.names
# + seqs.trim.unique.good.filter.unique.precluster.[SAMPLE_ID].map
mothur "#pre.cluster(fasta=${reads_name}.trim.unique.good.filter.unique.fasta, name=${reads_name}.trim.unique.good.filter.names, group=${reads_name}.trim.good.groups, diffs=2)" > step_8_precluster.log
rm *.map

# Identify chimeric sequences
# Output File Names: 
#  + seqs.trim.unique.good.filter.unique.precluster.uchime.chimeras
#  + seqs.trim.unique.good.filter.unique.precluster.uchime.accnos
mothur "#chimera.uchime(fasta=${reads_name}.trim.unique.good.filter.unique.precluster.fasta, name=${reads_name}.trim.unique.good.filter.unique.precluster.names, group=${reads_name}.trim.good.groups)" > step_9_uchime.log

# Remove chimeric sequences
# Output File Names: 
#  + seqs.trim.unique.good.filter.unique.precluster.pick.names
#  + seqs.trim.unique.good.filter.unique.precluster.pick.fasta
#  + seqs.trim.good.pick.groups
mothur "#remove.seqs(fasta=${reads_name}.trim.unique.good.filter.unique.precluster.fasta, accnos=${reads_name}.trim.unique.good.filter.unique.precluster.uchime.accnos, name=${reads_name}.trim.unique.good.filter.unique.precluster.names, group=${reads_name}.trim.good.groups, dups=true)" > step_10_remove_seqs.log

# Assign taxonomy to clean sequences
# (this commands runs continuously with group option)
# Output File Names: 
# + seqs.trim.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy
# + seqs.trim.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
mothur "#classify.seqs(fasta=${reads_name}.trim.unique.good.filter.unique.precluster.pick.fasta, reference=${trainset_fasta}, taxonomy=${trainset_tax}, cutoff=80, processors=$threads)" > step_11_classify_seqs.log

# Remove sequences with unexpected lineages
# Output File Names: 
#  + seqs.trim.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
#  + seqs.trim.unique.good.filter.unique.precluster.pick.pick.names
#  + seqs.trim.unique.good.filter.unique.precluster.pick.pick.fasta
#  + seqs.trim.good.pick.pick.groups
mothur "#remove.lineage(fasta=${reads_name}.trim.unique.good.filter.unique.precluster.pick.fasta, name=${reads_name}.trim.unique.good.filter.unique.precluster.pick.names, group=${reads_name}.trim.good.pick.groups, taxonomy=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon='${taxon_to_remove}')" > step_11_remove_lineage.log

# Cluster
# Output File Names: 
# + seqs.trim.unique.good.filter.unique.precluster.pick.pick.an.sabund
# + seqs.trim.unique.good.filter.unique.precluster.pick.pick.an.rabund
# + seqs.trim.unique.good.filter.unique.precluster.pick.pick.an.list
mothur "#cluster.split(fasta=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.fasta, taxonomy=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, name=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.names, taxlevel=3, method=$algorithm, processors=$threads)" > step_12_cluster_split.log

# Compute number of sequences in each OTU
# Output File Names: 
# seqs.trim.unique.good.filter.unique.precluster.pick.pick.an.shared
# seqs.trim.unique.good.filter.unique.precluster.pick.pick.an.[SAMPLE_ID].rabund
mothur "#make.shared(list=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.list, group=${reads_name}.trim.good.pick.pick.groups, label=0.03)" > step_13_make_shared.log

# Create BIOM table
# Output File Names: 
#   + seqs.trim.unique.good.filter.unique.precluster.pick.pick.[SAMPLE_ID].0.03.biom
mothur "#make.biom(shared=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.shared)" > step_14_make_biom.log
mv ${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.0.03.biom ${output_dir}/${reads_name}.biom

# Create representative FASTA file for OTUs
# Output File Names:
#  + seqs.trim.unique.good.filter.unique.precluster.pick.pick.an.0.03.rep.names
#  + seqs.trim.unique.good.filter.unique.precluster.pick.pick.an.0.03.rep.fasta
mothur "#get.oturep(fasta=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.fasta, list=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.list, group=${reads_name}.trim.good.pick.pick.groups, name=${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.names, method=abundance, label=0.03)" > step_15_get_oturep.log

# Convert FASTA alignment to FASTA file
sed 's|-||g' ${reads_name}.trim.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.0.03.rep.fasta > ${output_dir}/${reads_name}_rep.fa

# Assign taxonomy
echo "parallel_assign_taxonomy_rdp.py"
parallel_assign_taxonomy_rdp.py -i $output_dir/${reads_name}_rep.fa -o $output_dir/rdp_assigned_taxonomy \
-T --jobs_to_start $jobs_to_start --reference_seqs_fp $reference_fasta --id_to_taxonomy_fp $reference_taxonomy

# Add taxonomy to BIOM table
echo "biom add-metadata"
biom add-metadata -i ${output_dir}/${reads_name}.biom -o ${output_dir}/${reads_name}_wtax.biom --observation-metadata-fp $output_dir/rdp_assigned_taxonomy/${reads_name}_rep_tax_assignments.txt --observation-header OTUID,taxonomy,confidence
mv ${output_dir}/${reads_name}_wtax.biom ${output_dir}/${reads_name}.biom

# Align sequences command 
echo "parallel_align_seqs_pynast.py"
parallel_align_seqs_pynast.py -i $output_dir/${reads_name}_rep.fa -o $output_dir/pynast_aligned_seqs \
-T --jobs_to_start $jobs_to_start --template_fp $template_str

# Filter alignment command 
echo "filter_alignment.py"
filter_alignment.py -o $output_dir/pynast_aligned_seqs -i $output_dir/pynast_aligned_seqs/${reads_name}_rep_aligned.fasta 

# Build phylogenetic tree command 
echo "make_phylogeny.py"
make_phylogeny.py -i $output_dir/pynast_aligned_seqs/${reads_name}_rep_aligned_pfiltered.fasta -o $output_dir/rep_set.tre