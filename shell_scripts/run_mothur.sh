# usage run_mothur.sh

reads=$1                      # from split_libraries_fastq.py (demultiplexed and barcode-primer removed FASTQ file)
output_dir=$2
jobs_to_start=$3              # for parallel_assign_taxonomy_rdp.py
reference_fasta=$4            # for parallel_assign_taxonomy_rdp.py
reference_fasta_aligned=$5    # for mothur alignment
reference_taxonomy=$6         # for parallel_assign_taxonomy_rdp.py
study=$7
trainset_tax=$8               # for mothur classify.seqs
trainset_fasta=$9             # for mothur classify.seqs
declare -A mothur_trimlen=( ["1688"]="150" ["1686"]="150" ["1685"]="200" ["449"]="200" ["even"]="150" ["staggered"]="150" ["nematodes"]="200" ["632"]="100" ["2107"]="150")

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
mkdir $othur_output
touch $mothur_output/"mothur.file"
# Reverse-complement FASTQ reads
seqtk seq -r $reads > $output_dir/${reads_name}_rc.fastq
# Split FASTA files on sample IDs
# Forward (rename output files to include 'R1')
split_sequence_file_on_sample_ids.py -i $reads -o $output_dir/split --file_type fastq
for f in $output_dir/split/*.fastq
do
	mv $f ${f%.fastq}_R1.fastq
	# Create entry for sample in mothur mapping file for make.contigs()
	file=$(basename $f .fastq)
	printf "$file ${file}_R1.fastq\n" >> $mothur_output/"mothur.file"
done
# Reverse (rename output files to include 'R2')
split_sequence_file_on_sample_ids.py -i $output_dir/${reads_name}_rc.fastq -o $output_dir/split_rc --file_type fastq
for f in $output_dir/split_rc/*.fastq
do
	mv $f ${f%.fastq}_R2.fastq
done

# Move preprocessed sequences to output directory
mv $output_dir/split/*.fastq $mothur_output
mv $output_dir/split_rc/*.fastq $mothur_output

# Overlap forward and reverse reads
# Output File Names: 
#  + mothur.trim.contigs.fasta
#  + mothur.contigs.qual
#  + mothur.contigs.report
#  + mothur.scrap.contigs.fasta
#  + mothur.scrap.contigs.qual
#  + mothur.contigs.groups
mothur "#make.contigs(file=${mothur_output}/mothur.file, processors=$threads)"

# Filter ambiguous sequences and sequences longer than ${mothur_trimlen["$study"]}
# Output File Names: 
#  + mothur.trim.contigs.good.fasta
#  + mothur.trim.contigs.bad.accnos
#  + mothur.contigs.good.groups
mothur "#screen.seqs(fasta=${mothur_output}/mothur.trim.contigs.fasta, group=${mothur_output}/mothur.contigs.groups, maxambig=0, maxlength=${mothur_trimlen["$study"]})"

# Dereplication
# Output File Names: 
#  + mothur.trim.contigs.good.names
#  + mothur.trim.contigs.good.unique.fasta
mothur "#unique.seqs(fasta=${mothur_output}/mothur.trim.contigs.good.fasta)"

# Calculate the frequencies of each sequence in each sample
# Output File Names: 
#  + mothur.trim.contigs.good.count_table
mothur "#count.seqs(name=${mothur_output}/mothur.trim.contigs.good.names, group=${mothur_output}/mothur.contigs.good.groups)"

# Align sequences
# Output File Names:
#  + mothur.trim.contigs.good.unique.align
#  + mothur.trim.contigs.good.unique.align.report
mothur "#align.seqs(fasta=${mothur_output}/mothur.trim.contigs.good.unique.fasta, reference=${reference_aligned})"

# Run summary.seqs and use output file mothur.trim.contigs.good.unique.summary to
# determine the start and end positions of alignments
mothur "#summary.seqs(fasta=${mothur_output}/mothur.trim.contigs.good.unique.align, count=${mothur_output}/mothur.trim.contigs.good.count_table)" > ${mothur_output}/summary.txt

start_pos=$(awk '{if ($1=="Median:") printf $2}' ${mothur_output}/summary.txt)
end_pos=$(awk '{if ($1=="Median:") printf $3}' ${mothur_output}/summary.txt)

# Filter sequences that fail to overlap from start_pos to end_pos
# Output File Names: 
#  + mothur.trim.contigs.good.unique.good.summary
#  + mothur.trim.contigs.good.unique.good.align
#  + mothur.trim.contigs.good.unique.bad.accnos
#  + mothur.trim.contigs.good.good.count_table
mothur "#screen.seqs(fasta=${mothur_output}/mothur.trim.contigs.good.unique.align, count=${mothur_output}/mothur.trim.contigs.good.count_table, summary=${mothur_output}/mothur.trim.contigs.good.unique.summary, start=${start_pos}, end=${end_pos}, maxhomop=8)"
# Output File Names: 
#  + mothur.filter
#  + mothur.trim.contigs.good.unique.good.filter.fasta
mothur "#filter.seqs(fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.align, vertical=T, trump=.)"

# Dereplicate aligned sequences
# Output File Names: 
# + mothur.trim.contigs.good.unique.good.filter.count_table
# + mothur.trim.contigs.good.unique.good.filter.unique.fasta
mothur "#unique.seqs(fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.fasta, count=${mothur_output}/mothur.trim.contigs.good.good.count_table)"

# Pre-cluster the aligned sequences
# Output File Names: 
# + mothur.trim.contigs.good.unique.good.filter.unique.precluster.fasta
# + mothur.trim.contigs.good.unique.good.filter.unique.precluster.count_table
# + mothur.trim.contigs.good.unique.good.filter.unique.precluster.s1.map
mothur "#pre.cluster(fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.fasta, count=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.count_table, diffs=2)"

# Identify chimeric sequences and remove them
# Output File Names: 
# + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
mothur "#chimera.uchime(fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t); remove.seqs(fasta=current, accnos=current)"

# Assign taxonomy to clean sequences, final quality control step to remove sequences
# hitting taxons other than the expected
# Output File Names: 
# + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy
# + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
mothur "#classify.seqs(fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, reference=${trainset_fasta}, taxonomy=${trainset_tax}, cutoff=80)"

# Remove unexpected lineages
# Output File Names: 
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table
mothur "#remove.lineage(fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, taxonomy=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

# Cluster (generate mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list)
# METHOD: AVERAGE LINKAGE
# Output File Names: 
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list
mothur "#cluster.split(fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.15, method=average)"

# Compute number of sequences in each OTU
# Output File Names: 
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared
#  + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.s1.rabund
mothur "#make.shared(list=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, label=0.03)"

# Create BIOM table
# Output File Names: 
#   + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.biom
mothur "#make.biom(shared=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared)"
mv mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.biom ${output_dir}/${reads_name}.biom

# Create representative FASTA file for OTUs
# Output File Names: 
#   + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist
mothur "#dist.seqs(fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, processors=$threads)"
# Output File Names:
#   + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.names
#   + mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.fasta
mothur "#get.oturep(column=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist, list=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, name=${mothur_output}/mothur.trim.contigs.good.names, fasta=${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, label=0.03)"

# Convert FASTA alignment to FASTA file
sed 's|-||g' ${mothur_output}/mothur.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.fasta > ${output_dir}/${reads_name}_rep.fa

# Assign taxonomy
echo "parallel_assign_taxonomy_rdp.py"
parallel_assign_taxonomy_rdp.py -i $output_dir/${reads_name}_rep.fa -o $output_dir/rdp_assigned_taxonomy \
-T --jobs_to_start $jobs_to_start --reference_seqs_fp $reference_fasta --id_to_taxonomy_fp $reference_taxonomy

# Add taxonomy to BIOM table
echo "biom add-metadata"
biom add-metadata -i ${output_dir}/${reads_name}.biom -o ${output_dir}/${reads_name}_wtax.biom --observation-metadata-fp $output_dir/rdp_assigned_taxonomy/otus_tax_assignments.txt --observation-header OTUID,taxonomy,confidence

# Align sequences command 
echo "parallel_align_seqs_pynast.py"
parallel_align_seqs_pynast.py -i $output_dir/${reads_name}_rep.fa -o $output_dir/pynast_aligned_seqs \
-T --jobs_to_start $jobs_to_start --template_fp $template_str

# Filter alignment command 
echo "filter_alignment.py"
filter_alignment.py -o $output_dir/pynast_aligned_seqs -i $output_dir/pynast_aligned_seqs/otus_aligned.fasta 

# Build phylogenetic tree command 
echo "make_phylogeny.py"
make_phylogeny.py -i $output_dir/pynast_aligned_seqs/otus_aligned_pfiltered.fasta -o $output_dir/rep_set.tre