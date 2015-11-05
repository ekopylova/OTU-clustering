# usage run_micca.sh

reads=$1              # from split_libraries_fastq.py (not quality filtered FASTQ file)
output_dir=$2
jobs_to_start=$3      # for parallel_assign_taxonomy_rdp.py
reference_fasta=$4    # for parallel_assign_taxonomy_rdp.py
reference_taxonomy=$5 # for parallel_assign_taxonomy_rdp.py
study=$6
declare -A micca_trimlen=( ["1688"]="150" ["1686"]="150" ["1685"]="200" ["449"]="200" ["even"]="150" ["staggered"]="150" ["nematodes"]="200" ["632"]="100" ["2107"]="150")

mkdir $output_dir


# Filter reads (note: barcodes and primers assumed already removed and
# sequence IDs in QIIME's format SAMPLENAME_SEQID)
micca-preproc $reads -q 20 -l ${micca_trimlen["$study"]} -o $output_dir/pre

# Replace QIIME's SAMPLENAME and SEQID delimiter '_' by OTUCLUST's '||'
reads_file=$(basename $reads)
sed -i 's/_/||/' $output_dir/pre/${reads_file}

# De novo clustering (includes de-replication, sorting by abundance, de novo
# chimera removal, clustering using a greedy approach)
reads_name="${reads_file%.*}"
# Compute number of sequences in FASTQ file, run --derep-fast in OTUCLUST for
# datasets with 200000+ sequences
num_lines=$(wc -l $output_dir/pre/$reads | awk '{print $1/4}' | bc)
if [ $num_lines > 200000 ]
then
    otuclust -f fastq -s 0.97 -c --derep-fast \
        --out-clust $output_dir/"${reads_name}_otus.txt" \
        --out-rep $output_dir/"${reads_name}_rep.fa" $output_dir/pre/${reads_file}
else
    otuclust -f fastq -s 0.97 -c \
        --out-clust $output_dir/"${reads_name}_otus.txt" \
        --out-rep $output_dir/"${reads_name}_rep.fa" $output_dir/pre/${reads_file}
fi

# Replace OTUCLUST's SAMPLENAME and SEQID delimiter '||' by QIIME's '_'  
sed -i 's/||/_/g' $output_dir/"${reads_name}_otus.txt"

# Add OTU label to list of sequences in OTU
awk '{printf $1"\t"$0"\n"}' $output_dir/"${reads_name}_otus.txt" > $output_dir/"${reads_name}_otus2.txt"
mv $output_dir/"${reads_name}_otus2.txt" $output_dir/"${reads_name}_otus.txt"

# Replace OTUCLUST's SAMPLENAME and SEQID delimiter '||' by '_' in
# the representative FASTA file
sed -i 's/_/||/' $output_dir/"${reads_name}_rep.fa"

# Assign taxonomy
echo "parallel_assign_taxonomy_rdp.py"
parallel_assign_taxonomy_rdp.py -i $output_dir/"${reads_name}_rep.fa" -o $output_dir/rdp_assigned_taxonomy \
-T --jobs_to_start $jobs_to_start --reference_seqs_fp $reference_fasta --id_to_taxonomy_fp $reference_taxonomy

# Create BIOM table 
echo "make_otu_table.py"
make_otu_table.py -i $output_dir/"${reads_name}_otus.txt" -o $output_dir/otu_table.biom \
-t $output_dir/rdp_assigned_taxonomy/otus_tax_assignments.txt

# Align sequences command 
echo "parallel_align_seqs_pynast.py"
parallel_align_seqs_pynast.py -i $output_dir/"${reads_name}_rep.fa" -o $output_dir/pynast_aligned_seqs \
-T --jobs_to_start $jobs_to_start --template_fp $template_str

# Filter alignment command 
echo "filter_alignment.py"
filter_alignment.py -o $output_dir/pynast_aligned_seqs -i $output_dir/pynast_aligned_seqs/otus_aligned.fasta 

# Build phylogenetic tree command 
echo "make_phylogeny.py"
make_phylogeny.py -i $output_dir/pynast_aligned_seqs/otus_aligned_pfiltered.fasta -o $output_dir/rep_set.tre