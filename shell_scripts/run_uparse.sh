# usage run_usearch7_denovo.sh reads.fna output_dir chimera_db.fasta num_threads num_jobs_to_start reference_db.fasta reference_db_taxonomy.txt study_id pynast_template_fp q_value
# for 16S, db.fasta should be cs_gold.fa (downloaded from http://www.drive5.com/usearch/manual/uchime_ref.html)

u="usearch70"         # name of usearch executable
reads=$1              # from split_libraries_fastq.py
output_dir=$2
chimera_db_path=$3    # must be cs_gold.fa for 16S
threads=$4            # for usearch_global
jobs_to_start=$5      # for parallel_assign_taxonomy_rdp.py
reference_fasta=$6    # for parallel_assign_taxonomy_rdp.py
reference_taxonomy=$7 # for parallel_assign_taxonomy_rdp.py
study=$8
template_str=$9
q=${10}
declare -A uparse_trimlen=( ["1688"]="150" ["1686"]="150" ["1685"]="250" ["449"]="250" ["even"]="150" ["staggered"]="150" ["nematodes"]="250" ["632"]="100" ["2107"]="150" ["1919"]="250")

mkdir $output_dir

# filter the reads
fastq_ascii=""
if [ "$study" == "632" ]; then
  fastq_ascii="-fastq_ascii 64"
fi
echo $u -fastq_filter $reads -fastq_truncqual ${q} -fastq_trunclen ${uparse_trimlen["$study"]} -fastaout $output_dir/reads.fa ${fastq_ascii}
$u -fastq_filter $reads -fastq_truncqual ${q} -fastq_trunclen ${uparse_trimlen["$study"]} -fastaout $output_dir/reads.fa ${fastq_ascii}

# Dereplication
$u -derep_fulllength $output_dir/reads.fa -output $output_dir/derep.fa -sizeout

# Abundance sort and discard singletons
$u -sortbysize $output_dir/derep.fa -output $output_dir/sorted.fa -minsize 2

# OTU clustering
$u -cluster_otus $output_dir/sorted.fa -otus $output_dir/otus1.fa

# Chimera filtering using reference database
$u -uchime_ref $output_dir/otus1.fa -db $chimera_db_path -strand plus \
-nonchimeras $output_dir/otus.fa

# Map reads (including singletons) back to OTUs
$u -usearch_global $output_dir/reads.fa -db $output_dir/otus.fa -strand both \
-id 0.97 -threads $threads -uc $output_dir/map.uc

# Create OTU map
python /home/evko1434/scripts/get_otu_map_usearch_ref.py -i $output_dir/map.uc > $output_dir/seqs_otus.txt

# Assign taxonomy
echo "parallel_assign_taxonomy_rdp.py"
parallel_assign_taxonomy_rdp.py -i $output_dir/otus.fa -o $output_dir/rdp_assigned_taxonomy \
-T --jobs_to_start $jobs_to_start --reference_seqs_fp $reference_fasta --id_to_taxonomy_fp $reference_taxonomy

# Create BIOM table 
echo "make_otu_table.py"
make_otu_table.py -i $output_dir/seqs_otus.txt -o $output_dir/otu_table.biom \
-t $output_dir/rdp_assigned_taxonomy/otus_tax_assignments.txt

# Align sequences command 
echo "parallel_align_seqs_pynast.py"
parallel_align_seqs_pynast.py -i $output_dir/otus.fa -o $output_dir/pynast_aligned_seqs \
-T --jobs_to_start $jobs_to_start --template_fp $template_str

# Filter alignment command 
echo "filter_alignment.py"
filter_alignment.py -o $output_dir/pynast_aligned_seqs -i $output_dir/pynast_aligned_seqs/otus_aligned.fasta 

# Build phylogenetic tree command 
echo "make_phylogeny.py"
make_phylogeny.py -i $output_dir/pynast_aligned_seqs/otus_aligned_pfiltered.fasta -o $output_dir/rep_set.tre