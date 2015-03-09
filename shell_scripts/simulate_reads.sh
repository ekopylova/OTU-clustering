#!/bin/bash

# usage       : simulate_reads.sh
# purpose     : this scripts uses PrimerProspector, QIIME and ART to create an even
#               and staggered simulated community of reads, their OTU maps and BIOM tables
# dependencies: PrimerProspector (1.0.1), python_scripts (available on github), ART
#               (VanillaIceCream-03-11-2014)
# date        : 24 Nov 2014
# author      : Evguenia Kopylova (jenya.kopylov@gmail.com)

# directory to store all output
output_dir=$1

# Greengenes 97% OTU aligned database
gg_db=$2

# Greengenes 97% OTU taxonomy
gg_taxonomy=$3

# root directory "OTU-clustering" containing "python_scripts, shell_scripts and param_files"
otu_clustering=$4

python_scripts_dir=$otu_clustering/python_scripts
mapping_files=$otu_clustering/mapping_files

# Slice out V4 region from 16S rRNA (Greengenes 97% OTUs database)
# using primers 515F/806R (Caporaso JG et al., ISME J., 2012)
slice_aligned_region.py -p "515f" -x "GTGCCAGCMGCCGCGGTAA" -r "806r" \
    -y "GGACTACHVGGGTWTCTAAT" -o $output_dir/pp/ -f $gg_db

# Subsample 0.011% sequences from resulting V4 region database
# (want ~1000 species) (QIIME's script)
subsample_fasta.py -i $output_dir/pp/97_otus_alignment_region.fasta -p 0.011 \
    -o $output_dir/subsampled_species.fna

# Simulate reads even abundance using ART simulator
art_illumina -amp -sam -na -i $output_dir/subsampled_species.fna \
    -l 150 -f 100 -o $output_dir/art_even

# Quality filter using QIIME's split_library_fastq
split_libraries_fastq.py -i $output_dir/art_even.fq \
    -m $mapping_files/even_uneven_mapping_file.txt -o $output_dir/split_libraries_even \
    --barcode_type="not-barcoded" --sample_id s1

# Split libraries for UPARSE (even)
split_libraries_fastq.py -i $output_dir/art_even.fq \
    -m $mapping_files/even_uneven_mapping_file.txt --store_demultiplexed_fastq -r 1000 \
    -n 1000 -p 0.0001 --barcode_type="not-barcoded" -o $output_dir/split_libraries_even_for_uparse \
    --sample_ids s1

# Generate ground-truth OTU-map from even FASTA file
# (split libraries ran on art_even.fq for UPARSE was to put
# the file into QIIME's format, not for filtering)
python $python_scripts_dir/art_fasta_qiime_labels_to_otumap.py \
    $output_dir/split_libraries_even_for_uparse/seqs.fna $output_dir/art_even_otu_map.txt

# Generate ground-truth BIOM table (using QIIME's script)
make_otu_table.py -i $output_dir/art_even_otu_map.txt -t $gg_taxonomy \
    -o $output_dir/art_even_table.biom 

# Summarize taxonomy for even BIOM table
summarize_taxa.py -i $output_dir/art_even_table.biom -o $output_dir/summarize_taxa_even

# Create an uneven distribution FASTA file of subsampled
# species to be amplified with ART
python $python_scripts_dir/generate_abundance_file.py $output_dir/art_even_otu_map.txt \
    1 100 107600 $output_dir/subsampled_species.fna $output_dir/abundance_file.txt \
    $output_dir/subsampled_species_staggered.fasta

# Simulate uneven reads with ART
art_illumina -amp -sam -na -i $output_dir/subsampled_species_staggered.fasta \
    -l 150 -f 1 -o $output_dir/art_staggered

# Quality filter staggered reads using QIIME's split_library_fastq
split_libraries_fastq.py -i $output_dir/art_staggered.fq \
    -m $mapping_files/even_uneven_mapping_file.txt \
    -o $output_dir/split_libraries_staggered --barcode_type="not-barcoded" \
    --sample_id s1

# Split libraries for UPARSE (staggered)
split_libraries_fastq.py -i $output_dir/art_staggered.fq \
    -m $mapping_files/even_uneven_mapping_file.txt --store_demultiplexed_fastq \
    -r 1000 -n 1000 -p 0.0001 --barcode_type="not-barcoded" \
    -o $output_dir/split_libraries_staggered_for_uparse --sample_ids s1

# Generate ground-truth OTU-map from staggered FASTA file
python $python_scripts_dir/art_fasta_qiime_labels_to_otumap.py \
    $output_dir/split_libraries_staggered_for_uparse/seqs.fna \
    $output_dir/art_staggered_otu_map.txt

# Generate ground-truth BIOM table (using QIIME's script)
make_otu_table.py -i $output_dir/art_staggered_otu_map.txt -t $gg_taxonomy \
    -o $output_dir/art_staggered_table.biom

# Summarize taxonomy for staggered BIOM table
summarize_taxa.py -i $output_dir/art_staggered_table.biom \
    -o $output_dir/summarize_taxa_staggered