#!/usr/bin/env python
"""
Unit tests for generate_abundance_file.py
=========================================
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkstemp, mkdtemp
from os.path import join
from os import close
from shutil import rmtree

from skbio.util import remove_files
from skbio.parse.sequences import parse_fasta

from generate_abundance_file import create_abundance, output_fasta_abundance_files

# Test class and cases
class GenerateAbundance(TestCase):
    """ Tests for generate_abundance_file.py functionality """

    def setUp(self):
        """ Create temporary seqs file
        """
        self.root_dir = mkdtemp()

        # create temporary OTU map
        f, self.otu_map = mkstemp(prefix='otu_map_',
                                  suffix='.uc')
        close(f)

        # write OTU map to tmp file
        with open(self.otu_map, 'w') as tmp:
            tmp.write(otu_map)

        # create temporary representative sequence file
        f, self.rep_set_fp = mkstemp(prefix='rep_set_',
                                     suffix='.fasta')
        close(f)

        # write OTU map to tmp file
        with open(self.rep_set_fp, 'w') as tmp:
            tmp.write(rep_set)

        # create temporary staggered representative sequence file
        f, self.rep_set_staggered_fp = mkstemp(
            prefix='rep_set_staggered_',
            suffix='.fasta')
        close(f)

        # write staggered rep set to tmp file
        with open(self.rep_set_staggered_fp, 'w') as tmp:
            tmp.write(generated_staggered_set)

        self.files_to_remove = [self.otu_map,
                                self.rep_set_fp,
                                self.rep_set_staggered_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.root_dir)
        pass

    def test_create_abundance(self):
        """ Test create_abundance() method
            Expected to create a file with random
            abundances for each species
        """
        abund = create_abundance(otu_map_fp=self.otu_map,
                                 minimum=1,
                                 maximum=5)

        expected_ids = ["1639369", "1131243", "4332615", "4353300"]
        actual_ids = []
        total_abund = 0.0
        for key in abund:
            self.assertTrue(key in expected_ids)
            actual_ids.append(key)
            total_abund += abund[key]

        self.assertEqual(len(expected_ids), len(actual_ids))
        # total abundance should be very close to 1
        self.assertTrue(total_abund <= 1.0)
        self.assertTrue(total_abund >= 0.99)

    def test_output_fasta_abundance_files(self):
        """ Test output_fasta_abundance_files() method
            Expected to use a dictionary of random abundances
            and the total number of sequences to generate
            a FASTA file with number of sequences corresponding
            to the specified abundance for each
        """
        abund = {'1639369': 0.18, '1131243': 0.45, '4332615': 0.09, '4353300': 0.27}

        self.output_fasta_fp = join(self.root_dir, "output_staggered_species.fna")
        self.output_abundance_fp = join(self.root_dir, "output_staggered_abundances.txt")

        output_fasta_abundance_files(abund=abund,
            num_reads_to_generate=10,
            rep_set_fp=self.rep_set_fp,
            output_fasta_fp=self.output_fasta_fp,
            output_abundance_fp=self.output_abundance_fp)

        # check the output abundance file is correct
        actual_abund = {}
        with open(self.output_abundance_fp, "U") as abundance_fp:
            for line in abundance_fp:
                _id = line.strip().split('\t')[0]
                _abund = float(line.strip().split('\t')[1])
                if _id not in actual_abund:
                    actual_abund[_id] = float("{0:.2f}".format(_abund))
                else:
                    raise ValueError("Duplicate reference IDs are not supported: %s" % _id)

        for key in abund:
            self.assertTrue(key in actual_abund)
            self.assertEqual(abund[key], actual_abund[key])

        expected_staggered_set = {}
        with open(self.rep_set_staggered_fp, "U") as expect_fasta_f:
            for label, seq in parse_fasta(expect_fasta_f):
                if label not in expected_staggered_set:
                    expected_staggered_set[label] = [seq]
                else:
                    expected_staggered_set[label].append(seq)

        actual_staggered_set = {}
        # check the output FASTA file is correct
        with open(self.output_fasta_fp, "U") as out_fasta_f:
            for label, seq in parse_fasta(out_fasta_f):
                if label not in actual_staggered_set:
                    actual_staggered_set[label] = [seq]
                else:
                    actual_staggered_set[label].append(seq)

        self.assertEqual(len(expected_staggered_set), len(actual_staggered_set))

        for key in actual_staggered_set:
            self.assertTrue(key in expected_staggered_set)
            actual_staggered_set[key].sort()
            expected_staggered_set[key].sort()
            self.assertEqual(expected_staggered_set[key], actual_staggered_set[key])


otu_map = """1639369	s1_66100	s1_66101	s1_66102
1131243	s1_81000	s1_81001
4332615	s1_90500	s1_90501	s1_90502	s1_90503
4353300	s1_93000	s1_93001
"""

# V4 regions of selected Greengenes sequences (gg 13.8)
rep_set = """
>1639369
TACGGAGGGTGCAAGCGTTATCCGGATTCACTGGGTTTAAAGGGAGCGCAGGTGGGTAGGTAAGTCAGTGGTGAAA
TCTCCGAGCTTAACTCGGAAACTGCCGTTGATACTATCAGTCTTGAATATTGTGGAGGTTAGCGGAATATGTCATG
TAGCGGTGAAATGCATAGATATGACATAGAACACCAATTGCGAAGGCAGCTGGCTACACATATTGACACTGAGGCT
CGAAAGCGTGGGGATCAAACAGG
>1131243
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGTGTAGGTGGTTTTATAAGTTTGATGTGAAA
GCCCTGGGCTTAACCTGGGAATGGCGTTGAAAACTGTAAGACTGGAGTCGGACAGAGGGTGGTGGAATTTCCGGTG
TAGCAGTGAAATGCGTAGAGATCGGAAAGAACACCCGTGGCGAAGGCGGCCACCTGGGTCCACTGACACTGAGGCG
CGAAAGCGTGGGGAGCAAACAGG
>4332615
GACGGAGGGCGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTCTTAAGTCCTTTGTGAAA
GCCCGGGGCTCAACTCCGGAACTGCATCGGATACTGGGAGACTTGAGACAGGTAGAGGCTAGCGGAATTCCTGGTG
TAGCGGTGGAATGCGTAGATATCAGGAAGAACACCTGTGGCGAAGGCGGCTAGCTGGGCCTTCTGACGCTCATGTG
CGAAAGCGTGGGGAGCAAACAGG
>4353300
TACGGAGGGCGCGAGCGTTGTTCGGATTTATTGGGCGTAAAGGGCGCGTAGGTGGGCGAGTAAGTCTGGTGTGAAA
TCTTCCCGCTTAACGGGAAGAGGTCACTGGATACTGCTCGTCTTGAGGCCATTAGAGGAAGGCGGAATTCCTGGTG
TAGCGGTGGAATGCGTAGATATCAGGAAGAACACCGATTGCGAAGGCGGCCTTCTGGGATGCCTGACACTGAGGCG
CGAAAGCCAGGGGAGCGAACGGG
"""

generated_staggered_set = """>1639369
TACGGAGGGTGCAAGCGTTATCCGGATTCACTGGGTTTAAAGGGAGCGCAGGTGGGTAGGTAAGTCAGTGGTGAAA
TCTCCGAGCTTAACTCGGAAACTGCCGTTGATACTATCAGTCTTGAATATTGTGGAGGTTAGCGGAATATGTCATG
TAGCGGTGAAATGCATAGATATGACATAGAACACCAATTGCGAAGGCAGCTGGCTACACATATTGACACTGAGGCT
CGAAAGCGTGGGGATCAAACAGG
>1131243
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGTGTAGGTGGTTTTATAAGTTTGATGTGAAA
GCCCTGGGCTTAACCTGGGAATGGCGTTGAAAACTGTAAGACTGGAGTCGGACAGAGGGTGGTGGAATTTCCGGTG
TAGCAGTGAAATGCGTAGAGATCGGAAAGAACACCCGTGGCGAAGGCGGCCACCTGGGTCCACTGACACTGAGGCG
CGAAAGCGTGGGGAGCAAACAGG
>1131243
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGTGTAGGTGGTTTTATAAGTTTGATGTGAAA
GCCCTGGGCTTAACCTGGGAATGGCGTTGAAAACTGTAAGACTGGAGTCGGACAGAGGGTGGTGGAATTTCCGGTG
TAGCAGTGAAATGCGTAGAGATCGGAAAGAACACCCGTGGCGAAGGCGGCCACCTGGGTCCACTGACACTGAGGCG
CGAAAGCGTGGGGAGCAAACAGG
>1131243
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGTGTAGGTGGTTTTATAAGTTTGATGTGAAA
GCCCTGGGCTTAACCTGGGAATGGCGTTGAAAACTGTAAGACTGGAGTCGGACAGAGGGTGGTGGAATTTCCGGTG
TAGCAGTGAAATGCGTAGAGATCGGAAAGAACACCCGTGGCGAAGGCGGCCACCTGGGTCCACTGACACTGAGGCG
CGAAAGCGTGGGGAGCAAACAGG
>1131243
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGTGTAGGTGGTTTTATAAGTTTGATGTGAAA
GCCCTGGGCTTAACCTGGGAATGGCGTTGAAAACTGTAAGACTGGAGTCGGACAGAGGGTGGTGGAATTTCCGGTG
TAGCAGTGAAATGCGTAGAGATCGGAAAGAACACCCGTGGCGAAGGCGGCCACCTGGGTCCACTGACACTGAGGCG
CGAAAGCGTGGGGAGCAAACAGG
>4332615
GACGGAGGGCGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTCTTAAGTCCTTTGTGAAA
GCCCGGGGCTCAACTCCGGAACTGCATCGGATACTGGGAGACTTGAGACAGGTAGAGGCTAGCGGAATTCCTGGTG
TAGCGGTGGAATGCGTAGATATCAGGAAGAACACCTGTGGCGAAGGCGGCTAGCTGGGCCTTCTGACGCTCATGTG
CGAAAGCGTGGGGAGCAAACAGG
>4353300
TACGGAGGGCGCGAGCGTTGTTCGGATTTATTGGGCGTAAAGGGCGCGTAGGTGGGCGAGTAAGTCTGGTGTGAAA
TCTTCCCGCTTAACGGGAAGAGGTCACTGGATACTGCTCGTCTTGAGGCCATTAGAGGAAGGCGGAATTCCTGGTG
TAGCGGTGGAATGCGTAGATATCAGGAAGAACACCGATTGCGAAGGCGGCCTTCTGGGATGCCTGACACTGAGGCG
CGAAAGCCAGGGGAGCGAACGGG
>4353300
TACGGAGGGCGCGAGCGTTGTTCGGATTTATTGGGCGTAAAGGGCGCGTAGGTGGGCGAGTAAGTCTGGTGTGAAA
TCTTCCCGCTTAACGGGAAGAGGTCACTGGATACTGCTCGTCTTGAGGCCATTAGAGGAAGGCGGAATTCCTGGTG
TAGCGGTGGAATGCGTAGATATCAGGAAGAACACCGATTGCGAAGGCGGCCTTCTGGGATGCCTGACACTGAGGCG
CGAAAGCCAGGGGAGCGAACGGG
"""

if __name__ == '__main__':
    main()