#!/usr/bin/env python
"""
Unit tests for art_fasta_qiime_labels_to_otumap.py
==================================================
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

from art_fasta_qiime_labels_to_otumap import parse_fasta_to_otumap

# Test class and cases
class ParseQiimeLabels(TestCase):
    """ Tests for art_fasta_qiime_labels_to_otumap.py functionality """

    def setUp(self):
        """ Create temporary seqs file
        """
        self.root_dir = mkdtemp()

        # create temporary FASTA file with QIIME label format
        f, self.seqs = mkstemp(prefix='seqs_',
                               suffix='.txt')
        close(f)

        # write read sequences to tmp file
        with open(self.seqs, 'w') as tmp:
            tmp.write(seqs)

        self.files_to_remove = [self.seqs]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.root_dir)

    def test_parse_fasta_to_otumap(self):
        """ Test functionality of parse_fasta_to_otumap() method,
            given an QIIME formatted & ART simulated FASTA file,
            this function should output a ground-truth OTU map
        """
        output_otumap_fp = join(self.root_dir, "otu_map.txt")
        self.files_to_remove.append(output_otumap_fp)

        parse_fasta_to_otumap(self.seqs, output_otumap_fp)

        expected_otu_map = {"4483458": ["s1_107021", "s1_107022", "s1_107023"],
                            "813079": ["s1_0", "s1_1", "s1_2"],	
                            "630311": ["s1_246", "s1_247", "s1_248"]}

        actual_otu_map = {}

    	with open(output_otumap_fp, 'U') as output_fp:
    		for line in output_fp:
    			entry = line.strip().split('\t')[0]
    			if entry not in actual_otu_map:
    				actual_otu_map[entry] = line.strip().split('\t')[1:]
    			else:
    				raise ValueError("Cannot have duplicate OTU ids: %s" % entry)

    	self.assertEqual(len(expected_otu_map), len(actual_otu_map))

    	for key in actual_otu_map:
    		self.assertTrue(key in expected_otu_map)
    		actual_otu_map[key].sort()
    		expected_otu_map[key].sort()
    		self.assertEqual(actual_otu_map[key], expected_otu_map[key])



seqs = """>s1_0 813079-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGTCGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_1 813079-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
TACCCGCAGCTCACGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGTCGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_2 813079-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGTCGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_246 630311-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGACAGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_247 630311-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGACAGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_248 630311-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTCCCGGGAAATCTGACAGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_107021 4483458-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATGTGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTGTCGAAGAGGAAGGTGGAATTTCCGG
>s1_107022 4483458-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATGTGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s1_107023 4483458-1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATGTGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
"""


if __name__ == '__main__':
    main()
