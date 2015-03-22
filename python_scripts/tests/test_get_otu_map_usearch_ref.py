#!/usr/bin/env python
"""
Unit tests for get_otu_map_usearch_ref.py
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

from get_otu_map_usearch_ref import clusters_from_uc_file

# Test class and cases
class ParseUsearchUcFile(TestCase):
    """ Tests for get_otu_map_usearch_ref.py functionality """

    def setUp(self):
        """ Create temporary seqs file
        """
        self.root_dir = mkdtemp()

        # create temporary USEARCH61 .uc file
        f, self.usearch_file = mkstemp(prefix='usearch_clusters_',
                                       suffix='.uc')
        close(f)

        # write read sequences to tmp file
        with open(self.usearch_file, 'w') as tmp:
            tmp.write(usearch_uc)

        self.files_to_remove = [self.usearch_file]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.root_dir)

    def test_clusters_from_uc_file(self):
		""" Test the functionality of clusters_from_uc_file() method
		"""
		with open(self.usearch_file, "U") as usearch_f:
			clusters, unassigned_seqs = clusters_from_uc_file(
				uc_lines=usearch_f,
				otu_id_field=9)

		expected_clusters = {"299228": ["s1_32798;size=150;"],
							 "4364652": ["s1_67100;size=201;", "s1_83977;size=1;", "s1_83957;size=1;"],
							 "4348382": ["s1_3756;size=153;"],
							 "4383559": ["s1_95719;size=151;"],
							 "4372973": ["s1_47350;size=219;", "s1_47396;size=1;", "s1_47400;size=1;", "s1_47605;size=1;", "s1_47427;size=1;"],
							 "4404830": ["s1_33608;size=157;"]}

		expected_failures = ["s1_104394;size=135;", "s1_104394;size=135;", "s1_102755;size=129;"]

		self.assertEqual(len(expected_clusters), len(clusters))

		for key in clusters:
			self.assertTrue(key in expected_clusters)
			clusters[key].sort()
			expected_clusters[key].sort()
			self.assertEqual(clusters[key], expected_clusters[key])

		expected_failures.sort()
		unassigned_seqs.sort()

		self.assertEqual(expected_failures, unassigned_seqs)


usearch_uc = """H	27514	150	100.0	+	0	0	500I150M726I	s1_32798;size=150;	299228
H	85466	150	100.0	+	0	0	533I150M827I	s1_67100;size=201;	4364652
H	85466	150	98.7	+	0	0	533I150M827I	s1_83977;size=1;	4364652
H	85466	150	99.3	+	0	0	533I150M827I	s1_83957;size=1;	4364652
H	83967	150	99.3	+	0	0	452I150M831I	s1_3756;size=153;	4348382
H	87267	150	100.0	+	0	0	507I150M825I	s1_95719;size=151;	4383559
H	86256	150	100.0	+	0	0	485I150M806I	s1_47350;size=219;	4372973
H	86256	150	99.3	+	0	0	485I150M806I	s1_47396;size=1;	4372973
H	86256	150	99.3	+	0	0	485I150M806I	s1_47400;size=1;	4372973
H	86256	150	99.3	+	0	0	485I150M806I	s1_47605;size=1;	4372973
H	86256	150	98.0	+	0	0	485I150M806I	s1_47427;size=1;	4372973
H	89284	150	99.3	+	0	0	437I150M784I	s1_33608;size=157;	4404830
N	*	*	*	.	*	*	*	s1_104394;size=135;	*
N	*	*	*	.	*	*	*	s1_104394;size=135;	*
N	*	*	*	.	*	*	*	s1_102755;size=129;	*
"""

if __name__ == '__main__':
    main()