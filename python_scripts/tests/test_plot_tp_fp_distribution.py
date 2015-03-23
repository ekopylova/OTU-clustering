#!/usr/bin/env python
"""
Unit tests for plot_tp_fp_distribution.py
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
from tempfile import mkdtemp
from os.path import join, isfile
from os import close, makedirs, rename
from shutil import rmtree

from skbio.util import remove_files
from skbio.parse.sequences import parse_fasta

from plot_tp_fp_distribution import plot_tp_fp_distribution

# Test class and cases
class PlotTruePositiveFalsePositiveDistributions(TestCase):
    """ Tests for plot_tp_fp_distribution.py functionality """

    def setUp(self):
        """ Create temporary working directory
        """
        self.root_dir = mkdtemp()
        makedirs(join(self.root_dir, "16S", "de_novo"))

        self.taxonomy_mean_fp = join(self.root_dir, "16S", "de_novo", "test_taxonomy_mean.txt")
        self.taxonomy_stddev_fp = join(self.root_dir, "16S", "de_novo", "test_taxonomy_stdev.txt")

        with open(self.taxonomy_mean_fp, 'w') as test_mean_f:
            test_mean_f.write(taxonomy_mean)

        with open(self.taxonomy_stddev_fp, 'w') as test_std_f:
            test_std_f.write(taxonomy_stddev)

        self.files_to_remove = [self.taxonomy_mean_fp,
                                self.taxonomy_stddev_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.root_dir)

    def test_plot_tp_fp_distribution(self):
        """
        """
        tools = {"de_novo": ["swarm", "sumaclust", "uclust"]}
        studies = {"16S": ["test"]}
        datatypes = ["16S"]
        methods = ["de_novo"]
        plot_tp_fp_distribution(
			taxonomy_summary_dir=self.root_dir,
			datatypes=datatypes,
			methods=methods,
			tools=tools,
			studies=studies,
			output_dir=self.root_dir)

        self.assertTrue(isfile(join(self.root_dir, "tp_plot_de_novo_test.pdf")))


taxonomy_mean = """sumaclust	8250.0	33.0	88.0	0.0	
swarm	7758.0	138.0	28.0	15.0	
uclust	8167.0	177.0	90.0	0.0
"""

taxonomy_stddev = """sumaclust	14516.0	49.0	137.0	0.0	
swarm	13727.0	263.0	54.0	12.0	
uclust	14409.0	313.0	136.0	0.0
"""

if __name__ == '__main__':
    main()