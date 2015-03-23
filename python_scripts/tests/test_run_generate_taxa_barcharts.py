#!/usr/bin/env python
"""
Unit tests for run_generate_taxa_barcharts.py
=============================================
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

from run_generate_taxa_barcharts import generate_taxa_barcharts


# Test class and cases
class GenerateTaxaBarcharts(TestCase):
    """ Tests for run_generate_taxa_barcharts.py functionality """

    def setUp(self):
        """ Create temporary working directory & summarized taxonomy
            files for swarm, uclust and sumaclust
        """
        self.root_dir = mkdtemp()
        makedirs(join(self.root_dir, "16S", "de_novo", "swarm_test"))
        makedirs(join(self.root_dir, "16S", "de_novo", "uclust_test"))
        makedirs(join(self.root_dir, "16S", "de_novo", "sumaclust_test"))

        self.taxonomy_swarm_fp = join(self.root_dir, "16S", "de_novo", "swarm_test", "otu_table_mc2_L6.txt")
        self.taxonomy_uclust_fp = join(self.root_dir, "16S", "de_novo", "uclust_test", "otu_table_mc2_L6.txt")
        self.taxonomy_sumaclust_fp = join(self.root_dir, "16S", "de_novo", "sumaclust_test", "otu_table_mc2_L6.txt")

        with open(self.taxonomy_swarm_fp, 'w') as test_swarm_f:
            test_swarm_f.write(swarm_L6_classification)

        with open(self.taxonomy_uclust_fp, 'w') as test_uclust_f:
            test_uclust_f.write(uclust_L6_classification)

        with open(self.taxonomy_sumaclust_fp, 'w') as test_sumaclust_f:
            test_sumaclust_f.write(sumaclust_L6_classification)

        self.files_to_remove = [self.taxonomy_swarm_fp,
                                self.taxonomy_uclust_fp,
                                self.taxonomy_sumaclust_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.root_dir)

    def test_generate_taxa_barcharts(self):
        """ Test the functionality of generate_taxa_barcharts() method
            Expected to produce a PDF layered bargraph of abundant
            taxonomy distributions
        """
        tools = {"de_novo": ["swarm", "sumaclust", "uclust"],
                 "closed_ref": [],
                 "open_ref": []}
        studies = {"16S": ["test"]}
        studies_bac_mock = ["test"]
        studies_euk_mock = []
        datatypes = ["16S"]
        methods = ["de_novo"]
        generate_taxa_barcharts(
            output_dir=self.root_dir,
            taxa_summary_dir=self.root_dir,
            datatypes=datatypes,
            studies=studies,
            studies_bac_mock=studies_bac_mock,
            studies_euk_mock=studies_euk_mock,
            methods=methods,
            tools=tools,
            top_N_taxa_mock=9,
            top_N_taxa_env=9)

        self.assertTrue(isfile(join(self.root_dir, "16S", "barchart_test_top_9.pdf")))


swarm_L6_classification = """# Constructed from biom file
#OTU ID MockMiSeq.even.673461
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Dorea\t0.073351071653
k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus\t0.0837221947382
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Ruminococcus]\t0.0932866253961
k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides\t0.275264849157
k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Verrucomicrobiaceae;g__Akkermansia\t0.0426165358577
k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Porphyromonadaceae;g__Parabacteroides\t0.0245573219594
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__[Tissierellaceae];g__Anaerococcus\t0.0252971614232
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnospira\t0.026616393238
"""

uclust_L6_classification = """# Constructed from biom file
#OTU ID MockMiSeq.even.673461
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Dorea\t0.0710505380396
k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus\t0.0829083448478
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia\t0.00429175564209
k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides\t0.272570447599
k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Verrucomicrobiaceae;g__Akkermansia\t0.0423313654064
k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Porphyromonadaceae;g__Parabacteroides\t0.0236193108069
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnospira\t0.0252857681196
"""

sumaclust_L6_classification = """# Constructed from biom file
#OTU ID MockMiSeq.even.673461
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Dorea\t0.0701917233456
k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus\t0.08314949093
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia\t0.00430399347325
k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides\t0.272683377594
k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Verrucomicrobiaceae;g__Akkermansia\t0.0423115025682
k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Porphyromonadaceae;g__Parabacteroides\t0.0236594766943
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnospira\t0.0255117756262
"""

if __name__ == '__main__':
    main()