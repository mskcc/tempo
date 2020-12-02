#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration tests for running Tempo pipeline with different configurations

Python unittest docs:
https://docs.python.org/3.7/library/unittest.html

Tempo & Nextflow docs:
https://ccstempo.netlify.app/juno-setup.html#test-your-environment
https://github.com/mskcc/cas-ops/blob/master/tempo/Makefile
https://www.nextflow.io/docs/latest/config.html#environment-variables
"""
import sys
import os
import subprocess
import unittest
import getpass
import shutil
import glob
import json
from tempfile import TemporaryDirectory, mkdtemp
from serializeDir import DirSerializer
from settings import settings, test_files
from nextflow import Nextflow

class TestWorkflow(unittest.TestCase):
    """
    Test cases for running Tempo pipeline
    """
    THIS_DIR = os.path.dirname(os.path.abspath(__file__))

    def setUp(self):
        """make a new tmpdir for each test case"""
        self.preserve = False
        self.preserve_rename = True
        self.tmpdir = mkdtemp(dir = self.THIS_DIR)

    def tearDown(self):
        """
        remove the tmpdir after each test case unless it was preserved
        """
        if not self.preserve:
            shutil.rmtree(self.tmpdir)
        else:
            if self.preserve_rename:
                old_path = self.tmpdir
                new_path = os.path.join(self.THIS_DIR, str(self._testMethodName) + '.' + os.path.basename(self.tmpdir))
                print(">>> preserving tmpdir; ", old_path, ' -> ', new_path)
                shutil.move(old_path, new_path)

    def dontrun_test_full_pipeline(self):
        """
        Full pipeline test case

        https://ccstempo.netlify.app/juno-setup.html#test-your-environment
        nextflow run pipeline.nf \
            --mapping test_inputs/local/full_test_mapping.tsv \
            --pairing test_inputs/local/full_test_pairing.tsv \
            -profile test_singularity \
            --outDir results
            --somatic --germline --QC --aggregate

        NOTE: This one takes a really long time to finish so its disabled currently
        """
        args = [
        '--mapping', os.path.join(TEST_INPUTS_DIR, 'local/small_test_mapping.tsv'),
        '--pairing', os.path.join(TEST_INPUTS_DIR, 'local/small_test_pairing.tsv'),
        '-profile', 'juno', '--somatic', '--germline', '--QC', '--aggregate'
        ]

        # NOTE: make sure that tmpdir is in a location accessible by all compute nodes if using LSF execution
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = [settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG']])
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run(print_stdout = True, print_stderr = True, print_command = True)
        self.assertEqual(returncode, 0)

    def test_one_tool(self):
        """
        Small test case with one sample pair running only MuTect2
        Use local config to avoid LSF overhead

        takes about 3 minutes to complete
        """
        args = [
        '--mapping', test_files['tiny_test_mapping'],
        '--pairing', test_files['small_test_pairing'],
        '-profile', 'juno', '--somatic',
        '--tools', 'mutect2'
        ]
        configs = [ settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG'] ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run()
        self.assertEqual(returncode, 0)
        # TODO: more assertion criteria to check output

    def test_make_bam_part(self):
        """
        Runs in about 2 minutes

        TODO: gives error about ;
        No such variable: bamsForSvABA
        need to figure out fix for that
        """
        self.maxDiff = None
        args = [
        '--mapping', test_files['test_make_bam_and_qc'],
        '-profile', 'juno'
        ]
        configs = [ settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG'] ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run(validate = True, testcase = self)

        # make sure bam dir was created
        bam_dir = os.path.join(output_dir, "bams")
        self.assertTrue(os.path.exists(bam_dir))

        # check the bam dir contents
        # the .bam and .bai do not have consistent exact sizes due to random alignment discrepancies so cannot include their size & md5 in the output data
        serializer = DirSerializer(bam_dir,
            exclude_dirs = ['fastp'],
            exclude_sizes = ['.bai', '.bam'],
            exclude_md5s = ['.bai', '.bam']
            )

        output = serializer.data
        expected_output = {
            "1234N/1234N.bam": {},
            "1234N/1234N.bam.bai": {},
            "1234T/1234T.bam": {},
            "1234T/1234T.bam.bai": {}
        }
        self.assertDictEqual(output, expected_output)

        # check num lines on bam mapping; file contains full paths which will not be consistent due to tmpdirs
        bamMapping = os.path.join(self.tmpdir, 'bamMapping.tsv')
        self.assertTrue(os.path.exists(bamMapping))
        serializer = DirSerializer(bamMapping, exclude_sizes = ['bamMapping.tsv'], exclude_md5s = ['bamMapping.tsv'])
        output = serializer.data
        expected_output = {
            bamMapping : {
            "lines": 3
            }
        }
        self.assertDictEqual(output, expected_output)

    def test_manta_strelka(self):
        """
        """
        args = [
        '--bamMapping', test_files['test_somatic'],
        '--pairing', test_files['test_make_bam_and_qc_pairing'],
        "--tools", "manta,strelka2",
        "--somatic", "--germline",
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [ settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG'], settings['REFERENCES_CONFIG'] ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run(validate = True, testcase = self)


    def test_sv(self):
        """
        """
        args = [
        '--bamMapping', test_files['test_somatic'],
        '--pairing', test_files['test_make_bam_and_qc_pairing'],
        "--tools", "delly,manta",
        "--somatic", "--germline", "--aggregate",
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [ settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG'], settings['REFERENCES_CONFIG'] ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run(validate = True, testcase = self)

    def test_msisensor(self):
        """
        """
        args = [
        '--bamMapping', test_files['test_somatic'],
        '--pairing', test_files['test_make_bam_and_qc_pairing'],
        "--tools", "msisensor",
        "--somatic",
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [ settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG'], settings['REFERENCES_CONFIG'] ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run(validate = True, testcase = self)


    def test_pairing_file_validation_pipeline(self):
        """
        Should give an error;

        ERROR: Duplicatd inputs found in tests/test_pairing_duplicate.tsv

        [NORMAL_ID:1234N, TUMOR_ID:1234T]
        [NORMAL_ID:1234N, TUMOR_ID:1234T]
        """
        args = [
        "--mapping", test_files['test_make_bam_and_qc'],
        '--pairing', test_files['test_pairing_duplicate'],
        "--somatic",
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [ settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG'], settings['REFERENCES_CONFIG'] ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run()
        self.assertEqual(returncode, 1)

    def test_nonUnique_sampleLane_validation_pipeline(self):
        """
        Should give an error

        ERROR: Duplicatd inputs found in duplicate_samplelane_makebamqc.tsv

        1234N   test-data/testdata/tiny/normal/tiny_n_L001_R1_xxx.fastq.gz
        1235N   test-data/testdata/tiny/normal/tiny_n_L001_R1_xxx.fastq.gz
        """
        args = [
        "--mapping", test_files['duplicate_samplelane_makebamqc'],
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [ settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG'], settings['REFERENCES_CONFIG'] ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run()
        self.assertEqual(returncode, 1)


    def test_SvABA(self):
        """
        Test case for running SvABA

        Takes about 2 minutes to run
        """
        self.maxDiff = None
        args = [
        '--mapping', test_files['tiny_test_mapping'],
        '--pairing', test_files['small_test_pairing'],
        '-profile', 'juno', '--somatic', '--germline',
        '--tools', 'svaba'
        ]
        configs = [ settings['INTEGRATION_TEST_CONFIG'], settings['LOCAL_CONFIG'] ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run()
        self.assertEqual(returncode, 0)

        # path to output dir for SvABA
        svaba_output = os.path.join(output_dir, 'svaba')

        self.assertTrue(os.path.exists(svaba_output), "SvABA output dir does not exist")

        # serialize the dir listing
        # dont count lines, size, md5 for some files that are subject to change with timestamps, etc..
        serializer = DirSerializer(svaba_output,
            exclude_md5s = ['.gz', '.log', '.bam'],
            exclude_sizes = ['.gz', '.log', '.bam'],
            exclude_lines = ['.log'],
            )

        output = serializer.data

        expected_output = {
            "DU874145-T__DU874145-N.svaba.somatic.sv.vcf.gz": {
                "md5x": "f50e45832f2ca18bf5c3dea3e256f5bc",
                "lines": 145
            },
            "DU874145-T__DU874145-N.svaba.unfiltered.somatic.indel.vcf.gz": {
                "md5x": "794426b6f9ad14c7e40d8f06b4fd70a6",
                "lines": 123
            },
            "DU874145-T__DU874145-N.svaba.unfiltered.germline.indel.vcf.gz": {
                "md5x": "794426b6f9ad14c7e40d8f06b4fd70a6",
                "lines": 123
            },
            "DU874145-T__DU874145-N.svaba.unfiltered.germline.sv.vcf.gz": {
                "md5x": "f50e45832f2ca18bf5c3dea3e256f5bc",
                "lines": 145
            },
            "DU874145-T__DU874145-N.svaba.somatic.indel.vcf.gz": {
                "md5x": "794426b6f9ad14c7e40d8f06b4fd70a6",
                "lines": 123
            },
            "DU874145-T__DU874145-N.alignments.txt.gz": {
                "md5x": "d41d8cd98f00b204e9800998ecf8427e",
                "lines": 0
            },
            "DU874145-T__DU874145-N.log": {},
            "DU874145-T__DU874145-N.svaba.unfiltered.somatic.sv.vcf.gz": {
                "md5x": "f50e45832f2ca18bf5c3dea3e256f5bc",
                "lines": 145
            },
            "DU874145-T__DU874145-N.discordant.txt.gz": {
                "md5x": "db889d0da915dd1a8c362d6fe311d543",
                "lines": 1
            },
            "DU874145-T__DU874145-N.svaba.germline.indel.vcf.gz": {
                "md5x": "794426b6f9ad14c7e40d8f06b4fd70a6",
                "lines": 123
            },
            "DU874145-T__DU874145-N.svaba.germline.sv.vcf.gz": {
                "md5x": "f50e45832f2ca18bf5c3dea3e256f5bc",
                "lines": 145
            },
            "DU874145-T__DU874145-N.bps.txt.gz": {
                "md5x": "1d44ec335e15760826412cfafe72a5c1",
                "lines": 1
            },
            "DU874145-T__DU874145-N.contigs.bam": {}
        }
        self.assertDictEqual(output, expected_output)

if __name__ == "__main__":
    unittest.main()
