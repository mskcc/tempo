#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration tests for running Tempo pipeline with different configurations

https://ccstempo.netlify.app/juno-setup.html#test-your-environment
https://github.com/mskcc/cas-ops/blob/master/tempo/Makefile
https://www.nextflow.io/docs/latest/config.html#environment-variables
"""
import os
import subprocess
import unittest
import getpass
import shutil
import glob
import json
from tempfile import TemporaryDirectory, mkdtemp

try:
    from .serializeDir import DirSerializer
except ModuleNotFoundError:
    from serializeDir import DirSerializer

username = getpass.getuser()
USER_SCRATCH = os.path.join('/scratch', username)

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONF_DIR = os.path.join(PROJECT_DIR, 'conf')

TEST_INPUTS_DIR = os.path.join(PROJECT_DIR, 'test_inputs')
TEST_DATA = os.path.join(THIS_DIR, "test-data") # needs to be downloaded; `make test-data`

# # config files
# in this dir;
LOCAL_CONFIG = os.path.join(THIS_DIR, "local.config")
INTEGRATION_TEST_CONFIG = os.path.join(THIS_DIR, "test.config")
NEXTFLOW_CONFIG = os.path.join(PROJECT_DIR, "nextflow.config")
# in the root conf dir
TRAVIS_TEST_CONFIG = os.path.join(CONF_DIR, "test.config")
CONTAINERS_CONFIG = os.path.join(CONF_DIR, "containers.config")
RESOURCES_CONFIG = os.path.join(CONF_DIR, "resources.config")
REFERENCES_CONFIG = os.path.join(CONF_DIR, "references.config")
EXOME_CONFIG = os.path.join(CONF_DIR, "exome.config")
GENOME_CONFIG = os.path.join(CONF_DIR, "genome.config")
SINGULARITY_CONFIG = os.path.join(CONF_DIR, "singularity.config")
JUNO_CONFIG = os.path.join(CONF_DIR, "juno.config")

defaults = {
    'PIPELINE_SCRIPT' : os.path.join(PROJECT_DIR, 'pipeline.nf'),
    'TMPDIR' : USER_SCRATCH, # tmpdir var used inside Tempo pipeline
    'NXF_SINGULARITY_CACHEDIR' : '/juno/work/taylorlab/cmopipeline/singularity_images',
    'NXF_ANSI_LOG' : 'false'
}

class Nextflow(object):
    """
    Class to run a Nextflow pipeline in a temp dir for testing
    """
    def __init__(self,
        tmpdir, # TemporaryDirectory path to run in
        args = (), # CLI args to pass to `nextflow run`
        configs = (), # extra config files to use
        defaults = defaults, # default attributes to apply to this instance
        **kwargs # pass extra keyword args to override the defaults if needed
        ):
        self.tmpdir = tmpdir
        self.args = args

        # set all the default args
        for key, value in defaults.items():
            setattr(self, key, value)

        # overwrite with any key word args that were passed
        for key, value in kwargs.items():
            if key in defaults:
                setattr(self, key, value)

        # collect the config files passed
        self.configs = []
        # self.configs.append(self.INTEGRATION_TEST_CONFIG)
        if configs:
            for config in configs:
                self.configs.append(config)

    def run(self,
        print_stdout = False,
        print_stderr = False,
        print_command = False,
        validate = False, # check that a non-zero returncode was returned; requires testcase !
        testcase = None, # the unittest.TestCase instance to use for assertions
        ):
        """
        Run the Nextflow pipeline via subprocess

        >>> proc_stdout, proc_stderr, returncode, output_dir = nxf.run(print_stdout = True, print_stderr = True, print_command = True)
        >>> proc_stdout, proc_stderr, returncode, output_dir = nxf.run(validate = True, testcase = self)

        """
        # locations for Nextflow output items to write to tmpdir
        NXF_LOG = os.path.join(self.tmpdir, "nextflow.log")
        STDOUT_LOG = os.path.join(self.tmpdir, "stdout.log")
        STDERR_LOG = os.path.join(self.tmpdir, "stderr.log")
        NXF_WORK = os.path.join(self.tmpdir, "work")
        OUTPUT_DIR = os.path.join(self.tmpdir, "output")
        NXF_REPORT = os.path.join(self.tmpdir, "nextflow.html")
        NXF_TIMELINE = os.path.join(self.tmpdir, "timeline.html")
        NXF_TRACE = os.path.join(self.tmpdir, "trace.txt")
        BamMapping = os.path.join(self.tmpdir, "bamMapping.tsv") # bamMapping.tsv # params.outname
        FileTracking = os.path.join(self.tmpdir, "fileTracking.tsv") # fileTracking.tsv # params.fileTracking

        # add the extra items to the CLI args
        args = list(self.args)
        args += [
            '-with-report', NXF_REPORT,
            '-with-timeline', NXF_TIMELINE,
            '-with-trace', NXF_TRACE,
            '--outDir', OUTPUT_DIR,
            '--outname', BamMapping,
            '--fileTracking', FileTracking,
            ]

        # build args to use extra config files
        config_args = []
        for config in self.configs:
            config_args.append('-c')
            config_args.append(config)

        # build the full list of args for the command to run
        command = [
        'nextflow',
        '-log', NXF_LOG,
        *config_args,
        'run', self.PIPELINE_SCRIPT,
        *args
        ]

        # set up a dict of extra environment variables that we need to have set
        env = dict(
            os.environ,
            NXF_ANSI_LOG = self.NXF_ANSI_LOG,
            NXF_WORK = NXF_WORK,
            TMPDIR = self.TMPDIR,
            NXF_SINGULARITY_CACHEDIR = self.NXF_SINGULARITY_CACHEDIR
            )

        if print_command:
            print(' '.join(command))

        # need to run the command and 'tee' the stdout & stderr to file in the tmpdir but also print it to console and also return it
        proc_stdout = []
        proc_stderr = []
        with subprocess.Popen(command, env = env, stdout=subprocess.PIPE, stderr = subprocess.PIPE, bufsize=1, universal_newlines=True) as proc, open(STDOUT_LOG, 'w') as stdout, open(STDERR_LOG, 'w') as stderr:
            for line in proc.stdout:
                proc_stdout.append(line)
                stdout.write(line)
                stdout.flush()
                if print_stdout:
                    print(line, end = '')
            for line in proc.stderr:
                proc_stderr.append(line)
                stderr.write(line)
                stderr.flush()
                if print_stderr:
                    print(line, end = '')
        returncode = proc.returncode

        if validate:
            if returncode != 0:
                print(''.join(proc_stdout))
                print(''.join(proc_stderr))
            testcase.assertEqual(returncode, 0)

        return(proc_stdout, proc_stderr, returncode, OUTPUT_DIR)


class TestWorkflow(unittest.TestCase):
    """
    Test cases for running Tempo pipeline
    """
    def setUp(self):
        """make a new tmpdir for each test case"""
        self.preserve = False
        self.preserve_rename = True
        self.tmpdir = mkdtemp(dir = THIS_DIR)

    def tearDown(self):
        """
        remove the tmpdir after each test case unless it was preserved
        """
        if not self.preserve:
            shutil.rmtree(self.tmpdir)
        else:
            if self.preserve_rename:
                old_path = self.tmpdir
                new_path = os.path.join(THIS_DIR, str(self._testMethodName) + '.' + os.path.basename(self.tmpdir))
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

        # NOTE: make sure that tmpdir is in a location accessible by all computer nodes if using LSF execution
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = [INTEGRATION_TEST_CONFIG, LOCAL_CONFIG])
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run(print_stdout = True, print_stderr = True, print_command = True)
        self.assertEqual(returncode, 0)

    def test_one_tool(self):
        """
        Small test case with one sample pair running only MuTect2
        Use local config to avoid LSF overhead

        takes about 3 minutes to complete
        """
        args = [
        '--mapping', os.path.join(TEST_INPUTS_DIR, 'local/tiny_test_mapping.tsv'),
        '--pairing', os.path.join(TEST_INPUTS_DIR, 'local/small_test_pairing.tsv'),
        '-profile', 'juno', '--somatic',
        '--tools', 'mutect2'
        ]

        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = [INTEGRATION_TEST_CONFIG, LOCAL_CONFIG])
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run()
        self.assertEqual(returncode, 0)
        # TODO: more assertion criteria to check output

    def test_make_bam_part(self):
        """
        Runs in about 2 minutes
        """
        self.maxDiff = None
        args = [
        '--mapping', os.path.join(THIS_DIR, 'test_make_bam_and_qc.tsv'),
        '-profile', 'juno'
        # '-profile', 'juno,test', '-without-docker', '--profile_check=false'
        ]
        configs = [
            INTEGRATION_TEST_CONFIG,
            LOCAL_CONFIG,
            ]
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
        '--bamMapping', os.path.join(THIS_DIR, 'test_somatic.tsv'),
        '--pairing', os.path.join(THIS_DIR, 'test_make_bam_and_qc_pairing.tsv'),
        "--tools", "manta,strelka2",
        "--somatic", "--germline",
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [
            INTEGRATION_TEST_CONFIG,
            LOCAL_CONFIG,
            REFERENCES_CONFIG
            ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run(validate = True, testcase = self)


    def test_sv(self):
        """
        """
        args = [
        '--bamMapping', os.path.join(THIS_DIR, 'test_somatic.tsv'),
        '--pairing', os.path.join(THIS_DIR, 'test_make_bam_and_qc_pairing.tsv'),
        "--tools", "delly,manta",
        "--somatic", "--germline", "--aggregate",
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [
            INTEGRATION_TEST_CONFIG,
            LOCAL_CONFIG,
            REFERENCES_CONFIG
            ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = configs)
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run(validate = True, testcase = self)

    def test_msisensor(self):
        """
        """
        args = [
        '--bamMapping', os.path.join(THIS_DIR, 'test_somatic.tsv'),
        '--pairing', os.path.join(THIS_DIR, 'test_make_bam_and_qc_pairing.tsv'),
        "--tools", "msisensor",
        "--somatic",
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [
            INTEGRATION_TEST_CONFIG,
            LOCAL_CONFIG,
            REFERENCES_CONFIG
            ]
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
        "--mapping", os.path.join(THIS_DIR, "test_make_bam_and_qc.tsv"),
        '--pairing', os.path.join(THIS_DIR, 'test_pairing_duplicate.tsv'),
        "--somatic",
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [
            INTEGRATION_TEST_CONFIG,
            LOCAL_CONFIG,
            REFERENCES_CONFIG
            ]
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
        "--mapping", os.path.join(THIS_DIR, "duplicate_samplelane_makebamqc.tsv"),
        '-profile', 'juno,test,singularity',
        '-without-docker',
        '--profile_check=false',
        '--genome', "smallGRCh37",
        "--reference_base", "test-data/reference",
        '--genome_base', "test-data/reference"
        ]
        configs = [
            INTEGRATION_TEST_CONFIG,
            LOCAL_CONFIG,
            REFERENCES_CONFIG
            ]
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
        '--mapping', os.path.join(TEST_INPUTS_DIR, 'local/tiny_test_mapping.tsv'),
        '--pairing', os.path.join(TEST_INPUTS_DIR, 'local/small_test_pairing.tsv'),
        '-profile', 'juno', '--somatic', '--germline',
        '--tools', 'svaba'
        ]
        nxf = Nextflow(tmpdir = self.tmpdir, args = args, configs = [INTEGRATION_TEST_CONFIG, LOCAL_CONFIG])
        proc_stdout, proc_stderr, returncode, output_dir = nxf.run()
        self.assertEqual(returncode, 0)

        # path to output dir for SvABA
        svaba_output = os.path.join(output_dir, 'svaba')

        self.assertTrue(os.path.exists(svaba_output))

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
