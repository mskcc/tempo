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
from tempfile import TemporaryDirectory, mkdtemp

username = getpass.getuser()

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_INPUTS_DIR = os.path.join(PROJECT_DIR, 'test_inputs')
USER_SCRATCH = os.path.join('/scratch', username)
LOCAL_CONFIG = os.path.join(THIS_DIR, "local.config")

defaults = {
    'PIPELINE_SCRIPT' : os.path.join(PROJECT_DIR, 'pipeline.nf'),
    'TMPDIR' : USER_SCRATCH, # tmpdir var used inside Tempo pipeline
    'NXF_SINGULARITY_CACHEDIR' : '/juno/work/taylorlab/cmopipeline/singularity_images',
    'NXF_ANSI_LOG' : 'false',
    'TEST_CONFIG': os.path.join(THIS_DIR, "test.config")
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
        self.configs.append(self.TEST_CONFIG)
        if configs:
            for config in configs:
                self.configs.append(config)

    def run(self, print_stdout = False, print_stderr = False, print_command = False):
        """
        Run the Nextflow pipeline via subprocess
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

        return(proc_stdout, proc_stderr, returncode, OUTPUT_DIR)

def preserve_output(input_dir, prefix = 'output'):
    """
    Helper function to copy contents of a dir to another dir
    use during test dev to keep TemporaryDirectory contents from being deleted

    >>> preserve_output(tmpdir, prefix = self._testMethodName + '.')
    """
    output_dir = mkdtemp(prefix = prefix, dir = '.')
    print(">>> copying output to: ", output_dir)
    pattern = os.path.join(input_dir, '*')
    for old_path in glob.glob(pattern):
        new_path = os.path.join(output_dir, os.path.basename(old_path))
        shutil.move(old_path, new_path)


class TestWorkflow(unittest.TestCase):
    """
    Test cases for running Tempo pipeline
    """
    def _0_test_full_pipeline(self):
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
        with TemporaryDirectory(dir = THIS_DIR) as tmpdir:
            nxf = Nextflow(tmpdir = tmpdir, args = args)
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

        with TemporaryDirectory(dir = THIS_DIR) as tmpdir:
            nxf = Nextflow(tmpdir = tmpdir, args = args, configs = [LOCAL_CONFIG])
            proc_stdout, proc_stderr, returncode, output_dir = nxf.run()
            self.assertEqual(returncode, 0)
            # TODO: more assertion criteria to check output

    def test_SvABA(self):
        """
        Test case for running SvABA
        """
        args = [
        '--mapping', os.path.join(TEST_INPUTS_DIR, 'local/tiny_test_mapping.tsv'),
        '--pairing', os.path.join(TEST_INPUTS_DIR, 'local/small_test_pairing.tsv'),
        '-profile', 'juno', '--somatic', '--germline', '--QC', '--aggregate',
        '--tools', 'svaba'
        ]
        with TemporaryDirectory(dir = THIS_DIR) as tmpdir:
            nxf = Nextflow(tmpdir = tmpdir, args = args, configs = [LOCAL_CONFIG])
            proc_stdout, proc_stderr, returncode, output_dir = nxf.run(print_stdout = True, print_stderr = True, print_command = True)
            self.assertEqual(returncode, 0)

            svaba_output = os.path.join(output_dir, 'svaba')
            self.assertTrue(os.path.exists(svaba_output))
            
            preserve_output(tmpdir, prefix = self._testMethodName + '.')


if __name__ == "__main__":
    unittest.main()
