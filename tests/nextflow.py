#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for the Nexflow class
"""
import os
import subprocess
from settings import nextflow_args

class Nextflow(object):
    """
    Class to run a Nextflow pipeline in a temp dir for testing

    >>> args = ['--foo', 'bar']
    >>> configs = ['containers.config', 'resources.config']
    >>> nxf = Nextflow(tmpdir = '.', args = args, configs = configs)
    >>> nxf.run()
    """
    def __init__(self,
        tmpdir, # TemporaryDirectory path to run in
        args = (), # CLI args to pass to `nextflow run`
        configs = (), # extra Nextflow config files to use
        settings = { k:v for k,v in nextflow_args.items() }, # extra nextflow_args from the settings module
        **kwargs # pass extra keyword args to override the defaults if needed
        ):
        self.tmpdir = tmpdir
        self.args = args
        self.settings = settings

        # collect the Nextflow config files passed
        self.configs = []
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
        'run', self.settings['PIPELINE_SCRIPT'],
        *args
        ]

        # set up a dict of extra environment variables that we need to have set
        env = dict(
            os.environ,
            NXF_ANSI_LOG = self.settings['NXF_ANSI_LOG'],
            NXF_WORK = NXF_WORK,
            TMPDIR = self.settings['TMPDIR'],
            NXF_SINGULARITY_CACHEDIR = self.settings['NXF_SINGULARITY_CACHEDIR']
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
