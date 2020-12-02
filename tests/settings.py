#!/usr/bin/env python3
"""
Module for settings to run the test pipelines

Use configparser so that we can get more flexible and dynamic usages

Keys from config_str which are present in the user environment will override values here; interpolated values will be updated as well

Example;

$ ( TEST_DIR=foo ./settings.py )
"""
import os
import getpass
from configparser import ConfigParser, ExtendedInterpolation
config = ConfigParser(interpolation=ExtendedInterpolation())
config.optionxform = str # preserve case for keys

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(TEST_DIR)
username = getpass.getuser()

config_str = f"""
[settings]
TEST_DIR = {TEST_DIR}
PROJECT_DIR = {PROJECT_DIR}

USER_SCRATCH = /scratch/{username}

CONF_DIR = ${{PROJECT_DIR}}/conf
TEST_INPUTS_DIR = ${{PROJECT_DIR}}/test_inputs

# needs to be downloaded; `make test-data`
TEST_DATA = ${{TEST_DIR}}/test-data

# # config files
# in this dir;
LOCAL_CONFIG = ${{TEST_DIR}}/local.config
INTEGRATION_TEST_CONFIG = ${{TEST_DIR}}/test.config

# in project dir;
NEXTFLOW_CONFIG = ${{PROJECT_DIR}}/nextflow.config

# in the root conf dir;
TRAVIS_TEST_CONFIG = ${{CONF_DIR}}/test.config
CONTAINERS_CONFIG = ${{CONF_DIR}}/containers.config
RESOURCES_CONFIG = ${{CONF_DIR}}/resources.config
REFERENCES_CONFIG = ${{CONF_DIR}}/references.config
EXOME_CONFIG = ${{CONF_DIR}}/exome.config
GENOME_CONFIG = ${{CONF_DIR}}/genome.config
SINGULARITY_CONFIG = ${{CONF_DIR}}/singularity.config
JUNO_CONFIG = ${{CONF_DIR}}/juno.config

[test_files]
# mapping and pairing files to use with testing
small_test_mapping = ${{settings:TEST_INPUTS_DIR}}/local/small_test_mapping.tsv
small_test_pairing = ${{settings:TEST_INPUTS_DIR}}/local/small_test_pairing.tsv
tiny_test_mapping = ${{settings:TEST_INPUTS_DIR}}/local/tiny_test_mapping.tsv
test_make_bam_and_qc = ${{settings:TEST_DIR}}/test_make_bam_and_qc.tsv
test_somatic = ${{settings:TEST_DIR}}/test_somatic.tsv
test_make_bam_and_qc_pairing = ${{settings:TEST_DIR}}/test_make_bam_and_qc_pairing.tsv
test_pairing_duplicate = ${{settings:TEST_DIR}}/test_pairing_duplicate.tsv
duplicate_samplelane_makebamqc = ${{settings:TEST_DIR}}/duplicate_samplelane_makebamqc.tsv

[nextflow_args]
# tmpdir var used inside Tempo pipeline
TMPDIR = ${{settings:USER_SCRATCH}}
PIPELINE_SCRIPT = ${{settings:PROJECT_DIR}}/pipeline.nf
NXF_SINGULARITY_CACHEDIR = /juno/work/taylorlab/cmopipeline/singularity_images
NXF_ANSI_LOG = false
"""
config.read_string(config_str)

# import these objects into the other test modules
settings = config['settings']
nextflow_args = config['nextflow_args']
test_files = config['test_files']

# update config if the key is present in the environment
for key, value in config['settings'].items():
    if key in os.environ:
        config['settings'][key] = os.environ[key]

for key, value in config['nextflow_args'].items():
    if key in os.environ:
        config['nextflow_args'][key] = os.environ[key]

if __name__ == '__main__':
    for key, value in config['settings'].items():
        print(key, value)
    for key, value in config['nextflow_args'].items():
        print(key, value)
