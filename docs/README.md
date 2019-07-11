# Vaporware documentation

Vaporware is a computational pipeline for analyzing data from exome and genome sequencing. Its components are containerized and the pipeline runs on the (Juno high-performance computing cluster)[http://hpc.mskcc.org/] at Memorial Sloan Kettering Cancer Cencer and on Amazon Web Services (AWS).

These pages contain instructions on how to run the Vaporware pipeline. It also contains documentation on the bioinformatic components in the pipeline, some motivation for various parameter choices, plus sources and processing of all reference materials use. 

If there are any questions or comments, you are welcome to [raise an issue](https://github.com/mskcc/vaporware/issues/new?title=[User%20question]).

## Table of contents
1. [Installation](installation.md)
    - [Set up on Juno](juno-setup.md)
    - [Set up on AWS](aws-setup.md)
2. [Usage](usage.md)
3. [Bioinformatic components](bioinformatic-components.md)
    - [Somatic analyses](somatic-analyses.md)
    - [Germline analyses](germline-analyses.md)
4. [Reference files](reference-files.md)
    - [Intervals](intervals.md)
5. [Variant annotation and filtering](variant-annotation-and-filtering.md)
    - [Panel of normals for exomes](wes-panel-of-normals.md)
6. ...
