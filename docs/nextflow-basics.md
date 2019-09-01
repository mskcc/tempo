# Nextflow Basics

[Nextflow](https://nextflow.io) is a workflow framework that creation of computational pipelines that work in any POSIX-based environment. Nextflow is written in the Java-based [Groovy](https://groovy-lang.org/) language. A Nextflow script is executed by running `nextflow run script.nf`, followed by optional arguments. This page contains some useful information for Tempo users, for more read the [Nextflow documentation](https://www.nextflow.io/docs/latest/basic.html).

* __Nextflow runs scripts__: The basic component of Nextflow is a _process_. Code inside a process is executed as a Bash script, and is thus written essentially as on the command line. Therefore you can easily see what the contents of the different steps in Tempo are by peaking at the source code.

* __The number of dashes matters__: Nextflow has a quirk where its executor-specific, built-in flags are initiated with a single dash, whereas parameters we define require two. For example, to define which run profile to use, we call Nextflow's built-in feature `-profile`. Similarly, resume pipeline executions at a certain step, we use `-resume`. However, other parameters require double dashes. For instance, in order to provide mapping or pairing input file paths, we call arguments `--mapping` and `--pairing` respectively.

* __Configuration files__: Upon running a Nextflow script, you can load a `*.config` file with various kinds of preconfigured parameters. The `-profile` argument loads the configuration files associated with a user-defined profile, which for Tempo exist for running the pipeline on Juno and AWS. 

* __View the Nextflow log__: You can access the Nextflow cache metadata by running `nextflow log`. It contains information like TIMESTAMP, DURATION, and RUN_NAME, and a STATUS indicating a failed or successful run, among others. The values under RUN_NAME can also be submitted following the `-resume` flag to resume previously-run Nextflow jobs.

* __View intermediate output__: As the pipeline runs, everything needed to execute each process in the pipeline is located in the `work` in the run directory. Thus, you can peek at input and output files for each step of the pipeline in real time.

* __Run or skip specific tools:__ `pipeline.nf` has the argument `--tools`, which allows users to run only certain bioinformatic tools and skip others. The `--somatic` and `--germline` flags already have a preset of tools to include during a run, but you can limit this further by providing the `--tools` flag, followed by a comma-delimited string. For example, to use only DELLY for your somatic/germline runs, do `--somatic --germline --tools delly`; to use MuTect2, Manta, and Strelka2, do `--somatic --tools mutect2,manta,strelka2`.
