# Installation

## Installing Nextflow
[Nextflow](https://www.nextflow.io) requires Java 8 or later. You can check the version on your system with the command `java -version`.

Install Nextflow in the current directory by running:
```shell
curl -s https://get.nextflow.io | bash
```
Put the `nextflow` executable in a directory in your `PATH`, if you want to access it from anywhere. For more details, check out the [documentation](https://www.nextflow.io/docs/latest/getstarted.html).

We recommend Nextflow version 19.07.0 or later.

## Installing Tempo
Clone the [Tempo repository](http://github.com/mskcc/tempo):
```shell
git clone http://github.com/mskcc/tempo.git
```

You're now good to go!

For specifics on running Tempo in different environments, check out the documentation on [Juno](juno-setup.md) and [AWS](aws-setup.md).
