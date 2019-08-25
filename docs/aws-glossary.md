# AWS Glossary

This page provides a basic reference to cloud computing concepts and AWS terminology, with links to the official AWS documentation provided.

* __S3 Buckets__: [Amazon Simple Storage Service (S3)](https://aws.amazon.com/s3/) is designed for data storage in the cloud. Users create an S3 bucket within a specific AWS Region via the [AWS Command Line Interface (CLI)](https://aws.amazon.com/cli/) or the [AWS Management Console](https://aws.amazon.com/console/) web interface. This location is then used much like a standard subdirectory to upload, store, and download data. 

* __EC2 Instances__: [Amazon Elastic Compute Cloud (Amazon EC2)](https://aws.amazon.com/ec2/) provides the computational resource for AWS users, which come in a variety of instance sizes, operating systems, number of cores, and prices. 

* __EBS Volumes__: [Amazon Elastic Block Store (Amazon EBS)](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AmazonEBS.html) is the storage volumes used for EC2 instances for computation. 

* __AWS Batch__: [Amazon Batch](https://aws.amazon.com/batch/) dynamically provisions the optimal quantity and type of compute resources (both Spot or On-Demand) necessary for analysis. It functions analogously to a relatively simple job scheduler used by HPC environments (e.g. LSF or SGE) to schedule and execute batch computing workloads. You can read this [blog post](https://aws.amazon.com/blogs/compute/building-high-throughput-genomics-batch-workflows-on-aws-introduction-part-1-of-4/) for more information related to using AWS Batch for running Genomics Workflows.

* __EBS Autoscaling__: [Autoscalling EBS](https://docs.opendata.aws/genomics-workflows/core-env/create-custom-compute-resources/) EC2 instances have a running process which monitors disk usage and add more EBS volumes on the fly to expand the free space based on the capacity threshold. This feature increass in volume size, adjustments to performance, or changes in the volume type while the volume is in use, e.g. the disk size will increase automatically while processes are running. 
