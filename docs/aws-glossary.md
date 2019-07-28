# AWS Glossary

This page provides a basic reference to cloud computing concepts and AWS terminology, with links to the official AWS documentation provided.

* __S3 Buckets__: [Amazon Simple Storage Service (S3)](https://aws.amazon.com/s3/) is designed for data storage in the cloud. Users create an S3 bucket within a specific AWS Region via the [AWS Command Line Interface (CLI)](https://aws.amazon.com/cli/) or the [AWS Management Console](https://aws.amazon.com/console/) web interface. This location is then used much like a standard subdirectory to upload, store, and download data. 

* __EC2 Instances__: [Amazon Elastic Compute Cloud (Amazon EC2)](https://aws.amazon.com/ec2/) provides the computational resource for AWS users, which come in a variety of instance sizes, operating systems, number of cores, and prices. 

* __EBS Volumes__: [Amazon Elastic Block Store (Amazon EBS)](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AmazonEBS.html) is the storage volumes used for EC2 instances for computation. 

* __AWS Batch__: [AWS Batch](https://aws.amazon.com/batch/) is a relatively new feature for AWS customers which schedules and executes batch computing workloads. It works as a relatively simple job scheduler used by HPC clusters (e.g. SGE or LSF), using both EC2 instances and Spot instances.

* __EBS Autoscaling__: Feature which allows ncrease volume size, adjust performance, or change the volume type while the volume is in use, e.g. the disk size will increase automatically while processes are running. 
