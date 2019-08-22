# AWS Setup

_Note: These instructions assume the user has a moderate knowledge of AWS concepts, including S3 buckets, EC2 instances, AMIs, EBS volumes, AWS Batch, and autoscaling. For a succinct summary of these concepts, please check the [AWS glossary](aws-glossary.md), the [AWS documentation](https://docs.aws.amazon.com/AmazonS3), and other recent AWS tutorials circa 2018._

Building the compute environment for AWS Batch consist of two steps.

_Prerequisite: [Install](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) and [configure](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html) AWS CLI._

## Creating the S3 Bucket

You must first create an `S3 bucket` which will be used as working directory and where your outputs will be stored, as described [here](https://docs.aws.amazon.com/AmazonS3/latest/user-guide/create-bucket.html).

Use `S3` bucket path for `<AWS-S3-WORKDIR>` value when creating `awsbatch.config`.

## Building the AMI

Vaporware's Amazon Machine Image (AMI) is implemented with [EBS autoscaling script](https://docs.opendata.aws/genomics-workflows/core-env/create-custom-compute-resources) which will automatically increase disk space when needed.

To prepare the AMI for your compute environment, run the command below from the Vaporware root directory.

```shell
aws cloudformation create-stack \
    --stack-name vaporwareAMI \
    --template-body file://aws_cf_scripts/AMICreate.yaml \
    --parameters ParameterKey=AMIName,ParameterValue=vaporware-ami \
    --capabilities CAPABILITY_IAM
```

_Note: You can specify any name instead of `vaporwareAMI` for stack-name._

This command submits CloudFormation Stack for building the [AMI](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html) for execution, which we will use for building Batch Compute Environment. When you submit the command you will get StackId ARN as a response similar to this below.

```json
{
    "StackId": "arn:aws:cloudformation:us-east-1:474622381158:stack/vaporwareAMI/f8d47a90-41c8-11e9-98cc-0eb85d5eff94"
}
```

Building the AMI will take a few minutes. You can check the status of the build process by running `aws cloudformation describe-stacks --stack-name vaporwareAMI`.

When build is complete the `StackStatus` field in the response JSON will have the value `CREATE_COMPLETE` and the `Outputs` section will look similar to:

```json
"Outputs": [
                {
                    "OutputKey": "AMI",
                    "OutputValue": "ami-0599a84d7553f2fb7"
                }
            ]
```
We will use `OutputValue` for the next step.

## Building the Compute Environment

Run the command below from the Vaporware repo root directory and set `<AMI-ID>` parameter to `OutputValue` from previous step. If you skipped previous section do not add the `--parameters` argument when running the command.

```shell
aws cloudformation create-stack \
    --stack-name vaporwareAWSBatchCE \
    --template-body file://aws_cf_scripts/AWSBatchCreate.yaml \
    --capabilities CAPABILITY_IAM \
    --parameters ParameterKey=AmiId,ParameterValue=<AMI-ID>
```

_Note: You can specify any name instead of `vaporwareAWSBatchCE` for stack-name._

Building the AWS Batch Compute Environment will take a few minutes. You can check the status of the build process by running `aws cloudformation describe-stacks --stack-name vaporwareAWSBatchCE`.

When build is complete the `StackStatus` field in response JSON will have the value `CREATE_COMPLETE` and the `Outputs` section will look similar to:

```json
"Outputs": [
                {
                    "OutputKey": "ComputeEnvironmentArn",
                    "OutputValue": "arn:aws:batch:us-east-1:474622381158:compute-environment/ComputeEnvironment-c5b2c9926a6fb42"conf/awsbatch.config
                },
                {
                    "OutputKey": "JobQueueArn",
                    "OutputValue": "arn:aws:batch:us-east-1:474622381158:job-queue/JobQueue-f2090879374b5cc"
                }
            ]
```

Use `OutputValue` of `JobQueueArn` for `<AWS-BATCH-QUEUE-ARN>` value when creating `awsbatch.config`.

## Create Configuration File

Create the `awsbatch.config` file in `conf` directory using `conf/awsbatch.config.template`

Replace the following values:

- `<AWS-REGION>` with aws region where you built your compute environment.
- `<STORAGE-ENCRYPTION>` storage encryption for your bucket (default: 'AES256').
- `<AWS-BATCH-QUEUE-ARN>` ARN of your AWS Batch Job Queue built in [Building the Compute Environment](#Building-the-Compute-Environment) section.
- `<AWS-S3-WORKDIR>` S3 bucket used as working directory created in [Create the S3 Bucket](#Create-the-S3-Bucket) section.
