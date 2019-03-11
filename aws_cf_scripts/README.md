# Building AWS Batch Compute Environment

Building the compute environment for AWS Batch consist of two steps.

***Prerequisite: [Install](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) and [configure](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html) AWS CLI.***

## Create the S3 Bucket

You must first create `S3 bucket` which will be used as working directory and where your outputs will be stored, as described [here](https://docs.aws.amazon.com/AmazonS3/latest/user-guide/create-bucket.html).

Use `S3` bucket path for `<AWS-S3-WORKDIR>` value when creating `awsbatch.config`.

## Building the AMI

***NOTE: You need this step only if default disk size of 200GB is not enough for processing your data. If 200GB is suficient please got to the `Building the Compute Environment` section.***

Run the command below from the vaporware repo root directory and set `<DiskSize>` parameter to a [proper](#Calculating-the-DiskSize-value) value based on the input file sizes (NOTE: minimal disk size is 200GB).

`aws cloudformation create-stack --stack-name vaporwareAMI --template-body file://aws_cf_scripts/AMICreate.yaml --parameters ParameterKey=AMIName,ParameterValue=vaporware-ami ParameterKey=DiskSize,ParameterValue=<DiskSize> --capabilities CAPABILITY_IAM`

This command submits CloudFormation Stack for building the [AMI](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html) for execution, which we will use for building Batch Compute Environment. When you submit the command you will get StackId ARN as a responce similar to this below.

`{
    "StackId": "arn:aws:cloudformation:us-east-1:474622381158:stack/vaporwareAMI/f8d47a90-41c8-11e9-98cc-0eb85d5eff94"
}`

Building the AMI will last for a few minutes. You can check the status of the build process by running the command.

`aws cloudformation describe-stacks --stack-name vaporwareAMI`

When build is complete the `StackStatus` field in response json will have the value `CREATE_COMPLETE` and the `Outputs` section will look similar to:

```
"Outputs": [
                {
                    "OutputKey": "AMI",
                    "OutputValue": "ami-0599a84d7553f2fb7"
                }
            ]
```
We will use `OutputValue` for the next step.

## Building the Compute Environment

Run the command below from the vaporware repo root directory and set `<AMI-ID>` parameter to `OutputValue` from previous step. If you skiped previous section do not add --parameters argument when running the command.

`aws cloudformation create-stack --stack-name vaporwareAWSBatchCE --template-body file://aws_cf_scripts/AWSBatchCreate.yaml --parameters ParameterKey=AmiId,ParameterValue=<AMI-ID> --capabilities CAPABILITY_IAM`

Building the AWS Batch Compute Environment will last for a few minutes. You can check the status of the build process by running the command.

***NOTE: By default EC2 Instance Role will have access to all your buckets. If you want to limit access of EC2 Instance Role just to the bucket you created you need to add `--parameters ParameterKey=WorkDirBucket,ParameterValue=<AWS-S3-WORKDIR>` in the previous command.***

`aws cloudformation describe-stacks --stack-name aws cloudformation describe-stacks --stack-name `

When build is complete the `StackStatus` field in response json will have the value `CREATE_COMPLETE` and the `Outputs` section will look similar to:

```
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

## Create Configuration file

Create `awsbatch.config` file in `conf` directory using `conf/awsbatch.config.template`

Replace:

`<AWS-REGION>` with aws region where you built your compute environment

`<STORAGE-ENCRYPTION>` storage encryption for your bucket (default: 'AES256')

`<AWS-BATCH-QUEUE-ARN>` ARN of your AWS Batch Job Queue built in [Building the Compute Environment](#Building-the-Compute-Environment) section

`<AWS-S3-WORKDIR>` S3 bucket used as working directory created in [Create the S3 Bucket](#Create-the-S3-Bucket) section

## Calculating the DiskSize value

### TODO