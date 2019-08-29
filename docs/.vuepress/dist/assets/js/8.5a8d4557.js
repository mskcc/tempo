(window.webpackJsonp=window.webpackJsonp||[]).push([[8],{190:function(e,t,a){"use strict";a.r(t);var s=a(0),n=Object(s.a)({},function(){var e=this,t=e.$createElement,a=e._self._c||t;return a("ContentSlotsDistributor",{attrs:{"slot-key":e.$parent.slotKey}},[a("h1",{attrs:{id:"aws-setup"}},[a("a",{staticClass:"header-anchor",attrs:{href:"#aws-setup","aria-hidden":"true"}},[e._v("#")]),e._v(" AWS Setup")]),e._v(" "),a("p",[a("em",[e._v("Note: These instructions assume the user has a moderate knowledge of AWS concepts, including S3 buckets, EC2 instances, AMIs, EBS volumes, AWS Batch, and autoscaling. For a succinct summary of these concepts, please check the "),a("router-link",{attrs:{to:"/aws-glossary.html"}},[e._v("AWS glossary")]),e._v(", the "),a("a",{attrs:{href:"https://docs.aws.amazon.com/AmazonS3",target:"_blank",rel:"noopener noreferrer"}},[e._v("AWS documentation"),a("OutboundLink")],1),e._v(", and other recent AWS tutorials circa 2018.")],1)]),e._v(" "),a("p",[e._v("Building the compute environment for AWS Batch consist of two steps.")]),e._v(" "),a("div",{staticClass:"tip custom-block"},[a("p",{staticClass:"custom-block-title"},[e._v("NOTE")]),e._v(" "),a("p",[e._v("You need to first "),a("a",{attrs:{href:"https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html",target:"_blank",rel:"noopener noreferrer"}},[e._v("install"),a("OutboundLink")],1),e._v(" and "),a("a",{attrs:{href:"https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html",target:"_blank",rel:"noopener noreferrer"}},[e._v("configure"),a("OutboundLink")],1),e._v(" AWS command-line interface (CLI).")])]),e._v(" "),a("h2",{attrs:{id:"creating-the-s3-bucket"}},[a("a",{staticClass:"header-anchor",attrs:{href:"#creating-the-s3-bucket","aria-hidden":"true"}},[e._v("#")]),e._v(" Creating the S3 Bucket")]),e._v(" "),a("p",[e._v("You must first create an "),a("code",[e._v("S3 bucket")]),e._v(" which will be used as working directory and where your outputs will be stored, as described "),a("a",{attrs:{href:"https://docs.aws.amazon.com/AmazonS3/latest/user-guide/create-bucket.html",target:"_blank",rel:"noopener noreferrer"}},[e._v("here"),a("OutboundLink")],1),e._v(".")]),e._v(" "),a("p",[e._v("Use "),a("code",[e._v("S3")]),e._v(" bucket path for "),a("code",[e._v("<AWS-S3-WORKDIR>")]),e._v(" value when creating "),a("code",[e._v("awsbatch.config")]),e._v(".")]),e._v(" "),a("h2",{attrs:{id:"building-the-ami"}},[a("a",{staticClass:"header-anchor",attrs:{href:"#building-the-ami","aria-hidden":"true"}},[e._v("#")]),e._v(" Building the AMI")]),e._v(" "),a("p",[e._v("Tempo's Amazon Machine Image (AMI) is implemented with "),a("a",{attrs:{href:"https://docs.opendata.aws/genomics-workflows/core-env/create-custom-compute-resources",target:"_blank",rel:"noopener noreferrer"}},[e._v("EBS autoscaling script"),a("OutboundLink")],1),e._v(" which will automatically increase disk space when needed.")]),e._v(" "),a("p",[e._v("To prepare the AMI for your compute environment, run the command below from the Tempo root directory.")]),e._v(" "),a("div",{staticClass:"language-shell extra-class"},[a("pre",{pre:!0,attrs:{class:"language-text"}},[a("code",[e._v("aws cloudformation create-stack \\\n    --stack-name vaporwareAMI \\\n    --template-body file://aws_cf_scripts/AMICreate.yaml \\\n    --parameters ParameterKey=AMIName,ParameterValue=vaporware-ami \\\n    --capabilities CAPABILITY_IAM\n")])])]),a("p",[a("em",[e._v("Note: You can specify any name instead of "),a("code",[e._v("vaporwareAMI")]),e._v(" for stack-name.")])]),e._v(" "),a("p",[e._v("This command submits CloudFormation Stack for building the "),a("a",{attrs:{href:"https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html",target:"_blank",rel:"noopener noreferrer"}},[e._v("AMI"),a("OutboundLink")],1),e._v(" for execution, which we will use for building Batch Compute Environment. When you submit the command you will get StackId ARN as a response similar to this below.")]),e._v(" "),a("div",{staticClass:"language-json extra-class"},[a("pre",{pre:!0,attrs:{class:"language-json"}},[a("code",[a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("{")]),e._v("\n    "),a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"StackId"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token string"}},[e._v('"arn:aws:cloudformation:us-east-1:474622381158:stack/vaporwareAMI/f8d47a90-41c8-11e9-98cc-0eb85d5eff94"')]),e._v("\n"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("}")]),e._v("\n")])])]),a("p",[e._v("Building the AMI will take a few minutes. You can check the status of the build process by running "),a("code",[e._v("aws cloudformation describe-stacks --stack-name vaporwareAMI")]),e._v(".")]),e._v(" "),a("p",[e._v("When build is complete the "),a("code",[e._v("StackStatus")]),e._v(" field in the response JSON will have the value "),a("code",[e._v("CREATE_COMPLETE")]),e._v(" and the "),a("code",[e._v("Outputs")]),e._v(" section will look similar to:")]),e._v(" "),a("div",{staticClass:"language-json extra-class"},[a("pre",{pre:!0,attrs:{class:"language-json"}},[a("code",[a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"Outputs"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("[")]),e._v("\n                "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("{")]),e._v("\n                    "),a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"OutputKey"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token string"}},[e._v('"AMI"')]),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v("\n                    "),a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"OutputValue"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token string"}},[e._v('"ami-0599a84d7553f2fb7"')]),e._v("\n                "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("}")]),e._v("\n            "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("]")]),e._v("\n")])])]),a("p",[e._v("We will use "),a("code",[e._v("OutputValue")]),e._v(" for the next step.")]),e._v(" "),a("h2",{attrs:{id:"building-the-compute-environment"}},[a("a",{staticClass:"header-anchor",attrs:{href:"#building-the-compute-environment","aria-hidden":"true"}},[e._v("#")]),e._v(" Building the Compute Environment")]),e._v(" "),a("p",[e._v("Run the command below from the Vaporware repo root directory and set "),a("code",[e._v("<AMI-ID>")]),e._v(" parameter to "),a("code",[e._v("OutputValue")]),e._v(" from previous step. If you skipped previous section do not add the "),a("code",[e._v("--parameters")]),e._v(" argument when running the command.")]),e._v(" "),a("div",{staticClass:"language-shell extra-class"},[a("pre",{pre:!0,attrs:{class:"language-text"}},[a("code",[e._v("aws cloudformation create-stack \\\n    --stack-name vaporwareAWSBatchCE \\\n    --template-body file://aws_cf_scripts/AWSBatchCreate.yaml \\\n    --capabilities CAPABILITY_IAM \\\n    --parameters ParameterKey=AmiId,ParameterValue=<AMI-ID>\n")])])]),a("p",[a("em",[e._v("Note: You can specify any name instead of "),a("code",[e._v("vaporwareAWSBatchCE")]),e._v(" for stack-name.")])]),e._v(" "),a("p",[e._v("Building the AWS Batch Compute Environment will take a few minutes. You can check the status of the build process by running "),a("code",[e._v("aws cloudformation describe-stacks --stack-name vaporwareAWSBatchCE")]),e._v(".")]),e._v(" "),a("p",[e._v("When build is complete the "),a("code",[e._v("StackStatus")]),e._v(" field in response JSON will have the value "),a("code",[e._v("CREATE_COMPLETE")]),e._v(" and the "),a("code",[e._v("Outputs")]),e._v(" section will look similar to:")]),e._v(" "),a("div",{staticClass:"language-json extra-class"},[a("pre",{pre:!0,attrs:{class:"language-json"}},[a("code",[a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"Outputs"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("[")]),e._v("\n                "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("{")]),e._v("\n                    "),a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"OutputKey"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token string"}},[e._v('"ComputeEnvironmentArn"')]),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v("\n                    "),a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"OutputValue"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token string"}},[e._v('"arn:aws:batch:us-east-1:474622381158:compute-environment/ComputeEnvironment-c5b2c9926a6fb42"')]),e._v("conf/awsbatch.config\n                "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("}")]),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v("\n                "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("{")]),e._v("\n                    "),a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"OutputKey"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token string"}},[e._v('"JobQueueArn"')]),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v("\n                    "),a("span",{pre:!0,attrs:{class:"token property"}},[e._v('"OutputValue"')]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v(":")]),e._v(" "),a("span",{pre:!0,attrs:{class:"token string"}},[e._v('"arn:aws:batch:us-east-1:474622381158:job-queue/JobQueue-f2090879374b5cc"')]),e._v("\n                "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("}")]),e._v("\n            "),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("]")]),e._v("\n")])])]),a("p",[e._v("Use "),a("code",[e._v("OutputValue")]),e._v(" of "),a("code",[e._v("JobQueueArn")]),e._v(" for "),a("code",[e._v("<AWS-BATCH-QUEUE-ARN>")]),e._v(" value when creating "),a("code",[e._v("awsbatch.config")]),e._v(".")]),e._v(" "),a("h2",{attrs:{id:"create-configuration-file"}},[a("a",{staticClass:"header-anchor",attrs:{href:"#create-configuration-file","aria-hidden":"true"}},[e._v("#")]),e._v(" Create Configuration File")]),e._v(" "),a("p",[e._v("Create the "),a("code",[e._v("awsbatch.config")]),e._v(" file in "),a("code",[e._v("conf")]),e._v(" directory using "),a("code",[e._v("conf/awsbatch.config.template")])]),e._v(" "),a("p",[e._v("Replace the following values:")]),e._v(" "),a("ul",[a("li",[a("code",[e._v("<AWS-REGION>")]),e._v(" with aws region where you built your compute environment.")]),e._v(" "),a("li",[a("code",[e._v("<STORAGE-ENCRYPTION>")]),e._v(" storage encryption for your bucket (default: 'AES256').")]),e._v(" "),a("li",[a("code",[e._v("<AWS-BATCH-QUEUE-ARN>")]),e._v(" ARN of your AWS Batch Job Queue built in "),a("a",{attrs:{href:"#Building-the-Compute-Environment"}},[e._v("Building the Compute Environment")]),e._v(" section.")]),e._v(" "),a("li",[a("code",[e._v("<AWS-S3-WORKDIR>")]),e._v(" S3 bucket used as working directory created in "),a("a",{attrs:{href:"#Create-the-S3-Bucket"}},[e._v("Create the S3 Bucket")]),e._v(" section.")])])])},[],!1,null,null,null);t.default=n.exports}}]);