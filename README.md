# boltz-metaflow
Scripts for running the [boltz-1](https://github.com/jwohlwend/boltz) protein structure prediction model on the cloud using [metaflow](https://metaflow.org/).

The examples here are focused on benchmarking accuracy on antibody/antigen complexes but the infrastructure should be readily adaptable to other uses.

## Set up
Install the package:
```bash
pip install -e .
```

Build the docker image:
```bash
docker build -t boltz-1 .
```

You will also probably need to push the docker image to e.g. the AWS ECR registry.
You will need to create the ECR repository first (here called "boltz-1"). You need to do this only once ever:

```bash
aws ecr create-repository --repository-name boltz-1 --region us-east-2
```

Now we can push the docker image to it:

```bash
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
IMAGE=$ACCOUNT.dkr.ecr.us-east-2.amazonaws.com/boltz-1:latest
docker tag boltz-1 $IMAGE
docker push $IMAGE
```

The above is only an example. You may need to change your region or other details.

## Running inference
You first need to configure what predictions will be run.

```bash
python examples/many_seeds_8uzc.py --num-seeds 5 --out /tmp/many_seeds_8uzc.yaml
```

Then you can run the predictions locally:
```bash
python flows/predict.py run --input /tmp/many_seeds_8uzc.yaml --max-workers 1
```

or on AWS:
```bash
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
IMAGE=$ACCOUNT.dkr.ecr.us-east-2.amazonaws.com/boltz-1:latest
QUEUE=metaflow-gpu-g5
python flows/predict.py run --input /tmp/many_seeds_8uzc.yaml --with batch:image=$IMAGE,queue=$QUEUE
```

or the benchmark script (limiting to one prediction for testing)

```bash
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
IMAGE=$ACCOUNT.dkr.ecr.us-east-2.amazonaws.com/boltz-1:latest
QUEUE=metaflow-gpu-g5
python flows/antibody_benchmark.py run --spec examples/antibodies.yaml --with batch:image=$IMAGE,queue=$QUEUE --num_predictions 1
```

In order to run the above, you will need to first create and configure the specified AWS Batch job queue.
