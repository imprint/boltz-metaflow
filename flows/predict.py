import tempfile

from metaflow import FlowSpec, step, batch, pypi, environment, resources, Parameter, IncludeFile
import yaml

from . import common

class Boltz1PredictionFlow(FlowSpec):
    # Parameters for dynamic configuration
    input_yaml = IncludeFile(
        "input",
        required=True,
        is_text=True,
        help="Input yaml of prediction tasks")

    @step
    def start(self):
        self.tasks = yaml.safe_load(self.input_yaml)
        self.task_ids = list(self.tasks)
        self.next(self.predict, foreach="task_ids")

    @resources(cpu=4, memory=32000, shared_memory=8000, gpu=1)
    @step
    def predict(self):
        """
        Run an antibody/antigen complex prediction.
        """
        task_id = self.input
        task = self.tasks[task_id]
        prediction = common.run_boltz(
            task["input"],
            msa_files=task["msas"],
            extra_args=task["extra_args"] + [])

        self.result = {
            "task_id": task_id,
            "task": task,
            "result": prediction,
        }
        self.next(self.join)

    @step
    def join(self, inputs):
        """
        Join the parallel predictions.
        """
        self.predictions = [inp.result for inp in inputs]
        self.next(self.end)

    @step
    def end(self):
        """
        End the flow.
        """
        print("Flow is complete!")



if __name__ == "__main__":
    Boltz1PredictionFlow()
