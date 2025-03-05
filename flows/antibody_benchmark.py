import tempfile

from metaflow import FlowSpec, step, batch, pypi, environment, resources, Parameter, IncludeFile
import subprocess
import yaml
import json
import tempfile
import collections
import os
import glob
import time
import copy
import pandas
import numpy as np

import common

class Boltz1AntibodyBenchmarkFlow(FlowSpec):
    # Parameters for dynamic configuration
    spec = IncludeFile(
        "spec",
        required=True,
        is_text=True,
        help="Input yaml of antibodies")

    boltz_extra_args = Parameter(
        "boltz_extra_args",
        default="--step_scale 1.0 --diffusion_samples 10",
        help="Space-separated list of extra arguments to pass to the boltz predict command")

    num_predictions = Parameter(
        "num_predictions",
        default=None,
        type=int,
        help="Number of predictions to run (to limit the number of predictions for testing)")

    def expand_tasks(self, base_tasks):
        """
        This is where we implement parameter sweeps.
        """
        new_task_dict = collections.OrderedDict()
        for task_id, task in base_tasks.items():
            for num_recycles in range(0, 51, 5):
                task = copy.deepcopy(task)
                task["extra_args"].extend(["--recycling_steps", str(num_recycles)])
                task["info"]["recycles"] = num_recycles
                new_task_dict[f"{task_id}_recycles-{num_recycles}"] = task
        print("Expanded tasks from n=", len(base_tasks), "to n=", len(new_task_dict))
        print("Expanded tasks: ", new_task_dict.keys())
        return new_task_dict

    @resources(cpu=4, memory=32000, shared_memory=8000)
    @step
    def start(self):
        import prody

        self.parsed_spec = yaml.safe_load(self.spec)
        self.parsed_extra_args = self.boltz_extra_args.split()
        print("Loaded spec: ", self.parsed_spec)
        print("Loaded extra args: ", self.parsed_extra_args)

        base_tasks = collections.OrderedDict()
        self.complexes = []
        for (i, complex) in enumerate(self.parsed_spec['complexes']):
            self.complexes.append(complex)
            atoms, header = prody.parsePDB(complex["pdb"], header=True)

            sequence_dict = dict((item.chid, item.sequence) for item in header["polymers"])

            sequences_to_predict = [
                {
                    "protein": {
                        "id": "A",
                        "sequence": sequence_dict[complex['antigen_chain']],
                    }
                },
                {
                    "protein": {
                        "id": "B",
                        "sequence": sequence_dict[complex['heavy_chain']],
                    }
                },
                {
                    "protein": {
                        "id": "C",
                        "sequence": sequence_dict[complex['light_chain']],
                    }
                },
            ]
            boltz_input = {
                "sequences": sequences_to_predict,
            }
            task = common.make_task(boltz_input, add_msas=True)
            task["extra_args"].extend(self.parsed_extra_args)
            task["info"] = dict(complex)
            task["info"]["complex_id"] = i
            base_tasks[complex["name"]] = task

        self.tasks = self.expand_tasks(base_tasks)
        if self.num_predictions is not None:
            self.tasks = dict(list(self.tasks.items())[:self.num_predictions])

        for (name, task) in self.tasks.items():
            task["task_id"] = name

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

        assert len(prediction["structure"]) > 0

        self.result = {
            "task_id": task_id,
            "task": task,
            "result": prediction,
        }
        self.next(self.join)

    @resources(cpu=4, memory=128000, shared_memory=8000)
    @step
    def join(self, inputs):
        """
        Join the parallel predictions.
        """
        self.predictions = [inp.result for inp in inputs]
        self.merge_artifacts(inputs, exclude=["result"])
        self.next(self.dockq_score)

    @resources(cpu=4, memory=128000, shared_memory=8000)
    @step
    def dockq_score(self):
        import prody

        # Group predictions by complex_id
        predictions_by_complex = collections.defaultdict(list)
        for prediction in self.predictions:
            complex_id = prediction["task"]["info"]["complex_id"]
            predictions_by_complex[complex_id].append(prediction)

        all_scores = []
        for complex_id, predictions in predictions_by_complex.items():
            assert len(predictions) > 0

            complex = self.complexes[complex_id]
            reference, header = prody.parsePDB(complex["pdb"], header=True)
            reference = reference.select("protein")

            predicted_structures = {}
            model_name_to_task = {}
            for prediction in predictions:
                task = prediction["task"]
                result = prediction["result"]
                pdbs = result["structure"]
                assert len(pdbs) > 0
                for pdb_name, structure in pdbs.items():
                    model_name = f"{complex['name']}_{task['task_id']}_{pdb_name.replace('.pdb', '')}"
                    predicted_structures[model_name] = structure
                    model_name_to_task[model_name] = task

            scores = common.run_dockq(reference, predicted_structures)
            for model_name, score_dict in scores.items():
                score_dict["model_name"] = model_name
                score_dict["complex"] = complex['name']
                score_dict.update(model_name_to_task[model_name])
                all_scores.append(score_dict)

        self.summary_df = pandas.DataFrame.from_records(all_scores)
        print(self.summary_df)
        self.next(self.end)

    @step
    def end(self):
        """
        End the flow.
        """
        print("Flow is complete!")



if __name__ == "__main__":
    Boltz1AntibodyBenchmarkFlow()
