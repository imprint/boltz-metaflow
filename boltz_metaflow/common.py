import tempfile
import os
import yaml
import collections
import time
import subprocess
import glob
import json
import numpy as np

# Copied from boltz
import tempfile
import contextlib
from pathlib import Path


def run_dockq(reference_structure, predicted_structures):
    """
    Run the DockQ command on the given reference and predicted structures.

    Parameters
    ----------
    reference_structure : prody.AtomGroup
        The reference structure.
    predicted_structures : dict
        A dictionary mapping structure names to prody structures.

    Returns
    -------
    dict
        A dictionary mapping structure names to DockQ scores.
    """
    import prody
    assert len(predicted_structures) > 0
    with tempfile.TemporaryDirectory() as temp_dir:
        with contextlib.chdir(temp_dir):
            # Write the reference structure to a temporary file.
            reference_file = "reference.pdb"
            prody.writePDB(reference_file, reference_structure)

            predicted_files = {}
            result = collections.defaultdict(dict)
            for name, predicted_structure in predicted_structures.items():
                filename = name + ".pdb"
                prody.writePDB(filename, predicted_structure)
                predicted_files[filename] = name

                # Run the DockQ command.
                command = ["DockQ", filename, reference_file, "--json", "out.json"]
                try:
                    r = subprocess.check_output(command, stderr=subprocess.STDOUT).decode()
                except subprocess.CalledProcessError as e:
                    print("DockQ failed with return code", e.returncode)
                    print("DockQ output:\n", e.output.decode())
                    raise
                print(r)
                with open("out.json") as f:
                    result[name] = json.load(f)

            print("DockQ results:")
            print(json.dumps(result, indent=4))
            return result


def make_task(boltz_input, add_msas=True):
    assert "sequences" in boltz_input, list(boltz_input)

    if add_msas:
        protein_sequences = {}
        for item in boltz_input["sequences"]:
            if list(item.keys()) == ["protein"]:
                item = item["protein"]
                msa_name = item["id"] + ".csv"
                item["msa"] = msa_name
                protein_sequences[msa_name] = item["sequence"]

        print("Computing msas")
        msas = compute_msa(protein_sequences)
        print("Done computing msas")
    else:
        msas = {}

    return {
        "input": boltz_input,
        "msas": msas,
        "extra_args": [],
        "info": {},
    }


def compute_msa(
        data: dict[str, str],
        msa_server_url: str = "https://api.colabfold.com",
        msa_pairing_strategy: str = "greedy",
) -> None:
    """Compute the MSA for the input data.

    Parameters
    ----------
    data : dict[str, str]
        The input protein sequences.
    msa_server_url : str
        The MSA server URL.
    msa_pairing_strategy : str
        The MSA pairing strategy.
    """
    from boltz.data import const
    from boltz.data.msa.mmseqs2 import run_mmseqs2

    target_id = "target"
    with tempfile.TemporaryDirectory() as msa_dir:
        msa_dir = Path(msa_dir)

        if len(data) > 1:
            paired_msas = run_mmseqs2(
                list(data.values()),
                msa_dir / f"{target_id}_paired_tmp",
                use_env=True,
                use_pairing=True,
                host_url=msa_server_url,
                pairing_strategy=msa_pairing_strategy,
            )
        else:
            paired_msas = [""] * len(data)

        unpaired_msa = run_mmseqs2(
            list(data.values()),
            msa_dir / f"{target_id}_unpaired_tmp",
            use_env=True,
            use_pairing=False,
            host_url=msa_server_url,
            pairing_strategy=msa_pairing_strategy,
        )

        csv_strings = {}
        for idx, name in enumerate(data):
            # Get paired sequences
            paired = paired_msas[idx].strip().splitlines()
            paired = paired[1::2]  # ignore headers
            paired = paired[: const.max_paired_seqs]

            # Set key per row and remove empty sequences
            keys = [idx for idx, s in enumerate(paired) if s != "-" * len(s)]
            paired = [s for s in paired if s != "-" * len(s)]

            # Combine paired-unpaired sequences
            unpaired = unpaired_msa[idx].strip().splitlines()
            unpaired = unpaired[1::2]
            unpaired = unpaired[: (const.max_msa_seqs - len(paired))]
            if paired:
                unpaired = unpaired[1:]  # ignore query is already present

            # Combine
            seqs = paired + unpaired
            keys = keys + [-1] * len(unpaired)

            # Dump MSA
            csv_str = ["key,sequence"] + [f"{key},{seq}" for key, seq in zip(keys, seqs)]
            csv_strings[name] = csv_str

    return csv_strings


def run_boltz(sequences_dict, msa_files={}, extra_args=[]):
    """
    Run the boltz predict command on the given sequences and MSA files.

    Parameters
    ----------
    sequences_dict : dict
        A dictionary mapping sequence names to sequences.
    msa_files : dict
        A dictionary mapping MSA filenames to their contents.
    extra_args : list
        Extra arguments to pass to the predict command.

    Returns
    -------
    dict
        A dictionary containing the following keys:
            predict_returncode : int
                The return code of the predict command.
            predict_seconds : float
                The number of seconds the predict command took to run.
            structure : dict
                A dictionary mapping structure filenames to prody structures.
            mmcif : dict
                A dictionary mapping mmcif filenames to their contents.
            json : dict
                A dictionary mapping json filenames to their contents.
            npz : dict
                A dictionary mapping npz filenames to their contents.
    """
    import prody

    # Check that the boltz weights have been downloaded.
    assert os.path.exists(os.path.join(os.environ["HOME"], ".boltz", "boltz1_conf.ckpt"))
    with tempfile.TemporaryDirectory() as temp_dir:
        old_cwd = os.getcwd()
        try:
            os.chdir(temp_dir)

            # Write the sequences to a temporary file.
            sequences_file = "sequences.yaml"
            with open(sequences_file, "w") as f:
                yaml.dump(sequences_dict, f)

            # Write the MSA files to temporary files.
            for filename, msa_contents in msa_files.items():
                with open(filename, "w") as f:
                    f.writelines(f"{line}\n" for line in msa_contents)

            # Run the command.
            command = [
              "boltz",
              "predict",
              sequences_file,
              "--override",
              "--output_format", "pdb",
            ] + extra_args

            result = collections.defaultdict(dict)
            start = time.time()
            env = os.environ.copy()
            env["TORCH_CPP_LOG_LEVEL"] = "info"

            # Run and capture the output while also printing it
            r = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
            result["predict_returncode"] = r.returncode
            result["predict_seconds"] = time.time() - start
            result["predict_output"] = r.stdout.decode()
            print("Predict output:\n****\n", result["predict_output"], "\n****")
            print("Predict took", result["predict_seconds"], "seconds and returned", r.returncode)

            # Save structures (pdb, parsed)
            for filename in sorted(glob.glob("boltz_results_sequences/predictions/sequences/*.pdb")):
                result["structure"][os.path.basename(filename)] = prody.parsePDB(filename)

            # Save structures (mmcif, unparsed)
            for filename in sorted(glob.glob("boltz_results_sequences/predictions/sequences/*.cif")):
                result["mmcif"][os.path.basename(filename)] = open(filename).read()

            if not result["structure"] and not result["mmcif"]:
                raise ValueError("No structures were generated")

            # Save the json
            for filename in sorted(glob.glob("boltz_results_sequences/predictions/sequences/*.json")):
                result["json"][os.path.basename(filename)] = json.load(open(filename))

            # Save the npz
            for filename in sorted(glob.glob("boltz_results_sequences/predictions/sequences/*.npz")):
                result["npz"][os.path.basename(filename)] = dict(np.load(filename))
        finally:
            os.chdir(old_cwd)

        print("Done with run_boltz() call. Result keys: ", result.keys())
        return result
