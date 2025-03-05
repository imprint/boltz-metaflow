"""
Run boltz across seeds on an antibody/antigen complex (PDB 8UZC).
"""
import argparse
import yaml
import copy
import collections
import boltz_metaflow.common

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("--out", metavar="out.yaml")
parser.add_argument("--num-seeds", type=int, default=10)

boltz_input = yaml.safe_load(
"""
sequences:
  - protein:
      id: A
      sequence: MKTIIALSYILCLVFAQKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFKNESFNWTGVTQNGTSSACIRGSSSSFFSRLNWLTHLNYTYPALNVTMPNKEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKHSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGRGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRVQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIESIRNETYDHNVYRDEALNNRFQIKGVELKSGYKDGSGYIPEAPRDGQAYVRKDGEWVLLSTFLGSGLNDIFEAQKIEWHEGSHHHHHH
  - protein:
      id: B
      sequence: QVQLVESGGGVVQPGRSLRLSCATSGFTFSSYGIHWVRQAPGKGLGWVAMISFDGSKTYYADSVRGRFTISRDNSKNTLSLQMNSLRTEDTAVYYCAKERDRDGYNEGIYDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCD
  - protein:
      id: C
      sequence: EIVMTQSPATLSLSPGERATLSCRASQSAGFYLAWYQQKPGQAPRLLIYDTSNRATGIPARFSGRGSGTDFTLTINSLEPEDFAVYYCQQRYNWPITFGQGTRLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC
""".strip())

base_args = [
    "--step_scale", "1.0",
    "--recycling_steps", "1",
]

def run():
    args = parser.parse_args()

    base_task = boltz_metaflow.common.make_task(boltz_input, add_msas=True)

    tasks = {}
    for seed in range(args.num_seeds):
        name = f"seed_{seed}"
        tasks[name] = copy.copy(base_task)
        tasks[name]["extra_args"] = base_args + ["--seed", str(seed)]

    result = yaml.safe_dump(tasks)

    if args.out:
        with open(args.out, "w") as f:
            f.write(result)
        print("Wrote: ", args.out)
    else:
        print(result)

if __name__ == "__main__":
    run()