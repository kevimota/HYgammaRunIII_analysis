import os
import getpass
import subprocess

# from pprint import pprint

data_path = f"/home/{getpass.getuser()}/cernbox/rare_decays/MuonEG"
MC_path = f"/home/{getpass.getuser()}/cernbox/rare_decays/MC"


def find_folder(path):
    folder = []
    with os.scandir(path) as d1:
        for f1 in d1:
            if f1.name == "log":
                continue
            if f1.is_dir():
                folder.extend(find_folder(f1.path))
            else:
                if f1.name.find("_hist") > 0:
                    continue
                if not f1.name.endswith(".root"):
                    continue
                cut = f1.path.rfind("/")
                folder.append(f1.path[:cut])
    return set(folder)


def create_config(folder, is_MC):
    datasets = []
    for f in folder:
        splited_f = f.split("/")
        for sf in splited_f:
            if sf.find("Run") > 0:
                cut = sf.find("Run")
                name = sf[: cut + 8]
                if sf.endswith("_v2-v2"):
                    name += "-v2"
            elif sf.find("_MiniAOD") > 0:
                cut = sf.find("_MiniAOD")
                name = sf[:cut]
        dataset = {"name": name, "path": f, "is_mc": is_MC}
        datasets.append(dataset)
    return datasets


if __name__ == "__main__":
    datasets = []
    data_folder = find_folder(data_path)
    MC_folder = find_folder(MC_path)
    datasets.extend(create_config(data_folder, False))
    datasets.extend(create_config(MC_folder, True))

    for d in datasets:
        subprocess.run(
            [
                "python",
                "run_analysis.py",
                d["path"],
                "--parallel",
                "--workers",
                "10",
                "--dataset",
                d["name"],
            ],
            check=True,
        )

    # pprint(datasets)

    # python run_analysis.py ~/cernbox/rare_decays/MuonEG/MuonEG_Run2024G-MINIv6NANOv15-v3/260501_045551/0000 --parallel --workers 8 --dataset MuonEG2024G
