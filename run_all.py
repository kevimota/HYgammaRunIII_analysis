import os
import getpass
from collections import defaultdict

import typer
from utils import get_files

app = typer.Typer()

DEFAULT_DATA_PATH = f"/home/{getpass.getuser()}/cernbox/rare_decays/MuonEG"
DEFAULT_MC_PATH = f"/home/{getpass.getuser()}/cernbox/rare_decays/MC"


def find_folder(path):
    """Recursively find leaf directories containing .root files (excluding _hist)."""
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


def extract_name(folder_path, is_mc):
    """Extract dataset name from a leaf folder path."""
    parts = folder_path.split("/")
    for part in parts:
        if part.find("Run") > 0:
            cut = part.find("Run")
            name = part[: cut + 8]
            if part.endswith("_v2-v2"):
                name += "-v2"
            return name
        elif part.find("_MiniAOD") > 0:
            cut = part.find("_MiniAOD")
            return part[:cut]
    return folder_path.split("/")[-1]


def collect_mc_by_name(folders):
    """Group MC leaf folders by dataset name."""
    groups = defaultdict(list)
    for f in folders:
        name = extract_name(f, is_mc=True)
        groups[name].append(f)
    return groups


@app.command()
def compare(
    data_path: str = typer.Option(DEFAULT_DATA_PATH, "--data-path", help="Data directory"),
    mc_path: str = typer.Option(DEFAULT_MC_PATH, "--mc-path", help="MC directory"),
    workers: int = typer.Option(12, "--workers", "-w", help="Workers per executor"),
    parallel: bool = typer.Option(True, "--parallel/--serial", help="Use parallel executor"),
):
    from run_analysis import _run_compare, _save_output

    data_folders = find_folder(data_path)
    mc_folders = find_folder(mc_path)

    for folder in data_folders:
        name = extract_name(folder, is_mc=False)
        print(f"\n{'='*60}")
        print(f"Processing data: {name}")
        print(f"  {folder}")
        result, metrics = _run_compare(
            paths=[folder],
            dataset=name,
            parallel=parallel,
            workers=workers,
            is_mc=False,
        )
        _save_output(result, metrics, name, is_mc=False)

    mc_groups = collect_mc_by_name(mc_folders)
    for name, folders in mc_groups.items():
        for folder in folders:
            print(f"\n{'='*60}")
            print(f"Processing MC: {name}")
            print(f"  {folder}")
            result, metrics = _run_compare(
                paths=[folder],
                dataset=name,
                parallel=parallel,
                workers=workers,
                is_mc=True,
            )
            _save_output(result, metrics, name, is_mc=True)


@app.command()
def skim(
    data_path: str = typer.Option(DEFAULT_DATA_PATH, "--data-path", help="Data directory"),
    mc_path: str = typer.Option(DEFAULT_MC_PATH, "--mc-path", help="MC directory"),
    workers: int = typer.Option(12, "--workers", "-w", help="Workers per executor"),
    parallel: bool = typer.Option(True, "--parallel/--serial", help="Use parallel executor"),
):
    from run_skim import _run_skim, _save_skim_output

    # --- Data: single ntuple for all 2024 data ---
    data_folders = find_folder(data_path)
    all_data_files = []
    for folder in data_folders:
        all_data_files.extend(get_files([folder], exclude="hist"))
    all_data_files.sort()

    print(f"\n{'='*60}")
    print(f"Skimming data: MuonEG_Run2024")
    print(f"  {len(data_folders)} leaf dirs, {len(all_data_files)} files total")
    result, metrics = _run_skim(
        paths=data_folders,  # pass leaf dirs; _run_skim calls get_files internally
        dataset="MuonEG_Run2024",
        parallel=parallel,
        workers=workers,
        is_mc=False,
    )
    _save_skim_output(result, metrics, "MuonEG_Run2024", is_mc=False)

    # --- MC: one skim per MC sample ---
    mc_folders = find_folder(mc_path)
    mc_groups = collect_mc_by_name(mc_folders)

    for name, folders in mc_groups.items():
        print(f"\n{'='*60}")
        print(f"Skimming MC: {name}")
        print(f"  {len(folders)} leaf dir(s)")
        result, metrics = _run_skim(
            paths=folders,
            dataset=name,
            parallel=parallel,
            workers=workers,
            is_mc=True,
        )
        _save_skim_output(result, metrics, name, is_mc=True)


if __name__ == "__main__":
    app()
