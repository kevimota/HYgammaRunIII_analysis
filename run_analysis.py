"""run_analysis.py: run McDataCompProcessor over ROOT files.

Usage
-----
    python run_analysis.py compare path_1 [path_1 ...] --dataset NAME [--is-mc]
    python run_analysis.py compare --parallel path_1 ... --workers 8

Output
-----
    output/mc_data_comp/dataset.pkl   : pickled accumulator (histograms + cutflow)
"""

from pprint import pprint
import os

import typer
from coffea.processor import IterativeExecutor, FuturesExecutor, Runner
from coffea.util import save
from schema import OniaNanoSchema
from processor import McDataCompProcessor
from utils import get_files
import uproot
import awkward as ak
import numpy as np

app = typer.Typer()
OniaNanoSchema.warn_missing_crossrefs = False


def _run_compare(
    paths: list[str],
    dataset: str,
    parallel: bool = False,
    workers: int = 4,
    chunksize: int = 100_000,
    maxchunks: int | None = None,
    treename: str = "Events",
    is_mc: bool = False,
) -> tuple[dict, dict]:
    fileset = {
        dataset: {
            "files": get_files(paths, exclude="hist"),
            "treename": treename,
            "metadata": {"is_mc": is_mc},
        }
    }

    if parallel:
        executor = FuturesExecutor(workers=workers)
    else:
        executor = IterativeExecutor()

    runner = Runner(
        executor=executor,
        chunksize=chunksize,
        maxchunks=maxchunks,
        schema=OniaNanoSchema,
        format="root",
        savemetrics=True,
    )

    n_files = len(fileset[dataset]["files"])
    print(f"\nProcessing {n_files} file(s) with chunksize={chunksize}")
    print(f"Executor: {executor.__class__.__name__}")

    result, metrics = runner(
        fileset=fileset,
        processor_instance=McDataCompProcessor(),
    )

    return result, metrics


def _save_output(result: dict, metrics: dict, dataset: str, is_mc: bool):
    print("\n--- Metrics ---")
    pprint(metrics)

    cutflow = result["cutflow"]
    n_events = cutflow.get("n_events", 0)
    trigger_sel = cutflow.get("trigger_sel", 0)
    good_muons_sel = cutflow.get("good_muons_sel", 0)
    good_photon_sel = cutflow.get("flow_good_photon_sel", 0)
    upsilon_sel = cutflow.get("flow_upsilon_sel", 0)
    x_sel = cutflow.get("flow_x_sel", 0)

    print("\n--- Cutflow ---")
    print(f"  Total events:         {n_events:,}")
    print(f"  Trigger selection:    {trigger_sel:,}  ({100*trigger_sel/n_events:.1f}%)")
    print(f"  Muons (good):         {good_muons_sel:,}  ({100*good_muons_sel/n_events:.1f}%)")
    print(f"  Photon (good):        {good_photon_sel:,}  ({100*good_photon_sel/n_events:.1f}%)")
    print(f"  Upsilon mass:         {upsilon_sel:,}  ({100*upsilon_sel/n_events:.1f}%)")
    print(f"  X mass:               {x_sel:,}  ({100*x_sel/n_events:.1f}%)")

    out_dir = "output/mc_data_comp"
    os.makedirs(out_dir, exist_ok=True)
    result = {**result, "is_mc": is_mc, "dataset": dataset}
    save(result, f"{out_dir}/{dataset}.pkl")

    from coffea.processor import column_accumulator

    events = {}
    for k in result:
        if not isinstance(result[k], column_accumulator):
            continue
        events[k] = ak.Array(result[k].value)
    events = ak.zip(events)

    with uproot.recreate(f"{out_dir}/{dataset}.root") as f:
        f["events"] = events
        f["metadata"] = ak.Array([{"is_mc": is_mc, "dataset": dataset}])
        for grp in result["hists"]:
            for var in result["hists"][grp]:
                for cat in result["hists"][grp][var].axes["cat"]:
                    path = f"{grp}/{var}/{cat}"
                    f[path] = result["hists"][grp][var][:, cat]

    print(f"\nOutput saved to {out_dir}/{dataset}.pkl and {out_dir}/{dataset}.root")


@app.command()
def compare(
    paths: list[str] = typer.Argument(..., help="Input ROOT files or directories"),
    dataset: str = typer.Option(..., "--dataset", "-d", help="dataset name"),
    parallel: bool = typer.Option(False, "--parallel", "-p", help="use parallel executor"),
    workers: int = typer.Option(4, "--workers", "-w", help="number of workers"),
    chunksize: int = typer.Option(100_000, "--chunksize", help="entries per chunk"),
    maxchunks: int | None = typer.Option(None, "--maxchunks", help="max chunks to process"),
    treename: str = typer.Option("Events", "--treename", help="TTree name"),
    is_mc: bool = typer.Option(False, "--is-mc", help="dataset is MC"),
):
    result, metrics = _run_compare(
        paths=paths,
        dataset=dataset,
        parallel=parallel,
        workers=workers,
        chunksize=chunksize,
        maxchunks=maxchunks,
        treename=treename,
        is_mc=is_mc,
    )
    _save_output(result, metrics, dataset, is_mc)


if __name__ == "__main__":
    app()
