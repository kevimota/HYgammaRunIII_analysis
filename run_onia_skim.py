"""run_onia_skim.py: run OniaSkimProcessor over ROOT files.

Processes only data with OniaSkimProcessor — selects on trigger, good muons,
and onia mass window (8-12 GeV). Keeps all onia candidates per event
(no best-candidate picking) and saves them as a flat ntuple.

Usage
-----
    python run_onia_skim.py path_1 [path_1 ...] --dataset NAME
    python run_onia_skim.py --parallel path_1 ... --workers 8

Output
------
    output/onia_skim/{dataset}.root   : flat TTree of onia + muon columns
"""

from pprint import pprint
import os

import typer
from coffea.processor import IterativeExecutor, FuturesExecutor, Runner, column_accumulator
from schema import OniaNanoSchema
from processor import OniaSkimProcessor
from utils import get_files
import awkward as ak
import uproot

app = typer.Typer()
OniaNanoSchema.warn_missing_crossrefs = False


def _run_onia_skim(
    paths: list[str],
    dataset: str,
    parallel: bool = False,
    workers: int = 4,
    chunksize: int = 100_000,
    maxchunks: int | None = None,
    treename: str = "Events",
) -> tuple[dict, dict]:
    fileset = {
        dataset: {
            "files": get_files(paths, exclude="hist"),
            "treename": treename,
            "metadata": {"is_mc": False},
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
        processor_instance=OniaSkimProcessor(),
    )

    return result, metrics


def _save_onia_skim_output(result: dict, metrics: dict, dataset: str):
    print("\n--- Metrics ---")
    pprint(metrics)

    n_raw = result["n_raw"]["count"]

    cols = {}
    for k in result:
        if not isinstance(result[k], column_accumulator):
            continue
        cols[k] = ak.Array(result[k].value)

    n_sel = len(next(iter(cols.values()))) if cols else 0

    out_dir = "output/onia_skim"
    os.makedirs(out_dir, exist_ok=True)
    fname = f"{out_dir}/{dataset}.root"

    with uproot.recreate(fname) as f:
        if n_sel > 0:
            f["events"] = cols

    print(f"\n--- Summary ---")
    print(f"  Raw onia candidates:  {n_raw:,}")
    print(f"  Selected candidates:  {n_sel:,}")
    print(f"\n  Saved {fname}")


@app.command()
def onia_skim(
    paths: list[str] = typer.Argument(..., help="Input ROOT files or directories"),
    dataset: str = typer.Option(..., "--dataset", "-d", help="dataset name"),
    parallel: bool = typer.Option(False, "--parallel", "-p", help="use parallel executor"),
    workers: int = typer.Option(4, "--workers", "-w", help="number of workers"),
    chunksize: int = typer.Option(100_000, "--chunksize", help="entries per chunk"),
    maxchunks: int | None = typer.Option(None, "--maxchunks", help="max chunks to process"),
    treename: str = typer.Option("Events", "--treename", help="TTree name"),
):
    result, metrics = _run_onia_skim(
        paths=paths,
        dataset=dataset,
        parallel=parallel,
        workers=workers,
        chunksize=chunksize,
        maxchunks=maxchunks,
        treename=treename,
    )
    _save_onia_skim_output(result, metrics, dataset)


if __name__ == "__main__":
    app()
