"""run_skim.py: run SkimProcessor over ROOT files and save flat ntuples.

Produces two ROOT files per dataset:
  - output/skim/{dataset}_presel.root   (trigger + good muons + good photon)
  - output/skim/{dataset}_fullsel.root  (preselection + mass windows)

Usage
-----
    python run_skim.py skim path_1 [path_1 ...] --dataset NAME [--is-mc]
    python run_skim.py skim --parallel path_1 ... --workers 8

Output
------
    output/skim/{dataset}_presel.root
    output/skim/{dataset}_fullsel.root
"""

from pprint import pprint
import os

import typer
from coffea.processor import IterativeExecutor, FuturesExecutor, Runner, column_accumulator
from schema import OniaNanoSchema
from processor import SkimProcessor
from utils import get_files
from utils.weights import get_xsec_br, lumi_2024
import awkward as ak
import uproot

app = typer.Typer()
OniaNanoSchema.warn_missing_crossrefs = False


def write_stage(out_dir, dataset, stage, columns, n_raw, is_mc, weight_val):
    """Write one selection stage to a ROOT file."""
    fname = f"{out_dir}/{dataset}_{stage}.root"
    n_sel = len(next(iter(columns.values()))) if columns else 0

    with uproot.recreate(fname) as f:
        if n_sel > 0:
            f["events"] = columns
        f["metadata"] = ak.Array([{
            "dataset": dataset,
            "is_mc": is_mc,
            "n_raw": n_raw,
            "n_sel": n_sel,
            "weight": weight_val,
        }])
    return fname, n_sel


def _run_skim(
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
        processor_instance=SkimProcessor(),
    )

    return result, metrics


def _save_skim_output(result: dict, metrics: dict, dataset: str, is_mc: bool):
    print("\n--- Metrics ---")
    pprint(metrics)

    n_raw = result["n_raw"]["count"]

    raw = {}
    for k in result:
        if not isinstance(result[k], column_accumulator):
            continue
        raw[k] = ak.Array(result[k].value)

    presel_cols = {k[len("presel_"):]: v for k, v in raw.items() if k.startswith("presel_")}
    fullsel_cols = {k[len("fullsel_"):]: v for k, v in raw.items() if k.startswith("fullsel_")}

    xsec = get_xsec_br(dataset)
    if xsec > 0 and n_raw > 0:
        weight_val = xsec * lumi_2024 / n_raw
    else:
        weight_val = 1.0

    out_dir = "output/skim"
    os.makedirs(out_dir, exist_ok=True)

    fname_p, n_p = write_stage(out_dir, dataset, "presel", presel_cols, n_raw, is_mc, weight_val)
    fname_f, n_f = write_stage(out_dir, dataset, "fullsel", fullsel_cols, n_raw, is_mc, weight_val)

    print(f"\n--- Summary ---")
    print(f"  Raw events:           {n_raw:,}")
    print(f"  Preselected events:   {n_p:,}")
    print(f"  Full selected events: {n_f:,}")
    print(f"  Weight per event:     {weight_val:.6e}")
    print(f"\n  Saved {fname_p}")
    print(f"  Saved {fname_f}")


@app.command()
def skim(
    paths: list[str] = typer.Argument(..., help="Input ROOT files or directories"),
    dataset: str = typer.Option(..., "--dataset", "-d", help="dataset name"),
    parallel: bool = typer.Option(False, "--parallel", "-p", help="use parallel executor"),
    workers: int = typer.Option(4, "--workers", "-w", help="number of workers"),
    chunksize: int = typer.Option(100_000, "--chunksize", help="entries per chunk"),
    maxchunks: int | None = typer.Option(None, "--maxchunks", help="max chunks to process"),
    treename: str = typer.Option("Events", "--treename", help="TTree name"),
    is_mc: bool = typer.Option(False, "--is-mc", help="dataset is MC"),
):
    result, metrics = _run_skim(
        paths=paths,
        dataset=dataset,
        parallel=parallel,
        workers=workers,
        chunksize=chunksize,
        maxchunks=maxchunks,
        treename=treename,
        is_mc=is_mc,
    )
    _save_skim_output(result, metrics, dataset, is_mc)


if __name__ == "__main__":
    app()
