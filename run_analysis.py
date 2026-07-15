"""run_analysis.py: run McDataCompProcessor over ROOT files.

Usage
-----
    python -m run_analysis path_1 [path_1 ...]
    python -m run_analysis --parallel path_1 ...
    python run_analysis.py path_1 ...   (repo root must be on PYTHONPATH)

Output
-----
    output/mc_data_comp/dataset.pkl   : pickled accumulator (histograms + cutflow)
"""

from pprint import pprint
import argparse
import os

from coffea.processor import IterativeExecutor, FuturesExecutor, Runner
from coffea.util import save
from schema import OniaNanoSchema
from processor import McDataCompProcessor
from utils import get_files
import uproot

import awkward as ak
import numpy as np

OniaNanoSchema.warn_missing_crossrefs = False


def main():
    parser = argparse.ArgumentParser(
        description="Run McDataCompProcessor over NanoAOD-like ROOT files."
    )
    parser.add_argument(
        "paths",
        nargs="+",
        help="Input ROOT files (or file patterns supported by uproot).",
    )
    parser.add_argument("--dataset", type=str, help="dataset name")
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Use FuturesExecutor for parallel processing (default: IterativeExecutor).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of workers for Futures Executor (default: 4)",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=100_000,
        help="Maximum entries per chunk (default: 100_000).",
    )
    parser.add_argument(
        "--maxchunks",
        type=int,
        default=None,
        help="Maximum number of chunks to process per dataset (default: all).",
    )
    parser.add_argument(
        "--treename",
        default="Events",
        help="TTree name in the ROOT file (default: Events).",
    )
    parser.add_argument("--is_mc", action="store_true", help="Dataset is MC")
    args = parser.parse_args()

    # -----------------------------------------------------------------
    # Build fileset: dict of dataset -> [file list]
    # Here we treat all input files as belonging to a single dataset.
    # -----------------------------------------------------------------
    fileset = {
        args.dataset: {
            "files": get_files(args.paths, exclude="hist"),
            "treename": args.treename,
            "metadata": {"is_mc": args.is_mc},
        }
    }

    # -----------------------------------------------------------------
    # Executor selection
    #
    # IterativeExecutor (default): sequential, single-threaded.
    #   Good for debugging; use for quick local tests.
    #
    # FuturesExecutor: parallel via Python's concurrent.futures.
    #   Use when you have multiple cores and files.
    #   To switch: executor = FuturesExecutor(compression=1)
    # -----------------------------------------------------------------
    if args.parallel:
        executor = FuturesExecutor(workers=args.workers)
    else:
        executor = IterativeExecutor()

    # -----------------------------------------------------------------
    # Runner: orchestrates chunking, schema application, and accumulation.
    # -----------------------------------------------------------------
    runner = Runner(
        executor=executor,
        chunksize=args.chunksize,
        maxchunks=args.maxchunks,
        schema=OniaNanoSchema,
        format="root",
        savemetrics=True,
    )

    print(
        f"\nProcessing {len(fileset[args.dataset]["files"])} file(s) with chunksize={args.chunksize}"
    )
    print(f"Executor: {executor.__class__.__name__}")

    result, metrics = runner(
        fileset=fileset,
        processor_instance=McDataCompProcessor(),
    )

    print("\n--- Metrics ---")
    pprint(metrics)

    # -----------------------------------------------------------------
    # Print cutflow summary
    # -----------------------------------------------------------------
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
    print(
        f"  Muons (good):         {good_muons_sel:,}  ({100*good_muons_sel/n_events:.1f}%)"
    )
    print(
        f"  Photon (good):        {good_photon_sel:,}  ({100*good_photon_sel/n_events:.1f}%)"
    )
    print(f"  Upsilon mass:         {upsilon_sel:,}  ({100*upsilon_sel/n_events:.1f}%)")
    print(f"  X mass:               {x_sel:,}  ({100*x_sel/n_events:.1f}%)")

    # -----------------------------------------------------------------
    # Saving files
    # -----------------------------------------------------------------

    out_dir = "output/mc_data_comp"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    result = {**result, "is_mc": args.is_mc, "dataset": args.dataset}
    save(result, f"{out_dir}/{args.dataset}.pkl")

    from coffea.processor import column_accumulator

    events = {}

    for k in result:
        if not isinstance(result[k], column_accumulator):
            continue
        events[k] = ak.Array(result[k].value)

    events = ak.zip(events)

    with uproot.recreate(f"{out_dir}/{args.dataset}.root") as f:
        f["events"] = events
        f["metadata"] = ak.Array([{"is_mc": args.is_mc, "dataset": args.dataset}])
        for grp in result["hists"]:
            for var in result["hists"][grp]:
                for cat in result["hists"][grp][var].axes["cat"]:
                    path = f"{grp}/{var}/{cat}"
                    f[path] = result["hists"][grp][var][:, cat]

    print(f"\nOutput saved to {out_dir}/{args.dataset}.pkl and {out_dir}/{args.dataset}.root")


if __name__ == "__main__":
    main()
