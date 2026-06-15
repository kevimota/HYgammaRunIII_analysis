import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import sys, os
from copy import deepcopy
from coffea.util import load

# Allow running as `python utils/figure.py` from the project root
if __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from utils.weights import get_xsec_br, lumi_2024
from concurrent.futures import ProcessPoolExecutor, as_completed

plt.style.use(mh.style.CMS)
plt.rcParams.update(
    {
        "font.size": 16,
        "axes.titlesize": 18,
        "axes.labelsize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
    }
)


def sum_group(files):
    """Load, weight, and merge histograms from all pickle files in a dataset group.

    MC histograms are scaled by xsec*BR*lumi/n_events to produce absolute
    cross-section predictions at the data luminosity.
    Data files are left unscaled (raw event counts).
    """
    summed = None
    for f in files:
        data = load(f)
        this = data["hists"]
        n_events = data["cutflow"]["n_events"]

        # Determine per-event weight
        # MC: weight = xsec_br * lumi / n_events  (lumi_2024 in pb^-1)
        # Data / unknown: no weight (scale = 1)
        xsec_br = get_xsec_br(f)
        if xsec_br > 0 and n_events > 0:
            scale = xsec_br * lumi_2024 / n_events
        else:
            scale = 1.0

        if summed is None:
            summed = deepcopy(this)
            if scale != 1.0:
                for grp in summed:
                    for var in summed[grp]:
                        summed[grp][var] *= scale
        else:
            if scale != 1.0:
                for grp in summed:
                    for var in summed[grp]:
                        summed[grp][var] += this[grp][var] * scale
            else:
                for grp in summed:
                    for var in summed[grp]:
                        summed[grp][var] += this[grp][var]
    return summed


def sum_cutflow(files):
    """Sum cutflows across files with the same MC scaling as sum_group."""
    summed = None
    for f in files:
        data = load(f)
        cf = dict(data["cutflow"])
        n_events = cf["n_events"]

        xsec_br = get_xsec_br(f)
        if xsec_br > 0 and n_events > 0:
            scale = xsec_br * lumi_2024 / n_events
        else:
            scale = 1.0

        if summed is None:
            summed = {k: v * scale for k, v in cf.items()}
        else:
            for k in summed:
                summed[k] += cf[k] * scale
    return summed


def plot_cutflow_bar(cutflow, stages, ax=None, label=None, color=None, log=True):
    """Draw a cutflow bar chart onto an axes."""
    if ax is None:
        ax = plt.gca()
    vals = np.array([cutflow[s] for s in stages], dtype=float)
    n = len(stages)
    labels = [s.replace("_", " ").replace("flow ", "").title() for s in stages]
    x = np.arange(n)
    ax.bar(x, vals, color=color, label=label, width=0.6, edgecolor="gray", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=11)
    ax.set_ylabel("Events")
    if log:
        ax.set_yscale("log")
        lo = max(0.1, vals.min() * 0.5)
        hi = vals.max() * 5
        if lo >= hi:
            lo = hi / 100
        ax.set_ylim(lo, hi)
    else:
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 3), useMathText=True)
    for i, v in enumerate(vals):
        ax.text(
            i, v * 1.3, f"{v:.2e}", ha="center", va="bottom", fontsize=9, rotation=90
        )
    return ax


cutflow_stages = [
    "n_events",
    "trigger_sel",
    "good_muons_sel",
    "flow_good_photon_sel",
    "flow_upsilon_sel",
    "flow_x_sel",
]


def blind_hist(h, windows):
    """Zero out bin contents of a 2D histogram in specified mass windows.

    Parameters
    ----------
    h : hist.Hist
        Histogram with a 'mass' axis (axis 0) and a 'cat' axis (axis 1, StrCategory).
    windows : list of (lo, hi)
        Mass range(s) to blind, in GeV.
    """
    from hist import Hist

    mass_axis = h.axes["mass"]
    view = h.view(flow=False)
    for lo, hi in windows:
        i_lo = max(mass_axis.index(lo), 0)
        i_hi = min(mass_axis.index(hi) + 1, view.shape[0])
        view[i_lo:i_hi, :] = 0.0
    return h


def create_plot1d(hist1d, labels=None, log=False, ax=None, lumi=None, **kwargs):
    """Plot a 1D histogram with optional stats box."""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

    if ax is None:
        ax = plt.gca()

    artists = mh.histplot(hist1d, ax=ax, **kwargs)
    stairs_artists = [artist.stairs for artist in artists]
    if labels is not None:
        if len(labels) != len(stairs_artists):
            print("len of labels does not match artists")
        else:
            ax.legend(stairs_artists, labels)

    if lumi is not None:
        ax.text(
            1.0,
            1.0,
            f"{lumi:.2f}" + r" fb$^{-1}$ (13.6 TeV)",
            fontsize=18,
            horizontalalignment="right",
            verticalalignment="bottom",
            transform=ax.transAxes,
        )
    if log:
        ax.set_yscale("log")
    else:
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 3), useMathText=True)

    values = hist1d.values()
    if len(values) > 0:
        centers = hist1d.axes.centers[0]
        mean = np.sum(values * centers) / np.sum(values)
        std = np.sqrt(np.sum(values * ((centers - mean) ** 2)) / np.sum(values))
        annotation = TextArea(
            f"Total: {np.sum(values):.2e}\nMean: {mean:.2e}\nStd: {std:.2e}",
            textprops=dict(size=14),
        )
        at = AnchoredOffsetbox("upper right", child=annotation)
        at.patch.set_facecolor("None")
        ax.add_artist(at)
    return ax


# ---------------------------------------------------------------
# Dataset definitions
# ---------------------------------------------------------------
cats = ["all", "trigger", "good_muons", "good_photon", "upsilon_mass", "x_mass"]

mc_datasets = ["ggH", "VBF", "WH", "ZH", "bbH", "ttH", "qqZ", "HtoMuMuG", "ZGTo2MuG"]

datasets = {
    "Data": [
        "output/MuonEG_Run2024C.pkl",
        "output/MuonEG_Run2024D.pkl",
        "output/MuonEG_Run2024E.pkl",
        "output/MuonEG_Run2024F.pkl",
        "output/MuonEG_Run2024G.pkl",
        "output/MuonEG_Run2024H.pkl",
        "output/MuonEG_Run2024I.pkl",
        "output/MuonEG_Run2024I-v2.pkl",
    ],
    "ggH": [
        "output/GluGluToH_HToUps1SG_Ups1SToMuMu.pkl",
        "output/GluGluToH_HToUps2SG_Ups2SToMuMu.pkl",
        "output/GluGluToH_HToUps3SG_Ups3SToMuMu.pkl",
    ],
    "VBF": [
        "output/VBFToH_HToUps1SG_Ups1SToMuMu.pkl",
        "output/VBFToH_HToUps2SG_Ups2SToMuMu.pkl",
        "output/VBFToH_HToUps3SG_Ups3SToMuMu.pkl",
    ],
    "WH": [
        "output/WH_HToUps1SG_Ups1SToMuMu.pkl",
        "output/WH_HToUps2SG_Ups2SToMuMu.pkl",
        "output/WH_HToUps3SG_Ups3SToMuMu.pkl",
    ],
    "ZH": [
        "output/ZH_HToUps1SG_Ups1SToMuMu.pkl",
        "output/ZH_HToUps2SG_Ups2SToMuMu.pkl",
        "output/ZH_HToUps3SG_Ups3SToMuMu.pkl",
    ],
    "bbH": [
        "output/bbH_HToUps1SG_Ups1SToMuMu.pkl",
        "output/bbH_HToUps2SG_Ups2SToMuMu.pkl",
        "output/bbH_HToUps3SG_Ups3SToMuMu.pkl",
    ],
    "ttH": [
        "output/ttH_HToUps1SG_Ups1SToMuMu.pkl",
        "output/ttH_HToUps2SG_Ups2SToMuMu.pkl",
        "output/ttH_HToUps3SG_Ups3SToMuMu.pkl",
    ],
    "qqZ": [
        "output/ZToUps1SG_Ups1SToMuMu.pkl",
        "output/ZToUps2SG_Ups2SToMuMu.pkl",
        "output/ZToUps3SG_Ups3SToMuMu.pkl",
    ],
    "HtoMuMuG": ["output/ggH125_012j_NLO_FXFX_HtoMuMuGamma.pkl"],
    "ZGTo2MuG": ["output/ZGTo2MuG_mll_2to15_LO.pkl"],
}

# Histogram groups and their variables
hist_groups = {
    "muon_lead": ["pt", "eta", "phi"],
    "muon_trail": ["pt", "eta", "phi"],
    "photon": ["pt", "eta", "phi"],
    "onia": ["mass", "pt", "rap", "phi"],
    "x": ["mass", "pt", "rap", "phi", "n"],
}

lumi = lumi_2024 * 1e-3

# Module-level dict shared with worker processes (set via ProcessPoolExecutor initializer)
_group_data = None


def _init_worker(group_data_):
    global _group_data
    _group_data = group_data_


# --- Parallel plot helpers ------------------------------------------------


def _save_combined_plot(grp, var, cat, mc_datasets, dirname, lumi):
    fig, ax = plt.subplots(figsize=(8, 6))

    mc_hists = []
    mc_labels = []
    for mc_name in mc_datasets:
        if mc_name not in _group_data:
            continue
        h = _group_data[mc_name][grp][var][:, cat]
        mc_hists.append(h)
        mc_labels.append(mc_name)

    if mc_hists:
        mc_colors = plt.cm.tab20(np.linspace(0, 1, len(mc_hists)))
        mh.histplot(
            mc_hists,
            ax=ax,
            label=mc_labels,
            histtype="fill",
            stack=True,
            color=mc_colors,
        )

    if "Data" in _group_data:
        h_data = _group_data["Data"][grp][var][:, cat]
        mh.histplot(
            h_data,
            ax=ax,
            label="Data",
            histtype="errorbar",
            color="black",
            marker="o",
            markersize=6,
        )

    ax.legend(fontsize=12)
    log = var in ("pt", "n")
    if log:
        ax.set_yscale("log")
    else:
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 3), useMathText=True)
    ax.text(
        1.0,
        1.0,
        f"{lumi:.2f} fb$^{{-1}}$ (13.6 TeV)",
        fontsize=18,
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )

    fname = f"{dirname}/{cat}/{grp}_{var}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fname


def _save_per_dataset_plot(ds_name, ds_data, grp, var, cat, dirname, lumi):
    fig, ax = plt.subplots(figsize=(8, 6))
    h = ds_data[grp][var][:, cat]
    create_plot1d(h, ax=ax, lumi=lumi)

    log = var in ("pt", "n")
    if log:
        ax.set_yscale("log")
    else:
        ax.ticklabel_format(
            axis="y",
            style="sci",
            scilimits=(0, 3),
            useMathText=True,
        )

    fname = f"{dirname}/{cat}/{grp}_{var}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fname


def _save_comparison_plot(ds_name, ds_data, grp, var, dirname, lumi):
    fig, ax = plt.subplots(figsize=(8, 6))

    h_all = ds_data[grp][var][:, "all"]
    h_pre = ds_data[grp][var][:, "good_photon"]
    h_sel = ds_data[grp][var][:, "x_mass"]

    mh.histplot(
        [h_all, h_pre, h_sel],
        ax=ax,
        label=["Directly from NTuple", "After preselection", "After selection"],
        histtype="step",
        color=["#1f77b4", "#d62728", "#60c71f"],
        linewidth=2,
    )

    ax.legend(fontsize=12)
    log = var in ("pt", "n")
    if log:
        ax.set_yscale("log")
    else:
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 3), useMathText=True)
    ax.text(
        1.0,
        1.0,
        f"{lumi:.2f} fb$^{{-1}}$ (13.6 TeV)",
        fontsize=18,
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )

    fname = f"{dirname}/{grp}_{var}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fname


def _save_cutflow_plot(ds_name, cf, dirname, lumi, stages):
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_cutflow_bar(cf, stages, ax=ax, label=ds_name, log=True)
    ax.set_title(ds_name)
    ax.text(
        1.0,
        1.0,
        f"{lumi:.2f} fb$^{{-1}}$ (13.6 TeV)",
        fontsize=18,
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )
    fname = f"{dirname}/{ds_name}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fname


def _save_combined_cutflow_plot(group_cutflow, mc_datasets, dirname, lumi, stages):
    fig, ax = plt.subplots(figsize=(10, 7))
    active_mc = [m for m in mc_datasets if m in group_cutflow]
    mc_colors = plt.cm.tab20(np.linspace(0, 1, len(active_mc)))
    for i, mc_name in enumerate(active_mc):
        cf = group_cutflow[mc_name]
        vals = [cf[s] for s in stages]
        ax.plot(
            range(len(stages)),
            vals,
            marker="o",
            label=mc_name,
            color=mc_colors[i],
            linewidth=2,
        )
    if "Data" in group_cutflow:
        cf = group_cutflow["Data"]
        vals = [cf[s] for s in stages]
        ax.plot(
            range(len(stages)),
            vals,
            marker="s",
            label="Data",
            color="black",
            linewidth=3,
            linestyle="--",
        )
    labels = [s.replace("_", " ").replace("flow ", "").title() for s in stages]
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=11)
    ax.set_ylabel("Events")
    ax.set_yscale("log")
    ax.legend(fontsize=10, loc="upper right")
    ax.text(
        1.0,
        1.0,
        f"{lumi:.2f} fb$^{{-1}}$ (13.6 TeV)",
        fontsize=18,
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )
    fname = f"{dirname}/combined.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fname


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Make plots from analysis output")
    parser.add_argument(
        "--workers", type=int, default=4, help="Number of worker threads (default: 4)"
    )
    parser.add_argument(
        "--validation", action="store_true", help="create NTuple validation plots"
    )
    parser.add_argument(
        "--save-per-dataset",
        action="store_true",
        help="Save per-dataset validation plots",
    )
    args = parser.parse_args()

    import os

    if not os.path.exists("hists"):
        os.mkdir("hists")

    if args.validation:
        # ---------------------------------------------------------------
        # Load and sum all dataset groups (sequential)
        # ---------------------------------------------------------------
        print("Loading and summing histograms per dataset group...")
        group_data = {}
        for name, files in datasets.items():
            ok = True
            for f in files:
                if not os.path.exists(f):
                    print(f"  WARNING: {f} not found, skipping group '{name}'")
                    ok = False
                    break
            if ok:
                group_data[name] = sum_group(files)
                print(f"  {name}: {len(files)} file(s) loaded")

        # Blind data x mass histograms in Z and H mass windows
        if "Data" in group_data:
            blind_hist(
                group_data["Data"]["x"]["mass"],
                windows=[(80, 100), (120, 130)],
            )
            print("  Blinded data x-mass in Z (80-100 GeV) and H (120-130 GeV) windows")
        print()

        # ---------------------------------------------------------------
        # Create output directories
        # ---------------------------------------------------------------
        combined_dir = "hists/combined"
        comp_dir = "hists/comparison"

        for cat in cats:
            os.makedirs(f"{combined_dir}/{cat}", exist_ok=True)

        if args.save_per_dataset:
            for name in group_data:
                for cat in cats:
                    os.makedirs(f"hists/{name}/{cat}", exist_ok=True)

        for name in group_data:
            os.makedirs(f"{comp_dir}/{name}", exist_ok=True)

        # ---------------------------------------------------------------
        # Combined stacked plots: all MC groups stacked + data as errorbars
        # ---------------------------------------------------------------
        print("Saving stacked plots")
        tasks = []
        for grp, variables in hist_groups.items():
            for var in variables:
                for cat in cats:
                    tasks.append(
                        (
                            _save_combined_plot,
                            grp,
                            var,
                            cat,
                            mc_datasets,
                            combined_dir,
                            lumi,
                        )
                    )

        with ProcessPoolExecutor(
            max_workers=args.workers, initializer=_init_worker, initargs=(group_data,)
        ) as executor:
            futures = [executor.submit(f, *a) for f, *a in tasks]
            for future in as_completed(futures):
                print(f"  Saved {future.result()}")

        # ---------------------------------------------------------------
        # Per-dataset plots for reference (parallel, optional)
        # ---------------------------------------------------------------
        if args.save_per_dataset:
            tasks = []
            for ds_name, ds_data in group_data.items():
                ds_dir = f"hists/{ds_name}"
                for grp, variables in hist_groups.items():
                    for var in variables:
                        for cat in cats:
                            tasks.append(
                                (
                                    _save_per_dataset_plot,
                                    ds_name,
                                    ds_data,
                                    grp,
                                    var,
                                    cat,
                                    ds_dir,
                                    lumi,
                                )
                            )

            print(f"Saving per-dataset plots ({len(tasks)} total)")
            with ProcessPoolExecutor(
                max_workers=args.workers,
                initializer=_init_worker,
                initargs=(group_data,),
            ) as executor:
                futures = [executor.submit(f, *a) for f, *a in tasks]
                for future in as_completed(futures):
                    print(f"  Saved {future.result()}")

        # ---------------------------------------------------------------
        # Preselection vs full-selection comparison
        # ---------------------------------------------------------------
        print("Saving comparison plots")
        tasks = []
        for ds_name, ds_data in group_data.items():
            ds_comp_dir = f"{comp_dir}/{ds_name}"
            for grp, variables in hist_groups.items():
                for var in variables:
                    tasks.append(
                        (
                            _save_comparison_plot,
                            ds_name,
                            ds_data,
                            grp,
                            var,
                            ds_comp_dir,
                            lumi,
                        )
                    )

        with ProcessPoolExecutor(
            max_workers=args.workers, initializer=_init_worker, initargs=(group_data,)
        ) as executor:
            futures = [executor.submit(f, *a) for f, *a in tasks]
            for future in as_completed(futures):
                print(f"  Saved {future.result()}")

        # ---------------------------------------------------------------
        # Cutflow bar charts
        # ---------------------------------------------------------------
        cf_dir = "hists/cutflow"
        os.makedirs(cf_dir, exist_ok=True)

        # Load cutflows for all dataset groups (sequential)
        print("Loading cutflows...")
        group_cutflow = {}
        for name, files in datasets.items():
            ok = True
            for f in files:
                if not os.path.exists(f):
                    ok = False
                    break
            if ok:
                group_cutflow[name] = sum_cutflow(files)
                print(f"  {name}: cutflow loaded")

        # Per-dataset cutflow bar charts (parallel)
        tasks = []
        for ds_name, cf in group_cutflow.items():
            tasks.append(
                (
                    _save_cutflow_plot,
                    ds_name,
                    cf,
                    cf_dir,
                    lumi,
                    cutflow_stages,
                )
            )

        with ProcessPoolExecutor(
            max_workers=args.workers, initializer=_init_worker, initargs=(group_data,)
        ) as executor:
            futures = [executor.submit(f, *a) for f, *a in tasks]
            for future in as_completed(futures):
                print(f"  Saved {future.result()}")

        # Combined cutflow: all MC lines + data (single plot)
        fname = _save_combined_cutflow_plot(
            group_cutflow,
            mc_datasets,
            cf_dir,
            lumi,
            cutflow_stages,
        )
        print(f"  Saved {fname}")
