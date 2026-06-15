"""TestProcessor: coffea processor for testing purposes.

Output
------
    Accumulator containing histograms and cutflow counters.
"""

from coffea.processor import ProcessorABC, defaultdict_accumulator
import hist
import awkward as ak
import numpy as np
from utils import remove_none, fill_kin_hists


def fill_kin_hists(obj, hists, cat=None):
    for h in hists:
        if cat:
            hists[h].fill(obj[h], cat=cat)
        else:
            hists[h].fill(obj[h])


class TestProcessor(ProcessorABC):

    def __init__(self):
        """Initialize the processor."""

    def process(self, events):
        dataset = events.metadata["dataset"]

        # --- Accumulators definition --------------------------------------------

        muon_hists = {
            "pt": hist.Hist(
                hist.axis.Regular(100, 0, 200, name="pt", label=r"p$_{T}^{\mu}$ [GeV]"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "eta": hist.Hist(
                hist.axis.Regular(100, -3, 3, name="eta", label=r"$\eta_{\mu}$"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "phi": hist.Hist(
                hist.axis.Regular(100, -3.5, 3.5, name="phi", label=r"$\phi_{\mu}$"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
        }
        photon_hists = {
            "pt": hist.Hist(
                hist.axis.Regular(
                    100, 0, 200, name="pt", label=r"E$_{T}^{\gamma}$ [GeV]"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "eta": hist.Hist(
                hist.axis.Regular(100, -3, 3, name="eta", label=r"$\eta_{\gamma}$"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "phi": hist.Hist(
                hist.axis.Regular(
                    100, -3.5, 3.5, name="phi", label=r"$\phi_{\gamma}]$"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
        }
        onia_hists = {
            "mass": hist.Hist(
                hist.axis.Regular(50, 8, 12, name="mass", label=r"m$_{\mu\mu}$ [GeV]"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "pt": hist.Hist(
                hist.axis.Regular(
                    100, 0, 200, name="pt", label=r"p$_{T}^{\mu\mu}$ [GeV]"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "rap": hist.Hist(
                hist.axis.Regular(100, -3, 3, name="rap", label=r"y$_{\mu\mu}$"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "phi": hist.Hist(
                hist.axis.Regular(100, -3.5, 3.5, name="phi", label=r"$\phi_{\mu\mu}$"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
        }
        konia_hists = {
            "mass": hist.Hist(
                hist.axis.Regular(50, 8, 12, name="mass", label=r"m$_{k\mu\mu}$ [GeV]"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "pt": hist.Hist(
                hist.axis.Regular(
                    100, 0, 200, name="pt", label=r"p$_{T}^{k\mu\mu}$ [GeV]"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "rap": hist.Hist(
                hist.axis.Regular(100, -3, 3, name="rap", label=r"y$_{k\mu\mu}$"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "phi": hist.Hist(
                hist.axis.Regular(
                    100, -3.5, 3.5, name="phi", label=r"$\phi_{k\mu\mu}$"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
        }
        x_hists = {
            "mass": hist.Hist(
                hist.axis.Regular(
                    50, 60, 150, name="mass", label=r"m$_{\mu\mu\gamma}$ [GeV]"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "pt": hist.Hist(
                hist.axis.Regular(
                    100, 0, 200, name="pt", label=r"p$_{T}^{\mu\mu\gamma}$ [GeV]"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "rap": hist.Hist(
                hist.axis.Regular(100, -3, 3, name="rap", label=r"y$_{\mu\mu\gamma}$"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "phi": hist.Hist(
                hist.axis.Regular(
                    100, -3.5, 3.5, name="phi", label=r"$\phi_{\mu\mu\gamma}$"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
        }
        kx_hists = {
            "mass": hist.Hist(
                hist.axis.Regular(
                    50, 60, 150, name="mass", label=r"m$_{k\mu\mu\gamma}$ [GeV]"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "pt": hist.Hist(
                hist.axis.Regular(
                    100, 0, 200, name="pt", label=r"p$_{T}^{k\mu\mu\gamma}$ [GeV]"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "rap": hist.Hist(
                hist.axis.Regular(100, -3, 3, name="rap", label=r"y$_{k\mu\mu\gamma}$"),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
            "phi": hist.Hist(
                hist.axis.Regular(
                    100, -3.5, 3.5, name="phi", label=r"$\phi_{k\mu\mu\gamma}$"
                ),
                hist.axis.StrCategory(["all", "cut", "selected"], name="cat"),
            ),
        }

        cutflow = defaultdict_accumulator(int)

        # --- Cutflow counters -------------------------------------------------
        # ak.count returns the total number of events (non-jagged, length of the chunk)
        n_events = int(len(events))  # same as len(events)

        # --- Event selection cuts ---------------------------------------------
        # ak.count returns the total number of events (non-jagged, length of the chunk)
        trigger_mask = events.HLT.Mu17_Photon30_IsoCaloId
        events = events[trigger_mask]

        muons = events.selectedMuon
        photons = events.selectedPhoton
        onia = events.Onia
        x = events.X
        konia = events.kOnia
        kx = events.kX

        onia["rap"] = onia.rapidity
        konia["rap"] = konia.rapidity
        x["rap"] = x.rapidity
        kx["rap"] = kx.rapidity

        fill_kin_hists(ak.flatten(muons), muon_hists, cat="all")
        fill_kin_hists(ak.flatten(photons), photon_hists, cat="all")
        fill_kin_hists(ak.flatten(onia), onia_hists, cat="all")
        fill_kin_hists(ak.flatten(x), x_hists, cat="all")

        good_muons_mask = (muons.mediumPromptId) & (muons.pfRelIso03_all < 0.15)
        muons = ak.mask(muons, good_muons_mask)
        onia = ak.mask(
            onia,
            (~ak.is_none(muons[onia.muon1Idx].mass, axis=-1))
            & (~ak.is_none(muons[onia.muon2Idx].mass, axis=-1)),
        )
        x = ak.mask(x, ~ak.is_none(onia[x.oniaIdx].mass, axis=-1))

        photons = events.selectedPhoton

        good_photon_mask = (
            (photons.pt > 32)
            & np.logical_xor(
                photons.isScEtaEB,
                photons.isScEtaEE,
            )
            & (photons.mvaID_WP80)
            & (~photons.pixelSeed)
        )

        photons = ak.mask(photons, good_photon_mask)
        x = ak.mask(x, ~ak.is_none(photons[x.photonIdx].pt, axis=-1))

        upsilon_mask = (onia.mass > 8) & (onia.mass < 12)
        onia = ak.mask(onia, upsilon_mask)
        x = ak.mask(x, ~ak.is_none(onia[x.oniaIdx].mass, axis=-1))
        # onia = remove_none(onia)
        upsilon_mask = ak.any(upsilon_mask, axis=-1)

        x_mask = (x.mass > 60) & (x.mass < 150)
        x = x[x_mask]
        x = remove_none(x)
        x_mask = ak.any(x_mask, axis=-1)

        # --- Fill histograms ---------------------------------------------------
        onia["rap"] = onia.rapidity
        konia["rap"] = konia.rapidity
        x["rap"] = x.rapidity
        kx["rap"] = kx.rapidity

        fill_kin_hists(ak.flatten(remove_none(konia)), konia_hists, cat="all")
        fill_kin_hists(ak.flatten(remove_none(kx)), kx_hists, cat="all")

        fill_kin_hists(ak.flatten(remove_none(muons)), muon_hists, cat="cut")
        fill_kin_hists(ak.flatten(remove_none(photons)), photon_hists, cat="cut")
        fill_kin_hists(ak.flatten(remove_none(onia)), onia_hists, cat="cut")
        fill_kin_hists(ak.flatten(remove_none(x)), x_hists, cat="cut")

        fill_kin_hists(
            ak.flatten(remove_none(muons[x.muon1Idx])), muon_hists, cat="selected"
        )
        fill_kin_hists(
            ak.flatten(remove_none(muons[x.muon2Idx])), muon_hists, cat="selected"
        )
        fill_kin_hists(
            ak.flatten(remove_none(photons[x.photonIdx])), photon_hists, cat="selected"
        )
        fill_kin_hists(
            ak.flatten(remove_none(onia[x.oniaIdx])), onia_hists, cat="selected"
        )
        fill_kin_hists(ak.flatten(remove_none(x)), x_hists, cat="selected")

        # --- Accumulate cutflow counters ----------------------------------------
        cutflow["n_events"] = n_events
        cutflow["trigger_sel"] = ak.sum(trigger_mask)
        cutflow["good_muons_sel"] = ak.sum(ak.num(good_muons_mask) > 0)
        cutflow["good_photon_sel"] = ak.sum(good_photon_mask)
        cutflow["upsilon_sel"] = ak.sum(upsilon_mask)
        cutflow["x_sel"] = ak.sum(x_mask)

        out = {
            dataset: {
                "muon_hists": muon_hists,
                "photon_hists": photon_hists,
                "onia_hists": onia_hists,
                "konia_hists": konia_hists,
                "x_hists": x_hists,
                "kx_hists": kx_hists,
                "cutflow": cutflow,
            }
        }

        return out

    def postprocess(self, accumulator):
        return accumulator
