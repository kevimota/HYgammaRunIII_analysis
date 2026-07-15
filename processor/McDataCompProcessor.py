"""McDataCompProcessor: coffea processor for MC vs data comparison plots.

Produces histograms, cutflow, and selected event values (x_mass, onia_mass)
for stacked plot comparisons between MC samples and data.

Output
------
    Accumulator containing histograms, cutflow counters, and column accumulators.
"""

from coffea.processor import ProcessorABC, defaultdict_accumulator, column_accumulator
from hist import Hist
import awkward as ak
import numpy as np
from utils import fill_kin_hists


def fill_hists(x, hists, cat):
    for h in hists:
        if h == "muon_lead":
            fill_kin_hists(x.Muons.i0, hists[h], cat=cat)
        elif h == "muon_trail":
            fill_kin_hists(x.Muons.i1, hists[h], cat=cat)
        elif h == "photon":
            fill_kin_hists(x.Photon, hists[h], cat=cat)
        elif h == "onia":
            fill_kin_hists(x.Onia, hists[h], cat=cat)
        elif h == "x":
            fill_kin_hists(x, hists[h], cat=cat)


class McDataCompProcessor(ProcessorABC):

    def __init__(self):
        """Initialize the processor."""

    def process(self, events):
        dataset = events.metadata["dataset"]

        muon_lead_hists = {
            "pt": Hist.new.Regular(100, 0, 200, name="pt", label=r"p$_{T}^{\mu}$ [GeV]")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "eta": Hist.new.Regular(100, -3, 3, name="eta", label=r"$\eta_{\mu}$")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "phi": Hist.new.Regular(100, -3.5, 3.5, name="phi", label=r"$\phi_{\mu}$")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
        }

        muon_trail_hists = {
            "pt": Hist.new.Regular(100, 0, 200, name="pt", label=r"p$_{T}^{\mu}$ [GeV]")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "eta": Hist.new.Regular(100, -3, 3, name="eta", label=r"$\eta_{\mu}$")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "phi": Hist.new.Regular(100, -3.5, 3.5, name="phi", label=r"$\phi_{\mu}$")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
        }

        photon_hists = {
            "pt": Hist.new.Regular(
                100, 0, 200, name="pt", label=r"E$_{T}^{\gamma}$ [GeV]"
            )
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "eta": Hist.new.Regular(100, -3, 3, name="eta", label=r"$\eta_{\gamma}$")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "phi": Hist.new.Regular(
                100, -3.5, 3.5, name="phi", label=r"$\phi_{\gamma}]$"
            )
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
        }

        onia_hists = {
            "mass": Hist.new.Regular(
                50, 8, 12, name="mass", label=r"m$_{\mu\mu}$ [GeV]"
            )
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "pt": Hist.new.Regular(
                100, 0, 200, name="pt", label=r"p$_{T}^{\mu\mu}$ [GeV]"
            )
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "rap": Hist.new.Regular(100, -3, 3, name="rap", label=r"y$_{\mu\mu}$")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "phi": Hist.new.Regular(
                100, -3.5, 3.5, name="phi", label=r"$\phi_{\mu\mu}$"
            )
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
        }

        x_hists = {
            "mass": Hist.new.Regular(
                50, 60, 150, name="mass", label=r"m$_{\mu\mu\gamma}$ [GeV]"
            )
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "pt": Hist.new.Regular(
                100, 0, 200, name="pt", label=r"p$_{T}^{\mu\mu\gamma}$ [GeV]"
            )
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "rap": Hist.new.Regular(100, -3, 3, name="rap", label=r"y$_{\mu\mu\gamma}$")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "phi": Hist.new.Regular(
                100, -3.5, 3.5, name="phi", label=r"$\phi_{\mu\mu\gamma}$"
            )
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
            "n": Hist.new.Regular(10, -0.5, 9.5, name="n", label=r"X multiplicity")
            .StrCategory(
                [
                    "all",
                    "trigger",
                    "good_muons",
                    "good_photon",
                    "upsilon_mass",
                    "x_mass",
                ],
                name="cat",
            )
            .Double(),
        }

        hists = {
            "muon_lead": muon_lead_hists,
            "muon_trail": muon_trail_hists,
            "photon": photon_hists,
            "onia": onia_hists,
            "x": x_hists,
        }

        # --- Accumulators definition --------------------------------------------

        cutflow = defaultdict_accumulator(int)

        # --- Cutflow counters -------------------------------------------------
        # ak.count returns the total number of events (non-jagged, length of the chunk)
        n_events = int(len(events))  # same as len(events)

        # --- Event selection cuts ---------------------------------------------
        # ak.count returns the total number of events (non-jagged, length of the chunk)

        x = events.X
        onia = events.Onia[x.oniaIdx]
        onia["rap"] = onia.rapidity
        x["Onia"] = onia
        x["Muons"] = ak.zip(
            {
                "i0": events.selectedMuon[x.muon1Idx],
                "i1": events.selectedMuon[x.muon2Idx],
            }
        )
        x["Photon"] = events.selectedPhoton[x.photonIdx]
        x["rap"] = x.rapidity

        fill_hists(x, hists, "all")

        trigger_mask = events.HLT.Mu17_Photon30_IsoCaloId
        x = x[trigger_mask]

        good_muons_mask = (
            (x.Muons.i0.mediumPromptId)
            & (x.Muons.i0.pfRelIso03_all < 0.15)
            & (x.Muons.i1.mediumPromptId)
            & (x.Muons.i1.pfRelIso03_all < 0.15)
        )

        good_photon_mask = (
            (x.Photon.pt > 32)
            & np.logical_xor(
                x.Photon.isScEtaEB,
                x.Photon.isScEtaEE,
            )
            & (x.Photon.mvaID_WP80)
            & (~x.Photon.pixelSeed)
        )

        upsilon_mask = (x.Onia.mass > 8) & (x.Onia.mass < 12)

        x_mask = (x.mass > 60) & (x.mass < 150)

        # --- Accumulate cutflow counters ----------------------------------------
        cutflow["n_events"] = n_events
        cutflow["trigger_sel"] = ak.sum(trigger_mask)
        cutflow["good_muons_sel"] = ak.sum(ak.any(good_muons_mask, axis=1))
        cutflow["good_photon_sel"] = ak.sum(ak.any(good_photon_mask, axis=1))
        cutflow["upsilon_sel"] = ak.sum(ak.any(upsilon_mask, axis=1))
        cutflow["x_sel"] = ak.sum(ak.any(x_mask, axis=1))
        cutflow["flow_good_photon_sel"] = ak.sum(
            ak.any(good_muons_mask & good_photon_mask, axis=1)
        )
        cutflow["flow_upsilon_sel"] = ak.sum(
            ak.any(good_muons_mask & good_photon_mask & upsilon_mask, axis=1)
        )
        cutflow["flow_x_sel"] = ak.sum(
            ak.any(good_muons_mask & good_photon_mask & upsilon_mask & x_mask, axis=1)
        )

        # --- Fill histograms ----------------------------------------

        fill_hists(x, hists, "trigger")
        fill_hists(x[good_muons_mask], hists, "good_muons")
        fill_hists(x[good_muons_mask & good_photon_mask], hists, "good_photon")
        fill_hists(
            x[good_muons_mask & good_photon_mask & upsilon_mask], hists, "upsilon_mass"
        )
        fill_hists(
            x[good_muons_mask & good_photon_mask & upsilon_mask & x_mask],
            hists,
            "x_mass",
        )

        # --- Apply all cuts to x ----------------------------------------

        x = x[good_muons_mask & good_photon_mask & upsilon_mask & x_mask]

        args = ak.argsort(x.Onia.vProb)
        x = x[args]
        x = x[ak.num(x) > 0]
        if len(x) > 0:
            x = x[:, 0]
            x_mass = column_accumulator(x.mass.to_numpy())
            onia_mass = column_accumulator(x.Onia.mass.to_numpy())
        else:
            x_mass = column_accumulator(np.array([]))
            onia_mass = column_accumulator(np.array([]))

        out = {
            "x_mass": x_mass,
            "onia_mass": onia_mass,
            "hists": hists,
            "cutflow": cutflow,
        }

        return out

    def postprocess(self, accumulator):
        return accumulator
