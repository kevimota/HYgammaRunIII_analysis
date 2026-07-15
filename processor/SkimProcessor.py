"""SkimProcessor: coffea processor for creating flat ntuples of selected events.

Applies two selection stages:
  1. Preselection  — trigger + good muons + good photon
  2. Full selection — preselection + upsilon mass + X mass

Each stage picks the best X candidate per event by highest pT and saves
per-object kinematics as column accumulators with a stage-specific prefix.

Output
------
    column_accumulator arrays for presel_* and fullsel_* columns, plus n_raw.
"""

from coffea.processor import ProcessorABC, defaultdict_accumulator, column_accumulator
import awkward as ak
import numpy as np


def _save_candidate(out, prefix, x):
    if len(x) == 0:
        empty = column_accumulator(np.array([]))
        cols = [
            "cand_boson_mass", "cand_boson_pt", "cand_boson_eta", "cand_boson_phi", "cand_boson_rap",
            "cand_meson_mass", "cand_meson_pt", "cand_meson_eta", "cand_meson_phi", "cand_meson_rap",
            "muon1_pt", "muon1_eta", "muon1_phi",
            "muon2_pt", "muon2_eta", "muon2_phi",
            "gamma_pt", "gamma_eta", "gamma_phi",
        ]
        for c in cols:
            out[f"{prefix}_{c}"] = empty
        return

    out[f"{prefix}_cand_boson_mass"] = column_accumulator(ak.to_numpy(x.mass))
    out[f"{prefix}_cand_boson_pt"] = column_accumulator(ak.to_numpy(x.pt))
    out[f"{prefix}_cand_boson_eta"] = column_accumulator(ak.to_numpy(x.eta))
    out[f"{prefix}_cand_boson_phi"] = column_accumulator(ak.to_numpy(x.phi))
    out[f"{prefix}_cand_boson_rap"] = column_accumulator(ak.to_numpy(x.rap))

    out[f"{prefix}_cand_meson_mass"] = column_accumulator(ak.to_numpy(x.Onia.mass))
    out[f"{prefix}_cand_meson_pt"] = column_accumulator(ak.to_numpy(x.Onia.pt))
    out[f"{prefix}_cand_meson_eta"] = column_accumulator(ak.to_numpy(x.Onia.eta))
    out[f"{prefix}_cand_meson_phi"] = column_accumulator(ak.to_numpy(x.Onia.phi))
    out[f"{prefix}_cand_meson_rap"] = column_accumulator(ak.to_numpy(x.Onia.rap))

    out[f"{prefix}_muon1_pt"] = column_accumulator(ak.to_numpy(x.Muons.i0.pt))
    out[f"{prefix}_muon1_eta"] = column_accumulator(ak.to_numpy(x.Muons.i0.eta))
    out[f"{prefix}_muon1_phi"] = column_accumulator(ak.to_numpy(x.Muons.i0.phi))
    out[f"{prefix}_muon2_pt"] = column_accumulator(ak.to_numpy(x.Muons.i1.pt))
    out[f"{prefix}_muon2_eta"] = column_accumulator(ak.to_numpy(x.Muons.i1.eta))
    out[f"{prefix}_muon2_phi"] = column_accumulator(ak.to_numpy(x.Muons.i1.phi))

    out[f"{prefix}_gamma_pt"] = column_accumulator(ak.to_numpy(x.Photon.pt))
    out[f"{prefix}_gamma_eta"] = column_accumulator(ak.to_numpy(x.Photon.eta))
    out[f"{prefix}_gamma_phi"] = column_accumulator(ak.to_numpy(x.Photon.phi))


def _pick_best(x):
    if len(x) == 0:
        return x
    args = ak.argsort(x.pt, ascending=False)
    x = x[args]
    x = x[ak.num(x) > 0]
    if len(x) > 0:
        x = x[:, 0]
    return x


class SkimProcessor(ProcessorABC):

    def __init__(self):
        """Initialize the processor."""

    def process(self, events):
        dataset = events.metadata["dataset"]

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

        n_raw = len(events)

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

        presel_cut = good_muons_mask & good_photon_mask
        presel = _pick_best(x[presel_cut])

        upsilon_mask = (x.Onia.mass > 8) & (x.Onia.mass < 12)
        x_mass_mask = (x.mass > 60) & (x.mass < 150)

        full_sel_cut = presel_cut & upsilon_mask & x_mass_mask
        full_sel = _pick_best(x[full_sel_cut])

        out = {"n_raw": defaultdict_accumulator(int)}
        out["n_raw"]["count"] += n_raw

        _save_candidate(out, "presel", presel)
        _save_candidate(out, "fullsel", full_sel)

        return out

    def postprocess(self, accumulator):
        return accumulator
