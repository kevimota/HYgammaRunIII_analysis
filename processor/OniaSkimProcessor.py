from coffea.processor import ProcessorABC, defaultdict_accumulator, column_accumulator
import awkward as ak
import numpy as np


class OniaSkimProcessor(ProcessorABC):

    def __init__(self):
        pass

    def process(self, events):
        onia = events.Onia
        onia["rap"] = onia.rapidity
        onia["Muons"] = ak.zip(
            {
                "i0": events.selectedMuon[onia.muon1Idx],
                "i1": events.selectedMuon[onia.muon2Idx],
            }
        )

        n_raw = int(ak.sum(ak.num(onia, axis=1)))

        trigger_mask = events.HLT.Mu17_Photon30_IsoCaloId
        onia = onia[trigger_mask]

        good_muons_mask = (
            (onia.Muons.i0.mediumPromptId)
            & (onia.Muons.i0.pfRelIso03_all < 0.15)
            & (onia.Muons.i1.mediumPromptId)
            & (onia.Muons.i1.pfRelIso03_all < 0.15)
        )
        onia = onia[good_muons_mask]

        onia_mass_mask = (onia.mass > 8) & (onia.mass < 12)
        onia = onia[onia_mass_mask]

        onia = ak.flatten(onia)

        out = {"n_raw": defaultdict_accumulator(int)}
        out["n_raw"]["count"] += n_raw

        if len(onia) == 0:
            empty = column_accumulator(np.array([]))
            cols = [
                "onia_mass", "onia_pt", "onia_eta", "onia_phi", "onia_rap",
                "muon1_pt", "muon1_eta", "muon1_phi",
                "muon2_pt", "muon2_eta", "muon2_phi",
            ]
            for c in cols:
                out[c] = empty
            return out

        out["onia_mass"] = column_accumulator(ak.to_numpy(onia.mass))
        out["onia_pt"] = column_accumulator(ak.to_numpy(onia.pt))
        out["onia_eta"] = column_accumulator(ak.to_numpy(onia.eta))
        out["onia_phi"] = column_accumulator(ak.to_numpy(onia.phi))
        out["onia_rap"] = column_accumulator(ak.to_numpy(onia.rap))

        out["muon1_pt"] = column_accumulator(ak.to_numpy(onia.Muons.i0.pt))
        out["muon1_eta"] = column_accumulator(ak.to_numpy(onia.Muons.i0.eta))
        out["muon1_phi"] = column_accumulator(ak.to_numpy(onia.Muons.i0.phi))
        out["muon2_pt"] = column_accumulator(ak.to_numpy(onia.Muons.i1.pt))
        out["muon2_eta"] = column_accumulator(ak.to_numpy(onia.Muons.i1.eta))
        out["muon2_phi"] = column_accumulator(ak.to_numpy(onia.Muons.i1.phi))

        return out

    def postprocess(self, accumulator):
        return accumulator
