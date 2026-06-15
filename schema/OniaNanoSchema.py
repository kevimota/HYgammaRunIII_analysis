"""OniaNanoSchema: custom NanoAOD-like schema for H/Z -> mu+ mu- gamma onia analysis.

This schema groups flat columnar ROOT branches into awkward array zip collections,
following the NanoAOD naming convention (nXxx counter + Xxx_variable branches).
Each collection is exposed as an attribute of the NanoEvents object with full
Lorentz vector and candidate method support.

Collections
-----------
selectedMuon  : custom-selected muons (pt>5, |eta|<2.4, PF+Global, relIso<0.20)
refittedMuon  : muons after kinematic mass-constrained vertex fit
selectedPhoton: custom-selected photons (pt>10, |eta|<2.4, relIso<0.20)
Onia         : raw unfitted dimuon candidates (2 < m_mumu < 15 GeV)
kOnia        : mass-constrained kinematically-fitted dimuon candidates
X            : unfitted Onia + photon candidate (H/Z -> mu mu gamma)
kX           : kinematically-fitted kOnia + photon candidate

Index linking (resolved via *_IdxG global indexers in coffea >= 0.7)
-----------
Onia_muon1Idx / Onia_muon2Idx      : indices into selectedMuon (original muons)
kOnia_muon1Idx / kOnia_muon2Idx    : indices into refittedMuon (refitted muons)
X_oniaIdx   : index into Onia for this X candidate
X_photonIdx : index into selectedPhoton for this X candidate
kX_*Idx     : same pattern but referencing kOnia collections
"""

from coffea.nanoevents.schemas.nanoaod import NanoAODSchema


class OniaNanoSchema(NanoAODSchema):
    """Custom NanoAOD-like schema for Onia + photon analysis.

    Inherits from NanoAODSchema to reuse the auto-collection grouping logic:
    any branch ``n{name}`` is treated as a counts array, and all branches
    matching ``{name}_*`` are zipped into a jagged collection.
    The base class automatically converts ``*_Idx`` branches to global
    offsets-aware indices (named ``*_IdxG``), so that the indirection
    survives event-level slicing.

    Mixin types (``mixins``) determine the behaviour for each collection.
    ``Muon`` / ``Photon`` provide `.p4`, `.pt`, `.eta`, `.phi`, `.mass`,
    plus particle-specific fields (isolation, dxy/dz, etc.).
    ``Candidate`` extends that with ``.charge`` for charge-aware operations.
    """

    #: Mixin types per collection.
    #: - ``Muon``    : PtEtaPhiM + charge + muon-specific fields (dxy, dz, iso)
    #: - ``Photon``  : PtEtaPhiM + charge + photon-specific fields
    #: - ``Candidate``: PtEtaPhiM + charge (generic dimuon / X candidate)
    mixins = {
        "selectedMuon": "Muon",
        "refittedMuon": "Muon",
        "selectedPhoton": "Photon",
        "Onia": "Candidate",
        "kOnia": "Candidate",
        "X": "Candidate",
        "kX": "Candidate",
    }

    #: Cross-reference indexers: ``{branch}`` -> ``{target_collection}``.
    #: The base NanoAODSchema converts these to global ``{branch}G`` arrays
    #: using ``local2global`` with the target's offsets array.
    all_cross_references = {
        "Onia_muon1Idx": "selectedMuon",
        "Onia_muon2Idx": "selectedMuon",
        "kOnia_muon1Idx": "selectedMuon",
        "kOnia_muon2Idx": "selectedMuon",
        "kOnia_refittedMuon1Idx": "refittedMuon",
        "kOnia_refittedMuon2Idx": "refittedMuon",
        "X_oniaIdx": "Onia",
        "X_muon1Idx": "selectedMuon",
        "X_muon2Idx": "selectedMuon",
        "X_photonIdx": "selectedPhoton",
        "kX_oniaIdx": "Onia",
        "kX_refittedMuon1Idx": "refittedMuon",
        "kX_refittedMuon2Idx": "refittedMuon",
        "kX_photonIdx": "selectedPhoton",
    }

    #: Nested indexers: fixed-length groups of ``Idx1, Idx2, ...`` branches
    #: are combined into a single jagged ``{name}IdxG`` array by the base class.
    #: The base class converts each ``Idx{N}`` branch to ``Idx{N}G`` first,
    #: then combines them.  Here we list the final nested field names we want
    #: (without the ``G`` suffix), and the base class resolves the component
    #: branches automatically.
    nested_items = {
        "Onia_muonIdxG": ["Onia_muon1Idx", "Onia_muon2Idx"],
        "kOnia_muonIdxG": ["kOnia_muon1Idx", "kOnia_muon2Idx"],
        "kOnia_refittedMuonIdxG": ["kOnia_refittedMuon1Idx", "kOnia_refittedMuon2Idx"],
        "X_muonIdxG": ["X_muon1Idx", "X_muon2Idx"],
        "kX_muonIdxG": ["kX_muon1Idx", "kX_muon2Idx"],
        "kX_refittedMuonIdxG": ["kX_refittedMuon1Idx", "kX_refittedMuon2Idx"],
    }

    @classmethod
    def behavior(cls):
        """Merge NanoAOD, vector, and candidate behaviours."""
        from coffea.nanoevents.methods import nanoaod
        from coffea.nanoevents.methods import vector
        from coffea.nanoevents.methods import candidate

        beh = dict(nanoaod.behavior)
        beh.update(vector.behavior)
        beh.update(candidate.behavior)
        return beh