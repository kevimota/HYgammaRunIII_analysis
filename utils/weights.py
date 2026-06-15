# Branching ratios for Upsilon -> mu+ mu-
BR_Y1S_mumu = 0.0248
BR_Y2S_mumu = 0.0193
BR_Y3S_mumu = 0.0218

# Branching ratios for H -> Y(nS) gamma
BR_H_Y1Sgamma = 5.2e-9
BR_H_Y2Sgamma = 1.4e-9
BR_H_Y3Sgamma = 0.9e-9

# Branching ratios for Z -> Y(nS) gamma
BR_Z_Y1Sgamma = 4.8e-8
BR_Z_Y2Sgamma = 2.4e-8
BR_Z_Y3Sgamma = 1.9e-8

# Production cross sections [pb]
sigma_ggH = 52.23  # HIG-19-002   (pb)
sigma_VBF = 4.078
sigma_WH = 1.457
sigma_ZH = 0.9439
sigma_ttH = 0.5700
sigma_bbH = 0.5266

sigma_DY = 2124.08  # Drell-Yan (Z/gamma*)

# Total integrated luminosity [pb^-1]
lumi_2024 = 109.82e3  # 109.82 fb^-1 = 109.82e3 pb^-1

# ---------------------------------------------------------------
# xsec * BR for each production mode and Y(nS) state
# Units: pb (cross section) — final histograms will be scaled by lumi / n_events
# ---------------------------------------------------------------
weights_by_file = {
    # ---- ggH -> H -> Y(nS) gamma -> mu+ mu- gamma ----
    "GluGluToH_HToUps1SG_Ups1SToMuMu": sigma_ggH * BR_H_Y1Sgamma * BR_Y1S_mumu,
    "GluGluToH_HToUps2SG_Ups2SToMuMu": sigma_ggH * BR_H_Y2Sgamma * BR_Y2S_mumu,
    "GluGluToH_HToUps3SG_Ups3SToMuMu": sigma_ggH * BR_H_Y3Sgamma * BR_Y3S_mumu,
    # ---- VBF -> H -> Y(nS) gamma -> mu+ mu- gamma ----
    "VBFToH_HToUps1SG_Ups1SToMuMu": sigma_VBF * BR_H_Y1Sgamma * BR_Y1S_mumu,
    "VBFToH_HToUps2SG_Ups2SToMuMu": sigma_VBF * BR_H_Y2Sgamma * BR_Y2S_mumu,
    "VBFToH_HToUps3SG_Ups3SToMuMu": sigma_VBF * BR_H_Y3Sgamma * BR_Y3S_mumu,
    # ---- WH -> H -> Y(nS) gamma -> mu+ mu- gamma ----
    "WH_HToUps1SG_Ups1SToMuMu": sigma_WH * BR_H_Y1Sgamma * BR_Y1S_mumu,
    "WH_HToUps2SG_Ups2SToMuMu": sigma_WH * BR_H_Y2Sgamma * BR_Y2S_mumu,
    "WH_HToUps3SG_Ups3SToMuMu": sigma_WH * BR_H_Y3Sgamma * BR_Y3S_mumu,
    # ---- ZH -> H -> Y(nS) gamma -> mu+ mu- gamma ----
    "ZH_HToUps1SG_Ups1SToMuMu": sigma_ZH * BR_H_Y1Sgamma * BR_Y1S_mumu,
    "ZH_HToUps2SG_Ups2SToMuMu": sigma_ZH * BR_H_Y2Sgamma * BR_Y2S_mumu,
    "ZH_HToUps3SG_Ups3SToMuMu": sigma_ZH * BR_H_Y3Sgamma * BR_Y3S_mumu,
    # ---- bbH -> H -> Y(nS) gamma -> mu+ mu- gamma ----
    "bbH_HToUps1SG_Ups1SToMuMu": sigma_bbH * BR_H_Y1Sgamma * BR_Y1S_mumu,
    "bbH_HToUps2SG_Ups2SToMuMu": sigma_bbH * BR_H_Y2Sgamma * BR_Y2S_mumu,
    "bbH_HToUps3SG_Ups3SToMuMu": sigma_bbH * BR_H_Y3Sgamma * BR_Y3S_mumu,
    # ---- ttH -> H -> Y(nS) gamma -> mu+ mu- gamma ----
    "ttH_HToUps1SG_Ups1SToMuMu": sigma_ttH * BR_H_Y1Sgamma * BR_Y1S_mumu,
    "ttH_HToUps2SG_Ups2SToMuMu": sigma_ttH * BR_H_Y2Sgamma * BR_Y2S_mumu,
    "ttH_HToUps3SG_Ups3SToMuMu": sigma_ttH * BR_H_Y3Sgamma * BR_Y3S_mumu,
    # ---- qqZ -> Z -> Y(nS) gamma -> mu+ mu- gamma ----
    "ZToUps1SG_Ups1SToMuMu": sigma_DY * BR_Z_Y1Sgamma * BR_Y1S_mumu,
    "ZToUps2SG_Ups2SToMuMu": sigma_DY * BR_Z_Y2Sgamma * BR_Y2S_mumu,
    "ZToUps3SG_Ups3SToMuMu": sigma_DY * BR_Z_Y3Sgamma * BR_Y3S_mumu,
    # ---- ggH -> H -> mu+ mu- gamma (direct, no intermediate Y) ----
    # TODO: update BR(H->mu+mu-gamma) with the correct value for this sample
    "ggH125_012j_NLO_FXFX_HtoMuMuGamma": 2.1337e-03,  # placeholder — adjust!
    # ---- Z/gamma* -> mu+ mu- gamma (m_ll in [2,15]) ----
    # TODO: update cross section for this specific phase space
    "ZGTo2MuG_mll_2to15_LO": 7.9260e-02,  # placeholder — adjust!
}


def get_xsec_br(filename):
    """Return xsec * BR [pb] for a given pickle filename (without path/extension).

    Parameters
    ----------
    filename : str
        The basename of the pickle file, e.g. 'GluGluToH_HToUps1SG_Ups1SToMuMu'.

    Returns
    -------
    float
        Cross section times branching ratio in pb.
    """
    # Strip directory prefix and .pkl extension
    name = filename
    if name.endswith(".pkl"):
        name = name[:-4]
    if "/" in name:
        name = name.rsplit("/", 1)[1]

    return weights_by_file.get(name, 0.0)
