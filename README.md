# HYgammaRunIII Analysis

CMS Run 3 analysis of **H/Z → μ⁺μ⁻γ** via an intermediate Υ(nS) state,  
using NanoAOD-like ROOT data processed with the [coffea](https://coffea-hep.readthedocs.io/) columnar framework.

## Project structure

```
├── run_analysis.py              # McDataCompProcessor: MC vs data comparison
├── run_skim.py                  # SkimProcessor: flat ntuple creation
├── run_onia_skim.py             # OniaSkimProcessor: upsilon mass fitting ntuples
├── run_all.py                   # Batch runner (compare + skim + onia_skim)
├── processor/
│   ├── McDataCompProcessor.py    # Histograms, cutflow, event selection
│   ├── SkimProcessor.py          # Flat ntuple of selected events
│   ├── OniaSkimProcessor.py      # Upsilon-only flat ntuples (all candidates)
│   └── TestProcessor.py          # Test/legacy processor
├── schema/
│   └── OniaNanoSchema.py         # Custom NanoAOD schema (collections, cross-refs, mixins)
├── utils/
│   ├── figure.py                 # Plotting (parallel, --workers N)
│   ├── weights.py                # Cross-section × BR weights for MC scaling
│   └── utils.py                  # Helpers (file discovery, histogram filling)
├── output/                       # Generated outputs (gitignored)
├── hists/                        # Generated plots (gitignored)
└── requirements.txt
```

## Dependencies

Install with pip (conda environment recommended):

```bash
pip install -r requirements.txt
```

Requires Python ≥ 3.9, coffea ≥ 0.7, uproot ≥ 5, awkward ≥ 2, hist ≥ 2.

Listed libraries:
| Package    | Purpose |
|------------|---------|
| `coffea`   | Columnar processing framework (Runner, ProcessorABC, executors) |
| `uproot`   | ROOT I/O |
| `awkward`  | Jagged/columnar array operations |
| `hist`     | Histogramming |
| `numpy`    | Numeric operations |
| `typer`    | CLI interface |

## Schema: `OniaNanoSchema`

The project uses a custom NanoAOD-like schema for the H/Z → μμγ analysis
([`schema/OniaNanoSchema.py`](schema/OniaNanoSchema.py)), derived from `NanoAODSchema`.

### Collections

| Collection      | Mixin       | Contents |
|-----------------|-------------|----------|
| `selectedMuon`  | `Muon`      | Custom-selected muons (pT > 5, \|η\| < 2.4, PF+Global, relIso < 0.20) |
| `refittedMuon`  | `Muon`      | Muons after kinematic mass-constrained vertex fit |
| `selectedPhoton`| `Photon`    | Custom-selected photons (pT > 10, \|η\| < 2.4, relIso < 0.20) |
| `Onia`          | `Candidate` | Raw unfitted dimuon candidates (2 < m<sub>μμ</sub> < 15 GeV) |
| `kOnia`         | `Candidate` | Mass-constrained kinematically-fitted dimuon candidates |
| `X`             | `Candidate` | Unfitted Onia + photon candidate (H/Z → μμγ) |
| `kX`            | `Candidate` | Kinematically-fitted kOnia + photon candidate |

### Index linking (global cross-references)

Global indexers link collections across jagged nesting levels,
resolved via `*_IdxG` arrays by coffea ≥ 0.7:

| Index branch              | Target collection  |
|---------------------------|--------------------|
| `Onia_muon1Idx`           | `selectedMuon`     |
| `Onia_muon2Idx`           | `selectedMuon`     |
| `X_oniaIdx`               | `Onia`             |
| `X_photonIdx`             | `selectedPhoton`   |
| `X_muon1Idx`              | `selectedMuon`     |
| `X_muon2Idx`              | `selectedMuon`     |
| `kX_oniaIdx`              | `Onia`             |
| `kX_photonIdx`            | `selectedPhoton`   |
| *(and corresponding kOnia/kX indices)* |                    |

### Navigating from X in a processor

```python
x = events.X
onia = events.Onia[x.oniaIdx]               # dimuon matched to this X
mu1 = events.selectedMuon[x.muon1Idx]       # leading muon
mu2 = events.selectedMuon[x.muon2Idx]       # subleading muon
pho = events.selectedPhoton[x.photonIdx]    # photon
```

## Processors

### `McDataCompProcessor` ([`processor/McDataCompProcessor.py`](processor/McDataCompProcessor.py))

Produces histograms, cutflow counters, and selected-candidate values
for MC vs data comparison.

**Selection chain:**
1. All events
2. Trigger: `HLT_Mu17_Photon30_IsoCaloId`
3. Good muons: `mediumPromptId` + `pfRelIso03_all < 0.15` on both muons
4. Good photon: pT > 32, barrel/endcap-only, `mvaID_WP80`, no `pixelSeed`
5. Upsilon mass window: 8 < m<sub>μμ</sub> < 12 GeV
6. X mass window: 60 < m<sub>μμγ</sub> < 150 GeV

**Best-X selection:** picks the candidate with highest `Onia.vProb`
(vertex probability).

**Output:** histograms + cutflow + `x_mass`/`onia_mass` column accumulators.

### `SkimProcessor` ([`processor/SkimProcessor.py`](processor/SkimProcessor.py))

Creates flat ntuples of selected events.

**Selection chain:** same as McDataCompProcessor.

**Two-stage output:**
- **Preselection** (trigger + good muons + good photon)
- **Full selection** (preselection + upsilon + X mass windows)

**Best-X selection:** picks the candidate with highest pT.

**Output columns:** `cand_boson_*` (mass, pt, eta, phi, rap),
`cand_meson_*`, `muon1_*`, `muon2_*`, `gamma_*` — all prefixed with
`presel_` or `fullsel_` depending on the stage.

### `OniaSkimProcessor` ([`processor/OniaSkimProcessor.py`](processor/OniaSkimProcessor.py))

Data-only processor for upsilon mass fitting. Keeps **all** onia candidates
per event (no best-candidate picking) for maximum statistics.

**Selection chain:**
1. Trigger: `HLT_Mu17_Photon30_IsoCaloId`
2. Good muons: `mediumPromptId` + `pfRelIso03_all < 0.15` on both muons
3. Onia mass window: 8 < m<sub>μμ</sub> < 12 GeV

**Output columns** (flat, one row per onia candidate):

| Column | Source |
|--------|--------|
| `onia_mass`, `onia_pt`, `onia_eta`, `onia_phi`, `onia_rap` | Onia kinematics |
| `muon1_pt`, `muon1_eta`, `muon1_phi` | Leading muon (from `onia.muon1Idx`) |
| `muon2_pt`, `muon2_eta`, `muon2_phi` | Subleading muon (from `onia.muon2Idx`) |

### `TestProcessor` ([`processor/TestProcessor.py`](processor/TestProcessor.py))

Legacy test processor with kOnia/kX collections. Not used in production.

## Utilities

### `utils/utils.py`

| Function         | Purpose |
|------------------|---------|
| `get_files(paths, pattern, exclude)` | List `.root` files in directories, sorted naturally, excluding `_hist` files |
| `fill_kin_hists(obj, hists, cat)`    | Fill histogram dict from an awkward array |
| `remove_none(arr)`                   | Remove `None` entries from a jagged array |
| `natural_keys(text)`                 | Human-friendly sort key |

### `utils/weights.py`

Cross-section × branching ratio lookup for MC scaling:

| Function / Variable      | Value |
|--------------------------|-------|
| `get_xsec_br(filename)`  | Returns σ×BR [pb] for a dataset name, or `0.0` for data |
| `lumi_2024`              | 109.82 fb⁻¹ (109.82 × 10³ pb⁻¹) |
| `sigma_ggH`              | 52.23 pb |

To add a new MC sample, add its σ×BR to `weights_by_file` dict
in [`utils/weights.py`](utils/weights.py).

### `utils/figure.py`

Parallel plotting: loads McDataCompProcessor output pickles, merges
MC samples grouped by production mode, and produces stacked plots,
comparison plots, and cutflow bar charts.

**Blinding:** X-mass histograms blinded in Z (80–100 GeV) and
H (120–130 GeV) windows.

## CLI usage

### `run_analysis.py` — MC vs data comparison

```bash
python run_analysis.py /path/to/files --dataset mydata
python run_analysis.py /path/to/files --dataset mymc --is-mc
python run_analysis.py /path/to/files --dataset mydata --parallel -w 8
```

**Options:**
| Flag | Default | Description |
|------|---------|-------------|
| `paths` (positional) | required | Input ROOT files or directories |
| `--dataset` / `-d` | required | Dataset name |
| `--is-mc` | `False` | Mark as MC |
| `--parallel` / `-p` | `False` | Use `FuturesExecutor` |
| `--workers` / `-w` | `4` | Number of parallel workers |
| `--chunksize` | `100000` | Events per chunk |
| `--maxchunks` | `None` | Max chunks to process (all) |
| `--treename` | `Events` | TTree name |

**Output:** `output/mc_data_comp/{dataset}.pkl` + `output/mc_data_comp/{dataset}.root`

### `run_skim.py` — Flat ntuple creation

```bash
python run_skim.py /path/to/files --dataset mydata
python run_skim.py /path/to/files --dataset mymc --is-mc --parallel -w 8
```

Same options as `run_analysis.py`.

**Output:** two ROOT files per dataset:

| File | Contents |
|------|----------|
| `output/skim/{dataset}_presel.root` | Trigger + good muons + good photon |
| `output/skim/{dataset}_fullsel.root` | Preselection + upsilon + X mass windows |

Each output ROOT file has two TTrees:
- `events/`: `cand_boson_*` (mass, pt, eta, phi, rap), `cand_meson_*`, `muon1_*`, `muon2_*`, `gamma_*`
- `metadata/`: `dataset`, `is_mc`, `n_raw`, `n_sel`, `weight` (one entry)

The per-event `weight` column is `σ×BR × L / N_raw` for MC, `1` for data.

### `run_onia_skim.py` — Upsilon mass fitting ntuples (data only)

```bash
python run_onia_skim.py /path/to/files --dataset MuonEG_Run2024
python run_onia_skim.py /path/to/files --dataset MuonEG_Run2024 --parallel -w 8
```

**Options:**
| Flag | Default | Description |
|------|---------|-------------|
| `paths` (positional) | required | Input ROOT files or directories |
| `--dataset` / `-d` | required | Dataset name |
| `--parallel` / `-p` | `False` | Use `FuturesExecutor` |
| `--workers` / `-w` | `4` | Number of parallel workers |
| `--chunksize` | `100000` | Events per chunk |
| `--maxchunks` | `None` | Max chunks to process (all) |
| `--treename` | `Events` | TTree name |

**Output:** `output/onia_skim/{dataset}.root` — single `events/` TTree with
onia and muon columns, no metadata.

### `run_all.py` — Batch processing

```bash
python run_all.py compare                    # all datasets, compare mode
python run_all.py compare --workers 16       # custom parallelism
python run_all.py compare --data-path /custom --mc-path /custom

python run_all.py skim                       # all datasets, skim mode
python run_all.py skim --workers 16

python run_all.py onia-skim                  # upsilon mass fitting ntuples
python run_all.py onia-skim --workers 16
```

**Options (shared by `compare` and `skim`):**
| Flag | Default | Description |
|------|---------|-------------|
| `--data-path` | `~/cernbox/rare_decays/MuonEG` | Data directory |
| `--mc-path` | `~/cernbox/rare_decays/MC` | MC directory |
| `--workers` / `-w` | `12` | Workers per executor |
| `--parallel` / `--serial` | `parallel` | Use parallel executor |

**`compare` mode:** processes each leaf directory individually (one dataset per call),
same behavior as the original subprocess-based runner.

**`skim` mode:**
- **Data:** all `MuonEG_Run2024*` leaf directories combined into a single
  ntuple named `MuonEG_Run2024` (one output pair for the full year).
- **MC:** one `_run_skim()` call per MC sample (grouped by dataset name).

**`onia-skim` mode:** data only — automatically detects the data-taking year
from folder names and groups files by year (e.g. all `Run2024*` folders →
`MuonEG_Run2024`). When future Run 3 years are added, they will automatically
produce separate outputs (`MuonEG_Run2025`, etc.). No MC processing.

| Flag | Default | Description |
|------|---------|-------------|
| `--data-path` | `~/cernbox/rare_decays/MuonEG` | Data directory (MC not used) |
| `--workers` / `-w` | `12` | Workers per executor |
| `--parallel` / `--serial` | `parallel` | Use parallel executor |

### `utils/figure.py` — Plotting

```bash
python utils/figure.py                            # all plots (default: 4 workers)
python utils/figure.py --workers 16               # custom parallelism
python utils/figure.py --save-per-dataset         # individual dataset plots
```

**Options:**
| Flag | Default | Description |
|------|---------|-------------|
| `--workers` / `-w` | `4` | Number of parallel processes |
| `--save-per-dataset` | `False` | Save per-dataset plots |

**Output directories under `hists/`:**
| Directory | Contents |
|-----------|----------|
| `combined/{cat}/` | MC-stacked + data errorbar plots (one per var) |
| `comparison/{ds}/` | Preselection vs full-selection comparison |
| `cutflow/` | Per-dataset and combined cutflow bar charts |

## Output formats

### `run_analysis.py` — pickled accumulator

The `.pkl` file in `output/mc_data_comp/` is a dict with keys:

| Key | Type | Description |
|-----|------|-------------|
| `hists` | `dict` | Nested `{group: {var: hist.Hist}}` for μ₁, μ₂, γ, Onia, X |
| `cutflow` | `dict` | Event counts at each selection stage |
| `x_mass` | `column_accumulator` | Best X candidate masses |
| `onia_mass` | `column_accumulator` | Best Onia candidate masses |
| `dataset` | `str` | Dataset name |
| `is_mc` | `bool` | `True` for MC |

Cutflow stages: `n_events`, `trigger_sel`, `good_muons_sel`, `good_photon_sel`,
`flow_good_photon_sel`, `flow_upsilon_sel`, `flow_x_sel`, `upsilon_sel`, `x_sel`.

Histogram groups: `muon_lead`, `muon_trail`, `photon`, `onia`, `x`.

Selection categories: `all`, `trigger`, `good_muons`, `good_photon`, `upsilon_mass`, `x_mass`.

The `.root` file contains:
- `events/` TTree: `x_mass`, `onia_mass`
- `metadata/` TTree: `is_mc`, `dataset`
- Directories per histogram group: `{group}/{var}/{cat}`

### `run_skim.py` — ROOT ntuples

Two files per dataset: `{dataset}_presel.root` and `{dataset}_fullsel.root`.

**`events/` TTree columns:**

| Column group | Fields |
|--------------|--------|
| `cand_boson_*` | `mass`, `pt`, `eta`, `phi`, `rap` |
| `cand_meson_*` | `mass`, `pt`, `eta`, `phi`, `rap` |
| `muon1_*` | `pt`, `eta`, `phi` |
| `muon2_*` | `pt`, `eta`, `phi` |
| `gamma_*` | `pt`, `eta`, `phi` |

**`metadata/`:** `dataset`, `is_mc`, `n_raw`, `n_sel`, `weight` (single entry).

### `run_onia_skim.py` — ROOT flat ntuple

Single file: `output/onia_skim/{dataset}.root`

**`events/` TTree** (one row per onia candidate, no metadata TTree):

| Column | Description |
|--------|-------------|
| `onia_mass` | Dimuon invariant mass [GeV] |
| `onia_pt` | Dimuon transverse momentum [GeV] |
| `onia_eta` | Dimuon pseudorapidity |
| `onia_phi` | Dimuon azimuthal angle |
| `onia_rap` | Dimuon rapidity |
| `muon1_pt` | Leading muon pT [GeV] |
| `muon1_eta` | Leading muon η |
| `muon1_phi` | Leading muon φ |
| `muon2_pt` | Subleading muon pT [GeV] |
| `muon2_eta` | Subleading muon η |
| `muon2_phi` | Subleading muon φ |

### Plot output

All figures saved as `.png` (150 dpi) under `hists/`:

```bash
hists/
├── combined/{cat}/         # Stacked MC + data, one PNG per (group, var)
├── comparison/{ds}/        # Preselection vs full-selection, one PNG per (group, var)
└── cutflow/                # Bar charts: {ds}.png + combined.png
```
