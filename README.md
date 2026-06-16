# HYgammaRunIII Analysis

CMS Run 3 analysis of **H/Z → μ⁺μ⁻γ** via an intermediate Υ(nS) state,  
using NanoAOD-like ROOT data processed with the [coffea](https://coffea-hep.readthedocs.io/) columnar framework.

## Project structure

```
├── run_analysis.py              # CLI entry point to process ROOT files
├── run_all.py                   # Batch runner
├── processor/
│   └── OniaAnalysisProcessor.py  # Event selection, cutflow, histogram production
├── schema/
│   └── OniaNanoSchema.py         # Custom NanoAOD schema (collections, mixins, cross-refs)
├── utils/
│   ├── figure.py                 # Plotting (thread-parallel, --workers N)
│   ├── weights.py                # Cross-section × BR weights for MC scaling
│   └── utils.py                  # Helpers (file discovery, histogram filling)
├── output/                       # Pickled accumulators (gitignored)
├── hists/                        # Generated plots (gitignored)
└── requirements.txt
```

## Dependencies

Install with conda / pip:

```bash
pip install -r requirements.txt
```

Requires Python ≥ 3.9, coffea ≥ 0.7, uproot ≥ 5, awkward ≥ 2, hist ≥ 2.

## Usage

### Process ROOT files

```bash
python run_analysis.py path/to/files --dataset mydata
python run_analysis.py path/to/files --dataset mymc --is_mc
python run_analysis.py path/to/files --dataset mydata --parallel --workers 8
```

Output: `output/<dataset>.pkl` — pickled accumulator with histograms and cutflow.

for processing all datasets use:
```bash
python run_all.py
```

altering the values of the variables `data_path` and `MC_path`

### Create plots

#### Validation plots:

```bash
python utils/figure.py --validation                      # 8 worker threads (default)
python utils/figure.py --workers 16 --validation         # custom parallelism
python utils/figure.py --validation --save-per-dataset   # individual dataset plots
```

Output: `hists/combined/`, `hists/comparison/`, `hists/cutflow/` — `.png` files.

## Output format

A pickled accumulator `{dataset}.pkl` is created as a dict with keys:

| Key        | Contents |
|------------|----------|
| `hists`    | Nested dict of `hist.Hist` objects (μ₁, μ₂, γ, Onia, X kinematics) |
| `cutflow`  | Event counts after each selection stage |
| `x_mass`   | `column_accumulator` of X (μμγ) candidate masses |
| `onia_mass` | `column_accumulator` of Onia (μμ) candidate masses |
| `dataset`  | Dataset name |
| `is_mc`    | `True` for MC, `False` for data |

And a ROOT file `{dataset}.root` is also created with the following structure:

```
events/                        # TTree with selected event values
├── x_mass                     # X (μμγ) candidate mass [GeV]
└── onia_mass                  # Onia (μμ) candidate mass [GeV]

metadata/                      # TTree with one entry (file-level metadata)
├── is_mc                      # 1 for MC, 0 for data
└── dataset                    # Dataset name (string)

muon_lead/                     # Leading muon histograms
muon_trail/                    # Trailing muon histograms
photon/                        # Photon histograms
onia/                          # Onia histograms (pₜ, η, φ, mass)
x/                             # X candidate histograms (pₜ, η, φ, mass, n)
└── {var}/{cat}                # Each variable × selection category
```

Selection categories: `all`, `trigger`, `good_muons`, `good_photon`, `upsilon_mass`, `x_mass`.
