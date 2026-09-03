# Environment and within-host priority effects jointly predict disease risk in a multi-pathogen system

This repository contains the data products and R code supporting a study of long-term foliar fungal disease dynamics in tall fescue at Widener Farm, North Carolina (2017–2024). The analyses combine random forests, Bayesian hierarchical models, Cox proportional-hazards models, and niche-overlap estimates to evaluate the relative contributions of environmental conditions and pathogen–pathogen priority effects to disease risk.

The three focal diseases are crown rust (*Puccinia coronata*), anthracnose (*Colletotrichum cereale*), and brown patch (*Rhizoctonia solani*). The primary processed analysis data contain 29,965 observations and span 26 May 2017 through 7 November 2024.

## Repository layout

```text
WidenerLongTermAnalysis/
├── data/                         Raw, cleaned, and intermediate data
├── scripts/                      Numbered analysis scripts and shared functions
├── results/                      Model objects and script-generated diagnostics
├── Figures/                      Assembled manuscript and supplementary figures
├── WidenerLongTermAnalysis.Rproj
└── README.md
```

Run scripts from the repository root, or open `WidenerLongTermAnalysis.Rproj`, because all paths are relative to that directory.

## Reproducing the analyses

The fastest route uses the supplied processed data and pre-fitted models:

```text
data/Env+DiseaseData.Rdata ──┬──> 02_RandomForests.R
                             ├──> 04_NicheOverlap.R
                             └──> 05_BayesianModels.R

data/weather.pcs.Rdata + raw 2018 survey/treatment data
                             └──> 03_SurvivalAnalyses.R

2024 on-site station + DURH weather exports
                             └──> 06_GroundtruthWeather.R
                                  (standalone validation only)
```

To rebuild the processed disease–environment data, run `scripts/01_DataProcessing.R` first. Scripts `02`–`05` can then be run separately. `06_GroundtruthWeather.R` is not part of the analysis pipeline: it only validates the DURH weather record against a weather station installed at the Widener field site in 2024.

Model fitting in `02_RandomForests.R` and `05_BayesianModels.R` is computationally intensive. Both scripts contain `load()` calls for supplied fitted objects, but users who want to skip fitting must run the relevant sections selectively rather than source either script uninterrupted from top to bottom.

## Required intermediate products

All intermediate files required by the downstream scripts are present and load successfully.

| Product | Created by | Used by | Purpose |
|---|---|---|---|
| `data/Env+DiseaseData.Rdata` (`rf.data`) | `01_DataProcessing.R` | `02`, `04`, `05` | Primary merged disease and 7-day environmental dataset |
| `data/weather.pcs.Rdata` (`pc.roll`) | `01_DataProcessing.R` | `03` | Daily environmental PC scores and 7-day rolling summaries |
| `results/crown_rust_rf_mtry4_upsample.Rdata` | `02_RandomForests.R` | `02` | Fitted crown-rust random forest |
| `results/bp_rf_mtry12_roll7_upsample.Rdata` | `02_RandomForests.R` | `02` | Fitted brown-patch random forest |
| `results/w7_anth_rf_mtry12_upsample.Rdata` | `02_RandomForests.R` | `02` | Fitted anthracnose random forest |
| `results/bigdata.month.brm.rust.Rdata` | `05_BayesianModels.R` | `05` | Fitted crown-rust Bayesian model |
| `results/bigdata.month.brm.anth.Rdata` | `05_BayesianModels.R` | `05` | Fitted anthracnose Bayesian model |
| `results/bigdata.month.brm.bp.Rdata` | `05_BayesianModels.R` | `05` | Fitted brown-patch Bayesian model |

`data/Env+DiseaseData.csv` is supplied for inspection outside R. It includes the row-number column written by base R; the `.Rdata` file is the authoritative analysis input.

## Data sources and processing

`01_DataProcessing.R` combines:

- cleaned long-term plot surveys from `data/Widener_2021_22/` and the top level of `data/Widener (2022-24)/`;
- 2017–2019 disease surveys and treatment assignments from `data/Rita_longitudinal/`;
- control-plot observations from `data/2024PlotSurveys/PlotSummary2024_CNTL_ONLY.csv`; and
- Durham Wastewater Treatment Plant (DURH) weather exports in `data/EnviroData_DURH/`.

The script harmonizes disease records at the plant level, retains control observations where treatments were applied, converts weather variables to metric units, calculates four environmental principal components, produces 7-day rolling environmental summaries, calculates lagged plot-level disease prevalence, and joins disease and weather data.

It writes `data/MasterDataset.csv`, `data/weather.pcs.Rdata`, `data/Env+DiseaseData.Rdata`, `data/Env+DiseaseData.csv`, `data/long2018weather.Rdata`, and weather-loading and prevalence figures.

The survival analysis independently joins `data/weather.pcs.Rdata` to the 2018 longitudinal survey and treatment table. It does not load `data/long2018weather.Rdata` or `data/byleaf2018.Rdata`.

## Analysis scripts

### `01_DataProcessing.R`

Cleans and harmonizes disease surveys, processes DURH weather observations, performs the environmental PCA, calculates rolling summaries and lagged disease prevalence, and creates the main analysis dataset.

### `02_RandomForests.R`

Fits one temporally structured random-forest classifier per disease. Chronological train/test splits and sliding resamples reduce temporal leakage, and upsampling addresses class imbalance. Outputs include fitted objects, variable-importance plots, predictor-correlation plots, and accumulated local effect (ALE) plots.

### `03_SurvivalAnalyses.R`

Fits clustered Cox proportional-hazards models to 2018 longitudinal observations using time-varying environmental PCs and prior infection by the other pathogens. It writes one forest plot per focal disease.

### `04_NicheOverlap.R`

Uses `nicheROVER` to estimate disease-specific niches in four-dimensional environmental PC space and posterior pairwise overlap probabilities.

### `05_BayesianModels.R`

Fits plot-level hierarchical logistic models with lagged disease state, month, year, environmental PCs, co-occurring pathogens, and pathogen-by-environment interactions. Outputs include model objects, posterior predictive checks, coefficient summaries, log-odds plots, and conditional-effect plots.

### `06_GroundtruthWeather.R`

Performs a standalone validation of the off-site DURH weather record using measurements from a weather station installed at the Widener field site in 2024. The comparison covers temperature, humidity, wind, and precipitation. This script is not otherwise relevant to the analyses: its data and output do not feed into the disease models or any other analysis script.

### `helper_functions.R`

Contains shared random-forest and posterior-processing utilities, including ALE plotting, posterior extraction, class balancing, and year-based train/test splitting.

## Main outputs

- `Figures/` contains assembled main-text and supplementary figures, including the current files and editable PowerPoint assemblies.
- `results/` contains fitted model objects, diagnostic and component plots, ALE plots, Bayesian conditional effects, survival plots, niche overlap, and weather-validation output.
- `results/*_Brms_Summmary.csv` contains the Bayesian fixed-effect summaries. `05_BayesianModels.R` currently writes these to the repository root, whereas the supplied copies have been organized under `results/`.

## Software requirements

```r
install.packages(c(
  "ape", "caret", "contsurvplot", "cowplot", "ggpubr", "iml",
  "janitor", "nicheROVER", "patchwork", "pheatmap", "posterior",
  "pracma", "randomForest", "RColorBrewer", "readxl", "survival",
  "survminer", "themis", "tidybayes", "tidymodels", "tidyverse",
  "vegan", "vip", "zoo"
))
```

Bayesian models additionally require `brms`, `cmdstanr`, and a working CmdStan installation. See the [CmdStanR installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html). Package versions are not currently pinned with `renv` or another lockfile.

