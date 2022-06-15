### Cohen _et al._, Sampling strategies and pre-pandemic surveillance gaps for bat coronaviruses

This repository contains data and code to reproduce all analyses of bat coronavirus data. A public-facing version of this dataset with simplified fields is available in the (datacov)[https://github.com/viralemergence/datacov] repository.

```
DIRECTORY

data.csv - the basic dataset

01_data processing
|
├─ code_to_generate_data.csv.R - creates data.csv using three other files in the directory
|
| # components of raw data that go into data.csv
├─ raw-airtable-preprocessing.csv - uncleaned data
├─ MamPhy_...._target.tre - phylogenetic tree for processing
└─ taxonomy_mamPhy_5911species.csv - mammal taxonomy metadata

02_dissolve data
|
├─ script_to_generate_datasets_for_analyses.R - breaks down data.csv into smaller files for sub-analyses
|
| # These are used for prevalence analyses in code_for_meta-analysis_models.R
├─ set_infection_prevalence.csv  - all coronavirus records
├─ set_infection_prevalence_alphaonly.csv  - alphacoronavirus only
├─ set_infection_prevalence_betaonly.csv -  betacoronavirus only
|
| # These are used in geophylo analysis.R and descriptive stats not dependent on proportion positive
├─ set_other.csv  - all coronavirus records
├─ set_other_alphaonly.csv - alphacoronavirus only
└─ set_other_betaonly.csv -  betacoronavirus only 

03_analysis
|
├─ code_for_meta-analysis_models.R - runs the meta-analysis
├─ code_for_figure3.R - generates visualizations from the meta-analysis
|
├─ geophylo analysis.R - runs phylofactor and generates the maps and tree figures
└─ georegion.csv - geospatial metadata used for regional analyses 

04_outputs
|
| # some additional files that become figures or supplement
├─ Figure 1.png
├─ Figure 2.png
├─ Figure 3.png
├─ Table S3.csv
├─ Table S4.csv
├─ Table S5.csv
├─ Table S8.csv
└─ Table S9.csv
```


