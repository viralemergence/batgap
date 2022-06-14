# batgap
Gap analysis of bat coronaviruses

To create datacov.csv:
1) Download raw-airtable-preprocessing.csv from directory 00-preprocessing_and_cleaning
2) Download code_to_generate_datacov, Upham phylogeny ("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"), and taxonomy ("taxonomy_mamPhy_5911species.csv") from directory 01-generate_datacov
3) Run code_to_generate_datacov *requires code_to_generate_datacov csv, Upham phylogeny, and taxonomy*

code_to_generate_datacov cleans the dataset by standardizing bat species names against phylogeny, converting available latitude and longitude information to decimals, and creating several new columns used in the analyses (tissue_simplified, summer/fall/winter/spring sampling, euthanasia, and method_specific_simplified).

datacov.csv is stored in a separate repository: https://github.com/viralemergence/datacov
It includes 111 studies, with pooled bat species in same genus and rows not used in either prevalence or binary analyses

Directory https://github.com/viralemergence/batgap/tree/master/02-generate_datasets_from_datacov contains the script (script_to_generate_datasets_for_analyses) to produce the following datasets from datacov.csv, which can all be found here https://github.com/viralemergence/batgap/tree/master/02a-datasets_for_analyses:
1) set_infection_prevalence.csv contains rows to be used in https://github.com/viralemergence/batgap/blob/master/03-generate_and_run_REML_models/code_for_meta-analysis-models (pooled-coronavirus genera prevalence analyses)
3) set_infection_prevalence_alphaonly.csv contains rows to be used in (https://github.com/viralemergence/batgap/blob/master/03-generate_and_run_REML_models/code_for_meta-analysis-models (alphacoronavirus-genus only prevalence analyses))
4) set_infection_prevalence_betaonly.csv contains rows to be used in (https://github.com/viralemergence/batgap/blob/master/03-generate_and_run_REML_models/code_for_meta-analysis-models (betacoronavirus-genus only prevalence analyses)
5) set_other.csv contains rows used in ***pooled-coronavirus genera analyses not based on prevalence proportion positive***
6) set_other_alphaonly.csv contains rows used in ***alphacoronavirus-only analyses not based on prevalence proportion positive***
7) set_other_betaonly.csv contains rows used in ***betacoronavirus-only analyses not based on prevalence proportion positive***

Code to run phylogenetically-controlled meta-analyses using the above datasets is here: https://github.com/viralemergence/batgap/blob/master/03-generate_and_run_REML_models/code_for_meta-analysis-models

https://github.com/viralemergence/batgap/blob/master/03a-generate-figure-3-coefficient-plot/code_for_figure_3 contains the script for generating the phylogenetically-controlled meta-analysis models (same code as above)in addition to the script to generate figure 3 of the manuscript (coefficient plot).
