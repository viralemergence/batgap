# batgap
Gap analysis of bat coronaviruses

To get datasets:
1) Download airtable (prevalence tab) as a csv file
2) Download Upham phylogeny ("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre") and taxonomy ("taxonomy_mamPhy_5911species.csv") from data repository
3) Run airtable_to_datasets.R (in scripts repository) *requires airtable csv, Upham phylogeny, and taxonomy*

airtable_to_datasets.R cleans the dataset by standardizing bat species names, converting available latitude and longitude information to decimals, and creating several new columns used in the analyses (tissue_simplified, summer/fall/winter/spring sampling, euthanasia, and method_specific_simplified).
The script produces the following datasets, which can all be found in the data repository:
1) data_for_reader.csv ***not used in analyses***, includes pooled bat species in same genus and rows not used in either prevalence or binary analyses
2) set_infection_prevalence.csv contains back transform, bat family, and rows to be used in ***pooled-coronavirus genera prevalence analyses***
3) set_infection_prevalence_alphaonly.csv contains back transform, bat family, and rows to be used in ***alphacoronavirus genus-only prevalence analyses***
4) set_infection_prevalence_betaonly.csv contains back transform, bat family, and rows to be used in ***betacoronavirus genus-only prevalence analyses***
5) set_other_alphaonly.csv contains rows used in ***alphacoronavirus-only analyses not based on prevalence proportion positive***
6) set_other_betaonly.csv contains rows used in ***betacoronavirus-only analyses not based on prevalence proportion positive***
