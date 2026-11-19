News
-----------------------------

**Manuscript pre-print live!** 

The manuscript describing the methods used to estimate the TB household costs (previously referred to as TB patient costs) is live on MedRxiv. Read it [here](https://www.medrxiv.org/content/10.64898/2026.10.19.26365552v1)!

Data sources
-----------------------------

**Out-of-pocket percentage of current health expenditures**

World Health Organization (WHO) reported measure of the share of revenues from out-of-pocket payments as percentage of current health expenditure indicates how much is funded directly by households out-of-pocket expenditure on health. High out-of-pocket payments are associated with catastrophic and impoverishing household spending. Reported on an 0-100 scale. 

Read more [here](https://www.who.int/data/gho/data/indicators/indicator-details/GHO/out-of-pocket-expenditure-as-percentage-of-current-health-expenditure-(che)-(-)).

**Median household annual income** 

Individual median annual income estimated by scaling of Gross National Income (GNI) using the GINI coefficient. These values are scaled to household level using United Nations' Population Division's Household Composition data. 

**Tuberculosis case detection rate**

WHO estimated measure of case detection calculated by dividing the annual reported TB notifications by the WHO estimate of annual incidence. 

Read more [here](https://data.who.int/indicators/i/3805B1E/9A706FD).

**Universal Healthcare Coverage Index**

WHO index measure representing the coverage of essential health services (defined as the average coverage of essential services based on tracer interventions that include reproductive, maternal, newborn and child health, infectious diseases, non-communicable diseases and service capacity and access, among the general and the most disadvantaged population). Reported on a 0-100 scale. 

Read more [here](https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/data).


Predictive model overview
-----------------------------

We constructed a Bayesian meta-regression model to estimate the costs associated with direct medical (DM), direct non-medical (DNM), and indirect costs for each phase of the TB care cascade: pre-treatment, intensive phase treatment, and continuation phase treatment. This created nine cost strata. For each stratum, we defined a regression equation with a log link function, predicting the log of mean costs as a function of several covariates: an overall intercept, fraction of cohort receiving a RR/MDR-TB regimen, survey year, country median income, country WHO TB case detection rate, and country out-of-pocket expenses as a percentage of current health expenditures (OOP percentage of CHE). 

We modelled the cost data using a gamma likelihood, to account for key features of the data (continuous, non-negative, right-skewed). While some cost observations represented a single cost stratum, others represented the sum of several strata (e.g., total indirect costs, summed across treatment phases). Therefore, for each observation we set the mean of the data likelihood equal to the sum of the strata included in that observation, and estimated a common dispersion parameter accounting for observation-level variance.

Using the fitted model we predicted per-episode costs, stratified by RR/MDR-TB status, for all LMICs as defined by the World Bank in 2025 (Democratic People’s Republic of Korea, Kosovo, and Palestine were omitted due missing covariate data), using 2024 (the most recently available) values for each covariate. 

Summary of underlying systematic review
-----------------------------

We undertook a systematic review of published TB patient cost studies ([PROSPERO registration CRD42025159307](https://www.crd.york.ac.uk/PROSPERO/view/CRD42025159307)), including published national and sub-national studies conducted in low and middle-income countries over January 1, 2015 and December 31, 2024. Our search identified 2,430 publications for potential inclusion, of which 208 had full text reviewed. We included 113 of these publications, which described 119 unique surveys, and added 4 additional patient costs surveys conducted with the support of WHO, which were not captured by our search strategy. In total, we extracted data for 122 unique surveys. These surveys captured costs collected in 42 countries across six WHO-defined world regions: Africa region (N=15), Americas region (N=3), Eastern Mediterranean region (N=4), Europe region (N=2), South-East Asia region (N=8), and Western Pacific region (N=10).