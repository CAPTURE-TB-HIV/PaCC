# TB PaCC
Shiny tool for TB patient costs

## About TB PaCC 

The TB Patient Cost Calculator (TB PaCC)  estimates country-level TB patient costs for 129 low- and middle income countries. We chose the name TB PaCC to represent the significant planning and resources required among persons undergoing TB treatment. Furthermore, the economic burden of TB treatment can place financial strain on households, like carrying a heavy pack throughout a journey.

The underlying model is informed from data collected from 104 studies of TB patient costs, representing more than 55,000 persons treated for TB. The model uses routinely reported country- and year- matched data to estimate the cost of TB, stratified by three cost types (direct medical, direct non-medical, and indirect) and across the phases of the TB care cascade (pre-treatment, intensive treatment, and continuation treatment). All estimates represent TB patient costs in 2024 in 2024 USD.

Please note: TB patient cost data varies substantially across settings and study design. As such, our estimates should be used as a guide point for analyses.

## Repository layout

- `app/`: the Shiny app (`app.R`, data in `utilities/`, Markdown pages in `inst/md/`, images in `www/`).
- `docs/`: the static website made by [shinylive](https://posit-dev.github.io/r-shinylive/). It runs the app in the browser with no server. Do not edit it by hand; re-export it instead.
- `scripts/`: helper scripts.

## Running the app locally

```r
shiny::runApp("app")
```

## Exporting the website

`scripts/export.R` rebuilds `docs/` from `app/`. Run it from the repo root after changing the app:

```sh
Rscript scripts/export.R
```

Or, in RStudio with `PaCC.Rproj` open:

```r
source("scripts/export.R")
```

The script:

1. Checks that it is run from the repo root, and that shinylive 0.5.0 or newer is installed (`install.packages("shinylive")`).
2. Stops if the app calls `webr::` (for example `webr::install()`). shinylive would treat that as a need for the CRAN package `webr`, which is unrelated to webR and pulls in more than 100 other packages. Load packages with `library()` instead.
3. Deletes `docs/`, so packages the app no longer uses are removed.
4. Runs `shinylive::export("app", "docs", assets_version = "0.10.16")`. These assets show loading text under the spinner, so a slow load does not look broken.
5. Prints the size of `app.json` and of the R packages. Visitors download both before the app runs, so keep them small. Every package the app loads, and every package those packages depend on, is downloaded.

Keep `docs/` outside `app/`. If the website folder is inside the app folder, shinylive can copy the whole website into `app.json`.


