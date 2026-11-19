# Export the Shiny app in `app/` to a static website in `docs/`.
#
# Run from the repo root:
#   Rscript scripts/export.R
# or, in RStudio with PaCC.Rproj open:
#   source("scripts/export.R")

app_dir <- "app"
out_dir <- "docs"

if (!file.exists(file.path(app_dir, "app.R"))) {
  stop(
    "Could not find `", file.path(app_dir, "app.R"), "`. ",
    "Run this script from the repo root (the folder with PaCC.Rproj)."
  )
}

if (
  !requireNamespace("shinylive", quietly = TRUE) ||
    utils::packageVersion("shinylive") < "0.5.0"
) {
  stop('Please install shinylive 0.5.0 or newer: install.packages("shinylive")')
}

# A `webr::` call (such as `webr::install()`) makes shinylive think the app
# needs the CRAN package `webr`. That package is unrelated to webR, and it
# pulls in more than 100 other packages that every visitor downloads.
if ("webr" %in% renv::dependencies(app_dir, quiet = TRUE)$Package) {
  stop(
    "Found a `webr::` call in `", app_dir, "/`. Please remove it.\n",
    "shinylive would download the unrelated CRAN package `webr` and more ",
    "than 100 other packages. Add packages with `library()` instead."
  )
}

# Start from an empty folder, so packages the app no longer uses are not kept
unlink(out_dir, recursive = TRUE)

# Assets 0.10.16 and newer show loading text under the spinner
# ("Downloading R...", "Loading packages and starting app..."),
# so a slow load does not look broken.
shinylive::export(app_dir, out_dir, assets_version = "0.10.16")

# Report what a visitor downloads before the app runs
size_mb <- function(paths) sum(file.size(paths)) / 1e6
pkg_files <- list.files(
  file.path(out_dir, "shinylive", "webr", "packages"),
  pattern = "\\.tgz$",
  recursive = TRUE,
  full.names = TRUE
)
message(sprintf(
  "\napp.json: %.1f MB\nR packages: %d packages, %.1f MB",
  size_mb(file.path(out_dir, "app.json")),
  length(pkg_files),
  size_mb(pkg_files)
))
