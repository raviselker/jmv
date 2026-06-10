# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

jmv is the suite of analyses bundled with jamovi (ANOVAs, t-tests, regression,
frequencies, factor analysis, etc.). It is both a CRAN R package and a jamovi
module, so every analysis must work when called from R directly and when driven
by the jamovi UI.

## Commands

Dependencies are managed with renv (activates via .Rprofile).

```sh
# Run a single test file (preferred during development)
Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/testanova.R")'

# Run the full test suite
Rscript -e 'devtools::load_all(); testthat::test_dir("tests/testthat")'
```

Compiling the jamovi module (regenerates `R/*.h.R` from the yaml definitions)
requires the jamovi-compiler: `jmc . --install-to jmv` (see package.json).

## Architecture

Each analysis `<name>` is defined across several files that belong together:

- `jamovi/<name>.a.yaml` — options (arguments) definition
- `jamovi/<name>.r.yaml` — results definition (tables, columns, images, visibility rules)
- `jamovi/<name>.u.yaml` + `jamovi/js/<name>.events.js` — jamovi UI layout and events
- `R/<name>.h.R` — **generated** Options/Base R6 classes; never edit these by hand
- `R/<name>.b.R` — the analysis body: an R6 class inheriting `<name>Base`

The body class implements two phases driven by the jmvcore framework:

- `.init()` — runs before data is available; builds the result skeleton
  (table rows/columns based on options and factor levels)
- `.run()` — validates data, computes, and populates the skeleton with
  `table$setRow(rowKey=..., values)`

Analyses reuse each other through inheritance (e.g. `anovaClass` inherits
`ancovaClass`, so ANOVA behaviour usually lives in `R/ancova.b.R`). Shared
helpers live in `R/utils.R`, `R/utilsanova.R`, `R/utilsreg.R`.

Conventions that span files:

- Factor levels are base64-encoded (`toB64`/`fromB64`) before model fitting so
  special characters survive R formulas; decode with `fromB64` for display.
- Data problems are reported with
  `jmvcore::reject(.("..."), code=exceptions$dataError, ...)`; rejections
  surface as regular R errors for R users (test with `expect_error`).
- Non-fatal warnings use `jmvcore::Notice` via `setAnalysisNotice()` in
  `R/utils.R`.
- All user-facing strings are wrapped in `.()` for translation. The `.()`
  lookup requires a `self` in scope, which is why package-level helpers take a
  `self` argument even when otherwise unused. Translations live in the
  `jamovi/i18n` submodule; catalogs are regenerated separately ("Updated
  i18n" commits), so flag new strings when adding them.

## Code style

There is no formatter; match the surrounding code. The R code uses 4-space
indentation, no spaces around `=` in named arguments (`setRow(rowNo=i, ...)`),
and a spaced negation style (`if ( ! is.null(x))`). Tests use
GIVEN/WHEN/THEN comments and the spaced `arg = value` style.

## Committing

- **Small logical commits:** Break changes into small, focused commits with a
  single purpose.
- **Commit title:** Single sentence in imperative mood, max 50 characters, no
  trailing dot. No type prefixes (e.g., `feat:`, `fix:`).
- **Description:** Optional. Only to clarify functional choices (the "what"
  and "why"). Do not explain the "how" or anything already evident in the
  code. Max 72 characters per line.
- **No AI mentions:** Never mention AI assistants or tools in commit messages.
- **Propose first:** Always propose a draft commit message before committing.
