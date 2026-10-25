# KZN Unemployment & Deprivation Explorer

**STAT7006A Spatial Statistics – group R Shiny assignment (Wits, 2026)**

**Live dashboard:** https://lintlelenepa.shinyapps.io/kzn-unemployment/

## Research question
Is unemployment linked to deprivation in education and basic services across KwaZulu-Natal, and is this relationship spatial?

## Aim
To assess the association between unemployment and deprivation in education and basic services across the 51 local municipalities of KwaZulu-Natal, and to test for spatial dependence.

## Objectives
1. **Association** – maps, distributions and correlations between unemployment and deprivation indicators.
2. **Patterns** – (a) global spatial autocorrelation (Moran's I, analytical and Monte Carlo); (b) local patterns (LISA and Getis-Ord Gi*, with Benjamini–Hochberg adjustment).
3. **Spatial effects** – OLS with VIF and residual diagnostics, Lagrange Multiplier tests, spatial lag (SLM) and spatial error (SEM) models compared by AIC and residual Moran's I.

## Data
- Stats SA Census 2011 Municipal Fact Sheet (`data/census2011_kzn_municipal.csv`)
- Municipal Demarcation Board 2011 local municipality boundaries, KZN only (`data/kzn_boundaries.geojson`)
- Joined dataset used by the app: `data/municipalities.rds`

## How to run
1. Open `Unemployment_Shiny_App.Rproj` in RStudio.
2. Source `00_install_packages.R` once to install missing packages.
3. (Optional) Source `01_prepare_data.R` to rebuild `data/municipalities.rds`.
4. Open `app.R` and click **Run App**.
5. To publish an update, source `deploy.R`.

## Files
| File | Purpose |
|---|---|
| `app.R` | The Shiny dashboard |
| `00_install_packages.R` | Installs required packages |
| `01_prepare_data.R` | Joins census table to boundaries |
| `deploy.R` | Publishes to shinyapps.io |
| `data/` | Input and prepared data |

## Packages
shiny, bslib, sf, dplyr, spdep, spatialreg, leaflet, plotly, DT, rsconnect

## Group members
- Lintle Lenepa
- Sheofatso Padima
- Fanisa Mabunda
- Congress Mafukele
