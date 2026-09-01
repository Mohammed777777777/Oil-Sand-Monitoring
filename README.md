# Oil Sand Monitoring

Google Earth Engine (GEE) analysis of land-cover disturbance and wetland vegetation trends in the Athabasca Oil Sands Region, Alberta, using long-term Landsat time series (1985–2025).

## What's in this repo

### `CCDC_Implementation.ipynb`
Runs **CCDC (Continuous Change Detection and Classification)** across a 667-tile HUC8 watershed grid to detect structural breaks in Landsat spectral/index time series.

- Builds a harmonized Landsat surface-reflectance collection (1985–2025), with cloud/snow masking.
- Tiles the watershed into a fixed grid (target ~650 tiles, EPSG:32612) so each GEE export task stays within memory/time limits.
- Runs `ee.Algorithms.TemporalSegmentation.Ccdc` per tile with tuned parameters (min observations, chi-square probability, lambda, etc.) over 10 spectral/index bands (surface bands + NDVI, NBR, NDMI, EVI).
- Includes task submission, progress tracking, and automatic retry of failed tile exports (a full pipeline for running large batch jobs on GEE reliably).
- Produces a 39-band per-pixel summary mosaic (first break year, change count, change type, segment reliability, etc.) used downstream in the R analysis.
- A separate "Figure 3a" section pulls two example pixels (one vegetation-clearing, one oil-sands industrial) and reconstructs their actual CCDC-fitted NDVI trajectories for a publication figure.

### `CCDC_Implementation_2YearsConsistency.R`
Post-processing/validation script (Paper 2, Table T3 — "detection gradient") that consumes the 39-band CCDC mosaic exported above.

- Defines a locked detection hierarchy: **any break → loss-type break → reliable loss → temporally-consistent loss** (CCDC break year within ±2 years of a mapped Human Footprint Inventory (HFI) disturbance year).
- Evaluates this gradient separately for three HFI disturbance groups: vegetation clearing, oil-sands industrial footprint, and well pads.
- Rasterizes each HFI polygon group, computes per-group pixel counts/percentages at each detection stage, and cross-checks results against previously validated numbers (self-check with a 0.5 percentage-point tolerance).
- Exports full, paper-ready, and self-check CSV tables.

### `TimeSeries_Wetland_OSM.ipynb`
Complementary **LandTrendr**-based analysis focused on wetland vs. watershed vegetation trends.

- Loads OSM (Oil Sands Region) wetland, watershed, and HUC8 boundary assets from a GEE project.
- Builds annual NDVI medoid composites from the growing season (June–August) for 1985–2025, clipped to the study area (cross-sensor harmonization deliberately skipped — not needed for this comparison).
- Exports annual composites as GEE assets and verifies all 41 expected years exist before analysis.
- Computes and plots spatial mean NDVI per year separately for wetlands, full watersheds, and the watershed-minus-wetland difference geometry (to isolate the wetland signal from the surrounding catchment).
- Runs pixel-wise `LandTrendr` temporal segmentation on the annual NDVI stack to extract fitted/denoised vegetation trajectories, exports the fitted result as a multi-band GEE asset, and re-derives mean NDVI trends from the fitted output for comparison against the raw composites.

## Pipeline summary

```
Landsat 1985–2025 (GEE)
        │
        ├── CCDC_Implementation.ipynb ──► 39-band CCDC change mosaic
        │                                        │
        │                                        ▼
        │                     CCDC_Implementation_2YearsConsistency.R
        │                     (validates CCDC-detected disturbance against HFI records)
        │
        └── TimeSeries_Wetland_OSM.ipynb ──► annual NDVI composites + LandTrendr
                                             (wetland vs. watershed vegetation trend)
```

## Requirements

- Google Earth Engine account/project with access to the referenced assets (`projects/osmgee/...`, `projects/annual-ccdc-new/...`).
- Python: `earthengine-api`, `geemap`, `pandas`, `numpy`, `matplotlib`.
- R: `terra`.
