# Jaguar Habitat Projection Under Climate Change (SSP5-8.5)

This project models the future distribution of *Panthera onca* (jaguar) in Mexico and Central America using bioclimatic variables and ecological niche modeling in R.

## Objectives

- Clean and prepare GBIF occurrence records.
- Use WorldClim bioclimatic variables (current and 2050, SSP5-8.5).
- Train and evaluate a MaxNet model.
- Project range shifts and quantify loss/gain/stability of habitat.

## Project Structure

```
Panthera_onca_Climate_Modeling/
├── jaguar_model.R
├── bioclim_data/
├── output/
│   └── habitat_change_maps.pdf
├── figures/
│   └── suitability_maps.png
```

## Requirements

```r
install.packages(c("rgbif", "sf", "terra", "dplyr", "ENMeval", "maxnet", "ggplot2"))
```

## How to Run

1. Run `Panthera_onca_habitat_modeling.R` from RStudio.
2. Visuals and maps will be saved in `Mpas/` and `Report/`.
