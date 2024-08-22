# Bassing_et_al_Predator-Prey_Movement
R code associated with the publication Bassing, S. B., L. Satterfield, T. R. Ganz, M. DeVivo, B. N. Kertson, T. Roussin, A. J. Wirsing, and B. Gardner. 2024. Predator-prey space-use and landscape features influence animal movement behaviors in a large-mammal community. Ecology.

This repository archives the code used for analyses described in the manuscript. Data already formatted for resource selection and movement analyses are archived with Dryad at https://doi.org/10.5061/dryad.kh1893292. Raw data (coordinates of GPS-collar relocations) are considered sensitive and not publicly available. Please contact the Wildlife Chief Scientist of the Washington Department of Fish and Wildlife at (360) 902-2515 if you are interested in these data. 

Shorthand used throughout scripts and data:
-------------------------------------------
Species codes:
-coug, COUG = Cougar
-elk, ELK = Elk
-md, MD = Mule deer
-wtd, WTD = White-tailed deer
-wolf, WOLF = Wolf

Study area codes:
-NE = Northeast
-OK = Okanogan

Season codes:
-Summer18 = summer 2018
-Winter1819 = winter 2018-2019
-Summer19 = summer 2019
-Winter1920 = winter 2019-2020
-Summer20 = summer 2020
-Winter2021 = winter 2020-2021

Year codes:
-Year1 = 2018-2019
-Year2 = 2019-2020
-Year3 = 2020-2021

Covariate codes:
-Elev: Elevation (m) of observation location
-Slope: Slope (degrees) of terrain at observation location
-RoadDen: Total road length/1 km-sq at observation location
-Dist2Water: Distance (m) to nearest water
-HumanMod: Percentage of human modification to the landscape
-CanopyCover: Percentage of tree cover
-Dist2Edge: Distance (m) to nearest forested to non-forested habitat edge
-PercForestMix: Percentage of forested habitat within 250 m of observation location
-PercXGrass: Percentage of xeric grassland habitat within 250 m of observation location
-PercXShrub: Percentage of xeric shrubland habitat within 250 m of observation location
-Landcover: Numerical value representing landcover classification from Cascadia Partner Forum TerrAdapt:Cascadia tool (30m resolution)
-Landcover_type: Landcover classification label
