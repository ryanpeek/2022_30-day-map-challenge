# spatial packages to play with


# Packages ----------------------------------------------------------------

# standard
library(terra)
library(sf)
library(dplyr)

# new
library(supercells) # https://jakubnowosad.com/supercells/
library(osmdata) # https://docs.ropensci.org/osmdata/
library(CAST) # https://hannameyer.github.io/CAST/
library(mlr3) # https://mlr3.mlr-org.com
library(gstat) # http://r-spatial.github.io/gstat/
library(landscapemetrics) # https://r-spatialecology.github.io/landscapemetrics/
library(landscapetools) # https://ropensci.github.io/landscapetools/
library(RStoolbox) # http://bleutner.github.io/RStoolbox/rstbx-docu/RStoolbox.html
library(sfnetworks) # https://luukvdmeer.github.io/sfnetworks/
library(spdep) # https://r-spatial.github.io/spdep/
library(rmapshaper)
library(motif) # https://jakubnowosad.com/motif/
# hypertidy: https://github.com/hypertidy
#remotes::install_github("mdsumner/lazysf")
library(lazysf)
library(ncmeta) # https://github.com/hypertidy/ncmeta
library(geodist) # https://hypertidy.github.io/geodist/
library(tabularaster) # https://github.com/hypertidy/tabularaster

# Examples ----------------------------------------------------------------

## motif -----------------------------------

library(motif)
library(stars)
library(sf)
landcover <- read_stars(system.file("raster/landcover2015.tif", package = "motif"))
plot(landcover)

# find signature of forest
landcover_comp <- lsp_signature(landcover,
                               type = "composition",
                               threshold = 1,
                               normalization = "none")
landcover_comp
landcover_comp$signature
# change cell size
# local landscape will consist of 200 by 200 cells. In this example, each cell has a resolution of 300 by 300 meters, therefore a local landscape will be 60000 by 60000 meters

# use a co-occurrence matrix ("coma"): looks at its value, looks at the values of its neighbors, counts how many neighbors of each class our central cell has.
landcover_coma <- lsp_signature(landcover, type = "coma", window = 200)
landcover_coma$signature[[1]]

# sf layer
ecoregions <- read_sf(system.file("vector/ecoregions.gpkg", package = "motif"))
plot(ecoregions["ECO_NAME"])
landcover_coma_e <- lsp_signature(landcover, type = "coma", window = ecoregions["id"])
View(landcover_coma_e)
ecoregions$ECO_NAME
landcover_coma_e$signature[[1]]

# compare change:
# https://jakubnowosad.com/motif/articles/v4_compare.html


## osmdata ---------------------------------------------------------------------

library(dplyr)
library(osmdata)
cnty <- tigris::counties("CA") %>%
  filter(NAME=="El Dorado")

x <- opq(bbox = c(-0.27, 51.47, -0.20, 51.50)) %>% # Chiswick Eyot in London, U.K.
  add_osm_feature(key = 'name', value = "Thames", value_exact = FALSE) %>%
  osmdata_sf()

rivs <- opq(st_bbox(cnty)) %>%
  add_osm_feature(key = 'river',
                  key_exact = FALSE) %>%
  osmdata_sf()


## landscapemetrics --------------------------------------------------------


library(landscapemetrics)
library(landscapetools)

# demo
show_landscape(landscape, discrete = TRUE)

# euclid nearest neighbor distance on patch level
lsm_p_enn(landscape)

# total area and class edge
lsm_l_ta(landscape)
lsm_c_te(landscape)

# calculate all metrics on patch level
metrics_df <- calculate_lsm(landscape, level = "patch")

# check_landscape(landscape)
table(metrics_df$metric)
show_cores(landscape)
show_lsm(landscape, what = "lsm_p_circle", class=c(2,3), labels=TRUE)

## lazysf ------------------------------------------------------------------

library(lazysf)
library(dplyr)

# gpkg
f <- system.file("gpkg/nc.gpkg", package = "sf", mustWork = TRUE)
lazysf(f)

# from file
f <- "~/Downloads/kontur_population_20220630.gpkg"
(fin <- lazysf(f))
fin %>% distinct(h3, .keep_all = T) %>% dim()

# from url
url <- "https://github.com/Nowosad/spData/raw/master/inst/shapes/NY8_bna_utm18.gpkg"
(x <- lazysf(url))

# geojson
geojson <- file.path("https://raw.githubusercontent.com/SymbolixAU",
                     "geojsonsf/master/inst/examples/geo_melbourne.geojson")
lazysf(geojson)

## supercells --------------------------------------------------------------

library(supercells)
library(terra)
library(sf)

s = rast(system.file("ex/logo.tif", package = "terra"))
sc = supercells(s, 500, compactness = 50, transform = "to_LAB")
plot(sc)

vol = rast(system.file("raster/volcano.tif", package = "supercells"))
plot(vol)

vol_slic1 = supercells(vol, k = 100, compactness = 3)
plot(vol)
plot(st_geometry(vol_slic1), add = TRUE, lwd = 0.2)
