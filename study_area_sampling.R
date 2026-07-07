# Sampling locations for study areas

library(terra)
library(samplekmeans)
library(tidyverse)
library(dplyr)
library(tidyterra)
library(vctrs)
library(rcartocolor)

dir_data <- "C:/Users/anbm/OneDrive - SEGES Innovation PS/UDKIK/Data/"

study_areas <- paste0(
  dir_data,
  "Study_areas/Study_areas.shp"
) |>
  vect()

values(study_areas)

sampling_input <- paste0(
  dir_data,
  c(
    "Sampling_input/dhm2015_terraen_10m.tif",
    "Sampling_input/peat_probability_2025_resample.tif",
    "Sampling_input/vdtochn.tif"
  )
) |>
  rast()

study_area_idx <- 5

plot(study_areas[study_area_idx, ])

sampling_input_field <- crop(
  sampling_input,
  study_areas[study_area_idx, ]
) |>
  mask(
    study_areas[study_area_idx, ],
    touches = FALSE
  )

sampling_input_field <- ifel(
  sum(is.na(sampling_input_field)) == 0,
  sampling_input_field,
  NA
)

plot(
  sampling_input_field,
  nr = 1
)

# Transform input to percentiles to create clusters of roughly the same size.

sampling_input_pctile <- sapp(
  sampling_input_field,
  function(x) {
    x_ecdf <- ecdf(values(x))
    app(x, x_ecdf)
  }
)

field_means <- global(
  sampling_input_pctile,
  "mean",
  na.rm = TRUE
  ) |>
  unlist()
field_sds <- global(
  sampling_input_pctile,
  "sd",
  na.rm = TRUE
  ) |>
  unlist()

sampling_input_pctile <- sampling_input_pctile |>
  clamp(
    lower = field_means - field_sds*3,
    upper = field_means + field_sds*3
  )

plot(
  sampling_input_pctile,
  nr = 1
)

weights_field <- ifel(
  is.na(sum(sampling_input_field)),
  NA,
  1
) |>
  focal(
    w = 3,
    "sum",
    na.rm = TRUE,
    na.policy = "omit"
    )

candidates_field <- ifel(
  is.na(sum(sampling_input_field)),
  NA,
  1
) |>
  focal(w = 3, "mean") |>
  mask(
    study_areas[study_area_idx, ] |> as.lines(),
    inverse = TRUE
  )

plot(candidates_field)
plot(study_areas[study_area_idx,], add = TRUE)

set.seed(124)

myclusters <- sample_kmeans(
  input = sampling_input_pctile,
  clusters = round(study_areas[study_area_idx, ]$Shape_Area / 10000),
  use_xy = TRUE,
  sp_pts = TRUE,
  xy_weight = 2,
  candidates = candidates_field,
  weights = weights_field
)

plot(
  as.factor(myclusters$clusters)
)
points(
  myclusters$points,
  pch = 21,
  bg = "white"
)
plot(study_areas[study_area_idx,], add = TRUE)
text(
  myclusters$points,
  myclusters$points$ID,
  cex = 0.7,
  col = "black",
  pos = 3,
  hc = "white",
  hw = 0.1,
  halo = TRUE
)

plot(myclusters$distances)
points(
  myclusters$points,
  pch = 21,
  bg = "white"
)

cluster_areas <- table(myclusters$clusters |> as.data.frame())

cluster_areas

grid_samples <- spatSample(
  sampling_input_pctile,
  size = 20,
  method = "regular",
  as.points = TRUE
)

grid_samples <- sampling_input_pctile |>
  crds() |>
  (`-`)(5) |>
  (`/`)(30) |>
  round() |>
  (`*`)(30) |>
  (`+`)(5) |>
  as.data.frame() |>
  distinct() |>
  as.matrix() |>
  vect() |>
  mask(
    as.polygons(candidates_field)
  )

grid_samples

plot(candidates_field)
plot(study_areas[study_area_idx,], add = TRUE)
plot(
  grid_samples,
  add = TRUE
)
points(
  myclusters$points,
  pch = 21,
  bg = "white",
  add = TRUE
)


# END
