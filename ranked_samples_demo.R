# Ranked samples for UDKIK

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

study_area_idx <- 3

plot(study_areas[study_area_idx, ])

sampling_input <- paste0(
  dir_data,
  c(
    "Sampling_input/dhm2015_terraen_10m.tif",
    "Sampling_input/peat_probability_2025_resample.tif",
    "Sampling_input/vdtochn.tif"
  )
) |>
  rast()


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

candidates_field <- ifel(
  is.na(sum(sampling_input_field)),
  NA,
  1
) %>%
  mask(
    study_areas[5, ] |> as.lines(),
    inverse = TRUE
  )

set.seed(123)

myclusters <- sample_kmeans(
  input = sampling_input_field,
  clusters = round(study_areas[5, ]$Shape_Area / 10000),
  use_xy = TRUE,
  sp_pts = TRUE,
  # layer_weights = feature_weights_field,
  xy_weight = 2,
  candidates = candidates_field,
  pca = TRUE
)

plot(myclusters$clusters)
points(
  myclusters$points,
  pch = 21,
  bg = "white"
)

plot(myclusters$distances)
points(
  myclusters$points,
  pch = 21,
  bg = "white"
)

cluster_areas <- table(myclusters$clusters |> as.data.frame())

extra_pts <- list()
list_index <- 1

for (i in seq_len(nrow(myclusters$points))) {
  if(
    cluster_areas[i] > 30
  ) {
    rast_i <- mask(
      sampling_input_field,
      mask = myclusters$clusters,
      maskvalue = i,
      inverse = TRUE
    ) |>
      trim()

    candidates_i <- focal(
      rast_i,
      fun = function(x) {
        if (sum(is.na(x) > 0)) {
          return(NA)
        } else {
          return(1)
        }
      },
    )

    set.seed(123)

    clusters_i <- sample_kmeans(
      input = rast_i,
      clusters = round(cluster_areas[i] / 20),
      use_xy = TRUE,
      # weights = weights_i,
      xy_weight = 2,
      sp_pts = TRUE,
      candidates = candidates_i
    )

    extra_pts[[list_index]] <- clusters_i$points %>%
      mutate(
        cluster = i
      )

    list_index <- list_index + 1
  }
}


plot(clusters_i$clusters)
plot(rast_i)
points(
  clusters_i$points,
  pch = 21,
  bg = "white"
)

extra_pts <- do.call(rbind, extra_pts)

plot(sampling_input_field[[1]])
plot(myclusters$points, pch = 21, bg = "red", add = TRUE)
plot(extra_pts, pch = 21, bg = "white", add = TRUE)

all_pts <- bind_spat_rows(
  myclusters$points,
  extra_pts
) %>%
  mutate(
    cluster = case_when(
      is.na(cluster) ~ ID,
      .default = cluster
    )
  )

all_pts <- terra::extract(
  myclusters$distances,
  all_pts,
  bind = TRUE,
  ID = FALSE
) %>%
  arrange(
    cluster, distance
  ) %>%
  group_by(
    cluster
  ) %>%
  mutate(
    rank = rank(distance, ties.method = "first"),
    label = paste0(cluster, letters[rank])
  ) %>%
  ungroup()

plot(sampling_input_field[[1]])
plot(all_pts, pch = 21, bg = "white", add = TRUE)
plot(
    as.polygons(
      myclusters$clusters
      ),
  add = TRUE,
  lwd = 1)
text(
  all_pts,
  all_pts$label,
  cex = 0.7,
  col = "black",
  pos = 3
)

# Implement candidates

# END
