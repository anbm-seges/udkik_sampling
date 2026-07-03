# Ranked samples for UDKIK

library(terra)
library(samplekmeans)
library(tidyverse)
library(dplyr)
library(tidyterra)
library(vctrs)

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

study_areas <- "C:/Users/anbm/OneDrive - SEGES Innovation PS/UDKIK/Data/Study_areas/Study_areas.shp" |> vect()

values(study_areas)

plot(study_areas[5, ])

sampling_input <- c(
  "C:/Users/anbm/OneDrive - SEGES Innovation PS/UDKIK/Data/Sampling_input/dhm2015_terraen_10m.tif",
  "C:/Users/anbm/OneDrive - SEGES Innovation PS/UDKIK/Data/Sampling_input/peat_probability_2025_smooth.tif"
) |>
  rast()

# global_sds <- global(sampling_input, fun = "sd", na.rm = TRUE)
# #                              sd
# # dhm2015_terraen_10m   24.602547
# # peat_probability_2025 16.321100
# # vdtochn                5.671995
#
# field_sds <- global(sampling_input_field, fun = "sd", na.rm = TRUE)
# #                              sd
# # dhm2015_terraen_10m   0.9847146
# # peat_probability_2025 7.7782040
# # vdtochn               0.2480967

# feature_weights_field <- unlist(field_sds / global_sds)
#
# xy_weights_field <- sqrt(study_areas[5, ]$Shape_Area) / sqrt(45000*10^6)

sampling_input_field <- crop(sampling_input, study_areas[5, ]) %>%
  mask(study_areas[5, ])

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

extra_pts <- list()

for (i in seq_len(nrow(myclusters$points))) {
  rast_i <- mask(
    r,
    mask = myclusters$clusters,
    maskvalue = i,
    inverse = TRUE
  )

  weights_i <- mask(
    myclusters$distances,
    mask = myclusters$clusters,
    maskvalue = i,
    inverse = TRUE
  ) |>
    app(
      function(x) {
        1 / (x + 1)
      }
    )

  set.seed(123)

  clusters_i <- sample_kmeans(
    input = rast_i,
    clusters = 5,
    use_xy = TRUE,
    # weights = weights_i,
    xy_weight = 2,
    sp_pts = TRUE
  )

  extra_pts[[i]] <- clusters_i$points %>%
    mutate(
      cluster = i
    )
}


plot(clusters_i$clusters)
plot(rast_i)
plot(weights_i)
points(
  clusters_i$points,
  pch = 21,
  bg = "white"
)


extra_pts <- do.call(rbind, extra_pts)

plot(r)
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

plot(r)
plot(all_pts, pch = 21, bg = "white", add = TRUE)
plot(
  buffer(
    as.polygons(
      myclusters$clusters
      ),
    width = -75),
  add = TRUE,
  border = c("red", "green", "blue"),
  lwd = 1.5)
text(
  all_pts,
  all_pts$label,
  cex = 0.7,
  col = "black",
  pos = 3
)

# Implement candidates

# END
