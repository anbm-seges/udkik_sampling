# Sampling locations for study areas

library(terra)
library(samplekmeans)
library(tidyverse)
library(dplyr)
library(tidyterra)
library(vctrs)
library(rcartocolor)

source("function_colors.R")

dir_data <- "C:/Users/anbm/OneDrive - SEGES Innovation PS/UDKIK/Data/"
dir_sensors <- paste0(dir_data, "/Sensors data/TIFF Files/") |>
  list.dirs(
    recursive = FALSE
  )

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

list_clusters_10m <- list()
list_clusters_sensor <- list()
list_grid_samples <- list()

for (i in seq_len(nrow(study_areas))) {
  study_area_idx <- i


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


  set.seed(124)

  list_clusters_10m[[study_area_idx]] <- sample_kmeans(
    input = sampling_input_pctile,
    clusters = round(study_areas[study_area_idx, ]$Shape_Area / 10000),
    use_xy = TRUE,
    sp_pts = TRUE,
    xy_weight = 2,
    candidates = candidates_field,
    weights = weights_field,
    min_cluster_size = 50
  )

  # Grid samples

  list_grid_samples[[study_area_idx]] <- sampling_input_pctile |>
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

  # Clusters based on sensors data

  field_dualem <- paste0(dir_sensors[study_area_idx], "/DUALEM/") |>
    list.files(
      pattern = "\\.tif$",
      full.names = TRUE
    ) |>
    lapply(
      function(x) {
        rast(x) |>
          trim()
      }
    )

  field_gamma <- paste0(dir_sensors[study_area_idx], "/GAMMA/") |>
    list.files(
      pattern = "\\.tif$",
      full.names = TRUE
    ) |>
    lapply(
      function(x) {
        rast(x) |>
          trim()
      }
    )

  input_sensors <- c(
    field_dualem[[4]],
    field_gamma[[1]]
  )

  candidates_sensor <- candidates_field |>
    as.polygons() |>
    rasterize(
      input_sensors[[1]]
    )

  weights_sensor <- weights_field |>
    resample(
      input_sensors[[1]]
    ) |>
    cover(
      y = input_sensors[[1]]*0 + 3
    )

  set.seed(124)

  list_clusters_sensor[[study_area_idx]] <- sample_kmeans(
    input = input_sensors,
    clusters = round(study_areas[study_area_idx, ]$Shape_Area / 10000),
    use_xy = TRUE,
    sp_pts = TRUE,
    xy_weight = 2,
    candidates = candidates_sensor,
    weights = weights_sensor,
    min_cluster_size = 1250
  )
}

# Plot 10 m clusters

dir_plots <- paste0(
  "C:/Users/anbm/OneDrive - SEGES Innovation PS/UDKIK/Figures/"
)

pdf(
  paste0(dir_plots, "/figure_clusters_10m.pdf")
)

lapply(
  seq_len(nrow(study_areas)),
  function(x) {
    study_area_idx <- x

    plot(
      as.factor(list_clusters_10m[[study_area_idx]]$clusters),
      main = paste0(
        "Study area ",
        study_area_idx,
        ": Clusters based on 10 m rasters"
      ),
      col = get_map_colors(
        L = 80,
        n = nrow(list_clusters_10m[[study_area_idx]]$points),
        minC = 25,
        maxC = 75
      ),
      ext = ext(study_areas[study_area_idx,]),
      buffer = TRUE
    )
    plot(
      as.polygons(
        list_clusters_10m[[study_area_idx]]$clusters
      ),
      1,
      add = TRUE,
      alpha = 0.25,
      col = NA,
      legend = FALSE,
      border = "black"
    )
    plot(
      list_clusters_10m[[study_area_idx]]$points,
      pch = 21,
      bg = "white",
      add = TRUE
    )
    plot(study_areas[study_area_idx,], add = TRUE)
    text(
      list_clusters_10m[[study_area_idx]]$points,
      list_clusters_10m[[study_area_idx]]$points$ID,
      cex = 0.7,
      col = "black",
      pos = 3,
      hc = "white",
      hw = 0.1,
      halo = TRUE
    )
  }
)

dev.off()

# Plot grid samples

pdf(
  paste0(dir_plots, "/figure_grid_samples.pdf")
)

lapply(
  seq_len(nrow(study_areas)),
  function(x) {
    study_area_idx <- x

    plot(
      as.factor(list_clusters_10m[[study_area_idx]]$clusters),
      main = paste0(
        "Study area ",
        study_area_idx,
        ": Grid samples"
      ),
      col = get_map_colors(
        L = 80,
        n = nrow(list_clusters_10m[[study_area_idx]]$points),
        minC = 25,
        maxC = 75
      ),
      ext = ext(study_areas[study_area_idx,]),
      alpha = 0.5,
      buffer = TRUE
    )
    plot(
      as.polygons(
        list_clusters_10m[[study_area_idx]]$clusters
      ),
      1,
      add = TRUE,
      alpha = 0.25,
      col = NA,
      legend = FALSE,
      border = "gray50"
    )
    plot(
      list_grid_samples[[study_area_idx]],
      pch = 24,
      bg = "yellow",
      add = TRUE,
      cex = 0.6
    )
    plot(
      list_clusters_10m[[study_area_idx]]$points,
      pch = 21,
      bg = "white",
      add = TRUE
    )
    plot(study_areas[study_area_idx,], add = TRUE)
    text(
      list_clusters_10m[[study_area_idx]]$points,
      list_clusters_10m[[study_area_idx]]$points$ID,
      cex = 0.7,
      col = "black",
      pos = 3,
      hc = "white",
      hw = 0.1,
      halo = TRUE
    )
  }
)

dev.off()


# Plot sensor clusters

pdf(
  paste0(dir_plots, "/figure_clusters_sensor.pdf")
)

lapply(
  seq_len(nrow(study_areas)),
  function(x) {
    study_area_idx <- x

    plot(
      trim(
        as.factor(list_clusters_sensor[[study_area_idx]]$clusters)
      ),
      main = paste0(
        "Study area ",
        study_area_idx,
        ": Clusters based on sensors"
      ),
      col = get_map_colors(
        L = 80,
        n = nrow(list_clusters_sensor[[study_area_idx]]$points),
        minC = 25,
        maxC = 75
      ),
      ext = ext(study_areas[study_area_idx,]),
      buffer = TRUE
    )
    plot(
      as.polygons(
        list_clusters_sensor[[study_area_idx]]$clusters
      ),
      1,
      add = TRUE,
      alpha = 0.25,
      col = NA,
      legend = FALSE,
      border = "black"
    )
    plot(
      list_clusters_sensor[[study_area_idx]]$points,
      pch = 21,
      bg = "white",
      add = TRUE
    )
    plot(study_areas[study_area_idx,], add = TRUE)
    text(
      list_clusters_sensor[[study_area_idx]]$points,
      list_clusters_sensor[[study_area_idx]]$points$ID,
      cex = 0.7,
      col = "black",
      pos = 3,
      hc = "white",
      hw = 0.1,
      halo = TRUE
    )
  }
)

dev.off()

# Summarise costs

n_samples_big <- lapply(
  list_clusters_10m,
  function(x) {
    nrow(x$points)
  }
) |>
  unlist()

n_samples_grid <- lapply(
  list_grid_samples,
  function(x) {
    nrow(x)
  }
) |>
  unlist()

sum(n_samples_big*2450) + sum(n_samples_grid*375)

# Eliminate duplicate grid points?

# END
