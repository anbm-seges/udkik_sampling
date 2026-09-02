# Sampling locations for study areas

library(terra)
library(samplekmeans)
library(tidyverse)
library(dplyr)
library(tidyterra)
library(vctrs)
library(rcartocolor)
library(openxlsx)

source("function_colors.R")

dir_data <- "C:/Users/anbm/OneDrive - SEGES Innovation PS/UDKIK/Data/"
dir_sensors <- paste0(dir_data, "/Sensors data/TIFF Files/") |>
  list.dirs(
    recursive = FALSE
  )

sampling_zones_ha <- 2
grid_spacing <- 28

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

pixels_per_ha_10m <- 10000 / (prod(res(sampling_input)))


# 2m DEM for study areas

dem_2m_list <- list.files(
  "C:/Users/anbm/OneDrive - SEGES Innovation PS/UDKIK/Data/Dem_2m/Dem_2m",
  pattern = "\\.tif$",
  full.names = TRUE
) |>
  lapply(FUN = rast)


# Prepare loop

list_clusters_10m <- list()
list_clusters_sensor <- list()
list_samples_10m <- list()
list_samples_sensor <- list()
list_samples_grid <- list()
list_sampling_areas_10m <- list()
list_sampling_areas_sensor <- list()

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

  candidates_field <- ifel(
    is.na(sum(sampling_input_field)),
    NA,
    1
  ) |>
    mask(
      study_areas[study_area_idx, ] |>
        buffer(width = -10),
      touches = FALSE
    )

  # Create clusters for 10 m

  n_clusters_10m_i <- round(
    study_areas[study_area_idx, ]$Shape_Area * sampling_zones_ha / 10000
  )

  run_again <- TRUE
  plus_clusters_10m <- 0
  try_seed <- 5082

  while (run_again) {
    print(
      paste0(
        format(Sys.time(), "%Y-%m-%d %X"),
        ": Study area ", i,
        " - 10 m input - seed ",
        try_seed,
        " - ",
        n_clusters_10m_i + plus_clusters_10m,
        " clusters requested - target = ",
        n_clusters_10m_i
      )
    )

    list_clusters_10m[[study_area_idx]] <- sample_kmeans(
      input = sampling_input_pctile,
      clusters = n_clusters_10m_i + plus_clusters_10m,
      use_xy = TRUE,
      sp_pts = TRUE,
      xy_weight = 2,
      candidates = candidates_field,
      min_cluster_size = pixels_per_ha_10m / (sampling_zones_ha * 2),
      seed = try_seed
    )

    if (
      nrow(list_clusters_10m[[study_area_idx]]$points) < n_clusters_10m_i
    ) {
      plus_clusters_10m <- plus_clusters_10m + 1
      if (plus_clusters_10m > (n_clusters_10m_i) * 0.5) {
        plus_clusters_10m <- plus_clusters_10m - 1
      }
      try_seed <- try_seed + 1
    } else {
      if (
        nrow(list_clusters_10m[[study_area_idx]]$points) > n_clusters_10m_i
      ) {
        plus_clusters_10m <- plus_clusters_10m - 1
        try_seed <- try_seed + 1
      } else {
        run_again <- FALSE
        print(
          paste0(
            format(Sys.time(), "%Y-%m-%d %X"),
            ": Study area ", i,
            " - 10 m input - seed ",
            try_seed,
            " - ",
            nrow(list_clusters_10m[[study_area_idx]]$points),
            " clusters created"
          )
        )
      }
    }
  }

  list_samples_10m[[study_area_idx]] <- list_clusters_10m[[
    study_area_idx
    ]]$points |>
    mutate(
      cluster = ID,
      ID = row_number(),
      type = "10m",
      Study_area = study_area_idx
    ) |>
    bind_spat_cols(
      geom(list_clusters_10m[[
        study_area_idx
      ]]$points) |>
        as.data.frame() |>
        select(x, y)
    )

  names(list_clusters_10m[[i]]$clusters) <- "cluster"
  names(list_clusters_10m[[i]]$distances) <- "distance"


  # Areas of similarity, 10 m input

  mean_dist <- zonal(
    list_clusters_10m[[i]]$distances,
    list_clusters_10m[[i]]$clusters,
    mean,
    as.raster = TRUE
  )

  list_sampling_areas_10m[[i]] <- terra::intersect(
    buffer(
      list_clusters_10m[[i]]$points,
      15
    ),
    as.polygons(
      list_clusters_10m[[i]]$clusters
    )
  ) |>
    filter(
      ID == cluster
    ) |>
    terra::intersect(
      as.polygons(
        ifel(
          list_clusters_10m[[i]]$distances < mean_dist,
          1,
          NA
        )
      )
    ) |>
    terra::intersect(
      buffer(
        study_areas[study_area_idx, ],
        -10
      )
    ) |>
    buffer(-2) |>
    buffer(2)


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

  field_dem2m <- dem_2m_list[[study_area_idx]] |>
    project(
      y = field_dualem[[4]],
      mask =  TRUE,
      method = "cubicspline"
    )

  input_sensors <- c(
    field_dualem[[4]],
    field_gamma[[1]],
    field_dem2m
  )

  # Transform sensor input to percentiles to create clusters of roughly the
  # same size.

  input_sensors_pctile <- sapp(
    input_sensors,
    function(x) {
      x_ecdf <- ecdf(values(x))
      app(x, x_ecdf)
    }
  )

  sensor_means <- global(
    input_sensors_pctile,
    "mean",
    na.rm = TRUE
  ) |>
    unlist()

  sensor_sds <- global(
    input_sensors_pctile,
    "sd",
    na.rm = TRUE
  ) |>
    unlist()

  input_sensors_pctile <- input_sensors_pctile |>
    clamp(
      lower = sensor_means - sensor_sds * 3,
      upper = sensor_means + sensor_sds * 3
    )

  candidates_sensor <- study_areas[study_area_idx, ] |>
      buffer(width = -10) |>
    rasterize(
      input_sensors[[1]]
    ) |>
    mask(
      buffer(list_clusters_10m[[study_area_idx]]$points, 4),
      inverse = TRUE,
      touches = FALSE
    )

  pixels_per_ha_sens <- 10000 / (prod(res(input_sensors)))

  n_clusters_sens_i <- round(
    study_areas[study_area_idx, ]$Shape_Area * sampling_zones_ha / 10000
  )

  run_again <- TRUE
  plus_clusters_sens <- 0
  try_seed <- 5082

  while (run_again) {
    print(
      paste0(
        format(Sys.time(), "%Y-%m-%d %X"),
        ": Study area ", i,
        " - sensor input - seed ",
        try_seed,
        " - ",
        n_clusters_sens_i + plus_clusters_sens,
        " clusters requested - target = ",
        n_clusters_sens_i
      )
    )

    list_clusters_sensor[[study_area_idx]] <- sample_kmeans(
      input = input_sensors_pctile,
      clusters = n_clusters_sens_i + plus_clusters_sens,
      use_xy = TRUE,
      sp_pts = TRUE,
      xy_weight = 2,
      candidates = candidates_sensor,
      min_cluster_size = pixels_per_ha_sens / (sampling_zones_ha * 2),
      seed = try_seed
    )

    if (
      nrow(list_clusters_sensor[[study_area_idx]]$points) < n_clusters_sens_i
      ) {
      plus_clusters_sens <- plus_clusters_sens + 1
      if (plus_clusters_sens > (n_clusters_sens_i) * 0.5) {
        plus_clusters_sens <- plus_clusters_sens - 1
      }
      try_seed <- try_seed + 1
    } else {
      if (
        nrow(list_clusters_sensor[[study_area_idx]]$points) > n_clusters_sens_i
      ) {
        plus_clusters_sens <- plus_clusters_sens - 1
        try_seed <- try_seed + 1
      } else {
        run_again <- FALSE
        print(
          paste0(
            format(Sys.time(), "%Y-%m-%d %X"),
            ": Study area ", i,
            " - sensor input - seed ",
            try_seed,
            " - ",
            nrow(list_clusters_sensor[[study_area_idx]]$points),
            " clusters created"
          )
        )
      }
    }
  }

  names(list_clusters_sensor[[i]]$clusters) <- "cluster"
  names(list_clusters_sensor[[i]]$distances) <- "distance"

  list_samples_sensor[[study_area_idx]] <- list_clusters_sensor[[
    study_area_idx
  ]]$points |>
    mask(
      list_samples_10m[[study_area_idx]],
      inverse = TRUE
    ) |>
    mutate(
      cluster = ID,
      ID = row_number() +
        nrow(list_samples_10m[[study_area_idx]]),
      type = "sensor",
      Study_area = study_area_idx
    ) |>
    bind_spat_cols(
      geom(list_clusters_sensor[[
        study_area_idx
      ]]$points) |>
        as.data.frame() |>
        select(x, y)
    )

  # Areas of similarity, sensors

  mean_dist_sensor <- zonal(
    list_clusters_sensor[[i]]$distances,
    list_clusters_sensor[[i]]$clusters,
    mean,
    as.raster = TRUE
  )

  list_sampling_areas_sensor[[i]] <- terra::intersect(
    buffer(
      list_clusters_sensor[[i]]$points,
      15
    ),
    as.polygons(
      list_clusters_sensor[[i]]$clusters
    )
  ) |>
    filter(
      ID == cluster
    ) |>
    terra::intersect(
      as.polygons(
        ifel(
          list_clusters_sensor[[i]]$distances < mean_dist_sensor,
          1,
          NA
        )
      )
    ) |>
    terra::intersect(
      buffer(
        study_areas[study_area_idx, ],
        -10
      )
    ) |>
    buffer(-2) |>
    buffer(2)
}


# New grid sampling locations

for (i in seq_len(nrow(study_areas))) {
  study_area_idx <- i

  ext_i <- study_areas[study_area_idx, ] |> ext()

  list_samples_grid[[study_area_idx]] <- expand.grid(
    x = seq(
      from = round(ext_i[1] / grid_spacing) * grid_spacing,
      to = ext_i[2],
      by = grid_spacing
    ),
    y = seq(
      from = round(ext_i[3] / grid_spacing) * grid_spacing,
      to = ext_i[4],
      by = grid_spacing
    )
  ) |>
    vect(
      crs = crs(study_areas),
      keepgeom = TRUE
    ) |>
    mask(
      buffer(study_areas[study_area_idx, ], -10)
    ) |>
    mask(
      buffer(list_clusters_10m[[study_area_idx]]$points, 4),
      inverse = TRUE
    ) |>
    mask(
      buffer(list_clusters_sensor[[study_area_idx]]$points, 4),
      inverse = TRUE
    ) |>
    mutate(
      ID = row_number() +
        nrow(list_samples_10m[[study_area_idx]]) +
        nrow(list_samples_sensor[[study_area_idx]]),
      type = "grid",
      study_area = study_area_idx,
      cluster = NA
    )
}



# Params for colors

my_L <- 80
my_maxC <- 30
my_minC <- 25

# Plot 10 m clusters

legend_placements <- c(
  "topright",
  "topleft",
  "bottomright",
  "topright",
  "bottomleft"
)

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
        n = nrow(list_clusters_10m[[study_area_idx]]$points),
        L = my_L,
        minC = my_minC,
        maxC = my_maxC
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
      list_sampling_areas_10m[[study_area_idx]],
      1,
      add = TRUE,
      alpha = 0.25,
      col = "gray70",
      legend = FALSE,
      border = "gray50"
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
        n = nrow(list_clusters_sensor[[study_area_idx]]$points),
        L = my_L,
        minC = my_minC,
        maxC = my_maxC
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
      list_sampling_areas_sensor[[study_area_idx]],
      1,
      add = TRUE,
      alpha = 0.25,
      col = "gray70",
      legend = FALSE,
      border = "gray50"
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


# Plot grid samples

pdf(
  paste0(dir_plots, "/figure_grid_samples.pdf")
)

lapply(
  seq_len(nrow(study_areas)),
  function(x) {
    study_area_idx <- x

    name_grid_samples <- paste0(
      "Grid samples (spacing = ",
      grid_spacing,
      " m)"
    )

    shapes_pts <- c(
      "k-means samples (10 m input)" = 21,
      "k-means samples (sensor input)" = 23,
      "Grid samples" = 24
    )
    names(shapes_pts) <- c(
      "k-means samples (10 m input)",
      "k-means samples (sensor input)",
      name_grid_samples
    )
    colors_pts <- c(
      "k-means samples (10 m input)" = "white",
      "k-means samples (sensor input)" = "gray",
      "Grid samples" = "yellow"
    )
    names(colors_pts) <- c(
      "k-means samples (10 m input)",
      "k-means samples (sensor input)",
      name_grid_samples
    )

    pts_all_i <- rbind(
      list_samples_10m[[study_area_idx]] |>
        mutate(type = "k-means samples (10 m input)"),
      list_samples_sensor[[study_area_idx]] |>
        mutate(type = "k-means samples (sensor input)"),
      list_samples_grid[[study_area_idx]] |>
        mutate(type = name_grid_samples)
    )

    tidyterra::autoplot(
      dem_2m_list[[study_area_idx]]
    ) +
      geom_spatvector(
        data = study_areas[study_area_idx,],
        fill = NA,
        color = "black"
      ) +
      labs(
        fill = "Elevation (m)"
      ) +
      ggnewscale::new_scale_fill() +
      geom_spatvector(
        data = pts_all_i,
        aes(
          shape = type,
          fill = type
        ),
        color = "black",
        size = 1.5,
        show.legend = TRUE
      ) +
      ggtitle(
        paste0(
          "Study area ",
          study_area_idx,
          ": Grid samples"
        )
      ) +
      scale_shape_manual(
        name = "Sample type",
        values = shapes_pts,
        guide = guide_legend(order = 1)
      ) +
      scale_fill_manual(
        name = "Sample type",
        values = colors_pts,
        guide = guide_legend(order = 1)
      ) +
      theme_bw()
  }
)

dev.off()

# Plot 10m and sensor points together

pdf(
  paste0(dir_plots, "/figure_clusters_comparison.pdf")
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
        ": Clusters based on 10 m rasters and sensors"
      ),
      col = get_map_colors(
        n = nrow(list_clusters_10m[[study_area_idx]]$points),
        L = my_L,
        minC = my_minC,
        maxC = my_maxC
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
    plot(study_areas[study_area_idx,], add = TRUE)
    plot(
      list_clusters_10m[[study_area_idx]]$points,
      pch = 21,
      bg = "white",
      add = TRUE
    )
    plot(
      list_clusters_sensor[[study_area_idx]]$points,
      pch = 23,
      bg = "gray",
      add = TRUE
    )
    add_legend(
      legend_placements[study_area_idx],
      legend = c(
        "k-means samples (10 m input)",
        "k-means samples (sensor input)"
      ),
      pch = c(21, 23),
      cex = 0.6,
      pt.cex = 1,
      pt.bg = c("white", "gray"),
      bty = "n"
    )
  }
)

dev.off()

# Summarise costs

cost_grid <- 375
# cost_grid <- 675
cost_10m <- 675
cost_sensor <- 675

n_samples_10m <- lapply(
  list_clusters_10m,
  function(x) {
    nrow(x$points)
  }
) |>
  unlist()

n_samples_sensor <- lapply(
  list_clusters_sensor,
  function(x) {
    nrow(x$points)
  }
) |>
  unlist()

n_samples_grid <- lapply(
  list_samples_grid,
  function(x) {
    nrow(x)
  }
) |>
  unlist()

sampling_costs <- data.frame(
  Study_area = as.character(seq_len(nrow(study_areas))),
  n_samples_10m = n_samples_10m,
  n_samples_sensor = n_samples_sensor,
  n_samples_grid = n_samples_grid
) |>
  mutate(
    cost_10m = n_samples_10m * cost_10m,
    cost_sensor = n_samples_sensor * cost_sensor,
    cost_grid = n_samples_grid * cost_grid,
    total_cost = cost_10m + cost_sensor + cost_grid
  )

sampling_costs

sum_sampling_costs <- sampling_costs |>
  summarise(
    n_samples_10m = sum(n_samples_10m),
    n_samples_sensor = sum(n_samples_sensor),
    n_samples_grid = sum(n_samples_grid),
    cost_10m = sum(cost_10m),
    cost_sensor = sum(cost_sensor),
    cost_grid = sum(cost_grid),
    total_cost = sum(total_cost)
  ) |>
  mutate(
    Study_area = "Sum"
  )

sum_sampling_costs

bind_rows(
  sampling_costs,
  sum_sampling_costs
) |>
  write.xlsx(
    paste0(dir_data, "Sampling_costs_", Sys.Date(), ".xlsx"),
    overwrite = TRUE
  )



# END
