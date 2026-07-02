# Ranked samples for UDKIK

library(terra)
library(samplekmeans)
library(tidyverse)
library(dplyr)
library(tidyterra)

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

set.seed(123)

myclusters <- sample_kmeans(
  input = r,
  clusters = 3,
  use_xy = TRUE,
  sp_pts = TRUE
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
