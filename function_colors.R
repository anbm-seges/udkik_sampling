# Function to create map colors

get_map_colors <- function(
    L = 50,
    n = NULL,
    minC = 50,
    maxC = 100
) {

  require(colorspace)

  lab_mat <- expand.grid(
    L = L,
    A = seq(-maxC, maxC, 2),
    B = seq(-maxC, maxC, 2)
  ) |>
    as.matrix() |>
    LAB()

  lab_mat_ok <- lab_mat |>
    as("RGB") |>
    coords() |>
    as.matrix() |>
    apply(1, function(x) {
      max(x) < 1 & min(x) > 0
    })

  lab_mat_ok2 <- lab_mat |>
    as("polarLAB") |>
    coords() |>
    as.matrix() |>
    apply(1, function(x) {
      x[2] < maxC & x[2] > minC
    }
    )

  lab_mat_ok <- as.logical(lab_mat_ok*lab_mat_ok2)

  lab_mat_selected <- lab_mat[lab_mat_ok, ]

  kmeans(
    coords(lab_mat_selected),
    n
  )$centers |>
    LAB() |>
    as("RGB") |>
    coords() |>
    rgb()
}


# END
