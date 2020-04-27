#' Gaussian process distribution families for bamlss
#'
#' These functions set up Gaussian process (GP) distribution families for
#' [bamlss][bamlss::bamlss-package].
#'
#' `gp_constant_bamlss` uses a constant mean function and a Matérn covariance
#' function with the smoothness parameter \eqn{\nu} = 1.5. The distributional
#' parameters of `gp_constant_bamlss` are
#' - the **mean**,
#' - the **stddev** (standard deviation),
#' - and the **range**.
#'
#' `gp_weibull_bamlss` uses a Weibull growth curve as a mean function and a
#' Matérn covariance function with the smoothness parameter \eqn{\nu} = 1.5.
#' The distributional parameters of `gp_weibull_bamlss` are
#' - the **limit** (of the Weibull growth curve),
#' - the **shape** (of the Weibull growth curve),
#' - the **scale** (of the Weibull growth curve),
#' - the **stddev** (standard deviation),
#' - and the **range**.
#'
#' Optionally, a linear mean function can be added to the Weibull growth curve
#' by passing a list of design matrices `Z` to `gp_weibull_bamlss`. Each column
#' of the design matrices adds one distributional parameter with the
#' corresponding column name to the family.
#'
#' `gp_sphere_bamlss` uses a mean function on a sphere and an exponential
#' covariance function. It can be used to model tree crowns. The distributional
#' parameters of `gp_sphere_bamlss` are
#' - the **radius** (of the tree crown),
#' - the **south**-bound expansion (of the tree crown),
#' - the **height** (of the tree crown),
#' - the **stddev** (standard deviation),
#' - and the **range**.
#'
#' @usage
#' gp_constant_bamlss(y, D, initial_par, ...)
#' gp_weibull_bamlss(y, s, Z = NULL, D, q = NULL, initial_par, ...)
#' gp_sphere_bamlss(y, s, D, initial_par, ...)
#'
#' @param y The GPs to be used as response observations. A list of numeric
#'          vectors.
#' @param s The coordinates at which the GPs were observed. A list of numeric
#'          vectors or matrices.
#' @param Z Optionally, the design matrices for an additive linear mean
#'          function. A list of numeric matrices.
#' @param D The distance matrices of the coordinates. A list of numeric
#'          matrices.
#' @param q Optionally, the scaling factors of the standard deviation of the
#'          GPs at the coordinates. A list of numeric vectors.
#' @param initial_par The initial distributional parameters. A data frame with
#'                    one column per distributional parameter.
#' @param ... Not used.
#'
#' @return
#' A [`family.bamlss`][bamlss::family.bamlss] object.
#'
#' @examples
#' \dontrun{
#' N <- 30
#' n <- 30
#'
#' mean <- 0
#' stddev <- exp(0)
#' range <- exp(0)
#'
#' s <- seq.int(0, 1, length.out = n)
#' D <- as.matrix(dist(s))
#'
#' mu <- rep.int(mean, n)
#' Sigma <- stddev^2 * (1 + D / range) * exp(-D / range)
#' y <- MASS::mvrnorm(N, mu, Sigma)
#'
#' y <- unlist(apply(y, 1, list), recursive = FALSE)
#' D <- rep.int(list(D), N)
#'
#' initial_par <- data.frame(
#'   mean = rep.int(mean, N),
#'   stddev = rep.int(stddev, N),
#'   range = rep.int(range, N)
#' )
#'
#' family <- gp_constant_bamlss(y, D, initial_par)
#' model <- bamlss(rep.int(0, N) ~ 1, family)
#' joint <- sample_bamlss_gp(model)
#' }
#'
#' @aliases gp_weibull_bamlss gp_sphere_bamlss
#' @importFrom purrr map map_dbl map_int map2 map2_dbl
#' @importFrom rlang rep_along
#' @export

gp_constant_bamlss <- function(y, D, initial_par, ...) {
  # initialize cache

  n <- map_int(y, length)
  two_n <- 2 * n

  y_mu <- NULL
  R_det <- R_inv <- dR_dphi <- R_inv_dR_dphi <- NULL
  R_inv_y_mu <- y_mu_R_inv_y_mu <- NULL

  previous_par <- rep_along(initial_par, list(NULL))

  # update cache

  update_cache <- function(par) {
    updated <- check_if_updated(par, previous_par)

    if (updated["mean"]) {
      y_mu <<- map2(y, par$mean, `-`)
    }

    if (updated["range"]) {
      R_chol <- map2(D, par$range, ~ chol(rho_fun_matern(.x, .y)))
      R_det <<- map_dbl(R_chol, ~ 2 * sum(log(diag(.))))
      R_inv <<- map(R_chol, chol2inv)

      dR_dphi <<- map2(D, par$range, drho_dphi_matern)
      R_inv_dR_dphi <<- map2(R_inv, dR_dphi, `%*%`)
    }

    if (sum(updated[c("mean", "range")])) {
      R_inv_y_mu <<- map2(R_inv, y_mu, `%*%`)
      y_mu_R_inv_y_mu <<- map2_dbl(y_mu, R_inv_y_mu, `%*%`)
    }

    previous_par <<- par
  }

  update_cache(initial_par)

  # assemble family object

  family <- list(
    family = "gp_constant",
    names = c("mean", "stddev", "range"),
    links = c(mean = "identity", stddev = "log", range = "log"),
    d = function(y, par, log = FALSE) {
      update_cache(par)
      -(two_n * log(par$stddev) + R_det + y_mu_R_inv_y_mu / par$stddev^2) / 2
    },
    score = list(
      mean = function(y, par, ...) {
        update_cache(par)
        out <- map_dbl(R_inv_y_mu, sum)
        out / par$stddev^2
      },
      stddev = function(y, par, ...) {
        update_cache(par)
        # derivative w.r.t. predictor, not distributional parameter
        -n + y_mu_R_inv_y_mu / par$stddev^2
      },
      range = function(y, par, ...) {
        update_cache(par)
        trace <- map_dbl(R_inv_dR_dphi, ~ sum(diag(.)))
        qform <- map2_dbl(R_inv_y_mu, dR_dphi, ~ t(.x) %*% .y %*% .x)
        -(trace - qform / par$stddev^2) / 2
      }
    ),
    hess = list(
      mean = function(y, par, ...) {
        update_cache(par)
        out <- map_dbl(R_inv, sum)
        out / par$stddev^2
      },
      stddev = function(y, par, ...) {
        # derivative w.r.t. preditor, not distributional parameter
        two_n
      },
      range = function(y, par, ...) {
        update_cache(par)
        out <- map_dbl(R_inv_dR_dphi, ~ sum(t(.) * .))
        out / 2
      }
    ),
    initialize = list(
      mean = function(y, ...) initial_par$mean,
      stddev = function(y, ...) initial_par$stddev,
      range = function(y, ...) initial_par$range
    )
  )

  class(family) <- "family.bamlss"
  family
}


#' @importFrom purrr array_branch map map_dbl map_int map2 map2_dbl pmap
#'                   transpose
#' @importFrom rlang rep_along set_names
#' @export

gp_weibull_bamlss <- function(y, s, Z = NULL, D, q = NULL, initial_par, ...) {
  # process arguments

  if (!is.null(q)) {
    qq <- map(q, tcrossprod)
  } else {
    qq <- rep_along(y, list(NULL))
  }

  if (!is.null(Z)) {
    linear_names <- colnames(Z[[1]])

    if (is.null(linear_names)) {
      linear_names <- paste0("linear", 1:ncol(Z[[1]]))
    }

    for (i in 1:length(Z)) {
      colnames(Z[[i]]) <- linear_names
    }
  } else {
    linear_names <- character()
  }

  # initialize cache

  N <- length(y)
  one_to_N <- 1:N
  n <- map_int(y, length)
  two_n <- 2 * n

  weibull_mu <- 0
  dmu_dl <- dmu_da <- dmu_db <- NULL

  linear_mu <- 0
  dmu_dlinear <- transpose(map(Z, array_branch, margin = 2))
  linear_par <- matrix(nrow = N, ncol = length(linear_names))
  colnames(linear_par) <- linear_names

  y_mu <- NULL
  R_det <- R_inv <- dR_dphi <- R_inv_dR_dphi <- NULL
  R_inv_y_mu <- y_mu_R_inv_y_mu <- NULL

  previous_par <- rep_along(initial_par, list(NULL))

  # update cache

  update_cache <- function(par) {
    updated <- check_if_updated(par, previous_par)

    if (sum(updated[c("limit", "shape", "scale")])) {
      weibull_mu <<- pmap(list(s, par$limit, par$shape, par$scale), m_fun)

      dmu_dl <<- weibull_mu
      dmu_da <<- pmap(list(s, par$limit, par$shape, par$scale), dm_da)
      dmu_db <<- pmap(list(s, par$limit, par$shape, par$scale), dm_db)
    }

    for (nm in linear_names) {
      if (updated[nm]) linear_par[, nm] <<- par[[nm]]
    }

    if (sum(updated[linear_names])) {
      linear_mu <<- map(one_to_N, ~ drop(Z[[.]] %*% linear_par[.,]))
    }

    if (sum(updated[c("limit", "shape", "scale", linear_names)])) {
      y_mu <<- pmap(list(y, weibull_mu, linear_mu), ~ ..1 - ..2 - ..3)
    }

    if (updated["range"]) {
      R_chol <- pmap(
        .l = list(D, par$range, qq),
        .f = ~ chol(rho_fun_matern(..1, ..2, ..3))
      )

      R_det <<- map_dbl(R_chol, ~ 2 * sum(log(diag(.))))
      R_inv <<- map(R_chol, chol2inv)

      dR_dphi <<- pmap(list(D, par$range, qq), drho_dphi_matern)
      R_inv_dR_dphi <<- map2(R_inv, dR_dphi, `%*%`)
    }

    if (sum(updated[c("limit", "shape", "scale", linear_names, "range")])) {
      R_inv_y_mu <<- map2(R_inv, y_mu, `%*%`)
      y_mu_R_inv_y_mu <<- map2_dbl(y_mu, R_inv_y_mu, `%*%`)
    }

    previous_par <<- par
  }

  update_cache(initial_par)

  # assemble family object

  family <- list(
    family = "gp_weibull",
    names = c("limit", "shape", "scale", "stddev", "range"),
    links = c(limit = "log", shape = "log", scale = "log", stddev = "log",
              range = "log"),
    d = function(y, par, log = FALSE) {
      update_cache(par)
      -(two_n * log(par$stddev) + R_det + y_mu_R_inv_y_mu / par$stddev^2) / 2
    },
    score = list(
      limit = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(dmu_dl, R_inv_y_mu, `%*%`)
        out / par$stddev^2
      },
      shape = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(dmu_da, R_inv_y_mu, `%*%`)
        out / par$stddev^2
      },
      scale = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(dmu_db, R_inv_y_mu, `%*%`)
        out / par$stddev^2
      },
      stddev = function(y, par, ...) {
        update_cache(par)
        # derivative w.r.t. predictor, not distributional parameter
        -n + y_mu_R_inv_y_mu / par$stddev^2
      },
      range = function(y, par, ...) {
        update_cache(par)
        trace <- map_dbl(R_inv_dR_dphi, ~ sum(diag(.)))
        qform <- map2_dbl(R_inv_y_mu, dR_dphi, ~ t(.x) %*% .y %*% .x)
        -(trace - qform / par$stddev^2) / 2
      }
    ),
    hess = list(
      limit = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(dmu_dl, R_inv, ~ .x %*% .y %*% .x)
        out / par$stddev^2
      },
      shape = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(dmu_da, R_inv, ~ .x %*% .y %*% .x)
        out / par$stddev^2
      },
      scale = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(dmu_db, R_inv, ~ .x %*% .y %*% .x)
        out / par$stddev^2
      },
      stddev = function(y, par, ...) {
        # derivative w.r.t. preditor, not distributional parameter
        two_n
      },
      range = function(y, par, ...) {
        update_cache(par)
        out <- map_dbl(R_inv_dR_dphi, ~ sum(t(.) * .))
        out / 2
      }
    ),
    initialize = list(
      limit = function(y, ...) initial_par$limit,
      shape = function(y, ...) initial_par$shape,
      scale = function(y, ...) initial_par$scale,
      stddev = function(y, ...) initial_par$stddev,
      range = function(y, ...) initial_par$range
    )
  )

  # add linear mean function to family object

  linear_links <- rep_along(linear_names, "identity")

  linear_score <- map(linear_names, ~ function(y, par, ...) {
    update_cache(par)
    out <- map2_dbl(dmu_dlinear[[.]], R_inv_y_mu, `%*%`)
    out / par$stddev^2
  })

  linear_hess <- map(linear_names, ~ function(y, par, ...) {
    update_cache(par)
    out <- map2_dbl(dmu_dlinear[[.]], R_inv, ~ .x %*% .y %*% .x)
    out / par$stddev^2
  })

  linear_initialize <- map(linear_names, ~ function(y, ...) {
    initial_par[[.]]
  })

  linear_links <- set_names(linear_links, linear_names)
  linear_score <- set_names(linear_score, linear_names)
  linear_hess <- set_names(linear_hess, linear_names)
  linear_initialize <- set_names(linear_initialize, linear_names)

  family$names <- append(family$names, linear_names, after = 3)
  family$links <- append(family$links, linear_links, after = 3)
  family$score <- append(family$score, linear_score, after = 3)
  family$hess <- append(family$hess, linear_hess, after = 3)
  family$initialize <- append(family$initialize, linear_initialize, after = 3)

  class(family) <- "family.bamlss"
  family
}


#' @importFrom purrr map map_dbl map_int map2 map2_dbl
#' @importFrom rlang rep_along
#' @export

gp_sphere_bamlss <- function(y, s, D, initial_par, ...) {
  # process arguments

  Z <- map(
    .x = s,
    .f = ~ cbind(
      radius = rep.int(1, nrow(.)),
      south = cos(.[, 2]) * (cos(.[, 1]) + 1),
      height = .[, 2] + pi / 2
    )
  )

  Z_par <- Z

  # initialize cache

  one_to_N <- 1:length(y)
  n <- map_int(y, length)
  two_n <- 2 * n

  y_mu <- NULL
  R_det <- R_inv <- dR_dphi <- R_inv_dR_dphi <- NULL
  R_inv_y_mu <- y_mu_R_inv_y_mu <- NULL

  previous_par <- rep_along(initial_par, list(NULL))

  # update cache

  update_cache <- function(par) {
    updated <- check_if_updated(par, previous_par)

    if (updated["radius"]) {
      for (i in one_to_N) {
        Z_par[[i]][, "radius"] <<- Z[[i]][, "radius"] * par$radius[i]
      }
    }

    if (updated["south"]) {
      for (i in one_to_N) {
        Z_par[[i]][, "south"] <<- Z[[i]][, "south"] * par$south[i]
      }
    }

    if (updated["height"]) {
      for (i in one_to_N) {
        Z_par[[i]][, "height"] <<- Z[[i]][, "height"] * par$height[i]
      }
    }

    if (sum(updated[c("radius", "south", "height")])) {
      y_mu <<- map2(y, Z_par, ~ .x - rowSums(.y))
    }

    if (updated["range"]) {
      R_chol <- map2(D, par$range, ~ chol(rho_fun_exp(.x, .y)))
      R_det <<- map_dbl(R_chol, ~ 2 * sum(log(diag(.))))
      R_inv <<- map(R_chol, chol2inv)

      dR_dphi <<- map2(D, par$range, drho_dphi_exp)
      R_inv_dR_dphi <<- map2(R_inv, dR_dphi, `%*%`)
    }

    if (sum(updated[c("radius", "south", "height", "range")])) {
      R_inv_y_mu <<- map2(R_inv, y_mu, `%*%`)
      y_mu_R_inv_y_mu <<- map2_dbl(y_mu, R_inv_y_mu, `%*%`)
    }

    previous_par <<- par
  }

  update_cache(initial_par)

  # assemble family object

  family <- list(
    family = "gp_sphere",
    names = c("radius", "south", "height", "stddev", "range"),
    links = c(radius = "log", south = "log", height = "log", stddev = "log",
              range = "log"),
    d = function(y, par, log = FALSE) {
      update_cache(par)
      -(two_n * log(par$stddev) + R_det + y_mu_R_inv_y_mu / par$stddev^2) / 2
    },
    score = list(
      radius = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(Z_par, R_inv_y_mu, ~ .x[, "radius"] %*% .y)
        out / par$stddev^2
      },
      south = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(Z_par, R_inv_y_mu, ~ .x[, "south"] %*% .y)
        out / par$stddev^2
      },
      height = function(y, par, ...) {
        update_cache(par)
        out <- map2_dbl(Z_par, R_inv_y_mu, ~ .x[, "height"] %*% .y)
        out / par$stddev^2
      },
      stddev = function(y, par, ...) {
        update_cache(par)
        # derivative w.r.t. predictor, not distributional parameter
        -n + y_mu_R_inv_y_mu / par$stddev^2
      },
      range = function(y, par, ...) {
        update_cache(par)
        trace <- map_dbl(R_inv_dR_dphi, ~ sum(diag(.)))
        qform <- map2_dbl(R_inv_y_mu, dR_dphi, ~ t(.x) %*% .y %*% .x)
        -(trace - qform / par$stddev^2) / 2
      }
    ),
    hess = list(
      radius = function(y, par, ...) {
        update_cache(par)

        out <- map2_dbl(
          .x = Z_par,
          .y = R_inv,
          .f = ~ .x[, "radius"] %*% .y %*% .x[, "radius"]
        )

        out / par$stddev^2
      },
      south = function(y, par, ...) {
        update_cache(par)

        out <- map2_dbl(
          .x = Z_par,
          .y = R_inv,
          .f = ~ .x[, "south"] %*% .y %*% .x[, "south"]
        )

        out / par$stddev^2
      },
      height = function(y, par, ...) {
        update_cache(par)

        out <- map2_dbl(
          .x = Z_par,
          .y = R_inv,
          .f = ~ .x[, "height"] %*% .y %*% .x[, "height"]
        )

        out / par$stddev^2
      },
      stddev = function(y, par, ...) {
        # derivative w.r.t. preditor, not distributional parameter
        two_n
      },
      range = function(y, par, ...) {
        update_cache(par)
        out <- map_dbl(R_inv_dR_dphi, ~ sum(t(.) * .))
        out / 2
      }
    ),
    initialize = list(
      radius = function(y, ...) initial_par$radius,
      south = function(y, ...) initial_par$south,
      height = function(y, ...) initial_par$height,
      stddev = function(y, ...) initial_par$stddev,
      range = function(y, ...) initial_par$range
    )
  )

  class(family) <- "family.bamlss"
  family
}
