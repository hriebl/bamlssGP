#' @importFrom bamlssAPI grad_logpost

joint_grad <- function(model) {
  grad_stddev <- grad_logpost(model, "stddev", "p")
  grad_range <- grad_logpost(model, "range", "p")
  c(grad_stddev, grad_range)
}

#' @importFrom bamlssAPI hess_logpost
#' @importFrom purrr map_dbl

joint_hess <- function(model) {
  hess_stddev <- hess_logpost(model, "stddev", "p")
  hess_range <- hess_logpost(model, "range", "p")

  env <- environment(model$family$d)
  hess_cross <- map_dbl(env$R_inv_dR_dphi, ~ sum(diag(.)))
  X_stddev <- model$x$stddev$smooth.construct$model.matrix$X
  X_range <- model$x$range$smooth.construct$model.matrix$X
  hess_cross <- crossprod(X_stddev * hess_cross, X_range)

  rbind(
    cbind(hess_stddev, hess_cross),
    cbind(t(hess_cross), hess_range)
  )
}

propose_fisher <- function(state, grad, hess) {
  chol <- chol(hess)

  # faster than: tmp <- solve(hess) %*% grad

  tmp <- backsolve(
    r = chol,
    x = forwardsolve(
      l = chol,
      x = grad,
      upper.tri = TRUE,
      transpose = TRUE
    )
  )

  mean <- state + tmp # / 2
  prop <- drop(rmvnorm_p(1, mean, chol = chol))
  forward <- dmvnorm_p(prop, mean, hess, chol, log = TRUE)
  list(prop = prop, forward = forward)
}

backward_fisher <- function(state, prop, grad, hess) {
  chol <- chol(hess)

  tmp <- backsolve(
    r = chol,
    x = forwardsolve(
      l = chol,
      x = grad,
      upper.tri = TRUE,
      transpose = TRUE
    )
  )

  mean <- prop + tmp # / 2
  dmvnorm_p(state, mean, hess, chol, log = TRUE)
}


#' Sample from the posterior of a bamlss model with a Gaussian process response
#'
#' This function provides an alternative MCMC sampler for
#' [bamlss][bamlss::bamlss-package] models with
#' [Gaussian process (GP) responses][gp_constant_bamlss].
#' It is typically more efficient than the default bamlss sampler by sampling
#' the parametric predictors for the standard deviation and the range jointly,
#' thus reducing the autocorrelation of the chains.
#'
#' @param model A bamlss model with a GP response.
#' @param n The sample size/length of the chain to be drawn.
#' @param only_shared If `TRUE`, only the parameters for those covariates which
#'                    appear in the parametric predictors for both the standard
#'                    deviation and the range are sampled jointly. Not helpful
#'                    in my experience.
#'
#' @return
#' A list with the elements `model` and `samples`.
#'
#' @inherit gp_constant_bamlss examples
#'
#' @importFrom bamlssAPI accept apify logpost parameters predictors propose
#'                       set_parameters smooths update_logpost
#' @importFrom purrr map
#' @importFrom stats runif
#' @export

sample_bamlss_gp <- function(model, n = 1000, only_shared = FALSE) {
  model <- apify(model)

  s_smt <- map(predictors(model), function(predictor) {
    out <- smooths(model, predictor)
    if (any(predictor == c("stddev", "range"))) out <- out[out != "p"]
    out
  })

  names(s_smt) <- predictors(model)

  names_sd <- names(parameters(model, "stddev", "p"))
  names_rng <- names(parameters(model, "range", "p"))

  if (only_shared) {
    j_sd <- which(names_sd %in% names_rng)
    j_rng <- which(names_rng %in% names_sd)
  } else {
    j_sd <- seq_along(names_sd)
    j_rng <- seq_along(names_rng)
  }

  j_all <- c(j_sd, length(names_sd) + j_rng)

  s_sd <- setdiff(seq_along(names_sd), j_sd)
  s_rng <- setdiff(seq_along(names_rng), j_rng)

  sd_j <- seq_along(j_sd)
  rng_j <- length(j_sd) + seq_along(j_rng)

  samples <- map(predictors(model), function(predictor) {
    out <- map(smooths(model, predictor), function(smooth) {
      par <- parameters(model, predictor, smooth)
      out <- matrix(nrow = n, ncol = length(par))
      colnames(out) <- names(par)
      out
    })

    names(out) <- smooths(model, predictor)
    out
  })

  names(samples) <- predictors(model)

  for (i in 1:n) {
    if (i %% 100 == 0) message("Iteration ", i, "...")

    for (predictor in predictors(model)) {
      for (smooth in s_smt[[predictor]]) {
        prop <- propose(model, predictor, smooth)

        if (log(runif(1)) <= prop$alpha) {
          model <- accept(model, predictor, smooth, prop)
          model <- update_logpost(model)
        }

        par <- parameters(model, predictor, smooth)
        samples[[predictor]][[smooth]][i,] <- par
      }
    }

    if (length(s_sd)) {
      state <- prop_sd <- parameters(model, "stddev", "p")
      grad <- joint_grad(model)
      hess <- joint_hess(model)

      prop <- propose_fisher(state[s_sd], grad[s_sd], hess[s_sd, s_sd])

      prop_sd[s_sd] <- prop$prop
      prop_model <- set_parameters(model, "stddev", "p", prop_sd)
      prop_model <- update_logpost(prop_model)

      grad <- joint_grad(prop_model)
      hess <- joint_hess(prop_model)

      backward <- backward_fisher(
        state = state[s_sd],
        prop = prop$prop,
        grad = grad[s_sd],
        hess = hess[s_sd, s_sd]
      )

      alpha <- logpost(prop_model) - logpost(model) + backward - prop$forward
      if (log(runif(1)) <= alpha) model <- prop_model
    }

    if (length(s_rng)) {
      state <- prop_rng <- parameters(model, "range", "p")
      grad <- joint_grad(model)
      hess <- joint_hess(model)

      prop <- propose_fisher(state[s_rng], grad[s_rng], hess[s_rng, s_rng])

      prop_rng[s_rng] <- prop$prop
      prop_model <- set_parameters(model, "range", "p", prop_rng)
      prop_model <- update_logpost(prop_model)

      grad <- joint_grad(prop_model)
      hess <- joint_hess(prop_model)

      backward <- backward_fisher(
        state = state[s_rng],
        prop = prop$prop,
        grad = grad[s_rng],
        hess = hess[s_rng, s_rng]
      )

      alpha <- logpost(prop_model) - logpost(model) + backward - prop$forward
      if (log(runif(1)) <= alpha) model <- prop_model
    }

    if (length(j_all)) {
      prop_sd <- parameters(model, "stddev", "p")
      prop_rng <- parameters(model, "range", "p")
      state <- c(prop_sd, prop_rng)
      grad <- joint_grad(model)
      hess <- joint_hess(model)

      prop <- propose_fisher(state[j_all], grad[j_all], hess[j_all, j_all])

      prop_sd[j_sd] <- prop$prop[sd_j]
      prop_rng[j_rng] <- prop$prop[rng_j]
      prop_model <- set_parameters(model, "stddev", "p", prop_sd)
      prop_model <- set_parameters(prop_model, "range", "p", prop_rng)
      prop_model <- update_logpost(prop_model)

      grad <- joint_grad(prop_model)
      hess <- joint_hess(prop_model)

      backward <- backward_fisher(
        state = state[j_all],
        prop = prop$prop,
        grad = grad[j_all],
        hess = hess[j_all, j_all]
      )

      alpha <- logpost(prop_model) - logpost(model) + backward - prop$forward
      if (log(runif(1)) <= alpha) model <- prop_model
    }

    samples$stddev$p[i,] <- parameters(model, "stddev", "p")
    samples$range$p[i,] <- parameters(model, "range", "p")
  }

  list(model = model, samples = samples)
}
