# multivariate normal distribution --------------------------------------------

LOG2PI <- log(2 * pi)

dmvnorm_p <- function(x, mean, precision, chol = NULL, log = FALSE) {
  if (is.null(chol)) chol <- base::chol(precision)

  x <- x - mean
  det <- 2 * sum(log(diag(chol)))
  out <- (det - drop(x %*% precision %*% x) - length(x) * LOG2PI) / 2
  if (!log) out <- exp(out)

  out
}

#' @importFrom stats rnorm

rmvnorm_p <- function(n, mean, precision, chol = NULL) {
  if (is.null(chol)) chol <- base::chol(precision)

  d <- nrow(chol)
  tmp <- matrix(rnorm(d * n), d, n)
  out <- backsolve(chol, tmp)
  out <- out + mean

  out
}

#' @importFrom stats rnorm

rmvnorm_c <- function(n, mean, covariance, chol = NULL) {
  if (is.null(chol)) chol <- base::chol(covariance)

  d <- nrow(chol)
  tmp <- matrix(rnorm(d * n), d, n)
  out <- crossprod(chol, tmp)
  out <- out + mean

  out
}


# weibull mean function -------------------------------------------------------

#' @importFrom stats pweibull

m_fun <- function(s, l, a, b) {
  l * pweibull(s, a, b)
}

# derivative of the weibull mean function w.r.t. the predictors
# (rather than the distributional parameters!), assuming log links

#' @importFrom stats dweibull

dm_da <- function(s, l, a, b) {
  tmp <- (s / b)^a
  l * a * (log(s) - log(b)) * tmp * dweibull(tmp, 1, 1)
}

#' @importFrom stats dweibull

dm_db <- function(s, l, a, b) {
  tmp <- s / b
  -l * tmp * dweibull(tmp, a, 1)
}


# exponential correlation function --------------------------------------------

rho_fun_exp <- function(D, phi, qq = NULL) {
  out <- exp(-(D / phi))
  if (!is.null(qq)) out <- qq * out
  out
}

c_fun_exp <- function(D, sigma, phi, qq = NULL) {
  sigma^2 * rho_fun_exp(D, phi, qq)
}

# derivative of the exponential correlation function w.r.t. the predictor
# (rather than the distributional parameter!), assuming a log link

drho_dphi_exp <- function(D, phi, qq = NULL) {
  tmp <- D / phi
  out <- tmp * exp(-tmp)
  if (!is.null(qq)) out <- qq * out
  out
}


# matern correlation function -------------------------------------------------

rho_fun_matern <- function(D, phi, qq = NULL) {
  tmp <- D / phi
  out <- (1 + tmp) * exp(-tmp)
  if (!is.null(qq)) out <- qq * out
  out
}

c_fun_matern <- function(D, sigma, phi, qq = NULL) {
  sigma^2 * rho_fun_matern(D, phi, qq)
}

# derivative of the matern correlation function w.r.t. the predictor
# (rather than the distributional parameter!), assuming a log link

drho_dphi_matern <- function(D, phi, qq = NULL) {
  tmp <- D / phi
  out <- tmp^2 * exp(-tmp)
  if (!is.null(qq)) out <- qq * out
  out
}


# gauss correlation function --------------------------------------------------

rho_fun_gauss <- function(D, phi, qq = NULL) {
  out <- exp(-(D / phi)^2)
  if (!is.null(qq)) out <- qq * out
  out
}

c_fun_gauss <- function(D, sigma, phi, qq = NULL) {
  sigma^2 * rho_fun_gauss(D, phi, qq)
}


# other helpers ---------------------------------------------------------------

#' @importFrom purrr map2_lgl
#' @importFrom rlang set_names

check_if_updated <- function(par, previous_par) {
  out <- map2_lgl(par, previous_par, ~ !isTRUE(all.equal(.x, .y)))
  set_names(out, names(out))
}
