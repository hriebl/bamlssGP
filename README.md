
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bamlssGP

> Gaussian process responses for the bamlss package

<!-- badges: start -->

<!-- badges: end -->

This R package implements different Gaussian processes as response
distribution families for the bamlss package, i.e. for Bayesian additive
models for location, scale, and shape. It also provides an efficient
MCMC sampler for these models.

Take a look at the help pages below to get a better idea of what the
package does.

## Installation

You can install the development version of bamlssGP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hriebl/bamlssGP")
```

## Help pages

### Gaussian process distribution families for bamlss

#### Description

These functions set up Gaussian process (GP) distribution families for
bamlss.

#### Usage

``` r
gp_constant_bamlss(y, D, initial_par, ...)
gp_weibull_bamlss(y, s, Z = NULL, D, q = NULL, initial_par, ...)
gp_sphere_bamlss(y, s, D, initial_par, ...)
```

#### Arguments

|               |                                                                                                                     |
| ------------- | ------------------------------------------------------------------------------------------------------------------- |
| `y`           | The GPs to be used as response observations. A list of numeric vectors.                                             |
| `D`           | The distance matrices of the coordinates. A list of numeric matrices.                                               |
| `initial_par` | The initial distributional parameters. A data frame with one column per distributional parameter.                   |
| `...`         | Not used.                                                                                                           |
| `s`           | The coordinates at which the GPs were observed. A list of numeric vectors or matrices.                              |
| `Z`           | Optionally, the design matrices for an additive linear mean function. A list of numeric matrices.                   |
| `q`           | Optionally, the scaling factors of the standard deviation of the GPs at the coordinates. A list of numeric vectors. |

#### Details

`gp_constant_bamlss` uses a constant mean function and a Matérn
covariance function with the smoothness parameter *ν* = 1.5. The
distributional parameters of `gp_constant_bamlss` are

  - the **mean**,
  - the **stddev** (standard deviation),
  - and the **range**.

`gp_weibull_bamlss` uses a Weibull growth curve as a mean function and a
Matérn covariance function with the smoothness parameter *ν* = 1.5. The
distributional parameters of `gp_weibull_bamlss` are

  - the **limit** (of the Weibull growth curve),
  - the **shape** (of the Weibull growth curve),
  - the **scale** (of the Weibull growth curve),
  - the **stddev** (standard deviation),
  - and the **range**.

Optionally, a linear mean function can be added to the Weibull growth
curve by passing a list of design matrices `Z` to `gp_weibull_bamlss`.
Each column of the design matrices adds one distributional parameter
with the corresponding column name to the family.

`gp_sphere_bamlss` uses a mean function on a sphere and an exponential
covariance function. It can be used to model tree crowns. The
distributional parameters of `gp_sphere_bamlss` are

  - the **radius** (of the tree crown),
  - the **south**-bound expansion (of the tree crown),
  - the **height** (of the tree crown),
  - the **stddev** (standard deviation),
  - and the **range**.

#### Value

A `family.bamlss` object.

#### Examples

``` r
## Not run: 
N <- 30
n <- 30

mean <- 0
stddev <- exp(0)
range <- exp(0)

s <- seq.int(0, 1, length.out = n)
D <- as.matrix(dist(s))

mu <- rep.int(mean, n)
Sigma <- stddev^2 * (1 + D / range) * exp(-D / range)
y <- MASS::mvrnorm(N, mu, Sigma)

y <- unlist(apply(y, 1, list), recursive = FALSE)
D <- rep.int(list(D), N)

initial_par <- data.frame(
  mean = rep.int(mean, N),
  stddev = rep.int(stddev, N),
  range = rep.int(range, N)
)

family <- gp_constant_bamlss(y, D, initial_par)
model <- bamlss(rep.int(0, N) ~ 1, family)
joint <- sample_bamlss_gp(model)

## End(Not run)
```

### Sample from the posterior of a bamlss model with a Gaussian process response

#### Description

This function provides an alternative MCMC sampler for bamlss models
with Gaussian process (GP) responses. It is typically more efficient
than the default bamlss sampler by sampling the parametric predictors
for the standard deviation and the range jointly, thus reducing the
autocorrelation of the chains.

#### Usage

``` r
sample_bamlss_gp(model, n = 1000, only_shared = FALSE)
```

#### Arguments

|               |                                                                                                                                                                                                |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `model`       | A bamlss model with a GP response.                                                                                                                                                             |
| `n`           | The sample size/length of the chain to be drawn.                                                                                                                                               |
| `only_shared` | If `TRUE`, only the parameters for those covariates which appear in the parametric predictors for both the standard deviation and the range are sampled jointly. Not helpful in my experience. |

#### Value

A list with the elements `model` and `samples`.

#### Examples

``` r
## Not run: 
N <- 30
n <- 30

mean <- 0
stddev <- exp(0)
range <- exp(0)

s <- seq.int(0, 1, length.out = n)
D <- as.matrix(dist(s))

mu <- rep.int(mean, n)
Sigma <- stddev^2 * (1 + D / range) * exp(-D / range)
y <- MASS::mvrnorm(N, mu, Sigma)

y <- unlist(apply(y, 1, list), recursive = FALSE)
D <- rep.int(list(D), N)

initial_par <- data.frame(
  mean = rep.int(mean, N),
  stddev = rep.int(stddev, N),
  range = rep.int(range, N)
)

family <- gp_constant_bamlss(y, D, initial_par)
model <- bamlss(rep.int(0, N) ~ 1, family)
joint <- sample_bamlss_gp(model)

## End(Not run)
```
