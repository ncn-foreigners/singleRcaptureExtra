#' @title Generating data from foreign objects
#'
#' @description
#' An S3 method for \code{stats::simulate} to handle \code{singleRforeign} and
#' \code{singleRfamily} classes.
#'
#' @param object an object representing a fitted model.
#' @param nsim a numeric scalar specifying:
#' \itemize{
#'    \item number of response vectors to simulate in \code{simulate.singleRStaticCountData}, defaults to \code{1L}.
#'    \item number of units to draw in \code{simulate.singleRfamily}, defaults to \code{NROW(eta)}.
#' }
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param truncated logical value indicating whether to sample from truncated or
#' full distribution.
#' @param ... additional optional arguments.
#' @return a \code{data.frame} with \code{n} rows and \code{nsim} columns.
#' @seealso [stats::simulate()] [singleRcapture::estimatePopsize()]
#' @importFrom stats simulate
#' @method simulate singleRforeign
#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM weights
#' @exportS3Method
#' @name simulate
#' @export
simulate.singleRforeign <- function(object,
                                    nsim = 1,
                                    seed = NULL,
                                    truncated = FALSE,
                                    ...) {
  simulateInternal(object$foreignObject, nsim, seed, truncated)
}

simulateInternal <- function(object, ...)
  UseMethod("simulateInternal")

#' @importFrom stats runif
simulateInternal.zerotrunc <- function(object,
                                       nsim,
                                       seed,
                                       truncated,
                                       ...) {
  n <- object$foreignObject$n
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    stats::runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  val <- countreg::rztnbinom(
    n  = n * nsim,
    mu = rep(exp(model.matrix(object$foreignObject) %*% object$foreignObject$coefficients), nsim),
    theta = object$foreignObject$theta
  )

  if (!truncated) {
    prob0 <- 1 - predict(object$foreignObject, type = "zero")
    val[as.logical(rbinom(n = n * nsim, size = 1, prob = rep(prob0, nsim)))] <- 0
  }

  dim(val) <- c(n, nsim)
  val <- as.data.frame(val)
  names(val) <- paste0("sim_", seq_len(nsim))
  attr(val, "seed") <- RNGstate
  val
}

#' @importFrom stats runif
simulateInternal.vglm <- function(object,
                                  nsim,
                                  seed,
                                  truncated,
                                  ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    stats::runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  eta <- predict(object)
  n   <- nrow(eta)
  val <- matrix(object@family@simslot(object, nsim), ncol = nsim)
  dim(val) <- c(n, nsim)
  val <- as.data.frame(val)
  colnames(val) <- paste0("sim_", seq_len(nsim))

  if (!truncated) {

    PW <- tryCatch(
      expr = {fittedvlm(object, type.fitted = "prob0")},
      error = function(e) {
        if (object@family@vfamily[1] == "oapospoisson") {
          links <- (strsplit(object@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
          lambda <- VGAM::eta2theta(
            eta[, c(FALSE, TRUE), drop = FALSE], links[2],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )

          return(exp(-lambda) / (1 - lambda * exp(-lambda)))
        } else if (object@family@vfamily[1] == "oipospoisson") {
          links <- (strsplit(object@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
          lambda <- VGAM::eta2theta(
            eta[, c(FALSE, TRUE), drop = FALSE], links[2],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          return(exp(-lambda))
        } else {
          stop("TODO")
        }
      }
    )

    val[as.logical(stats::rbinom(
      n = n * nsim, size = 1, prob = rep(PW, length.out = n * nsim)
    )), ] <- 0
  }

  attr(val, "seed") <- RNGstate
  val
}
