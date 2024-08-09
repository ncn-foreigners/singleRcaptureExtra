# TODO
simulate.singleRforeign <- function(object,
                                    nsim = 1,
                                    seed = NULL,
                                    truncated = FALSE,
                                    ...) {
  simulateInternal(object$foreignObject, nsim, seed, truncated)
}

simulateInternal <- function(object, ...)
  UseMethod("simulateInternal")


simulateInternal.zerotrunc <- function(object,
                                       nsim,
                                       seed,
                                       truncated,
                                       ...) {
  n <- object$foreignObject$n
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
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
  print(val)

  if (!truncated) {
    prob0 <- 1 - predict(object$foreignObject, type = "zero")
    val[rbinom(n = n * nsim, size = 1, prob = rep(prob0, nsim))] <- 0
  }

  dim(val) <- c(n, nsim)
  val <- as.data.frame(val)
  names(val) <- paste0("sim_", seq_len(nsim))
  val
}
