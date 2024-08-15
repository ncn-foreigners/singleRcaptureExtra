#' @importFrom VGAM vgam
bootVGAM <- function(object,
                     B = 500,
                     trace = FALSE,
                     N, visT = FALSE,
                     bootType = c("parametric",
                                  "semiparametric",
                                  "nonparametric"),
                     data = NULL,
                     ...) {
  if (visT) {
    plot(
      1, type = "n",
      xlab = "Bootstrap sample",
      ylab = "Value of population size estimator",
      main = expression(paste("Plot of values of ", hat(N), " obtained from bootstrap samples")),
      sub = "Points will be added in real time",
      xlim = c(0, B + 1), ylim = c(0, 2 * N)
    )
  }

  wt <- object@prior.weights
  n  <- nobs(object)
  if (length(wt) != n) {
    wt <- rep(1, n)
  }
  cll <- object@call
  strappedStatistic <- numeric(B)

  eta <- object@predictors

  # pre actions for all types of bootstrap
  cll$weights <- as.symbol("wtStrap")
  cll$data <- as.symbol("mf")
  cll$etastart <- as.symbol("etaStrap")

  if (!is.null(cll$offset)) {
    cll$offset <- as.symbol("offsetStrap")
    offset <- object@offset
  }

  if (bootType == "parametric") {
    y <- as.integer(model.response(model.frame(object)))

    # getting links
    links <- getLinksBlurb(object@family@blurb)
    # getting contribution
    contr <- getPw(object)
  }

  k <- 1
  while (k <= B) {
    strap <- switch (bootType,
      "nonparametric" = {
        strap <- sample.int(replace = TRUE, n = n)
        mf  <- data[strap, , drop = FALSE]
        wtStrap <- as.numeric(wt[strap])
        etaStrap <- eta[strap, , drop = FALSE]

        if (!is.null(cll$offset)) {
          offsetStrap <- offset[strap, , drop = FALSE]
        }
      },
      "semiparametric" = {
        strap <- sum(rbinom(size = 1, n = N, prob = n / N))
        strap <- sample.int(replace = TRUE, n = n, size = strap)

        mf  <- data[strap, , drop = FALSE]
        wtStrap <- as.numeric(wt[strap])
        etaStrap <- eta[strap, , drop = FALSE]

        if (!is.null(cll$offset)) {
          offsetStrap <- offset[strap, , drop = FALSE]
        }
      },
      "parametric" = {
        nn <- floor(N) + stats::rbinom(n = 1, size = 1, prob = N - floor(N))
        strap <- sample.int(replace = TRUE, n = length(y),
                            size = nn, prob = wt * contr)

        wtStrap <- as.numeric(wt[strap])
        etaStrap <- eta[as.numeric(strap), , drop = FALSE]
        mf <- data[strap, , drop = FALSE]

        if (!is.null(cll$offset)) {
          offsetStrap  <- offset[strap, , drop = FALSE]
        }

        yStrap <- object@extra$singleRcaptureSimulate(
          nn, eta = etaStrap + if (!is.null(cll$offset)) offset else 0,
          links = links
        )

        cond <- yStrap > 0

        yStrap   <- yStrap[cond]
        wtStrap  <- wtStrap[cond]
        etaStrap <- etaStrap[cond, , drop = FALSE]
        mf <- mf[cond, , drop = FALSE]

        if (!is.null(cll$offset)) {
          offsetStrap <- offsetStrap[cond, , drop = FALSE]
        }

        mf[, as.character(cll$formula[[2]])] <- yStrap
      }
    )

    est <- NULL
    try(
      suppressWarnings(
        est <- eval(cll, envir = environment())
      ),
      silent = TRUE
    )

    k <- k + 1
    #trace <- TRUE
    if (is.null(est)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      strappedStatistic[k - 1] <- sum(wtStrap * getPw(est))

      if (visT) graphics::points(k - 1, strappedStatistic[k - 1], pch = 1)
      if (isTRUE(trace)) {
        cat("Bootstrap iteration nr. ", k - 1,
            "\nEstimated population size: ", strappedStatistic[k - 1],
            ", Observed population: ", nrow(etaStrap),
            "\nSummary of additive predictors:",
            "\n", sep = "")
        print(summary(est@predictors))
        cat("-----------------------\n")
      }
      #if (visT) graphics::points(k - 1, est, pch = 1)
    }
  }

  strappedStatistic
}


multiCoreBootVGAM <- function(object,
                              B = 500,
                              trace = FALSE,
                              N, visT = FALSE,
                              bootType = c("parametric",
                                           "semiparametric",
                                           "nonparametric"),
                              data = NULL,
                              cores,
                              ...) {
  cl <- parallel::makeCluster(cores)
  links <- getLinksBlurb(object@family@blurb)

  wt <- object@prior.weights
  n  <- nobs(object)
  if (length(wt) != n) {
    wt <- rep(1, n)
  }
  cll <- object@call
  strappedStatistic <- numeric(B)

  eta <- object@predictors

  # pre actions for all types of bootstrap
  cll$weights <- as.symbol("wtStrap")
  cll$data <- as.symbol("mf")
  cll$etastart <- as.symbol("etaStrap")

  offset <- object@offset
  if (!is.null(cll$offset)) {
    cll$offset <- as.symbol("offsetStrap")
  }
  y <- as.integer(model.response(model.frame(object)))

  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  strappedStatistic <- foreach::`%dopar%`(
    obj = foreach::foreach(
      k = 1:B,
      .combine = c,
      .packages = c("VGAM", "singleRcapture"),
      .export = c(
        "getPw", "getPw.vglm", "links", "cll", "object", "y",
        "eta", "data", "wt", "n", "offset", "mf", "bootType"
      )
    ),
    ex  = {
      est <- NULL
      while (is.null(est)) {
        switch (bootType,
          "nonparametric" = {
            strap <- sample.int(replace = TRUE, n = n)
            mf  <- data[strap, , drop = FALSE]
            wtStrap <- as.numeric(wt[strap])
            etaStrap <- eta[strap, , drop = FALSE]

            if (!is.null(cll$offset)) {
              offsetStrap <- offset[strap, , drop = FALSE]
            }
          },
          "semiparametric" = {
            strap <- sum(rbinom(size = 1, n = N, prob = n / N))
            strap <- sample.int(replace = TRUE, n = n, size = strap)

            mf  <- data[strap, , drop = FALSE]
            wtStrap <- as.numeric(wt[strap])
            etaStrap <- eta[strap, , drop = FALSE]

            if (!is.null(cll$offset)) {
              offsetStrap <- offset[strap, , drop = FALSE]
            }
          },
          "parametric" = {
            contr <- getPw(object)
            N <- sum(contr * wt)

            nn <- floor(N) + stats::rbinom(n = 1, size = 1, prob = N - floor(N))
            strap <- sample.int(replace = TRUE, n = length(y),
                                size = nn, prob = wt * contr)

            wtStrap <- as.numeric(wt[strap])
            etaStrap <- eta[as.numeric(strap), , drop = FALSE]
            mf <- data[strap, , drop = FALSE]

            if (!is.null(cll$offset)) {
              offsetStrap  <- offset[strap, , drop = FALSE]
            }

            yStrap <- object@extra$singleRcaptureSimulate(
              nn, eta = etaStrap + if (!is.null(cll$offset)) offset else 0,
              links = links
            )

            cond <- yStrap > 0

            yStrap   <- yStrap[cond]
            wtStrap  <- wtStrap[cond]
            etaStrap <- etaStrap[cond, , drop = FALSE]
            mf <- mf[cond, , drop = FALSE]

            if (!is.null(cll$offset)) {
              offsetStrap <- offsetStrap[cond, , drop = FALSE]
            }

            mf[, as.character(cll$formula[[2]])] <- yStrap
          }
        )

        try(
          est <- eval(cll, envir = environment()),
          silent = FALSE
        )
      }

      sum(wtStrap * getPw(est))
    }
  )

  strappedStatistic
}
