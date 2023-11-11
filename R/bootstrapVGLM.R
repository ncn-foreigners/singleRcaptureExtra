# Add multicore
#' @importFrom stats as.formula rpois terms
#' @importFrom VGAM vglm.fit
#' @importFrom stats rbinom
bootVGLM <- function(object,
                     B = 500,
                     trace = FALSE,
                     N, visT = FALSE,
                     bootType = c("parametric",
                                  "semiparametric",
                                  "nonparametric"),
                     ...) {
  strappedStatistic <- vector("numeric", length = B)

  n   <- nobs(object)
  wt  <- object@prior.weights
  mf  <- model.frame(object)
  if (length(wt) != NROW(mf)) {
    wt <- rep(1, NROW(mf))
  }
  eta <- object@predictors
  extra <- object@extra
  extra$type.fitted <- if (object@family@vfamily[1] %in% c("oipospoisson", "oapospoisson")) extra$type.fitted else "prob0"
  # without names it takes less memory
  y <- as.vector(model.response(mf))

  # getting links
  # TODO:: this is beyond moronic
  links <- getLinksBlurb(object@family@blurb)

  offset <- object@offset
  if (any(dim(offset) != dim(object@predictors))) {
    offset <- matrix(0, nrow = NROW(object@predictors),
                        ncol = NCOL(object@predictors))
  }
  cf <- coef(object)

  control <- object@control
  control$noWarning <- TRUE

  if (visT) {
    plot(
      1, type = "n",
      xlab = "Bootstrap sample",
      ylab = "Value of population size estimator",
      main = expression(paste("Plot of values of ", hat(N), " obtained from bootstrap samples")),
      sub = "Points will be added in real time",
      xlim = c(0, B + 1), ylim = c(0, 2 * sum(N))
    )
  }
  k <- 1

  contr <- N
  N <- sum(N)
  while (k <= B) {
    switch (
      bootType,
      "nonparametric" = {
        strap <- sample.int(replace = TRUE, n = n)
        yStrap  <- as.numeric(y[strap])
        wtStrap <- as.numeric(wt[strap])

        offsetStrap  <- offset[strap, , drop = FALSE]
        mfStrap <- mf[strap, , drop = FALSE]

        if (length(object@contrasts)) {
          xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
        } else {
          xStrap <- model.matrix(object@terms, mfStrap)
        }
        etaStrap <- eta[as.numeric(strap), , drop = FALSE]
        oasgn <- attr(xStrap, "assign")
        attr(xStrap, "assign") <- VGAM:::attrassigndefault(xStrap, attr(mfStrap, "terms"))
        attr(xStrap, "orig.assign.lm") <- oasgn
      },
      "parametric" = {
        nn <- floor(N) + stats::rbinom(n = 1, size = 1, prob = N - floor(N))
        strap <- sample.int(replace = TRUE, n = length(y), size = nn, prob = contr)
        wtStrap <- as.numeric(wt[strap])
        etaStrap <- eta[as.numeric(strap), , drop = FALSE]

        offsetStrap  <- offset[strap, , drop = FALSE]
        mfStrap <- mf[strap, , drop = FALSE]

        if (length(object@contrasts)) {
          xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
        } else {
          xStrap <- model.matrix(object@terms, mfStrap)
        }

        # edit this for other fam funcs
        yStrap <- object@extra$singleRcaptureSimulate(
          nn, eta = xStrap %*% cf, links = links
        )

        wtStrap <- wtStrap[yStrap > 0]
        mfStrap <- mfStrap[yStrap > 0, , drop = FALSE]

        etaStrap     <- etaStrap[yStrap > 0, , drop = FALSE]
        offsetStrap  <- offsetStrap[yStrap > 0, , drop = FALSE]

        if (length(object@contrasts)) {
          xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
        } else {
          xStrap <- model.matrix(object@terms, mfStrap)
        }

        yStrap <- yStrap[yStrap > 0]

        oasgn <- attr(xStrap, "assign")
        attr(xStrap, "assign") <- VGAM:::attrassigndefault(xStrap, attr(mfStrap, "terms"))
        attr(xStrap, "orig.assign.lm") <- oasgn
      },
      "semiparametric" = {
        strap <- sum(rbinom(size = 1, n = N, prob = n / N))

        strap <- sample.int(replace = TRUE, n = n, size = strap)
        yStrap  <- as.numeric(y[strap])
        wtStrap <- as.numeric(wt[strap])

        offsetStrap  <- offset[strap, , drop = FALSE]
        mfStrap <- mf[strap, , drop = FALSE]

        if (length(object@contrasts)) {
          xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
        } else {
          xStrap <- model.matrix(object@terms, mfStrap)
        }
        etaStrap <- eta[as.numeric(strap), , drop = FALSE]
        oasgn <- attr(xStrap, "assign")
        attr(xStrap, "assign") <- VGAM:::attrassigndefault(xStrap, attr(mfStrap, "terms"))
        attr(xStrap, "orig.assign.lm") <- oasgn
      }
    )

    theta <- NULL
    suppressWarnings(
      try(
        theta <- vglm.fit(
          y = yStrap, x = xStrap,
          etastart = etaStrap,
          constraints = object@constraints,
          extra = extra, w = wtStrap, offset = offsetStrap,
          qr.arg = TRUE, Terms = object@terms,
          function.name = "vglm", family = object@family,
          control = control,
          Xm2 = object@Xm2, Ym2 = object@Ym2, mustart = NULL
        ),
        silent = FALSE
      )
    )
    k <- k + 1

    if (is.null(theta)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      theta <- theta$predictors
      if (isTRUE(trace)) {
        switch (
          bootType,
          "parametric" = {1},
          "semiparametric" = {1},
          "nonparametric" = {print(summary(theta))}
        )
      }

      if (object@family@vfamily[1] == "oapospoisson") {
        links <- (strsplit(object@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          theta[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )

        est <- sum(wtStrap / (1 - exp(-lambda) / (1 - lambda * exp(-lambda))))
      } else if (object@family@vfamily[1] == "oipospoisson") {
        links <- (strsplit(object@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          theta[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )
        est <- sum(wtStrap / (1 - exp(-lambda)))
      } else {
        est <- sum(wtStrap / (1 - object@family@linkinv(
          eta = theta,
          extra = list("type.fitted" = "prob0")
        )))
      }

      if (visT) graphics::points(k - 1, est, pch = 1)
      if (isTRUE(trace)) cat(" Estimated population size: ", est,"\n",sep = "")
      #if (visT) graphics::points(k - 1, est, pch = 1)

      strappedStatistic[k - 1] <- est
    }
  }

  strappedStatistic
}
