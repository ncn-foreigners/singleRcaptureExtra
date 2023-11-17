#' @importFrom countreg zerotrunc
#' @importFrom stats var
#' @importFrom stats nobs
#' @importFrom stats weights
#' @importFrom stats quantile
#' @importFrom stats coef
#' @importFrom stats predict
#' @importFrom stats model.frame
#' @importFrom stats model.response
#' @rdname foreignMethods
#' @export
estimatePopsize.zerotrunc <- function(formula,
                                      subset   = NULL,
                                      naAction = NULL,
                                      popVar   = c("analytic",
                                                   "bootstrap"),
                                      control = NULL,
                                      ...) {
  ### TODO adjust for weights being counts
  if (is.null(control)) control <- controlEstPopCountreg()
  if (missing(popVar)) popVar <- "analytic"

  sizeObserved <- stats::nobs(formula)
  eta <- cbind(as.numeric(log(predict(formula, type = "count"))))
  signlevel <- control$alpha
  sc <- qnorm(p = 1 - signlevel / 2)

  if (formula$dist == "negbin") {
    eta <- cbind(eta, log(formula$theta))
  }

  wg <- stats::weights(formula)
  if (is.null(wg)) wg <- rep(1, sizeObserved)
  wg <- as.numeric(wg)

  offset <- formula$offset
  if (is.null(offset)) {
    offset <- switch (formula$dist,
      "poisson"   = cbind(rep(0, sizeObserved)),
      "negbin"    = cbind(rep(0, sizeObserved), rep(0, sizeObserved)),
      "geometric" = cbind(rep(0, sizeObserved))
    )
  }

  trcount <- control$trcount

  family <- switch (formula$dist,
    "poisson"   = singleRcapture::ztpoisson(),
    "negbin"    = singleRcapture::ztnegbin(alphaLink = "neglog"),
    "geometric" = singleRcapture::ztgeom()
  )

  POP <- switch(popVar,
    "analytic" = {
      strappedStatistic <- "No bootstrap performed"
      N <- family$pointEst(pw = wg, eta = eta) + trcount
      X <- if (formula$dist != "negbin") model.matrix(formula) else
        singleRcapture:::singleRinternalGetXvlmMatrix(
          parNames = c("lambda", "theta"),
          formulas = list(formula$formula, ~ 1),
          X        = model.frame(formula)
        )

      variation <- as.numeric(family$popVar(
        eta = eta,
        pw  = wg,
        cov = if (formula$dist != "negbin") vcov(formula) else {
          solve(-family$makeMinusLogLike(
            y = model.response(model.frame(formula)),
            X = X,
            weight = wg,
            deriv = 2,
            offset = offset
          )(c(stats::coef(formula), log(formula$theta))))
        },
        Xvlm = X
      ))

      if (!is.finite(variation))
        stop("Computed variance is infinite/NaN/NULL")
      sd <- sqrt(variation)
      G <- exp(sc * sqrt(log(1 + variation / ((N - sizeObserved) ^ 2))))

      confidenceInterval <- data.frame(t(data.frame(
        "normal" = c(lowerBound = max(N - sc * sd, sizeObserved),
                     upperBound = N + sc * sd),
        "logNormal" = c(lowerBound = sizeObserved + max((N - sizeObserved) / G, 0),
                        upperBound = sizeObserved + (N - sizeObserved) * G)
      )))

      structure(
        list(
          pointEstimate = N,
          variance = variation,
          confidenceInterval = confidenceInterval,
          boot = NULL,
          control = control
        ),
        class = "popSizeEstResults"
      )
    },
    "bootstrap" = {
      N <- family$pointEst(pw = wg, eta = eta) + trcount
      ### Assigning necessary values
      hwm <- switch (formula$dist,
        "poisson"   = c(length(formula$coefficients)),
        "negbin"    = c(length(formula$coefficients), 1),
        "geometric" = c(length(formula$coefficients))
      )

      formulas <- list(formula$formula)
      if (formula$dist == "negbin") formulas <- append(formulas, ~ 1)

      modelFrame <- model.frame(formula)
      y <- as.vector(model.response(modelFrame))
      X <- model.matrix(formula)


      if (control$cores > 1) {
        funBoot <- switch(
          control$bootType,
          "parametric"     = singleRcapture:::parBootMultiCore,
          "semiparametric" = singleRcapture:::semparBootMultiCore,
          "nonparametric"  = singleRcapture:::noparBootMultiCore
        )
        strappedStatistic <- funBoot(
          y = y, X = X, hwm = hwm, eta = eta, N = N,
          family   = family,
          formulas = formulas,
          beta     = if (formula$dist != "negbin") stats::coef(formula) else
            c(stats::coef(formula), log(formula$theta)),
          weights  = wg,
          trcount  = trcount,
          numboot  = control$B,
          cores    = control$cores,
          method   = control$fittingMethod,
          controlBootstrapMethod = control$bootstrapFitcontrol,
          Xvlm = if (formula$dist != "negbin") X else
            singleRcapture:::singleRinternalGetXvlmMatrix(
              parNames = c("lambda", "theta"),
              formulas = singleRcapture:::singleRinternalMergeFormulas(
                list(formula$formula, ~ 1)
              ), X = modelFrame
            ),
          modelFrame = modelFrame,
          offset     = offset
        )
      } else {
        funBoot <- switch(
          control$bootType,
          "parametric"     = singleRcapture:::parBoot,
          "semiparametric" = singleRcapture:::semparBoot,
          "nonparametric"  = singleRcapture:::noparBoot
        )
        strappedStatistic <- funBoot(
          y = y, X = X, hwm = hwm, eta = eta, N = N,
          family   = family,
          formulas = formulas,
          beta     = if (formula$dist != "negbin") stats::coef(formula) else
            c(stats::coef(formula), log(formula$theta)),
          weights  = wg,
          trcount  = trcount,
          numboot  = control$B,
          method   = control$fittingMethod,
          controlBootstrapMethod = control$bootstrapFitcontrol,
          Xvlm = if (formula$dist != "negbin") X else
            singleRcapture:::singleRinternalGetXvlmMatrix(
              parNames = c("lambda", "theta"),
              formulas = singleRcapture:::singleRinternalMergeFormulas(
                list(formula$formula, ~ 1)
              ), X = modelFrame
            ),
          modelFrame = modelFrame,
          offset     = offset,
          visT       = control$bootstrapVisualTrace,
          trace      = control$traceBootstrapSize
        )
      }
      if (N < stats::quantile(strappedStatistic, .05)) {
        warning("bootstrap statistics unusually high, try higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)\n")
      } else if (N > stats::quantile(strappedStatistic, .95)) {
        warning("bootstrap statistics unusually low, try higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)\n")
      }
      if (max(strappedStatistic) > N ^ 1.5) {
        warning("Outlier(s) in statistics from bootstrap sampling detected, consider higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)\n")
      }

      variation <- stats::var(strappedStatistic)

      if (!is.finite(variation))
        stop("Computed variance is infinite/NaN/NULL")
      sd <- sqrt(variation)
      if (control$confType == "percentilic") {
        confidenceInterval <- stats::quantile(strappedStatistic,
                                              c(signlevel / 2,
                                                1 - signlevel / 2))
        names(confidenceInterval) <- c("lowerBound", "upperBound")
      } else if (control$confType == "normal") {
        G <- exp(sc * sqrt(log(1 + variation / ((N - sizeObserved) ^ 2))))

        confidenceInterval <- data.frame(t(data.frame(
          "normal" = c(lowerBound = max(N - sc * sd,
                                        sizeObserved),
                       upperBound = N + sc * sd),
          "logNormal" = c(lowerBound = max(sizeObserved + (N - sizeObserved) / G,
                                           sizeObserved),
                          upperBound = sizeObserved + (N - sizeObserved) * G)
        )))
      } else if (control$confType == "basic") {
        confidenceInterval <- 2 * N - stats::quantile(strappedStatistic,
                                                      c(1 - signlevel / 2,
                                                        signlevel / 2))
        names(confidenceInterval) <- c("lowerBound", "upperBound")
      }

      structure(
        list(
          pointEstimate = N,
          variance = variation,
          confidenceInterval = confidenceInterval,
          boot = if (isTRUE(control$keepbootStat)) strappedStatistic else NULL,
          control = control
        ),
        class = "popSizeEstResults"
      )
    }
  )

  structure(
    list(
      foreignObject   = formula,
      call            = match.call(),
      sizeObserved    = sizeObserved,
      populationSize  = POP,
      pacakgeInfo     = "countreg::zerotrunc"
    ),
    class = c("singleRforeign", "singleRStaticCountData", "singleR")
  )
}
