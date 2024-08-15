#' @importClassesFrom VGAM vgam
#' @importFrom VGAM dtheta.deta
#' @importFrom VGAM eta2theta
#' @importFrom VGAM fittedvlm
#' @importFrom stats weights
#' @importFrom stats model.matrix
#' @importFrom stats nobs
#' @importFrom stats qnorm
#' @importFrom stats vcov
#' @importFrom singleRcapture estimatePopsize
#' @importFrom VGAMdata oapospoisson
#' @importFrom VGAMdata oipospoisson
#' @importFrom VGAM pospoisson
#' @importFrom VGAM posnegbinomial
#' @rdname foreignMethods
#' @export
estimatePopsize.vgam <- function(formula,
                                 subset   = NULL,
                                 naAction = NULL,
                                 popVar   = c("analytic",
                                              "bootstrap",
                                              "no"),
                                 control = controlEstPopVgam(),
                                 derivFunc = NULL,
                                 ...) {
  if (missing(popVar)) popVar <- "analytic"

  if (popVar == "bootstrap" && is.null(control$data)) {
    control$data <- tryCatch(
      expr = {eval(formula@call$data, envir = parent.frame())},
      error = function (e) {
        stop("Bootstrap with vgam class needs a data argument in control.vgam to work. An attempt was made to manually evaluate data slot from vgam class object call but it failed.")
      }
    )
    if (is.null(control$data) | !is.data.frame(control$data)) {
      stop("Bootstrap with vgam class needs a data argument in control.vgam to work. An attempt was made to manually evaluate data slot from vgam class object call but it failed.")
    } else {
      message("Bootstrap with vgam class needs a data argument in control.vgam to work. This data argument was manually evaluated from call of vgam class object.")
    }
  }
  sizeObserved <- nobs(formula)
  signlevel <- control$alpha
  trcount <- control$trcount
  numboot <- control$B

  PW <- getPw(object = formula)
  wg <- formula@prior.weights
  if (is.null(wg) | !length(wg)) wg <- rep(1, sizeObserved)

  POP <- switch (
    popVar,
    "analytic" = {
      N <- sum(wg * PW)
      if (is.null(derivFunc)) {
        if (formula@family@vfamily[1] %in% c("oapoisson", "pospoisson",
                                             "oipospoisson", "posnegbinomial",
                                             "posbinomial", "oiposbinomial")) {
          derivFunc <- switch (formula@family@vfamily[1],
            "oapoisson" = function(eta) {
              links <- (strsplit(formula@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
              TFvec <- c(TRUE, FALSE)

              lambda <- VGAM::eta2theta(
                eta[, !TFvec, drop = FALSE], links[2],
                list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE)
              )
              pobs1 <- VGAM::eta2theta(
                eta[, TFvec, drop = FALSE], links[1],
                list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE)
              )

              dlambda.deta <- VGAM::dtheta.deta(
                lambda, links[2], earg = list(
                  theta = NULL, bvalue = NULL,
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE
                )
              )

              dpobs1.deta <- VGAM::dtheta.deta(
                pobs1, links[1], earg = list(
                  theta = NULL,
                  bvalue = NULL, inverse = FALSE,
                  deriv = 0, short = TRUE,
                  tag = FALSE
                )
              )

              as.vector(t(
                cbind(dpobs1.deta, dlambda.deta) *
                  cbind(rep(0, NROW(eta)), (exp(lambda) - 1) / (exp(lambda) - lambda - 1) ^ 2)
              ))
            },
            "pospoisson" = function(eta) {
              links <- (strsplit(formula@family@blurb[c(3)], split = "\\(") |> unlist())[c(1)]
              lambda <- eta2theta(
                eta, "loglink",
                list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE)
              )
              dlambda.deta <- VGAM::dtheta.deta(
                lambda, links[1], earg = list(
                  theta = NULL, bvalue = NULL,
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE
                )
              )
              (exp(lambda) / (exp(lambda) - 1) ^ 2) * dlambda.deta
            },
            "oipospoisson" = function(eta) {
              links <- (strsplit(formula@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
              TFvec <- c(TRUE, FALSE)

              lambda <- VGAM::eta2theta(
                eta[, !TFvec, drop = FALSE], links[2],
                list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE)
              )
              dlambda.deta <- VGAM::dtheta.deta(
                lambda, links[2], earg = list(
                  theta = NULL, bvalue = NULL,
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE
                )
              )

              as.vector(t(cbind(rep(0, NROW(eta)), (exp(lambda) / (exp(lambda) - 1) ^ 2) * dlambda.deta)))
            },
            "posnegbinomial" = function(eta) {
              links <- (strsplit(formula@family@blurb[c(3,5)], split = "\\(") |> unlist())[c(1,3)]
              TFvec <- c(TRUE, FALSE)

              lambda <- VGAM::eta2theta(
                eta[, TFvec, drop = FALSE], links[1],
                list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE)
              )
              size <- VGAM::eta2theta(
                eta[, !TFvec, drop = FALSE], links[2],
                list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE)
              )

              dlambda.deta <- VGAM::dtheta.deta(
                lambda, links[1], earg = list(
                  theta = NULL, bvalue = NULL,
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE
                )
              )

              dsize.deta <- VGAM::dtheta.deta(
                size, links[2], earg = list(
                  theta = NULL,
                  bvalue = NULL, inverse = FALSE,
                  deriv = 0, short = TRUE,
                  tag = FALSE
                )
              )
              tmp1 <- size / (lambda + size)
              as.vector(t(
                cbind(-size*tmp1^size/((lambda+size)*(tmp1^size-1)^2),
                      tmp1^size*((size + lambda)*log(tmp1)+lambda)/((lambda + size)*(tmp1^size-1)^2)) *
                  cbind(dlambda.deta, dsize.deta)
              ))
            }
          )
        } else {
          stop("family slot not recognised, please provide a derivFunc argument with function for computing derivatives of 1/(1-prob0) with respect to linear predictors (or use bootstrap)")
        }
      }

      dd <- t(rep(wg, formula@family@infos()$M1) * derivFunc(formula@predictors)) %*% model.matrix(formula, type = "vlm")

      variation <- sum(wg * (1 - PW ^ -1) * PW ^ 2)
      variation <- variation + dd %*% vcov(formula) %*% t(dd)

      sc <- qnorm(p = 1 - signlevel / 2)
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
          control = list(alpha = signlevel)
        ),
        class = "popSizeEstResults"
      )},
    "bootstrap" = {
      N <- sum(wg * PW)

      formula@extra$singleRcaptureSimulate <- switch(formula@family@vfamily[1],
        "oapospoisson" = function(n, eta, links) {
          lambda <- VGAM::eta2theta(
            eta[, c(FALSE, TRUE), drop = FALSE], links[2],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          pobs1 <- VGAM::eta2theta(
            eta[, !c(FALSE, TRUE), drop = FALSE], links[1],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          singleRcapture::ztHurdlepoisson(lambdaLink = "log",
                          piLink = "logit")$simulate(
            n = n, lower = -1,
            eta = cbind(log(lambda),
                        VGAM::logitlink(pobs1))
          )
        },
        "pospoisson" = function(n, eta, links) {
          lambda <- VGAM::eta2theta(
            eta, links,
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          stats::rpois(n = n, lambda = lambda)
        },
        "oipospoisson" = function(n, eta, links) {
          lambda <- VGAM::eta2theta(
            eta[, c(FALSE, TRUE), drop = FALSE], links[2],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          pstr1 <- VGAM::eta2theta(
            eta[, !c(FALSE, TRUE), drop = FALSE], links[1],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          singleRcapture::ztoipoisson(lambdaLink = "log",
                          omegaLink = "logit")$simulate(
            n = n, lower = -1,
            eta = cbind(log(lambda),
                        VGAM::logitlink(pstr1))
          )
        },
        "posnegbinomial" = function(n, eta, links) {
          lambda <- VGAM::eta2theta(
            eta[, !c(FALSE, TRUE), drop = FALSE], links[1],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          size <- VGAM::eta2theta(
            eta[, c(FALSE, TRUE), drop = FALSE], links[2],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          stats::rnbinom(n = n, size = size, mu = lambda)
        }
      )

      if (control$cores == 1) {
        strappedStatistic <- bootVGAM(
          formula,
          B = control$B,
          trace = control$traceBootstrapSize,
          N = N,
          visT = control$bootstrapVisualTrace,
          bootType = control$bootType,
          data = control$data
        )
      } else {
        strappedStatistic <- multiCoreBootVGAM(
          formula,
          B = control$B,
          trace = control$traceBootstrapSize,
          N = N,
          visT = control$bootstrapVisualTrace,
          bootType = control$bootType,
          cores = control$cores,
          data = control$data
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
      if (control$sd == "normalMVUE") {
        sd <- sd / (sqrt(2 / (sizeObserved - 1)) *
              exp(lgamma(sizeObserved / 2) - lgamma((sizeObserved- 1) / 2)))
      }

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
          boot = if (control$keepbootStat) strappedStatistic else NULL,
          control = control
        ),
        class = "popSizeEstResults"
      )
    },
    "no" = {
      N <- sum(wg * PW)
      structure(list(
        pointEstimate = N
      ), class = "popSizeEstResults")
    }
  )

  structure(
    list(
      foreignObject  = formula,
      call           = match.call(),
      sizeObserved   = sizeObserved,
      populationSize = POP,
      derivFunc      = derivFunc
    ),
    class = c("singleRadditive", "singleRforeign", "singleRStaticCountData", "singleR")
  )
}
