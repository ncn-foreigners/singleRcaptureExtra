#' @importClassesFrom VGAM vglm
#' @importFrom VGAM dtheta.deta eta2theta fittedvlm posnegbinomial pospoisson
#' @importFrom stats weights model.matrix nobs qnorm vcov quantile rnbinom
#' @importFrom singleRcapture estimatePopsize ztoipoisson ztHurdlepoisson
#' @importFrom VGAMdata oapospoisson oipospoisson
#' @rdname foreignMethods
#' @export
estimatePopsize.vglm <- function(formula,
                                 subset   = NULL,
                                 naAction = NULL,
                                 popVar   = c("analytic",
                                              "bootstrap"),
                                 control = controlEstPopVglm(),
                                 derivFunc = NULL,
                                 ...) {
  # Add posbinomial, oiposbinomial
  if (missing(popVar)) popVar <- "analytic"
  sizeObserved <- nobs(formula)
  signlevel <- control$alpha

  ### TODO:: This is terrible, don't use tryCatch if you don't have to
  PW <- tryCatch(
    expr = {fittedvlm(formula, type.fitted = "prob0")},
    error = function(e) {
      if (formula@family@vfamily[1] == "oapospoisson") {
        links <- (strsplit(formula@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          formula@predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )

        return(exp(-lambda) / (1 - lambda * exp(-lambda)))
      } else if (formula@family@vfamily[1] == "oipospoisson") {
        links <- (strsplit(formula@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          formula@predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )
        return(exp(-lambda))
      } else {
        stop("estimatePopsize.vglm method is only for objects with family slots with possibility of type.fitted = prob0 (and a few select ones)")
      }
    }
  )

  wg <- formula@prior.weights
  if (is.null(wg) | !length(wg)) wg <- rep(1, sizeObserved)

  POP <- switch (popVar,
    "analytic" = {
      N <- sum(wg / (1 - PW))
      if (is.null(derivFunc)) {
        if (formula@family@vfamily[1] %in% c("oapospoisson", "pospoisson",
                                             "oipospoisson", "posnegbinomial")) {
          derivFunc <- switch (formula@family@vfamily[1],
            "oapospoisson" = function(eta) {
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

              ### TODO:: this is a terrible way of ordering, fix it !!!!
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
              as.vector(t(cbind(-size*tmp1^size/((lambda+size)*(tmp1^size-1)^2),
                tmp1^size*((size + lambda)*log(tmp1)+lambda)/((lambda + size)*(tmp1^size-1)^2)) *
                cbind(dlambda.deta, dsize.deta)))
            }
          )
        } else if (!missing(derivFunc)) {
          derivFunc(formula)
        } else {
          stop("family slot not recognised, please provide a derivFunc argument with function for computing derivatives of 1/(1-prob0) with respect to linear predictors (or use bootstrap)")
        }
      }

      ### TODO:: vcov for oa/oi pospoisson is different in singleRcapture
      dd <- t(rep(wg, formula@family@infos()$M1) * derivFunc(formula@predictors)) %*% model.matrix(formula, type = "vlm")
      #print(dd)

      variation <- sum(wg * PW / (1 - PW)^2)
      vc <- vcov(formula)

      #print(dd %*% vc %*% t(dd))
      #print(variation)
      variation <- variation + dd %*% vc %*% t(dd)

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
          control = control
        ),
        class = "popSizeEstResults"
      )},
    "bootstrap" = {
      N <- wg / (1 - PW)

      # Asign sampling functions for each distribution
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

      # Bootstrap call
      if (control$cores == 1) {
        strappedStatistic <- bootVGLM(
          formula,
          B = control$B,
          trace = control$traceBootstrapSize,
          N = N,
          visT = control$bootstrapVisualTrace,
          bootType = control$bootType
        )
      } else {
        strappedStatistic <- multiCoreBootVGLM(
          formula,
          B = control$B,
          trace = control$traceBootstrapSize,
          N = N,
          visT = control$bootstrapVisualTrace,
          bootType = control$bootType,
          cores = control$cores
        )
      }
      N <- sum(N)


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
        sd <- sd / (sqrt(2 / (sizeObserved - 1)) * exp(lgamma(sizeObserved / 2) - lgamma((sizeObserved- 1) / 2)))
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
    }
  )

  structure(
    list(
      foreignObject  = formula,
      call           = match.call(),
      sizeObserved   = sizeObserved,
      populationSize = POP,
      pacakgeInfo    = "VGAM::vglm"
    ),
    class = c("singleRforeign", "singleRStaticCountData", "singleR")
  )
}
