#' @importClassesFrom VGAM vglm
#' @importFrom VGAM dtheta.deta
#' @importFrom VGAM eta2theta
#' @importFrom VGAM fittedvlm
#' @importFrom stats weights
#' @importFrom stats model.matrix
#' @importFrom stats nobs
#' @importFrom stats qnorm
#' @importFrom stats vcov
#' @importFrom singleRcapture estimatePopsize
#' @rdname foreignMethods
#' @export
estimatePopsize.vglm <- function(formula,
                                 subset   = NULL,
                                 naAction = NULL,
                                 popVar   = c("analytic",
                                              "bootstrap"),
                                 controlForeign = NULL,
                                 derivFunc = NULL,
                                 ...) {
  if (missing(popVar)) popVar <- "analytic"
  sizeObserved <- nobs(formula)
  # siglevel <- controlForeign$alpha
  # trcount <- controlForeign$trcount
  # numboot <- controlForeign$B
  siglevel <- .05
  trcount <- 0
  numboot <- 500

  ### Get try-catch for that
  tryCatch(
    expr = {PW <- fittedvlm(formula, type.fitted = "prob0")},
    error = function(e) {
      stop("estimatePopsize.vglm method is only for objects with family slots with possibility of type.fitted = prob0")
    }
  )
  wg <- weights(formula, type = "prior")

  POP <- switch (popVar,
    "analytic" = {
      N <- sum(wg / (1 - PW))
      if (is.null(derivFunc)) {
        if (formula@family@vfamily[1] %in% c("oapoisson", "pospoisson",
                                             "oipospoisson", "posnegbinomial",
                                             "posbinomial", "oiposbinomial")) {
          derivFunc <- switch (formula@family@vfamily[1],
            "oapoisson" = function(eta) {

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
              # lambda created by eval
              (exp(lambda) / (exp(lambda) - 1) ^ 2) * dlambda.deta
            },
            "oipospoisson" = function(eta) {

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
              c(-size*tmp1^size/((lambda+size)*(tmp1^size-1)^2),
                tmp1^size*((size + lambda)*log(tmp1)+lambda)/((lambda + size)*(tmp1^size-1)^2)) *
                c(dlambda.deta, dsize.deta)
            },
            "posbinomial" = function(eta) {

            },
            "oiposbinomial" = function(eta) {

            }
          )
        } else {
          stop("family slot not recognised, please provide a derivFunc argument with function for computing derivatives of 1/(1-prob0) with respect to linear predictors (or use bootstrap)")
        }
      }

      dd <- t(rep(wg, formula@family@infos()$M1) * derivFunc(formula@predictors)) %*% model.matrix(formula, type = "vlm")

      variation <- sum(wg * PW / (1 - PW)^2)
      variation <- variation + dd %*% vcov(formula) %*% t(dd)

      sc <- qnorm(p = 1 - siglevel / 2)
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
          control = controlForeign
        ),
        class = "popSizeEstResults"
      )},
    "bootstrap" = {
      N <- sum(wg / (1 - PW))

    }
  )

  structure(
    list(
      foreignObject  = formula,
      call           = match.call(),
      sizeObserved   = sizeObserved,
      populationSize = POP
    ),
    class = c("singleRforeign", "singleRStaticCountData", "singleR")
  )
}

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
#' @rdname foreignMethods
#' @export
estimatePopsize.vgam <- function(formula,
                                 subset   = NULL,
                                 naAction = NULL,
                                 popVar   = c("analytic",
                                              "bootstrap"),
                                 controlForeign = NULL,
                                 derivFunc = NULL,
                                 ...) {
  if (missing(popVar)) popVar <- "analytic"
  sizeObserved <- nobs(formula)
  # siglevel <- controlForeign$alpha
  # trcount <- controlForeign$trcount
  # numboot <- controlForeign$B
  siglevel <- .05
  trcount <- 0
  numboot <- 500

  ### Get try-catch for that
  tryCatch(
    expr = {PW <- fittedvlm(formula, type.fitted = "prob0")},
    error = function(e) {
      stop("estimatePopsize.vgam method is only for objects with family slots with possibility of type.fitted = prob0")
    }
  )
  wg <- formula@prior.weights
  if (is.null(wg) | !length(wg)) wg <- rep(1, sizeObserved)

  POP <- switch (
    popVar,
    "analytic" = {
      N <- sum(wg / (1 - PW))
      if (is.null(derivFunc)) {
        if (formula@family@vfamily[1] %in% c("oapoisson", "pospoisson",
                                             "oipospoisson", "posnegbinomial",
                                             "posbinomial", "oiposbinomial")) {
          derivFunc <- switch (formula@family@vfamily[1],
            "oapoisson" = function(eta) {

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
              # lambda created by eval
              (exp(lambda) / (exp(lambda) - 1) ^ 2) * dlambda.deta
            },
            "oipospoisson" = function(eta) {

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
              c(-size*tmp1^size/((lambda+size)*(tmp1^size-1)^2),
                tmp1^size*((size + lambda)*log(tmp1)+lambda)/((lambda + size)*(tmp1^size-1)^2)) *
                c(dlambda.deta, dsize.deta)
            },
            "posbinomial" = function(eta) {

            },
            "oiposbinomial" = function(eta) {

            }
          )
        } else {
          stop("family slot not recognised, please provide a derivFunc argument with function for computing derivatives of 1/(1-prob0) with respect to linear predictors (or use bootstrap)")
        }
      }

      dd <- t(rep(wg, formula@family@infos()$M1) * derivFunc(formula@predictors)) %*% model.matrix(formula, type = "vlm")

      variation <- sum(wg * PW / (1 - PW)^2)
      variation <- variation + dd %*% vcov(formula) %*% t(dd)

      sc <- qnorm(p = 1 - siglevel / 2)
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
          control = controlForeign
        ),
        class = "popSizeEstResults"
      )},
    "bootstrap" = {
      N <- sum(wg / (1 - PW))

    }
  )

  structure(
    list(
      foreignObject  = formula,
      call           = match.call(),
      sizeObserved   = sizeObserved,
      populationSize = POP
    ),
    class = c("singleRforeign", "singleRStaticCountData", "singleR")
  )
}
