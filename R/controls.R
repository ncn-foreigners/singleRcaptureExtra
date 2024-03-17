#' @rdname controls
#' @export
controlEstPopCountreg <- function(trcount = 0,
                                  cores   = 1L,
                                  alpha   = .05,
                                  B       = 500,
                                  fittingMethod = c("IRLS", "optim"),
                                  bootstrapFitcontrol = singleRcapture::controlMethod(
                                    epsilon     = 1e-3,
                                    maxiter     = 20,
                                    optimMethod = "Nelder-Mead",
                                    silent      = TRUE
                                  ),
                                  traceBootstrapSize   = FALSE,
                                  bootstrapVisualTrace = FALSE,
                                  confType = c("percentilic",
                                               "normal",
                                               "basic"),
                                  bootType = c("parametric",
                                               "semiparametric",
                                               "nonparametric"),
                                  keepbootStat = TRUE) {
  if (missing(fittingMethod)) fittingMethod <- "IRLS"
  if (missing(bootType)) bootType <- "parametric"
  if (missing(confType)) confType <- "percentilic"

  if (!isTRUE(is.finite(alpha)) || isTRUE(0 > alpha || 1 < alpha) || isTRUE(length(alpha) > 1))
    stop("Argument alpha should be a numeric value between 0 and 1 (of length 1).")

  if (!isTRUE(is.infinite(trcount)) && !missing(trcount) || isTRUE(length(trcount) > 1))
    stop("Argument trcount should be a proper numeric value (of length 1).")

  if (!missing(bootType) && !isTRUE(bootType %in% c("parametric",
                                                    "semiparametric",
                                                    "nonparametric")))
    stop("Argument bootType should be a character with value either in c(parametric, semiparametric, nonparametric).")

  if (!isTRUE(is.finite(B)) || isTRUE(B < 0) || isTRUE(length(B) > 1))
    stop("Argument B should be a numeric value greater than 0 (of length 1).")

  if (!isTRUE(all(c(is.logical(traceBootstrapSize), is.logical(bootstrapVisualTrace)))) ||
      isTRUE(length(bootstrapVisualTrace) > 1) || isTRUE(length(traceBootstrapSize) > 1))
    stop("Arguments traceBootstrapSize and bootstrapVisualTrace should be logical values (of length 1).")

  list(
    bootstrapVisualTrace = bootstrapVisualTrace,
    bootstrapFitcontrol  = bootstrapFitcontrol,
    traceBootstrapSize   = traceBootstrapSize,
    fittingMethod        = fittingMethod,
    keepbootStat         = keepbootStat,
    bootType             = bootType,
    confType             = confType,
    trcount              = trcount,
    cores                = cores,
    alpha                = alpha,
    B                    = B
  )

}

#' @rdname controls
#' @export
controlEstPopVglm <- function(alpha = 0.05,
                              bootType = c("parametric",
                                           "semiparametric",
                                           "nonparametric"),
                              B = 500,
                              confType = c("percentilic",
                                           "normal",
                                           "basic"),
                              keepbootStat = TRUE,
                              traceBootstrapSize = FALSE,
                              bootstrapVisualTrace = FALSE,
                              bootstrapFitcontrol = VGAM::vglm.control(
                                epsilon = 1e-3,
                                maxit   = 10,
                                noWarning = TRUE
                              ),
                              sd = c("sqrtVar",
                                     "normalMVUE"),
                              cores = 1L) {

  if (missing(bootType)) bootType <- "parametric"
  if (missing(confType)) confType <- "percentilic"
  if (missing(sd))       sd <- "sqrtVar"

  # I'm using !isTRUE instead of isFALSE because I do not wan't nulls to pass the tests
  if (!isTRUE(is.finite(alpha)) || isTRUE(0 > alpha || 1 < alpha) || isTRUE(length(alpha) > 1))
    stop("Argument alpha should be a numeric value between 0 and 1 (of length 1).")

  if (!isTRUE(sd %in% c("sqrtVar", "normalMVUE")) && !missing(sd))
    stop("Argument sd should be a character with value either sqrtVar or normalMVUE.")

  if (!missing(bootType) && !isTRUE(bootType %in% c("parametric",
                                                    "semiparametric",
                                                    "nonparametric")))
    stop("Argument bootType should be a character with value either in c(parametric, semiparametric, nonparametric).")

  if (!isTRUE(is.finite(B)) || isTRUE(B < 0) || isTRUE(length(B) > 1))
    stop("Argument B should be a numeric value greater than 0 (of length 1).")

  if (!isTRUE(all(c(is.logical(traceBootstrapSize), is.logical(bootstrapVisualTrace)))) ||
      isTRUE(length(bootstrapVisualTrace) > 1) || isTRUE(length(traceBootstrapSize) > 1))
    stop("Arguments traceBootstrapSize and bootstrapVisualTrace should be logical values (of length 1).")


  list(
    bootstrapVisualTrace = bootstrapVisualTrace,
    bootstrapFitcontrol  = bootstrapFitcontrol,
    traceBootstrapSize   = traceBootstrapSize,
    keepbootStat         = keepbootStat,
    bootType             = bootType,
    confType             = confType,
    cores                = cores,
    alpha                = alpha,
    B                    = B,
    sd                   = sd
  )
}

#' @rdname controls
#' @export
controlEstPopVglm <- function(alpha = 0.05,
                              bootType = c("parametric",
                                           "semiparametric",
                                           "nonparametric"),
                              B = 500,
                              confType = c("percentilic",
                                           "normal",
                                           "basic"),
                              keepbootStat = TRUE,
                              traceBootstrapSize = FALSE,
                              bootstrapVisualTrace = FALSE,
                              bootstrapFitcontrol = VGAM::vgam.control(
                                epsilon = 1e-3,
                                maxit   = 10,
                                noWarning = TRUE
                              ),
                              sd = c("sqrtVar",
                                     "normalMVUE"),
                              cores = 1L,
                              data = NULL) {

  if (missing(bootType)) bootType <- "parametric"
  if (missing(confType)) confType <- "percentilic"
  if (missing(sd))       sd <- "sqrtVar"

  # I'm using !isTRUE instead of isFALSE because I do not wan't nulls to pass the tests
  if (!isTRUE(is.finite(alpha)) || isTRUE(0 > alpha || 1 < alpha) || isTRUE(length(alpha) > 1))
    stop("Argument alpha should be a numeric value between 0 and 1 (of length 1).")

  if (!isTRUE(sd %in% c("sqrtVar", "normalMVUE")) && !missing(sd))
    stop("Argument sd should be a character with value either sqrtVar or normalMVUE.")

  if (!missing(bootType) && !isTRUE(bootType %in% c("parametric",
                                                    "semiparametric",
                                                    "nonparametric")))
    stop("Argument bootType should be a character with value either in c(parametric, semiparametric, nonparametric).")

  if (!isTRUE(is.finite(B)) || isTRUE(B < 0) || isTRUE(length(B) > 1))
    stop("Argument B should be a numeric value greater than 0 (of length 1).")

  if (!isTRUE(all(c(is.logical(traceBootstrapSize), is.logical(bootstrapVisualTrace)))) ||
      isTRUE(length(bootstrapVisualTrace) > 1) || isTRUE(length(traceBootstrapSize) > 1))
    stop("Arguments traceBootstrapSize and bootstrapVisualTrace should be logical values (of length 1).")


  list(
    bootstrapVisualTrace = bootstrapVisualTrace,
    bootstrapFitcontrol  = bootstrapFitcontrol,
    traceBootstrapSize   = traceBootstrapSize,
    keepbootStat         = keepbootStat,
    bootType             = bootType,
    confType             = confType,
    cores                = cores,
    alpha                = alpha,
    B                    = B,
    sd                   = sd
  )
}
