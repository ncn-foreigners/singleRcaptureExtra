#' Controls for zerotrunc class
#'
#' @param trcount a
#' @param cores a
#' @param alpha a
#' @param B a
#' @param fittingMethod a
#' @param bootstrapFitcontrol a
#' @param traceBootstrapSize a
#' @param bootstrapVisualTrace a
#' @param confType a
#' @param bootType a
#' @param keepbootStat a
#'
#' @return a
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


#' @title Control parameters for \code{vglm} class
#'
#' @param alpha Numeric value indicating significance level. By default \code{.05}.
#' @param bootType Type of bootstrap by default \code{parametric}. For more detail
#' see [singleRcapture::estimatePopsize()].
#' @param B Number of bootstrap samples to be drawn. By default \code{500}.
#' @param confType Type of confidence interval to use for bootstrap variance
#' estimation.
#' @param keepbootStat Boolean value indicating whether to keep statistics from
#' bootstrap samples.
#' @param traceBootstrapSize Boolean value indicating whether to print sample
#' sizes (and estimate population sizes) on each iteration of bootstrap
#' algorithm. Only matters when \code{cores = 1}.
#' @param bootstrapVisualTrace Logical value indicating whether to make an
#' in real time plot of estimated statistics for \code{cores = 1} or whether
#' to make progress bar when \code{cores > 1}.
#' @param bootstrapFitcontrol A [VGAM::vglm.control()] type object controlling
#' behaviour of \code{vglm.fit} function used for fitting.
#' @param sd Type of standard error estimate.
#' @param cores Number of cores to use be default \code{1}.
#'
#' @seealso
#' [singleRcapture::estimatePopsize()] for details on bootstrap
#' [estimatePopsize.vglm()] for details on computation
#' [singleRcapture::controlPopVar()] for comparison with \code{singleRcapture}
#' [VGAM::vglm.control()] for fitting controls in bootstrap
#' [VGAM::vglm.fit()] for information on fitting algorithm
#' @return A list with desired control specifications.
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
