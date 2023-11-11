#' Title
#'
#' @param alpha a
#' @param bootType a
#' @param B a
#' @param confType a
#' @param keepbootStat a
#' @param traceBootstrapSize a
#' @param bootstrapVisualTrace a
#' @param bootstrapFitcontrol a
#' @param sd a
#' @param cores a
#'
#' @return a
#' @export
controlPopVarVglm <- function(alpha = 0.05,
                              bootType = c("parametric", "semiparametric", "nonparametric"),
                              B = 500,
                              confType = c("percentilic", "normal", "basic"),
                              keepbootStat = TRUE,
                              traceBootstrapSize = FALSE,
                              bootstrapVisualTrace = FALSE,
                              bootstrapFitcontrol = NULL,
                              sd = c("sqrtVar", "normalMVUE"),
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
