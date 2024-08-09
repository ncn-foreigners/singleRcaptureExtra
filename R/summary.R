#' @title AAA
#' @description
#' A summary method for \code{singleRforeign} class
#' @details
#' Computes summary statistics for population size estimation and calls
#' summary on a foreign object. The results can be printed either as (default)
#' \code{print(summary(object))} where both the summary of population size
#' estimation and of the foreign object will be printed or alternatively the
#' later statistics can be omitted in printing with the call
#' \code{print(summary(object), summaryForeign = FALSE)}.
#'
#' @param object object of \code{singleRforeign} class
#' @param popSizeEst object of \code{popSizeEstResults} class such as the ones
#' created by [singleRcapture::redoPopEstimation()]. If not specified
#' population size estimation results will be drawn from the \code{object}.
#' results
#' @param ... passed to \code{summary} method on foreign object provided
#' on call to \code{estimatePopsize}.
#'
#' @importFrom stats sd
#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM summary
#' @method summary singleRforeign
#' @export
summary.singleRforeign <- function(object,
                                   popSizeEst,
                                   ...) {
  if (is.numeric(object$populationSize$boot)) {
    n <- length(object$populationSize$boot)
    m <- sum((object$populationSize$boot - mean(object$populationSize$boot)) ^ 3) / n
    s <- sd(object$populationSize$boot)
    skew <- m / (s ^ 3)
  } else {
    skew <- NULL
  }

  structure(
    list(
      summaryFreignObject = summary(object$foreignObject, ...),
      call = object$call,
      populationSize = if (missing(popSizeEst)) object$populationSize else popSizeEst,
      sizeObserved = object$sizeObserved,
      skew = skew,
      control = object$controlForeign
    ),
    class = "summarysingleRforeign"
  )
}
