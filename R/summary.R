#' @title AAA
#' @description
#' AA
#' @details
#' aa
#'
#' @param object a
#' @param popSizeEst a
#' @param ... a
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
