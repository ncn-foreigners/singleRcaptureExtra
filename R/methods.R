#' @importFrom stats sd
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

#' @export
print.summarysingleRforeign <- function(x,
                                        summaryForeign = TRUE,
                                        ...) {
  cat(
    "\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    "-----------------------",
    "\nPopulation size estimation results: ",
    "\nPoint estimate ", x$populationSize$pointEstimate,
    "\nObserved proportion: ", round(100 * x$sizeObserved / x$populationSize$pointEstimate,
                                     digits = 1),
    "% (N obs = ", x$sizeObserved, ")",
    if (!is.null(x$skew)) "\nBoostrap sample skewness: ",
    if (!is.null(x$skew)) x$skew,
    if (!is.null(x$skew)) "\n0 skewness is expected for normally distributed variable\n---",
    if (isTRUE(x$call$popVar == "bootstrap")) "\nBootstrap Std. Error " else "\nStd. Error ",
    sqrt(x$populationSize$variance), "\n",
    (1 - x$control$alpha) * 100,
    "% CI for the population size:\n",
    sep = ""
  )

  print(x$populationSize$confidenceInterval)

  cat((1 - x$populationSize$control$alpha) * 100,
      "% CI for the share of observed population:\n",
      sep = "")

  dd <- as.data.frame(x$populationSize$confidenceInterval)

  if (ncol(dd) == 1) {
    vctpop <- sort(x$populationSize$confidenceInterval, decreasing = TRUE)
    names(vctpop) <- rev(names(vctpop))
    print(100 * x$sizeObserved / vctpop)
  } else {
    print(data.frame(
      lowerBound = 100 * x$sizeObserved / dd[, 2],
      upperBound = 100 * x$sizeObserved / dd[, 1],
      row.names = rownames(dd)
    ))
  }

  if (isTRUE(summaryForeign)) {
    cat("\n-------------------------------", sep = "",
        "\n-- Summary of foreign object --",
        "\n-------------------------------\n")
    print(x$summaryFreignObject)
  }

  invisible(x)
}

#' @method print singleRforeign
#' @importFrom stats df.residual
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @export
print.singleRforeign <- function(x, ...) {
  cat("Call: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coef(x$foreignObject))
  cat("\nDegrees of Freedom:", x$foreignObject$df.null,
      "Total (i.e NULL);", df.residual(x$foreignObject),
      "Residual")

  cat("\nAIC: ", signif(AIC(x$foreignObject)),
      "\nBIC: ", signif(BIC(x$foreignObject)))
  cat("\n-----------------------------------\n",
      "Population size estimation results:\n",
      sep = "")
  print(x$populationSize)

  invisible(x)
}
