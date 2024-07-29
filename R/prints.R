#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
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
    (1 - x$populationSize$control$alpha) * 100,
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
#' @importFrom stats nobs
#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @export
print.singleRforeign <- function(x, ...) {
  cat("Call: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coef(x$foreignObject))
  cat("\nDegrees of Freedom:", stats::nobs(x$foreignObject) - 1,
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
