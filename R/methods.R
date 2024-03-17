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

#### TODO::
## - redoPopEstimation
## - stratifyPopsize
## - simulate (with truncared and usual)
## - dfpopsize
## - marginalFreq (make it a method in singleRcapture)


#' AA
#'
#' @param model a
#' @param cores a
#' @param trace asd
#' @param ... b
#' @importFrom singleRcapture dfpopsize
#' @method dfpopsize singleRadditive
#' @return a
#' @export
dfpopsize.singleRadditive <- function(model,
                                      cores = 1L,
                                      trace = FALSE,
                                      ...) {
  # add cores
  tryCatch(
    mf <- eval(model$foreignObject@call$data),
    error = function(e) {
      stop("Call to data object in dfpopsize failed")
    }
  )

  cll <- model$foreignObject@call

  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))

    res <- foreach::`%dopar%`(
      obj = foreach::foreach(
        k = 1:NROW(mf),
        .combine = c,
        ## TODO:: figure out something that requires less maintenance
        .packages = c("VGAM", "singleRcapture"),
        .export = c("estimatePopsize.vgam")
      ),
      ex = {
        dd <- mf[-k, , drop = FALSE]
        cll$data <- as.symbol("dd")
        est <- eval(cll, envir = environment())
        estimatePopsize(est)$populationSize$pointEstimate
      }
    )
  } else {
    res <- vector(mode = "numeric", length = as.integer(NROW(mf)))

    for (k in 1:NROW(mf)) {
      if (trace) cat("Current itteration: ", k, "\n", sep = "")
      dd <- mf[-k, , drop = FALSE]
      cll$data <- as.symbol("dd")
      est <- eval(cll, envir = environment())
      res[k] <- estimatePopsize(est)$populationSize$pointEstimate
    }
  }

  model$populationSize$pointEstimate - res
}

#' AAA
#'
#' @param model a
#' @param dfbeta a
#' @param ... a
#'
#' @importFrom singleRcapture dfpopsize
#' @method dfpopsize singleRforeign
#' @return a
#' @export
dfpopsize.singleRforeign <- function(model,
                                     dfbeta = NULL,
                                     ...) {
  if (is.null(dfbeta)) dfbeta <- dfbeta(model, ...)

  X <- model.frame(model)
  y <- model.response(X)

  for (variable in vector) {

  }
}

# TODO make them all methods

#' AA
#'
#' @param object A
#' @param includezeros A
#' @param range A
#' @param ... A
#'
#' @return A
#' @importFrom VGAM dgaitdpois dgaitdnbinom
#' @importFrom VGAMdata doipospois doapospois
#' @export
marginalFreqVglm <- function(object,
                             includezeros = TRUE,
                             range, ...) {
  y <- if (is.null(object$y)) as.numeric(stats::model.response(model.frame(object))) else object$y
  if (missing(range)) {range <- (min(y):max(y))}
  y <- table(y)[names(table(y)) %in% as.character(range)]
  y <- y[!is.na(y)]

  # PMF for truncated distributions:

  links <- getLinksBlurb(object$foreignObject@family@blurb)
  probFun <- switch(object$foreignObject@family@vfamily[1],
    "oapospoisson" = function(x, eta) {
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
      VGAMdata::doapospois(x, lambda = lambda, pobs1 = pobs1)
    },
    "pospoisson" = function(x, eta) {
      lambda <- VGAM::eta2theta(
        eta, links,
        list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
             short = TRUE, tag = FALSE)
      )
      VGAM::dgaitdpois(x = x, lambda.p = lambda, truncate = 0)
    },
    "oipospoisson" = function(x, eta) {
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
      VGAMdata::doipospois(x = x, lambda = lambda, pstr1 = pstr1)
    },
    "posnegbinomial" = function(x, eta) {
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
      VGAM::dgaitdnbinom(x = x, size.p = size, munb.p = lambda)
    }
  )

  eta <- predict(object$foreignObject)
  res <- apply(
    eta, 1,
    FUN = function(x) {probFun(
      x = range,
      eta = x
    )}
  )
  res <- rowSums(res)
  names(res) <- as.character(range)

  if(isTRUE(includezeros)) {
    res <- c(object$populationSize$pointEstimate - length(object$y), res)
    names(res)[1] <- "0"
  }

  structure(list(
    table = res, y = y,
    df = length(y) - length(coef(object$foreignObject)) - 1,
    name = "Foreign object"
  ), class = c("singleRmargin"))
}

