#' @title AAA
#' @description
#' AA
#' @details
#' aa
#'
#' @param object a
#' @param stratas a
#' @param alpha description
#' @param cov a
#' @param ... a
#'
#' @return A
#' @seealso [singleRcapture::stratifyPopsize()] [estimatePopsize()]
#' @importFrom singleRcapture stratifyPopsize
#' @method stratifyPopsize singleRforeign
#' @export
stratifyPopsize.singleRforeign <- function(object,
                                           stratas,
                                           alpha,
                                           cov = NULL,
                                           ...) {
  # TODO:: a method for
  # If stratas is unspecified get all levels of factors in modelFrame
  terms <- terms(object)
  if (missing(stratas)) {
    stratas <- names(which(attr(terms, "dataClasses") == "factor"))
    stratas <- stratas[stratas %in% attr(terms, "term.labels")]
    if (!length(stratas)) {
      stratas <- names(which(attr(terms, "dataClasses") == "character"))
      stratas <- stratas[stratas %in% attr(terms, "term.labels")]
    }
    if (!length(stratas)) {
      stop("No stratas argument was provided and no factors or character columns are present in model.frame.")
    }
  }
  # If there are no factors or characters and no stratas was provided throw error
  # if significance level is unspecified set it to 5%
  if (missing(alpha)) alpha <- .05

  # convert stratas to list for all viable types of specifying the argument
  if (inherits(stratas, "formula")) {
    mf <- model.frame(stratas, model.frame(object))
    mmf <- model.matrix(
      stratas, data = mf,
      contrasts.arg = lapply(
        subset(mf, select = sapply(mf, is.factor)), # this makes it so that all levels of factors are encoded
        contrasts, contrasts = FALSE
      )
    )
    trm <- attr(mf, "terms")
    stratas <- list()
    for (k in attr(trm, "term.labels")) {
      if (k %in% colnames(mf)) {
        if (is.integer(mf[, k]) | is.character(mf[, k])) {
          for (t in unique(mf[,k])) {
            stratas[[paste0(k, "==", t)]] <- mf[,k] == t
          }
        } else if (is.factor(mf[, k])) {
          for (t in levels(mf[, k])) {
            stratas[[paste0(k, "==", t)]] <- mf[,k] == t
          }
        }
      } else {
        tLevs <- colnames(mmf)[attr(mmf, "assign") == which(attr(trm, "term.labels") == k)]
        for (t in tLevs) {
          stratas[[as.character(t)]] <- mmf[, t] == 1
        }
      }
    }
  } else if (is.list(stratas)) {
    if (!all(sapply(stratas, is.logical)))
      stop("Invalid way of specifying subpopulations in stratas. If stratas argument is a list it is expected that this list contains logical vectors only.")

    if (length(stratas[[1]]) != object$sizeObserved)
      stop("Elements of stratas object should have length equal to number of observed units.")

  } else if (is.logical(stratas)) {
    if (length(stratas) != object$sizeObserved)
      stop("Stratas object should have length equal to number of observed units.")

    stratas <- list(strata = stratas)
  } else if (is.character(stratas)) {
    modelFrame <- model.frame(object)
    out <- list()
    for (k in stratas) {
      if (!(k %in% colnames(modelFrame)))
        stop("Variable specified in stratas is not present in model frame.")

      #if (!(is.factor(modelFrame[, k])) & !(is.character(modelFrame[, k])))
      #  stop("Variable specified in stratas is not a factor or a character vector.")

      if (is.factor(modelFrame[, k])) {
        # this makes a difference on factor that is not present
        for (t in levels(modelFrame[, k])) {
          out[[paste0(as.character(k), "==", t)]] <- (modelFrame[, k] == t)
        }
      } else {
        for (t in unique(modelFrame[, k])) {
          out[[paste0(as.character(k), "==", t)]] <- (modelFrame[, k] == t)
        }
      }
    }
    stratas <- out
  } else {
    # a formula, a list with logical vectors specifying different sub
    # populations or a single logical vector or a vector
    # with names of factor variables.
    errorMessage <- paste0(
      "Invalid way of specifying subpopulations in stratas.\n",
      "Please provide either:\n",
      "(1) - a list with logical vectors specifying different sub populations\n",
      "(2) - a single logical vector\n",
      "(3) - a formula\n",
      "(4) - a vector with names of variables by which stratas will be created\n"
    )
    stop(errorMessage)
  }

  internalStratPop(object$foreignObject, stratas, alpha, cov, object$derivFunc, ...)
}
