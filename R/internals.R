# Internals, not documented
# TODO:: this is beyond moronic
getLinksBlurb <- function(x) {
  x <- x[which(x == "Links:    "):length(x)]
  x <- x[!(1:length(x) %% 2)]
  unlist(strsplit(x, split = "\\("))[1:(2*length(x)) %% 2]
}

#' @importFrom stats terms
internalGetXvlmMatrixFixed <- function(X, formulas, parNames, contrasts = NULL) {
  if (length(formulas[[1]]) == 3) {
    formulas[[1]][[2]] <- NULL
  }
  terms <- attr(X, "terms")

  if (attr(terms, "response") != 0) {
    #X <- X[, colnames(X)[-attr(terms, "response")], drop = FALSE]
  }
  nPar <- length(parNames)
  Xses <- list()

  for (k in 1:nPar) {
    # TODO:: Add contrasts here
    if (length(attr(terms(formulas[[k]], data = X), "term.labels")) != 0) {
      if (attr(terms(formulas[[k]], data = X), "intercept") == 0) {
        Xses[[k]] <- model.matrix(
          ~ . - 1,
          data = X[, intersect(attr(terms(formulas[[k]], data = X), "term.labels"),
                               colnames(X)), drop = FALSE]
        )
      } else {
        Xses[[k]] <- model.matrix(
          terms,
          data = X
        )
        print(Xses[[k]] |> head())
        stop("abc")
      }
    } else {
      Xses[[k]] <- model.matrix(
        ~ 1,
        X[, intersect(attr(terms(formulas[[k]], data = X), "term.labels"),
                      colnames(X)), drop = FALSE]
      )
      if (attr(terms(formulas[[k]], data = X), "intercept") == 0)
        warning(paste0(
          "One of formulas for the model has no variable ",
          "and no intercept and will be coerced to ~ 1."
        ))
    }
    if (k != 1) {
      colnames(Xses[[k]]) <- paste0(colnames(Xses[[k]]), ":", parNames[k])
    }
  }
  hwm <- sapply(Xses, ncol)

  # TODO:: Low priority but this could be much better
  # (without the need to allocate memmory for each X in Xses)
  # and for Xvlm which just stacks them
  Xvlm <- matrix(0, nrow = nPar * nrow(X), ncol = sum(hwm))
  colnames(Xvlm) <- unlist(sapply(X = Xses, FUN = colnames))
  row <- 0
  col <- 0
  for (k in Xses) {
    Xvlm[(row + 1):(row + nrow(k)), (col + 1):(col + ncol(k))] <- as.matrix(k)
    row <- row + nrow(k)
    col <- col + ncol(k)
  }
  attr(Xvlm, "hwm") <- hwm
  Xvlm
}
