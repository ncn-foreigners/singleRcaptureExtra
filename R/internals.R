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
        # print(Xses[[k]] |> head())
        # stop("abc")
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
#
internalStratPop <- function(object, stratas, alpha, cov, derivFunc, ...) {
  UseMethod("internalStratPop")
}
#
internalStratPop.vglm <- function(object, stratas, alpha, cov, derivFunc, ...) {
  priorWeights <- VGAM::weights(object)
  eta <- object@predictors
  Xvlm <- model.matrix(object, "vlm")

  # get covariance matrix
  if (is.function(cov)) cov <- cov(object, ...)
  if (is.null(cov)) cov <- vcov(object, ...)

  obs <- vector(mode = "numeric", length = length(stratas))
  est <- vector(mode = "numeric", length = length(stratas))
  stdErr <- vector(mode = "numeric", length = length(stratas))

  cnfNormal <- matrix(nrow = length(stratas), ncol = 2)
  cnfLogNormal <- matrix(nrow = length(stratas), ncol = 2)
  sc <- qnorm(p = 1 - alpha / 2)

  if (length(sc) != length(stratas)) {
    sc <- rep(sc, length.out = length(stratas))
  }

  PW <- tryCatch(
    expr = {fittedvlm(object, type.fitted = "prob0")},
    error = function(e) {
      if (object@family@vfamily[1] == "oapospoisson") {
        links <- (strsplit(object@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          object@predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )

        return(exp(-lambda) / (1 - lambda * exp(-lambda)))
      } else if (object@family@vfamily[1] == "oipospoisson") {
        links <- (strsplit(object@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          object@predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )
        return(exp(-lambda))
      } else {
        stop("estimatePopsize.vglm method is only for objects with family slots with possibility of type.fitted = prob0 (and a few select ones)")
      }
    }
  )

  dd <- derivFunc(eta)

  for (k in 1:length(stratas)) {
    cond <- stratas[[k]]
    obs[k] <- sum(cond)
    xx <- rep(cond, object@family@infos()$M1)

    if (obs[k] > 0) {
      est[k] <- sum(1 / (1 - PW[cond]))

      ddd <- t((rep(priorWeights, object@family@infos()$M1) * dd)[xx]) %*%
        Xvlm[rep(cond, each = object@family@infos()$M1), , drop = FALSE]

      variation <- priorWeights * PW / (1 - PW) ^ 2
      variation <- sum(variation[cond])
      variation <- variation + ddd %*% cov %*% t(ddd)

      G <- exp(sc[k] * sqrt(log(1 + variation / ((est[k] - obs[k]) ^ 2))))

      stdErr[k] <- sqrt(variation)
      cnfNormal[k, ] <- est[k] + c(-sc[k] * stdErr[k], sc[k] * stdErr[k])
      cnfLogNormal[k, ] <- obs[k] + c((est[k] - obs[k]) / G, (est[k] - obs[k]) * G)
    } else {
      est[k] <- 0
      stdErr[k] <- 0
      cnfNormal[k, ] <- c(0, 0)
      cnfLogNormal[k, ] <- c(0, 0)
    }
  }

  result <- data.frame(
    obs, est, 100 * obs / est, stdErr,
    cnfNormal[, 1], cnfNormal[, 2],
    cnfLogNormal[, 1], cnfLogNormal[, 2],
    names(stratas), alpha
  )

  bounds <- c("LowerBound", "UpperBound")

  colnames(result) <- c(
    "Observed", "Estimated",
    "ObservedPercentage", "StdError",
    paste0("normal", bounds),
    paste0("logNormal", bounds),
    "name", "confLevel"
  )

  result
}

internalStratPop.vgam <- function(object, stratas, alpha, cov, derivFunc, ...) {
  priorWeights <- VGAM::weights(object)
  eta <- object@predictors
  Xvlm <- model.matrix(object, "vlm")

  # get covariance matrix
  if (is.function(cov)) cov <- cov(object, ...)
  if (is.null(cov)) cov <- vcov(object, ...)

  obs <- vector(mode = "numeric", length = length(stratas))
  est <- vector(mode = "numeric", length = length(stratas))
  stdErr <- vector(mode = "numeric", length = length(stratas))

  cnfNormal <- matrix(nrow = length(stratas), ncol = 2)
  cnfLogNormal <- matrix(nrow = length(stratas), ncol = 2)
  sc <- qnorm(p = 1 - alpha / 2)

  if (length(sc) != length(stratas)) {
    sc <- rep(sc, length.out = length(stratas))
  }

  PW <- tryCatch(
    expr = {fittedvlm(object, type.fitted = "prob0")},
    error = function(e) {
      if (object@family@vfamily[1] == "oapospoisson") {
        links <- (strsplit(object@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          object@predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )

        return(exp(-lambda) / (1 - lambda * exp(-lambda)))
      } else if (object@family@vfamily[1] == "oipospoisson") {
        links <- (strsplit(object@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          object@predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )
        return(exp(-lambda))
      } else {
        stop("estimatePopsize.vglm method is only for objects with family slots with possibility of type.fitted = prob0 (and a few select ones)")
      }
    }
  )

  dd <- derivFunc(eta)

  for (k in 1:length(stratas)) {
    cond <- stratas[[k]]
    obs[k] <- sum(cond)
    xx <- rep(cond, object@family@infos()$M1)

    if (obs[k] > 0) {
      est[k] <- sum(1 / (1 - PW[cond]))

      ddd <- t((rep(priorWeights, object@family@infos()$M1) * dd)[xx]) %*%
        Xvlm[rep(cond, each = object@family@infos()$M1), , drop = FALSE]

      variation <- priorWeights * PW / (1 - PW) ^ 2
      variation <- sum(variation[cond])
      variation <- variation + ddd %*% cov %*% t(ddd)

      G <- exp(sc[k] * sqrt(log(1 + variation / ((est[k] - obs[k]) ^ 2))))

      stdErr[k] <- sqrt(variation)
      cnfNormal[k, ] <- est[k] + c(-sc[k] * stdErr[k], sc[k] * stdErr[k])
      cnfLogNormal[k, ] <- obs[k] + c((est[k] - obs[k]) / G, (est[k] - obs[k]) * G)
    } else {
      est[k] <- 0
      stdErr[k] <- 0
      cnfNormal[k, ] <- c(0, 0)
      cnfLogNormal[k, ] <- c(0, 0)
    }
  }

  result <- data.frame(
    obs, est, 100 * obs / est, stdErr,
    cnfNormal[, 1], cnfNormal[, 2],
    cnfLogNormal[, 1], cnfLogNormal[, 2],
    names(stratas), alpha
  )

  bounds <- c("LowerBound", "UpperBound")

  colnames(result) <- c(
    "Observed", "Estimated",
    "ObservedPercentage", "StdError",
    paste0("normal", bounds),
    paste0("logNormal", bounds),
    "name", "confLevel"
  )

  result
}

internalStratPop.zerotrunc <- function(object, stratas, alpha, cov, derivFunc, ...) {
  family <- switch (object$dist,
    "poisson"   = singleRcapture::ztpoisson(),
    "negbin"    = singleRcapture::ztnegbin(alphaLink = "neglog"),
    "geometric" = singleRcapture::ztgeom()
  )

  Xvlm <- model.matrix(object)
  eta <- Xvlm %*% coef(object)

  if (object$dist == "negbin") {
    eta <- cbind(eta, log(object$theta))
    Xvlm <- internalGetXvlmMatrixFixed(
      parNames = c("lambda", "theta"),
      formulas = list(object$formula, ~ 1),
      X        = model.frame(object)
    )
  }
  priorWeights <- object$weights
  if (is.null(priorWeights)) priorWeights <- rep(1, nrow(eta))

  # get covariance matrix
  if (is.function(cov)) cov <- cov(object, ...)
  if (is.null(cov)) {
    if (object$dist != "negbin") {
      cov <- vcov(object)
    } else {
      cov <- solve(-family$makeMinusLogLike(
        y = model.response(model.frame(object)),
        X = Xvlm,
        weight = priorWeights,
        deriv = 2,
        offset = c(object$offset, rep(0, nrow(eta)))
      )(c(stats::coef(object), log(object$theta))))
    }
  }

  obs <- vector(mode = "numeric", length = length(stratas))
  est <- vector(mode = "numeric", length = length(stratas))
  stdErr <- vector(mode = "numeric", length = length(stratas))

  cnfNormal <- matrix(nrow = length(stratas), ncol = 2)
  cnfLogNormal <- matrix(nrow = length(stratas), ncol = 2)
  sc <- qnorm(p = 1 - alpha / 2)

  if (length(sc) != length(stratas)) {
    sc <- rep(sc, length.out = length(stratas))
  }

  PW <- family$pointEst(pw = priorWeights, eta = eta, contr = TRUE)

  for (k in 1:length(stratas)) {
    cond <- stratas[[k]]
    obs[k] <- sum(cond)

    if (obs[k] > 0) {
      est[k] <- sum(PW[cond])

      variation <- as.numeric(family$popVar(
        eta = eta[cond, , drop = FALSE],
        pw  = priorWeights[cond],
        cov = cov,
        Xvlm = Xvlm[rep(cond, ncol(eta)), , drop = FALSE]
      ))

      G <- exp(sc[k] * sqrt(log(1 + variation / ((est[k] - obs[k]) ^ 2))))

      stdErr[k] <- sqrt(variation)
      cnfNormal[k, ] <- est[k] + c(-sc[k] * stdErr[k], sc[k] * stdErr[k])
      cnfLogNormal[k, ] <- obs[k] + c((est[k] - obs[k]) / G, (est[k] - obs[k]) * G)
    } else {
      est[k] <- 0
      stdErr[k] <- 0
      cnfNormal[k, ] <- c(0, 0)
      cnfLogNormal[k, ] <- c(0, 0)
    }
  }

  result <- data.frame(
    obs, est, 100 * obs / est, stdErr,
    cnfNormal[, 1], cnfNormal[, 2],
    cnfLogNormal[, 1], cnfLogNormal[, 2],
    names(stratas), alpha
  )

  bounds <- c("LowerBound", "UpperBound")

  colnames(result) <- c(
    "Observed", "Estimated",
    "ObservedPercentage", "StdError",
    paste0("normal", bounds),
    paste0("logNormal", bounds),
    "name", "confLevel"
  )

  result
}
