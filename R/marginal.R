#' @rdname marginals
#' @importFrom VGAM dgaitdpois dgaitdnbinom
#' @importFrom VGAMdata doipospois doapospois
#' @export
marginalFreqVglm <- function(object,
                             includezeros = TRUE,
                             range, ...) {
  y <- if (is.null(object$foreignObject@y)) {
    as.numeric(stats::model.response(model.frame(object)))
  } else {
    as.numeric(object$foreignObject@y)
  }

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
    res <- c(object$populationSize$pointEstimate - length(eta), res)
    names(res)[1] <- "0"
  }

  structure(list(
    table = res, y = y,
    df = length(y) - length(coef(object$foreignObject)) - 1,
    name = "Foreign object"
  ), class = c("singleRmargin"))
}


#' @rdname marginals
#' @export
marginalFreqCountreg <- function(object,
                                 includezeros = TRUE,
                                 range, ...) {
  y <- if (is.null(object$foreignObject$y)){
    as.numeric(stats::model.response(model.frame(object)))
  } else {
    object$foreignObject$y
  }

  if (missing(range)) {range <- (min(y):max(y))}

  res <- predict(object$foreignObject, type = "prob")

  res <- res[, colnames(res) %in% as.character(range), drop = FALSE]
  res <- colSums(res)

  if(isTRUE(includezeros)) {
    res <- c(object$populationSize$pointEstimate - length(y), res)
    names(res)[1] <- "0"
  }

  y <- table(y)[names(table(y)) %in% as.character(range)]
  y <- y[!is.na(y)]

  structure(list(
    table = res, y = y,
    df = length(y) - length(coef(object$foreignObject)) - 1,
    name = "Foreign object"
  ), class = c("singleRmargin"))
}
