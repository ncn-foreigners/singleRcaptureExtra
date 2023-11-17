# Add multicore
#' @importFrom VGAM vgam.fit
bootVGAM <- function(object,
                     B = 500,
                     trace = FALSE,
                     N, visT = FALSE,
                     bootType = c("parametric",
                                  "semiparametric",
                                  "nonparametric"),
                     ...) {
  strappedStatistic <- vector("numeric", length = B)

  n   <- nobs(object)
  wt  <- object@prior.weights

  #X  <- model.matrix(object, "vlm")
  X  <- model.matrix(object, "lm")
  #lmMatAsgn <- attr(object, "assign")
  eta <- object@predictors
  NPRED <- NCOL(eta)
  y <- model.response(model.frame(object))
  constraints <- constraints(object, type = "term")

  if (length(wt) != NROW(eta)) {
    wt <- rep(1, NROW(eta))
  }

  extra <- object@extra
  extra$type.fitted <- if (object@family@vfamily[1] %in% c("oipospoisson", "oapospoisson")) extra$type.fitted else "prob0"

  links <- getLinksBlurb(object@family@blurb)

  offset <- object@offset
  if (any(dim(offset) != dim(object@predictors))) {
    offset <- matrix(0, nrow = NROW(object@predictors),
                     ncol = NCOL(object@predictors))
  }
  cf <- coef(object)

  control <- object@control
  control$noWarning <- TRUE

  # Mirroring VGAM::vgam fitting setup
  mtsave <- terms(as.formula(object@call$formula), specials = c("s", "sm.os", "sm.ps"),
                  data = model.frame(object))
  aa <- attributes(mtsave)
  smoothers <- aa$specials
  mgcv.sm.os <- length(smoothers$sm.os) > 0
  mgcv.sm.ps <- length(smoothers$sm.ps) > 0
  mgcv.sm.PS <- length(smoothers$sm.PS) > 0
  any.sm.os.terms <- mgcv.sm.os
  any.sm.ps.terms <- mgcv.sm.ps || mgcv.sm.PS
  mgcv.s <- length(smoothers$s) > 0
  nonparametric <- length(smoothers$s) > 0
  if (nonparametric) {
    ff <- apply(aa$factors[smoothers[["s"]], , drop = FALSE],
                2, any)
    smoothers[["s"]] <- if (any(ff))
      seq(along = ff)[aa$order == 1 & ff]
    else NULL
    smooth.labels <- aa$term.labels[unlist(smoothers)]
  } else {
    function.name <- "vglm"
  }
  are.sm.os.terms <- length(smoothers$sm.os) > 0
  are.sm.ps.terms <- (length(smoothers$sm.ps) + length(smoothers$sm.PS)) > 0
  if (are.sm.os.terms || are.sm.ps.terms) {
    if (length(smoothers$sm.os) > 0) {
      ff.sm.os <- apply(aa$factors[smoothers[["sm.os"]],
                                   , drop = FALSE], 2, any)
      smoothers[["sm.os"]] <- if (any(ff.sm.os))
        seq(along = ff.sm.os)[aa$order == 1 & ff.sm.os]
      else NULL
      smooth.labels <- aa$term.labels[unlist(smoothers)]
    }
    if (length(smoothers$sm.ps) > 0) {
      ff.sm.ps <- apply(aa$factors[smoothers[["sm.ps"]],
                                   , drop = FALSE], 2, any)
      smoothers[["sm.ps"]] <- if (any(ff.sm.ps))
        seq(along = ff.sm.ps)[aa$order == 1 & ff.sm.ps]
      else NULL
      smooth.labels <- aa$term.labels[unlist(smoothers)]
    }
  }

  if (visT) {
    plot(
      1, type = "n",
      xlab = "Bootstrap sample",
      ylab = "Value of population size estimator",
      main = expression(paste("Plot of values of ", hat(N), " obtained from bootstrap samples")),
      sub = "Points will be added in real time",
      xlim = c(0, B + 1), ylim = c(0, 2 * N)
    )
  }
  k <- 1

  while (k <= B) {
    strap <- switch (
      bootType,
      "nonparametric" = sample.int(replace = TRUE, n = n)
    )

    ystrap  <- as.numeric(y[strap])
    wtStrap <- as.numeric(wt[strap])

    offsetStrap  <- offset[strap, , drop = FALSE]
    mfStrap <- mf[strap, , drop = FALSE]
    xStrap <- model.matrix(object@terms, mfStrap)
    attr(xStrap, "assign") <- VGAM:::attrassigndefault(xStrap, attr(mfStrap, "terms"))

    if (are.sm.os.terms || are.sm.ps.terms) {
      assignx <- attr(xStrap, "assign")
      which.X.sm.osps <- assignx[smooth.labels]
      Data <- mf[, names(which.X.sm.osps), drop = FALSE]
      attr(Data, "class") <- NULL
      S.arg <- lapply(Data, attr, "S.arg")
      sparlist <- lapply(Data, attr, "spar")
      ridge.adj <- lapply(Data, attr, "ridge.adj")
      fixspar <- lapply(Data, attr, "fixspar")
      ps.int <- lapply(Data, attr, "ps.int")
      knots <- lapply(Data, attr, "knots")
      term.labels <- aa$term.labels
    }

    sm.osps.list <- if (any.sm.os.terms || any.sm.ps.terms) {
      list(indexterms = if (any.sm.os.terms) ff.sm.os else ff.sm.ps,
           intercept = aa$intercept, which.X.sm.osps = which.X.sm.osps,
           S.arg = S.arg, sparlist = sparlist, ridge.adj = ridge.adj,
           term.labels = term.labels, fixspar = fixspar, orig.fixspar = fixspar,
           ps.int = ps.int, knots = knots, assignx = assignx)
    }

    theta <- NULL
    try(
      theta <- VGAM::vgam.fit(
        y = y, x = xStrap, mf = mf, w = wtStrap, offset = offsetStrap,
        Xm2 = object@Xm2, Ym2 = object@Ym2,
        etastart = eta[as.numeric(strap), , drop = FALSE],
        constraints = object@constraints,
        extra = extra, qr.arg = TRUE, Terms = mtsave, family = object@family,
        control = object@control, mustart = NULL,
        nonparametric = nonparametric,
        smooth.labels = smooth.labels,
        function.name = function.name,
        sm.osps.list = sm.osps.list
      ),
      silent = FALSE
    )
    k <- k + 1

    if (is.null(theta)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      if (isTRUE(trace)) {print(summary(theta$predictors))}

      if (object@family@vfamily[1] == "oapospoisson") {
        links <- (strsplit(object@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          theta$predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )

        est <- sum(wt / (1 - exp(-lambda) / (1 - lambda * exp(-lambda))))
      } else if (object@family@vfamily[1] == "oipospoisson") {
        links <- (strsplit(object@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          theta$predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )
        est <- sum(wt / (1 - exp(-lambda)))
      } else {
        est <- sum(wt / (1 - theta$fitted.values))
      }

      if (visT) graphics::points(k - 1, est, pch = 1)
      print(k)
      if (isTRUE(trace)) cat(" Estimated population size: ", est,"\n",sep = "")
      #if (visT) graphics::points(k - 1, est, pch = 1)

      strappedStatistic[k - 1] <- est
    }
  }

  strappedStatistic
}
