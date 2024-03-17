# Add multicore
#' @importFrom stats as.formula rpois terms
#' @importFrom VGAM vglm.fit lm2vlm.model.matrix predictvglm
#' @importFrom stats rbinom
bootVGLM <- function(object,
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

  # if (length(object@contrasts)) {
  #   lmMatAsgn  <- model.matrix(object@terms, model.frame(object), object@contrasts)
  # } else {
  #   lmMatAsgn  <- model.matrix(object@terms, model.frame(object))
  # }
  # lmMatAsgn <- attr(lmMatAsgn, "assign")
  # print(lmMatAsgn)
  # print(attr(model.matrix(object, "lm"), "assign"))
  # stop("abc")

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

  # getting links
  links <- getLinksBlurb(object@family@blurb)

  offset <- object@offset
  if (any(dim(offset) != dim(object@predictors))) {
    offset <- matrix(0, nrow = NROW(object@predictors),
                        ncol = NCOL(object@predictors))
  }
  cf <- coef(object)

  control <- object@control
  control$noWarning <- TRUE

  if (visT) {
    plot(
      1, type = "n",
      xlab = "Bootstrap sample",
      ylab = "Value of population size estimator",
      main = expression(paste("Plot of values of ", hat(N), " obtained from bootstrap samples")),
      sub = "Points will be added in real time",
      xlim = c(0, B + 1), ylim = c(0, 2 * sum(N))
    )
  }
  k <- 1

  contr <- N
  N <- sum(N)
  while (k <= B) {
    switch (
      bootType,
      "nonparametric" = {
        strap <- sample.int(replace = TRUE, n = n)
        yStrap  <- as.numeric(y[strap])
        wtStrap <- as.numeric(wt[strap])
        etaStrap <- eta[as.numeric(strap), , drop = FALSE]

        offsetStrap  <- offset[strap, , drop = FALSE]
        xStrap <- X[strap, , drop = FALSE]

        # Figure out if non treatment constrasts break this
        # if (length(object@contrasts)) {
        #   xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
        # } else {
        #   xStrap <- model.matrix(object@terms, mfStrap)
        # }

        attr(xStrap, "contrasts") <- object@contrasts
        attr(xStrap, "assign") <- attr(X, "assign")
        attr(xStrap, "orig.assign.lm") <- attr(X, "orig.assign.lm")
      },
      "parametric" = {
        nn <- floor(N) + stats::rbinom(n = 1, size = 1, prob = N - floor(N))
        strap <- sample.int(replace = TRUE, n = length(y), size = nn, prob = contr)
        wtStrap <- as.numeric(wt[strap])
        etaStrap <- eta[as.numeric(strap), , drop = FALSE]

        offsetStrap  <- offset[strap, , drop = FALSE]
        xStrap <- X[strap, , drop = FALSE]

        # if (length(object@contrasts)) {
        #   xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
        # } else {
        #   xStrap <- model.matrix(object@terms, mfStrap)
        # }

        yStrap <- object@extra$singleRcaptureSimulate(
          nn, eta = VGAM::predictvglm(
            type = "link", object,
            newdata = as.data.frame(model.frame(object)[strap, , drop = FALSE])
          ), links = links
        )

        wtStrap <- wtStrap[yStrap > 0]
        xStrap  <- xStrap[yStrap > 0, , drop = FALSE]

        etaStrap     <- etaStrap[yStrap > 0, , drop = FALSE]
        offsetStrap  <- offsetStrap[yStrap > 0, , drop = FALSE]

        yStrap <- yStrap[yStrap > 0]

        attr(xStrap, "contrasts") <- object@contrasts
        attr(xStrap, "assign") <- attr(X, "assign")
        attr(xStrap, "orig.assign.lm") <- attr(X, "orig.assign.lm")
      },
      "semiparametric" = {
        strap <- sum(rbinom(size = 1, n = N, prob = n / N))

        strap <- sample.int(replace = TRUE, n = n, size = strap)
        yStrap  <- as.numeric(y[strap])
        wtStrap <- as.numeric(wt[strap])

        offsetStrap  <- offset[strap, , drop = FALSE]
        xStrap <- X[strap, , drop = FALSE]

        # if (length(object@contrasts)) {
        #   xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
        # } else {
        #   xStrap <- model.matrix(object@terms, mfStrap)
        # }
        etaStrap <- eta[as.numeric(strap), , drop = FALSE]

        attr(xStrap, "contrasts") <- object@contrasts
        attr(xStrap, "assign") <- attr(X, "assign")
        attr(xStrap, "orig.assign.lm") <- attr(X, "orig.assign.lm")
      }
    )

    if (isTRUE(trace)) cat("Iteration number: ", k,
                           " sample size: ", length(yStrap),
                           sep = "")
    theta <- NULL
    suppressWarnings(
      try(
        theta <- vglm.fit(
          y = yStrap, x = xStrap,
          etastart = etaStrap,
          constraints = constraints,
          extra = extra, w = wtStrap, offset = offsetStrap,
          qr.arg = TRUE, Terms = object@terms,
          function.name = "vglm", family = object@family,
          control = control,
          Xm2 = object@Xm2, Ym2 = object@Ym2, mustart = NULL
        ),
        silent = TRUE
      )
    )
    k <- k + 1

    if (is.null(theta)) {
      k <- k - 1
      if (isTRUE(trace)) cat("\n",sep = "")
    } else {
      theta <- theta$predictors

      if (object@family@vfamily[1] == "oapospoisson") {
        links <- (strsplit(object@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          theta[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )

        est <- sum(wtStrap / (1 - exp(-lambda) / (1 - lambda * exp(-lambda))))
      } else if (object@family@vfamily[1] == "oipospoisson") {
        links <- (strsplit(object@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          theta[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )
        est <- sum(wtStrap / (1 - exp(-lambda)))
      } else {
        est <- sum(wtStrap / (1 - object@family@linkinv(
          eta = theta,
          extra = list("type.fitted" = "prob0")
        )))
      }

      if (visT) graphics::points(k - 1, est, pch = 1)
      if (isTRUE(trace)) cat(" Estimated population size: ", est,"\n",sep = "")
      #if (visT) graphics::points(k - 1, est, pch = 1)

      strappedStatistic[k - 1] <- est
    }
  }

  strappedStatistic
}

#' @importFrom stats as.formula rpois terms
#' @importFrom VGAM vglm.fit lm2vlm.model.matrix predictvglm
#' @importFrom stats rbinom
multiCoreBootVGLM <- function(object,
                              B = 500,
                              trace = FALSE,
                              N, visT = FALSE,
                              bootType = c("parametric",
                                           "semiparametric",
                                           "nonparametric"),
                              cores,
                              ...) {
  cl <- parallel::makeCluster(cores)

  n   <- nobs(object)
  wt  <- object@prior.weights

  X  <- model.matrix(object, "lm")
  eta <- object@predictors
  NPRED <- NCOL(eta)
  y <- model.response(model.frame(object))
  constraints <- constraints(object, type = "term")

  if (length(wt) != NROW(eta)) {
    wt <- rep(1, NROW(eta))
  }

  extra <- object@extra
  extra$type.fitted <- if (object@family@vfamily[1] %in% c("oipospoisson", "oapospoisson")) extra$type.fitted else "prob0"

  # getting links
  links <- getLinksBlurb(object@family@blurb)

  offset <- object@offset
  if (any(dim(offset) != dim(object@predictors))) {
    offset <- matrix(0, nrow = NROW(object@predictors),
                     ncol = NCOL(object@predictors))
  }
  cf <- coef(object)

  control <- object@control
  control$noWarning <- TRUE

  contr <- N
  N <- sum(N)

  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  strappedStatistic <- foreach::`%dopar%`(
    obj = foreach::foreach(
      k = 1:B,
      .combine = c,
      .packages = c("VGAM", "singleRcapture")
    ),
    ex  = {
      theta <- NULL
      while (is.null(theta)) {
        switch (bootType,
          "nonparametric" = {
            strap <- sample.int(replace = TRUE, n = n)
            yStrap  <- as.numeric(y[strap])
            wtStrap <- as.numeric(wt[strap])
            etaStrap <- eta[as.numeric(strap), , drop = FALSE]

            offsetStrap  <- offset[strap, , drop = FALSE]
            xStrap <- X[strap, , drop = FALSE]

            # Figure out if non treatment constrasts break this
            # if (length(object@contrasts)) {
            #   xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
            # } else {
            #   xStrap <- model.matrix(object@terms, mfStrap)
            # }

            attr(xStrap, "contrasts") <- object@contrasts
            attr(xStrap, "assign") <- attr(X, "assign")
            attr(xStrap, "orig.assign.lm") <- attr(X, "orig.assign.lm")
          },
          "parametric" = {
            nn <- floor(N) + stats::rbinom(n = 1, size = 1, prob = N - floor(N))
            strap <- sample.int(replace = TRUE, n = length(y), size = nn, prob = contr)
            wtStrap <- as.numeric(wt[strap])
            etaStrap <- eta[as.numeric(strap), , drop = FALSE]

            offsetStrap  <- offset[strap, , drop = FALSE]
            xStrap <- X[strap, , drop = FALSE]

            # if (length(object@contrasts)) {
            #   xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
            # } else {
            #   xStrap <- model.matrix(object@terms, mfStrap)
            # }
            #dd <- VGAM::model.framevlm(object)[strap, , drop = FALSE]

            yStrap <- object@extra$singleRcaptureSimulate(
              nn, eta = etaStrap, links = links
            )

            wtStrap <- wtStrap[yStrap > 0]
            xStrap  <- xStrap[yStrap > 0, , drop = FALSE]

            etaStrap     <- etaStrap[yStrap > 0, , drop = FALSE]
            offsetStrap  <- offsetStrap[yStrap > 0, , drop = FALSE]

            yStrap <- yStrap[yStrap > 0]

            attr(xStrap, "contrasts") <- object@contrasts
            attr(xStrap, "assign") <- attr(X, "assign")
            attr(xStrap, "orig.assign.lm") <- attr(X, "orig.assign.lm")
          },
          "semiparametric" = {
            strap <- sum(rbinom(size = 1, n = N, prob = n / N))

            strap <- sample.int(replace = TRUE, n = n, size = strap)
            yStrap  <- as.numeric(y[strap])
            wtStrap <- as.numeric(wt[strap])

            offsetStrap  <- offset[strap, , drop = FALSE]
            xStrap <- X[strap, , drop = FALSE]

            # if (length(object@contrasts)) {
            #   xStrap <- model.matrix(object@terms, mfStrap, object@contrasts)
            # } else {
            #   xStrap <- model.matrix(object@terms, mfStrap)
            # }
            etaStrap <- eta[as.numeric(strap), , drop = FALSE]

            attr(xStrap, "contrasts") <- object@contrasts
            attr(xStrap, "assign") <- attr(X, "assign")
            attr(xStrap, "orig.assign.lm") <- attr(X, "orig.assign.lm")
          }
        )


        suppressWarnings(
          try(
            theta <- vglm.fit(
              y = yStrap, x = xStrap,
              etastart = etaStrap,
              constraints = constraints,
              extra = extra, w = wtStrap, offset = offsetStrap,
              qr.arg = TRUE, Terms = object@terms,
              function.name = "vglm", family = object@family,
              control = control,
              Xm2 = object@Xm2, Ym2 = object@Ym2, mustart = NULL
            ),
            silent = TRUE
          )
        )
      }


      theta <- theta$predictors

      if (object@family@vfamily[1] == "oapospoisson") {
        links <- (strsplit(object@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          theta[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )

        est <- sum(wtStrap / (1 - exp(-lambda) / (1 - lambda * exp(-lambda))))
      } else if (object@family@vfamily[1] == "oipospoisson") {
        links <- (strsplit(object@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
        lambda <- VGAM::eta2theta(
          theta[, c(FALSE, TRUE), drop = FALSE], links[2],
          list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
               short = TRUE, tag = FALSE)
        )
        est <- sum(wtStrap / (1 - exp(-lambda)))
      } else {
        est <- sum(wtStrap / (1 - object@family@linkinv(
          eta = theta,
          extra = list("type.fitted" = "prob0")
        )))
      }

      est
    }
  )

  strappedStatistic
}
