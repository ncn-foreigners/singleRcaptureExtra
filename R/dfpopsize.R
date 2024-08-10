#' @importFrom singleRcapture dfpopsize
#' @importFrom VGAM dgaitdpois dgaitdnbinom
#' @method dfpopsize singleRadditive
#' @rdname dfpop
#' @export
dfpopsize.singleRadditive <- function(model,
                                      data = NULL,
                                      cores = 1L,
                                      trace = FALSE,
                                      ...) {
  # add cores
  if (is.null(data)) {
    tryCatch(
      mf <- eval(model$foreignObject@call$data, envir = parent.frame()),
      error = function(e) {
        stop(paste0(
          "Call to data object in dfpopsize failed, ",
          "please provite the data argument in ",
          "dfpopsize.singleRadditive method."
        ))
      }
    )
  }

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
        .packages = c("VGAM", "singleRcapture", "VGAMdata"),
        .export = c("estimatePopsize.vgam", "vgam",
                    "dgaitdpois", "dgaitdnbinom")
      ),
      ex = {
        dd <- mf[-k, , drop = FALSE]
        cll$data <- as.symbol("dd")
        est <- eval(cll)
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

#' @importFrom singleRcapture dfpopsize
#' @method dfpopsize singleRforeign
#' @rdname dfpop
#' @export
dfpopsize.singleRforeign <- function(model,
                                     dfbeta = NULL,
                                     ...) {
  if (is.null(dfbeta)) dfbeta <- dfbeta(model, ...)

  X <- model.frame(model)
  y <- model.response(X)
  XXX <- model.matrix(model)
  N <- model$populationSize$pointEstimate
  wg <- model$foreignObject@prior.weights
  if (is.null(wg) | !length(wg)) wg <- rep(1, nobs(model))

  cf <- coef(model)
  res <- vector(mode = "numeric", length = nobs(model))

  for (kkk in 1:nrow(dfbeta)) {
    eta <- matrix(XXX %*% (cf - dfbeta[kkk, ]),
                  nrow = nobs(model), byrow = TRUE)[-kkk, , drop = FALSE]
    AA <- model$foreignObject
    AA@predictors <- eta

    PW <- tryCatch(
      expr = {fittedvlm(AA, type.fitted = "prob0")},
      error = function(e) {
        if (AA@family@vfamily[1] == "oapospoisson") {
          links <- (strsplit(AA@family@blurb[c(5, 7)], split = "\\(") |> unlist())[c(1,3)]
          lambda <- VGAM::eta2theta(
            AA@predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )

          return(exp(-lambda) / (1 - lambda * exp(-lambda)))
        } else if (AA@family@vfamily[1] == "oipospoisson") {
          links <- (strsplit(AA@family@blurb[c(3, 5)], split = "\\(") |> unlist())[c(1,3)]
          lambda <- VGAM::eta2theta(
            AA@predictors[, c(FALSE, TRUE), drop = FALSE], links[2],
            list(theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                 short = TRUE, tag = FALSE)
          )
          return(exp(-lambda))
        } else {
          stop("estimatePopsize.vglm method is only for objects with family slots with possibility of type.fitted = prob0 (and a few select ones)")
        }
      }
    )

    res[kkk] <- sum(wg[-kkk] / (1 - PW))
  }

  N - res
}
