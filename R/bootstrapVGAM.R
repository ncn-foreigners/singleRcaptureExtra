# Add multicore
#' @importFrom VGAM vgam.fit
bootVGAM <- function(object,
                     B = 500,
                     trace = FALSE,
                     N, visT = FALSE,
                     bootType = c("parametric",
                                  "semiparametric",
                                  "nonparametric"),
                     data = NULL,
                     ...) {
  stop("Doesn't work yet")
  # data here
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

    try(
      theta <- vgam(
        # todo
      ),
      silent = FALSE
    )

    k <- k + 1

    if (is.null(theta)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      if (isTRUE(trace)) {print(summary(theta$predictors))}

      # change to values from bootstrap
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
