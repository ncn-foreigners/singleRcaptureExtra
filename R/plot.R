#' Title
#'
#' @param x A
#' @param plotType A
#' @param confIntStrata A
#' @param histKernels A
#' @param dfpop A
#' @param ... A
#'
#' @return A
#' @importFrom graphics abline barplot hist lines matplot legend boxplot panel.smooth axis text arrows par points
#' @importFrom stats density dlnorm dnorm
#' @export
plot.singleRforeign <- function(x,
                                plotType = c("bootHist",
                                             "marginal",
                                             "rootogram",
                                             "dfpopContr",
                                             "dfpopBox",
                                             "strata"),
                                confIntStrata = c("normal", "logNormal"),
                                histKernels = TRUE,
                                dfpop, ...) {
  if (missing(plotType)) stop("Argument plotType must be provided")
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  if ((plotType == "bootHist") && (!is.numeric(x$populationSize$boot))) {
    stop("Trying to plot bootstrap results with no bootstrap performed")
  }

  plotType <- match.arg(plotType)

  switch(plotType,
    marginal = {
      M <- marginalFreqVglm(x);
      FF <- M$table;
      FF[names(M$y)] <- M$y;
      FF[setdiff(names(M$table), names(M$y))] <- 0;
      graphics::matplot(
        y = cbind(M$table, FF),
        x = 0:max(x$y), type = "o",
        col = 1:2, pch = 21:22, lty = 2,
        main = "Plot of observed and fitted marginal frequencies",
        ylab = "Frequency",
        xlab = "Counts",
        ...);
      legend("topright",
             legend = c(x$model$family,
                        "Observed"),
             col = 1:2, pch = 21:22)
    },
    bootHist = {
      h <- graphics::hist(
        x$populationSize$boot,
        ylab = "Number of bootstrap samples",
        xlab = expression(hat(N)),
        main = "Bootstrap of population size estimates",
        ...
      );
      if (isTRUE(histKernels)) {
        xxx <- density(x$populationSize$boot, kernel = "epanechnikov");
        graphics::lines(
          x = xxx$x,
          y = dlnorm(
            x       = xxx$x,
            meanlog = log(mean(x$populationSize$boot) /
                            sqrt(1 + var(x$populationSize$boot) /
                                   mean(x$populationSize$boot) ^ 2)),
            sdlog   = sqrt(log(1 + var(x$populationSize$boot) / mean(x$populationSize$boot) ^ 2))
          ) * length(x$populationSize$boot) * diff(h$breaks)[1],
          lty = 2,
          col = 8
        )
        graphics::lines(
          x = xxx$x,
          y = dnorm(
            x    = xxx$x,
            mean = mean(x$populationSize$boot),
            sd   = sd(x$populationSize$boot)
          ) * length(x$populationSize$boot) * diff(h$breaks)[1],
          lty = 3,
          col = 9
        )
        graphics::lines(
          x = xxx$x,
          y = xxx$y * length(x$populationSize$boot) * diff(h$breaks)[1],
          lty = 4,
          col = 10
        )
        graphics::legend(
          "topright",
          c("log-normal density", "normal density", "Epanechnikov kernel"),
          lty = 2:4,
          col = 8:10
        )
      }
    },
    rootogram = {
      M <- marginalFreqVglm(x);
      FF <- M$table;
      FF[names(M$y)] <- M$y;
      FF[setdiff(names(M$table), names(M$y))] <- 0;
      FF <- FF[-1];
      bp <- graphics::barplot(
        sqrt(FF),
        offset = sqrt(M$table[-1]) - sqrt(FF),
        ylab = expression(sqrt("Frequency")), # This looks just ever so slightly fancier
        xlab = "captures",
        ylim = c(min(sqrt(M$table[-1]) - sqrt(FF)) - 1, max(sqrt(M$table[-1]) + 1)),
        ...);
      graphics::lines(bp, sqrt(M$table[-1]),
                      type = "o",
                      pch = 19,
                      lwd = 2,
                      col = 2);
      graphics::abline(h = 0,
                       lty = 2)
    },
    dfpopContr = {
      if (missing(dfpop)) dfpop <- dfpopsize(x, ...);
      # TODO:: use predict method here
      contr <- x$model$pointEst(
        pw = x$priorWeights,
        eta = x$linearPredictors,
        contr = TRUE
      );
      plot(x = dfpop, y = contr,
           main = paste0("Observation deletion effect on point estimate of",
                         "\npopulation size estimate vs observation contribution"),
           xlab = "Deletion effect", ylab = "Observation contribution",
           ...);
      abline(a = 0, b = 1, col = "red")
    },
    dfpopBox = {
      if (missing(dfpop)) dfpop <- dfpopsize(x, ...);
      graphics::boxplot(
        dfpop,
        ylab = "Deletion effect",
        main = paste0("Boxplot of observation deletion effect on",
                      "\npoint estimate of population size estimate"),
        ...
      )
    },
    strata = {
      if (missing(confIntStrata)) confIntStrata <- "logNormal"
      result <- stratifyPopsize(x, ...)
      est <- result[, 2]
      obs <- result[, 1]
      nm <- result[, 9]
      if (confIntStrata == "logNormal") cnf <- result[, 7:8]
      else cnf <- result[, 5:6]
      # OY ####
      # plot(x = 1:NROW(result), est,
      #      ylim = range(cnf),
      #      xlab = "", ylab="Sub population size estimate",
      #      main="Confidence intervals and point estimates for specified sub populations\nObserved population sizes are presented as navy coloured points",
      #      xaxt = "n", pch = 19
      # )
      # points(x = 1:NROW(result), obs, col = "navy", pch = 19)
      # axis(side = 1, at = 1:NROW(result), labels = FALSE)
      # text(x = 1:NROW(result), y=graphics::par("usr", no.readonly = TRUE)[3] - (range(cnf)[2] - range(cnf)[1]) / 20, adj = 1,
      #      nm, srt = 30, cex = .75,
      #      xpd = TRUE)
      # arrows(1:NROW(result), cnf[ ,1], 1:NROW(result), cnf[ ,2],
      #        length=0.05, angle=90, code=3)
      # OX ####

      tilt <- 0 # maybe add to parameters??
      plot(y = 1:NROW(result), x = est,
           xlim = range(cnf),
           xlab = "Sub population size estimate", ylab="",
           main = paste0(
             "Confidence intervals and point estimates for specified sub populations\n",
             "Observed population sizes are presented as navy coloured points"
           ),
           yaxt = "n", pch = 19
      )
      points(y = 1:NROW(result), x = obs, col = "navy", pch = 19)
      axis(side = 2, at = 1:NROW(result), labels = FALSE)
      text(
        y = 1:NROW(result),
        x = graphics::par("usr", no.readonly = TRUE)[3] - (range(cnf)[2] - range(cnf)[1]) / 20,
        adj = 1,
        nm,
        srt = tilt,
        cex = .6,
        xpd = TRUE
      )
      arrows(cnf[ ,1], 1:NROW(result),
             cnf[ ,2], 1:NROW(result),
             length = 0.05,
             angle  = 90,
             code   = 3)
    }
  )

  invisible()
}
