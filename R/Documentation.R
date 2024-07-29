#' @import mathjaxr
NULL
#' \loadmathjax
#'
#' @title \code{estimatePopsize} methods for objects from other regression packages
#'
#' @description \code{estimatePopsize} methods currently implemented for the
#' \code{"zerotrunc", "vlgm", "vgam"} classes from \code{"countreg", "VGAM"}
#' packages. It is assumed that the response vector (i.e. the dependent variable)
#' corresponds to the number of times a given unit was observed in the data.
#'
#' Population size is then usually estimated by Horvitz-Thompson type estimator:
#'
#' \mjsdeqn{\hat{N} = \sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)} =
#' \sum_{k=1}^{N_{obs}}\frac{1}{1-\mathbb{P}(Y_{k}=0)}}
#'
#' where \mjseqn{I_{k}=I_{Y_{k} > 0}} are indicator variables, with value 1 if
#' kth unit was observed at least once and 0 otherwise. The full estimate is
#' obtained after plugging model based estimate for
#' \mjseqn{\dfrac{1}{\mathbb{P}(Y_{k}>0)}}.
#'
#' @param formula A fitted object containing some zero truncated regression
#' based on which population size estimation is to be done.
#' @param subset Not yet implemented.
#' @param naAction Not yet implemented.
#' @param popVar A character value indicating which variance estimation method
#' to use for population size estimation, either \code{"analytic"} or
#' \code{"bootstrap"}. The type of bootstrap can be specified in control argument
#' and the whole theory with bootstrap is available in [singleRcapture::estimatePopsize()].
#' @param control A list with control information for population size estimation
#' based on \code{formula} object and for bootstrap. Can (and should) be created
#' via functions written specifically to create these objects i.e. by
#' \code{controlEstPopCountreg} or \code{controlEstPopVglm}.
#' @param derivFunc A function of \code{formula} object that computes the derivative
#' of \mjseqn{(\hat{N}|I_{1},\dots I_{N})} to compute the:
#'  \mjsdeqn{\mathbb{E}\left(\text{var}
#' \left(\hat{N}|I_{1},\dots,I_{n}\right)\right)=
#' \left.\left(\frac{\partial(N|I_1,\dots,I_N)}{\partial\boldsymbol{\beta}}\right)^{T}
#' \text{cov}\left(\boldsymbol{\beta}\right)
#' \left(\frac{\partial(N|I_1,\dots,I_N)}{\partial\boldsymbol{\beta}}\right)
#' \right|_{\boldsymbol{\beta}=\hat{\boldsymbol{\beta}}}}
#' term full derivation can be found at [singleRcapture::estimatePopsize()]
#' Only used for non-standard \code{vglmff} class object i.e.
#' for anything other than \code{c("posbinomial", "oiposbinomial", "oapospoisson",
#' "pospoisson", "oipospoisson", "posnegbinomial")} this will be expanded in
#' future versions.
#' @param ... A
#' @details
#'
#' Most of the theory on single source capture-recapture (SSCR) is contained in
#' details of help file [singleRcapture::estimatePopsize()].
#'
#' As for the implementation itself starting with \code{zerotrunc} class i.e.
#' objects created by [countreg::zerotrunc()] there are only 3 possible
#' distributions those are zero truncated \code{c("poisson", "negbin", "geometric")}
#' all of which are also implemented in \code{singleRcapture}
#' (see [singleRcapture::ztpoisson()] for full list). Since for every model that
#' can be fitted in \code{countreg::zerotrunc} is also possible to make an equivalent
#' object in \code{singleRcapture} implementation of this method is most just
#' code via converting the \code{zerotrunc} class object information to that
#' used in \code{estimatePopsize.default} and then calling internals from \code{singleRcapture}.
#'
#' Additionally, instead of calling \code{countreg::zerotrunc} for fitting on
#' data from bootstrap samples (if \code{popVar} argument was set to \code{"bootstrap"})
#' \code{estimatePopsize.fit} is called to fit these objects (it is planned to
#' add additional control argument to choosing the fitting engine an option in controls).
#'
#'
#' @return An object with following \code{S3} classes:
#' \code{"singleRforeign", "singleRStaticCountData", "singleR"}
#' containing:
#' \itemize{
#'  \item{\code{foreignObject} - the \code{formula} object provided on call. This is necessary for proper functionality of methods.}
#'  \item{\code{call} - call made to create object.}
#'  \item{\code{sizeObserved} - number of observed units.}
#'  \item{\code{populationSize} - a \code{popSizeEstResults} class object containing population size estimates. Can be extracted via [singleRcapture::popSizeEst()]}
#'  \item{\code{pacakgeInfo} - character specifying package and function which created \code{formula} object.}
#' }
#' @seealso [singleRcapture::estimatePopsize()] - For original implementation and more theory.
#' [controlEstPopCountreg()] and [controlEstPopVglm()] - For controls provided to each \code{estimatePopsize} method.
#' [countreg::zerotrunc()] - for \code{zerotrunc} class.
#' @name foreignMethods
NULL


#' @title Control parameters for \code{vgam}, \code{vglm}, \code{countreg} class
#'
#' @param alpha Numeric value indicating significance level. By default \code{.05}.
#' @param bootType Type of bootstrap by default \code{parametric}. For more detail
#' see [singleRcapture::estimatePopsize()].
#' @param B Number of bootstrap samples to be drawn. By default \code{500}.
#' @param confType Type of confidence interval to use for bootstrap variance
#' estimation.
#' @param keepbootStat Boolean value indicating whether to keep statistics from
#' bootstrap samples.
#' @param traceBootstrapSize Boolean value indicating whether to print sample
#' sizes (and estimate population sizes) on each iteration of bootstrap
#' algorithm. Only matters when \code{cores = 1}.
#' @param bootstrapVisualTrace Logical value indicating whether to make an
#' in real time plot of estimated statistics for \code{cores = 1} or whether
#' to make progress bar when \code{cores > 1}.
#' @param bootstrapFitcontrol A [VGAM::vglm.control()] type object controlling
#' behaviour of \code{vglm.fit} function used for fitting.
#' @param sd Type of standard error estimate.
#' @param cores Number of cores to use be default \code{1}.
#' @param fittingMethod For \code{countreg} class which fitting algorithm from
#' \code{singleRcapture} to use (bootstrap for \code{countreg} class calls
#' bootstrap from \code{singleRcapture}).
#' @param trcount todo
#' @param data for \code{vgam} class only, a original \code{data.frame} object
#' on which model was fitted. Only required for bootstrap.
#'
#' @seealso
#' [singleRcapture::estimatePopsize()] for details on bootstrap and fitting algorithms
#' [estimatePopsize.vglm()] for details on computation
#' [singleRcapture::controlPopVar()] for comparison with \code{singleRcapture}
#' [VGAM::vglm.control()] for fitting controls in bootstrap
#' [VGAM::vglm.fit()] for information on fitting algorithm
#' @return A list with desired control specifications.
#' @name controls
NULL

#' @title Observed and fitted marginal Frequencies
#' @author Piotr Chlebicki
#'
#' @param object object of \code{singleR} class.
#' @param includezeros logical value indicating whether to include one counts in the zero-one truncated models.
#' @param range optional argument specifying range of selected Y values.
#' @param ... currently does nothing.
#'
#' @return A list with observed name of the fitted model family degrees of freedom and observed and fitted marginal frequencies.
#' @seealso [singleRcapture::marginalFreq()] -- for original \code{singleRcapture} implementation. [estimatePopsize.vgam()]
#' @name marginals
NULL

#' AA
#'
#' @param model a
#' @param cores a
#' @param trace asd
#' @param dfbeta a
#' @param ... b
#'
#' @return a
#' @name dfpop
NULL
