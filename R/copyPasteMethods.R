# All methods defined in this file are just basically calling original methods
# for foreign objects they contain. No need for doccumentation
#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importFrom stats predict
#' @method predict singleRforeign
#' @importFrom utils methods
#' @exportS3Method
predict.singleRforeign <- function(object, ...) predict(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM hatvalues
#' @importFrom stats vcov
#' @method vcov singleRforeign
#' @importFrom utils methods
#' @exportS3Method
vcov.singleRforeign <- function(object, ...) vcov(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM hatvalues
#' @importFrom stats hatvalues
#' @method hatvalues singleRforeign
#' @importFrom utils methods
#' @exportS3Method
hatvalues.singleRforeign <- function(model, ...) {
  stopifnot("No hatvalues method for foreign object or package needs to be loaded" =
              length(suppressWarnings(intersect(methods(class = class(model$foreignObject)), methods("hatvalues")))) > 0)

  hatvalues(model$foreignObject, ...)
}

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM dfbeta
#' @importFrom stats dfbeta
#' @method dfbeta singleRforeign
#' @importFrom utils methods
#' @exportS3Method
dfbeta.singleRforeign <- function(model, ...) {
  stopifnot("No dfbeta method for foreign object or package needs to be loaded" =
              length(suppressWarnings(intersect(methods(class = class(model$foreignObject)), methods("dfbeta")))) > 0)

  dfbeta(model$foreignObject, ...)
}

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM residuals
#' @importFrom stats residuals
#' @method residuals singleRforeign
#' @exportS3Method
residuals.singleRforeign <- function(object, ...) residuals(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importFrom stats cooks.distance
#' @method cooks.distance singleRforeign
#' @importFrom utils methods
#' @exportS3Method
cooks.distance.singleRforeign <- function(model, ...) {
  stopifnot("No cooks.distance method for foreign object or package needs to be loaded" =
              length(suppressWarnings(intersect(methods(class = class(model$foreignObject)), methods("cooks.distance")))) > 0)

  cooks.distance(model$foreignObject, ...)
}

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importFrom stats family
#' @method family singleRforeign
#' @importFrom utils methods
#' @exportS3Method
family.singleRforeign <- function(object, ...) {
  stopifnot("No family method for foreign object or package needs to be loaded" =
              length(suppressWarnings(intersect(methods(class = class(object$foreignObject)), methods("family")))) > 0)

  family(object$foreignObject, ...)
}

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM AIC
#' @importFrom stats AIC
#' @method AIC singleRforeign
#' @importFrom utils methods
#' @exportS3Method
AIC.singleRforeign <- function(object, ...) AIC(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM BIC
#' @importFrom stats BIC
#' @method BIC singleRforeign
#' @importFrom utils methods
#' @exportS3Method
BIC.singleRforeign <- function(object, ...) BIC(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM logLik
#' @importFrom stats logLik
#' @method logLik singleRforeign
#' @exportS3Method
logLik.singleRforeign <- function(object, ...) logLik(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM coef
#' @importFrom VGAM coef
#' @method coef singleRforeign
#' @exportS3Method
coef.singleRforeign <- function(object, ...) coef(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM Coef
#' @importFrom VGAM Coef
#' @method Coef singleRforeign
#' @exportS3Method
Coef.singleRforeign <- function(object, ...) Coef(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM model.matrix
#' @importFrom stats model.matrix
#' @importFrom utils methods
#' @method model.matrix singleRforeign
#' @exportS3Method
model.matrix.singleRforeign <- function(object, ...) {
  stopifnot("No model.matrix method for foreign object or package needs to be loaded" =
              length(suppressWarnings(intersect(methods(class = class(object$foreignObject)), methods("model.matrix")))) > 0)

  model.matrix(object$foreignObject, ...)
}

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM model.frame
#' @importFrom stats model.frame
#' @importFrom utils methods
#' @method model.frame singleRforeign
#' @exportS3Method
model.frame.singleRforeign <- function(formula, ...) {
  stopifnot("No model.frame method for foreign object or package needs to be loaded" =
              length(suppressWarnings(intersect(methods(class = class(formula$foreignObject)), methods("model.frame")))) > 0)

  model.frame(formula$foreignObject, ...)
}

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM fitted
#' @importFrom stats fitted
#' @importFrom utils methods
#' @method fitted singleRforeign
#' @exportS3Method
fitted.singleRforeign <- function(object, ...) {
  stopifnot("No fitted method for foreign object or package needs to be loaded" =
            length(suppressWarnings(intersect(methods(class = class(object$foreignObject)), methods("fitted")))) > 0)

  fitted(object$foreignObject, ...)
}

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM nobs
#' @importFrom stats nobs
#' @method nobs singleRforeign
#' @exportS3Method
nobs.singleRforeign <- function(object, ...) nobs(object$foreignObject, ...)

#' @importClassesFrom VGAM vglm
#' @importClassesFrom VGAM vgam
#' @importMethodsFrom VGAM df.residual
#' @importFrom stats df.residual
#' @method df.residual singleRforeign
#' @importFrom utils methods
#' @exportS3Method
df.residual.singleRforeign <- function(object, ...) {
  stopifnot("No df.residual method for foreign object or package needs to be loaded" =
              length(suppressWarnings(intersect(methods(class = class(object$foreignObject)), methods("df.residual")))) > 0)

  df.residual(object$foreignObject, ...)
}
