#detach("package:VGAM", unload = TRUE)
library(VGAM)
library(VGAMdata)

# TODO use generated data

obj <- vgam(TOTAL_SUB ~ s(log_size) + s(log_distance) + C_TYPE,
            family = pospoisson(), data = farmsubmission)

obj1 <- vglm(TOTAL_SUB ~ .,
             family = pospoisson(), data = farmsubmission)

obj2 <- vglm(TOTAL_SUB ~ .,
             family = oipospoisson(), data = farmsubmission)

obj3 <- vglm(TOTAL_SUB ~ .,
             family = oapospoisson(), data = farmsubmission)

expect_silent(AA <- estimatePopsize(obj))
expect_silent(AA1 <- estimatePopsize(obj1))

expect_silent(
  AA1.1 <- estimatePopsize(obj1, popVar = "bootstrap",
                           control = controlEstPopVglm(
                             B = 50,
                             bootstrapFitcontrol = vglm.control(epsilon = .01)
                           ))
)

expect_silent(
  AA1.2 <- estimatePopsize(
    obj1, popVar = "bootstrap",
    control = controlEstPopVglm(bootType = "semiparametric",
                                B = 50,
                                bootstrapFitcontrol = vglm.control(epsilon = .01))
  )
)


expect_silent(
  AA1.3 <- estimatePopsize(
    obj1, popVar = "bootstrap",
    control = controlEstPopVglm(bootType = "parametric",
                                B = 50,
                                bootstrapFitcontrol = vglm.control(epsilon = .01))
  )
)

expect_silent(AA2 <- estimatePopsize(obj2))

expect_silent(
  AA2.1 <- estimatePopsize(obj2, popVar = "bootstrap",
                           control = controlEstPopVglm(
                             B = 50,
                             bootstrapFitcontrol = vglm.control(epsilon = .01)
                           ))
)

expect_silent(
  AA2.2 <- estimatePopsize(
    obj2, popVar = "bootstrap",
    control = controlEstPopVglm(bootType = "semiparametric",
                                B = 50,
                                bootstrapFitcontrol = vglm.control(epsilon = .01))
  )
)


expect_silent(
  AA2.3 <- estimatePopsize(
    obj2, popVar = "bootstrap",
    control = controlEstPopVglm(bootType = "parametric",
                                B = 50,
                                bootstrapFitcontrol = vglm.control(epsilon = .01))
  )
)

expect_silent(AA3 <- estimatePopsize(obj3))

expect_silent(
  AA3.1 <- estimatePopsize(obj3, popVar = "bootstrap",
                           control = controlEstPopVglm(
                             B = 50,
                             bootstrapFitcontrol = vglm.control(epsilon = .01)
                           ))
)

expect_silent(
  AA3.2 <- estimatePopsize(
    obj3, popVar = "bootstrap",
    control = controlEstPopVglm(bootType = "semiparametric",
                                B = 50,
                                bootstrapFitcontrol = vglm.control(epsilon = .01))
  )
)


expect_silent(
  AA3.3 <- estimatePopsize(
    obj3, popVar = "bootstrap",
    control = controlEstPopVglm(bootType = "parametric",
                                B = 50,
                                bootstrapFitcontrol = vglm.control(epsilon = .01))
  )
)

expect_silent(fitted(AA))
expect_silent(fitted(AA1))

expect_silent(fitted(obj))
expect_silent(fitted(obj1))

expect_silent(df.residual(AA))
expect_silent(df.residual(AA1))

expect_silent(df.residual(obj))
expect_silent(df.residual(obj1))

expect_silent(model.frame(AA))
expect_silent(model.frame(AA1))

expect_silent(model.frame(obj))
expect_silent(model.frame(obj1))

expect_silent(model.matrix(AA))
expect_silent(model.matrix(AA1))

expect_silent(model.matrix(obj))
expect_silent(model.matrix(obj1))

expect_error(dfbeta(AA))

expect_silent(plot(AA3.3, plotType = "bootHist", ylim =c(0, 50)))

expect_error(dfbeta(obj))

expect_silent(model.matrix(obj))

if (isTRUE(tolower(Sys.getenv("TEST_SINGLERCAPTURE_MULTICORE_DEVELOPER")) == "true")) {
  obj4 <- vglm(TOTAL_SUB ~ .,
               family = posnegbinomial(), data = farmsubmission)

  expect_silent(AA4 <- estimatePopsize(obj4))
}
