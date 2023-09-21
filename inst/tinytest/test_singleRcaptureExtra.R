#detach("package:VGAM", unload = TRUE)
library(VGAM)

obj <- vgam(TOTAL_SUB ~ s(log_size) + s(log_distance) + C_TYPE,
            family = pospoisson(), data = farmsubmission)

obj1 <- vglm(TOTAL_SUB ~ .,
             family = pospoisson(), data = farmsubmission)

AA <- estimatePopsize(obj)
AA1 <- estimatePopsize(obj1)

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

expect_error(dfbeta(obj))
