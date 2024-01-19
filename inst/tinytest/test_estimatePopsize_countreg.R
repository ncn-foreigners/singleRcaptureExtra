library(countreg)

# generate data
# set.seed(123)
# N <- 2400
# x1 <- rbinom(n = N, prob = .5, size = 1)
# x2 <- runif(n = N)
#
# y1 <-   rpois(n = N, lambda = exp(-1 + x1 * x2 + .3 * x2 - .1 * x1))
# y2 <-   rgeom(n = N, prob = 1 / (exp(-1 + x1 * x2 + .3 * x2 - .1 * x1) + 1))
# y3 <- rnbinom(n = N, mu = exp(-1 + x1 * x2 + .3 * x2 - .1 * x1),
#               size = .3) # log(1.62482361) ~~ 1.2
# df <- data.frame(y1 = y1, y2 = y2, y3 = y3, x1 = x1, x2 = x2)
# write.csv(x = df, file = "inst/tinytest/countreg_test_dataframe.csv",
#           row.names = FALSE)

#df <- read.csv("inst/tinytest/countreg_test_dataframe.csv", header = TRUE)
df <- read.csv("countreg_test_dataframe.csv", header = TRUE)

mm1 <- zerotrunc(
  formula = y1 ~ x1 * x2, dist = "poisson",
  data    = df[df$y1>0,]
)

mm2 <- zerotrunc(
  formula = y2 ~ x1 * x2, dist = "geometric",
  data    = df[df$y2>0,]
)

mm3 <- zerotrunc(
  formula = y3 ~ x1 * x2, dist = "negbin",
  data    = df[df$y3>0,]
)

# test for parameter estimation (not actually needed)

# expect_true(
#   all(confint(mm1)[,1] < c(-1, -.1, .3, 1)) &
#   all(confint(mm1)[,2] > c(-1, -.1, .3, 1))
# )
#
# expect_true(
#   all(confint(mm2)[,1] < c(-1, -.1, .3, 1)) &
#   all(confint(mm2)[,2] > c(-1, -.1, .3, 1))
# )
#
# expect_true(
#   all(confint(mm3)[,1] < c(-1, -.1, .3, 1)) &
#   all(confint(mm3)[,2] > c(-1, -.1, .3, 1))
# )

# estimate popsize test

expect_silent(est1 <- estimatePopsize(mm1))
expect_silent(est2 <- estimatePopsize(mm2))
expect_silent(est3 <- estimatePopsize(mm3))

# conf int coverage
expect_true(
  all(c(est1$populationSize$confidenceInterval[,1] < 2400,
        est2$populationSize$confidenceInterval[,1] < 2400,
        est3$populationSize$confidenceInterval[,1] < 2400)) &
  all(c(est1$populationSize$confidenceInterval[,2] > 2400,
        est2$populationSize$confidenceInterval[,2] > 2400,
        est3$populationSize$confidenceInterval[,2] > 2400))
)

expect_silent(est1.1 <- estimatePopsize(mm1, popVar = "bootstrap"))
expect_silent(est2.1 <- estimatePopsize(mm2))
