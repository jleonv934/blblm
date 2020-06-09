library(pryr)
test_that("blblm class works", {
  fit3 <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  expect_equal(otype(fit3), "S3")
})
