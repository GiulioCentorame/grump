set.seed(102)

data <- lm(mpg ~ cyl, data = mtcars)
beta <- coef(summary(data))[, "Estimate"]["cyl"]
se <- coef(summary(data))[, "Std. Error"]["cyl"]

test_that("the Wald test obtains the same p value from known beta and SE", {
  expect_equal(wald_test_p(beta, se), 4.675579e-19)
})
