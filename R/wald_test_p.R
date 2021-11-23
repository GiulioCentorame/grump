#' Simple calculation of p-value from a Wald test
#'
#' @importFrom stats pnorm
#'
#' @param beta Regression coefficient
#' @param se Standard error of the model
#'
#' @return P-value
#' @export
#'
#' @examples
#' data <- lm(mpg ~ cyl, data = mtcars)
#' beta <- coef(summary(data))[, "Estimate"]["cyl"]
#' se <- coef(summary(data))[, "Std. Error"]["cyl"]
#'
#' wald_test_p(beta, se)
wald_test_p <- function(beta, se){
  2*pnorm(q = abs(beta/se),
          mean = 0,
          sd = 1,
          lower.tail = FALSE,
          log.p = FALSE)
}
