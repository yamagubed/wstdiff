#' Welch-Satterthwaite Approximation for Univariate t-Differences
#'
#' @description
#' Approximates the distribution of the difference between two independent
#' non-standardized t-distributed random variables using the Welch-Satterthwaite
#' method.
#'
#' @details
#' For two independent non-standardized t-distributed random variables:
#' \itemize{
#'   \item X1 ~ t(mu1, sigma1^2, nu1)
#'   \item X2 ~ t(mu2, sigma2^2, nu2)
#' }
#'
#' The difference Z = X1 - X2 is approximated as:
#' Z ~ t(mu1 - mu2, sigma_star^2, nu_star)
#'
#' where the effective parameters are computed through moment matching:
#' \itemize{
#'   \item sigma_star matches the variance of Z
#'   \item nu_star is derived from fourth moment matching
#' }
#'
#' The method requires nu1 > 4 and nu2 > 4 for the existence of fourth moments.
#' The approximation quality improves as degrees of freedom increase and
#' approaches exactness as nu -> infinity (normal limit).
#'
#' @param mu1 Location parameter of first distribution
#' @param sigma1 Scale parameter of first distribution (must be > 0)
#' @param nu1 Degrees of freedom of first distribution (must be > 4)
#' @param mu2 Location parameter of second distribution
#' @param sigma2 Scale parameter of second distribution (must be > 0)
#' @param nu2 Degrees of freedom of second distribution (must be > 4)
#'
#' @return An S3 object of class "ws_tdiff_univariate" containing:
#'   \item{mu_diff}{Location parameter of difference (mu1 - mu2)}
#'   \item{sigma_star}{Effective scale parameter (Equation 1 from paper)}
#'   \item{nu_star}{Effective degrees of freedom (Equation 2 from paper)}
#'   \item{input_params}{List of input parameters for reference}
#'   \item{method}{Character string "univariate"}
#'
#' @references
#' Yamaguchi, Y., Homma, G., Maruo, K., & Takeda, K.
#' Welch-Satterthwaite Approximation for Difference of Non-Standardized
#' t-Distributed Variables. (unpublished).
#'
#' @examples
#' # Example 1: Different scale parameters
#' result <- ws_tdiff_univariate(
#'   mu1 = 0, sigma1 = 1, nu1 = 10,
#'   mu2 = 0, sigma2 = 1.5, nu2 = 15
#' )
#' print(result)
#'
#' # Example 2: Equal parameters (special case)
#' result_equal <- ws_tdiff_univariate(
#'   mu1 = 5, sigma1 = 2, nu1 = 20,
#'   mu2 = 3, sigma2 = 2, nu2 = 20
#' )
#' # Should match ws_tdiff_equal_params(5-3, 2, 20)
#'
#' @seealso
#' \code{\link{ws_tdiff_equal_params}} for the special case of equal parameters
#' \code{\link{dtdiff}}, \code{\link{ptdiff}}, \code{\link{qtdiff}}, \code{\link{rtdiff}}
#' for distribution functions
#'
#' @export
#' @importFrom stats dt pt qt rt
ws_tdiff_univariate <- function(mu1, sigma1, nu1, mu2, sigma2, nu2) {
  # Input validation
  if (nu1 <= 4 || nu2 <= 4) {
    stop("Both nu1 and nu2 must be greater than 4 for fourth moment to exist")
  }
  if (sigma1 <= 0 || sigma2 <= 0) {
    stop("Scale parameters must be positive")
  }
  if (!is.numeric(c(mu1, sigma1, nu1, mu2, sigma2, nu2))) {
    stop("All parameters must be numeric")
  }

  # Compute variance components (Var(X) = sigma^2 * nu/(nu-2) for t-distribution)
  var1 <- sigma1^2 * nu1 / (nu1 - 2)
  var2 <- sigma2^2 * nu2 / (nu2 - 2)

  # Effective scale parameter (Equation 1 from paper)
  sigma_star <- sqrt(var1 + var2)

  # Effective degrees of freedom (Equation 2 from paper)
  # Based on Welch-Satterthwaite principle applied to fourth moments
  numerator <- (var1 + var2)^2
  term1 <- (sigma1^2)^2 * nu1^2 / ((nu1 - 2)^2 * (nu1 - 4))
  term2 <- (sigma2^2)^2 * nu2^2 / ((nu2 - 2)^2 * (nu2 - 4))
  denominator <- term1 + term2
  nu_star <- numerator / denominator

  # Prepare result object
  result <- list(
    mu_diff = mu1 - mu2,
    sigma_star = sigma_star,
    nu_star = nu_star,
    input_params = list(
      mu1 = mu1, sigma1 = sigma1, nu1 = nu1,
      mu2 = mu2, sigma2 = sigma2, nu2 = nu2
    ),
    method = "univariate"
  )
  class(result) <- c("ws_tdiff_univariate", "list")
  return(result)
}
