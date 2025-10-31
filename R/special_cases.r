#' Equal Parameters Special Case
#'
#' @description
#' Computes the Welch-Satterthwaite approximation for the special case where
#' both distributions have identical parameters.
#'
#' @details
#' When X1 ~ t(mu, sigma^2, nu) and X2 ~ t(mu, sigma^2, nu) are independent,
#' the difference Z = X1 - X2 simplifies to:
#' \itemize{
#'   \item Location: mu_diff = 0
#'   \item Scale: sigma_star = sigma * sqrt(2*nu/(nu-2))
#'   \item Degrees of freedom: nu_star = 2*(nu - 4)
#' }
#'
#' This special case provides validation for the general formulas and
#' computational efficiency when parameters are known to be equal.
#'
#' @param mu Common location parameter
#' @param sigma Common scale parameter (must be > 0)
#' @param nu Common degrees of freedom (must be > 4)
#'
#' @return An S3 object of class "ws_tdiff_univariate" with the simplified parameters
#'
#' @examples
#' # Equal parameters case
#' result <- ws_tdiff_equal_params(mu = 0, sigma = 1, nu = 10)
#' print(result)
#' # nu_star should be 2*(10-4) = 12
#'
#' # Verify against general formula
#' general <- ws_tdiff_univariate(0, 1, 10, 0, 1, 10)
#' all.equal(result$nu_star, general$nu_star)
#'
#' @export
ws_tdiff_equal_params <- function(mu, sigma, nu) {
  if (nu <= 4) stop("Degrees of freedom must be greater than 4")
  if (sigma <= 0) stop("Scale parameter must be positive")

  # Simplified formulas from Proposition 1
  sigma_star <- sigma * sqrt(2 * nu / (nu - 2))
  nu_star <- 2 * (nu - 4)

  result <- list(
    mu_diff = 0,  # Difference of equal means is zero
    sigma_star = sigma_star,
    nu_star = nu_star,
    input_params = list(mu = mu, sigma = sigma, nu = nu),
    method = "equal_params"
  )
  class(result) <- c("ws_tdiff_univariate", "list")
  return(result)
}
