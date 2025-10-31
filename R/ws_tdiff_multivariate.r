#' Welch-Satterthwaite Approximation for Multivariate t-Differences (Independent)
#'
#' @description
#' Approximates the distribution of differences between two independent
#' p-dimensional vectors with independent t-distributed components.
#'
#' @details
#' This function applies the univariate Welch-Satterthwaite approximation
#' component-wise when all components are mutually independent. Each
#' component difference Zj = X1j - X2j is approximated independently using
#' the univariate method.
#'
#' This approach is optimal for:
#' \itemize{
#'   \item Marginal inference on specific components
#'   \item Cases where components have different tail behaviors
#'   \item Maintaining computational efficiency in high dimensions
#' }
#'
#' @param mu1 Location vector of first distribution (length p)
#' @param sigma1 Scale vector of first distribution (length p, all > 0)
#' @param nu1 Degrees of freedom vector of first distribution (length p, all > 4)
#' @param mu2 Location vector of second distribution (length p)
#' @param sigma2 Scale vector of second distribution (length p, all > 0)
#' @param nu2 Degrees of freedom vector of second distribution (length p, all > 4)
#'
#' @return An S3 object of class "ws_tdiff_multivariate_independent" containing:
#'   \item{mu_diff}{Location vector of difference}
#'   \item{sigma_star}{Vector of effective scale parameters}
#'   \item{nu_star}{Vector of effective degrees of freedom}
#'   \item{p}{Dimension of the vectors}
#'   \item{method}{Character string "multivariate_independent"}
#'
#' @examples
#' result <- ws_tdiff_multivariate_independent(
#'   mu1 = c(0, 1), sigma1 = c(1, 1.5), nu1 = c(10, 12),
#'   mu2 = c(0, 0), sigma2 = c(1.2, 1), nu2 = c(15, 20)
#' )
#' print(result)
#'
#' @seealso \code{\link{ws_tdiff_multivariate_general}} for correlated components
#'
#' @export
ws_tdiff_multivariate_independent <- function(mu1, sigma1, nu1, mu2, sigma2, nu2) {
  p <- length(mu1)

  # Input validation
  if (length(sigma1) != p || length(nu1) != p ||
      length(mu2) != p || length(sigma2) != p || length(nu2) != p) {
    stop("All input vectors must have the same length")
  }
  if (any(nu1 <= 4) || any(nu2 <= 4)) {
    stop("All degrees of freedom must be greater than 4")
  }
  if (any(sigma1 <= 0) || any(sigma2 <= 0)) {
    stop("All scale parameters must be positive")
  }

  # Apply univariate approximation to each component
  sigma_star <- numeric(p)
  nu_star <- numeric(p)

  for (j in 1:p) {
    result <- ws_tdiff_univariate(mu1[j], sigma1[j], nu1[j],
                                  mu2[j], sigma2[j], nu2[j])
    sigma_star[j] <- result$sigma_star
    nu_star[j] <- result$nu_star
  }

  result <- list(
    mu_diff = mu1 - mu2,
    sigma_star = sigma_star,
    nu_star = nu_star,
    p = p,
    method = "multivariate_independent"
  )
  class(result) <- c("ws_tdiff_multivariate_independent", "list")
  return(result)
}

#' Welch-Satterthwaite Approximation for General Multivariate t-Differences
#'
#' @description
#' Approximates the distribution of differences between two independent
#' multivariate t-distributed random vectors with arbitrary covariance
#' structure. This implements Theorem 3 from Yamaguchi et al. (2025).
#'
#' @details
#' This function handles the general case where components may be correlated
#' within each multivariate t-distribution. The approximation uses a single
#' scalar degrees of freedom parameter to capture the overall tail behavior.
#'
#' The iterative algorithm (Section 4.3 of the paper):
#' \enumerate{
#'   \item Initialize with sum of covariance matrices
#'   \item Compute effective degrees of freedom using trace formulas
#'   \item Update scale matrix
#'   \item Iterate until convergence
#' }
#'
#' Note: For high dimensions with heterogeneous component behaviors,
#' consider using \code{\link{ws_tdiff_multivariate_independent}} instead.
#'
#' @param mu1 Location vector of first distribution (length p)
#' @param Sigma1 Scale matrix of first distribution (p x p, positive definite)
#' @param nu1 Degrees of freedom of first distribution (must be > 4)
#' @param mu2 Location vector of second distribution (length p)
#' @param Sigma2 Scale matrix of second distribution (p x p, positive definite)
#' @param nu2 Degrees of freedom of second distribution (must be > 4)
#' @param max_iter Maximum iterations for convergence (default: 10)
#' @param tol Convergence tolerance (default: 1e-6)
#'
#' @return An S3 object of class "ws_tdiff_multivariate_general" containing:
#'   \item{mu_diff}{Location vector of difference}
#'   \item{Sigma_star}{Effective scale matrix}
#'   \item{nu_star}{Effective degrees of freedom (scalar)}
#'   \item{converged}{Logical indicating convergence}
#'   \item{iterations}{Number of iterations performed}
#'   \item{method}{Character string "multivariate_general"}
#'
#' @examples
#' Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
#' result <- ws_tdiff_multivariate_general(
#'   mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
#'   mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
#' )
#' print(result)
#'
#' @export
ws_tdiff_multivariate_general <- function(mu1, Sigma1, nu1, mu2, Sigma2, nu2,
                                          max_iter = 10, tol = 1e-6) {
  # Input validation
  if (nu1 <= 4 || nu2 <= 4) {
    stop("Both nu1 and nu2 must be greater than 4")
  }
  if (!is.matrix(Sigma1) || !is.matrix(Sigma2)) {
    stop("Sigma1 and Sigma2 must be matrices")
  }
  if (nrow(Sigma1) != ncol(Sigma1) || nrow(Sigma2) != ncol(Sigma2)) {
    stop("Sigma matrices must be square")
  }

  p <- length(mu1)
  if (nrow(Sigma1) != p || nrow(Sigma2) != p || length(mu2) != p) {
    stop("Dimensions of mu and Sigma must be consistent")
  }
  if (!isSymmetric(Sigma1) || !isSymmetric(Sigma2)) {
    warning("Sigma matrices should be symmetric; using symmetrized version")
    Sigma1 <- (Sigma1 + t(Sigma1)) / 2
    Sigma2 <- (Sigma2 + t(Sigma2)) / 2
  }

  # Check positive definiteness
  ev1 <- eigen(Sigma1, only.values = TRUE)$values
  ev2 <- eigen(Sigma2, only.values = TRUE)$values
  if (any(ev1 <= 0) || any(ev2 <= 0)) {
    stop("Sigma matrices must be positive definite")
  }

  # Initialize with sum of covariance matrices (Section 4.3 step 1)
  Sigma_init <- Sigma1 * nu1 / (nu1 - 2) + Sigma2 * nu2 / (nu2 - 2)
  Sigma_star <- Sigma_init
  converged <- FALSE

  # Iterative procedure (Section 4.3)
  for (iter in 1:max_iter) {
    # Compute traces for current Sigma_star
    trace_Sigma2 <- sum(diag(Sigma_star %*% Sigma_star))
    trace_Sigma <- sum(diag(Sigma_star))

    # Numerator of Equation (6)
    numerator <- 2 * (trace_Sigma2 + trace_Sigma^2)

    # Q4 computation (denominator of Equation 6)
    Q4_term1 <- nu1^2 / ((nu1 - 2)^2 * (nu1 - 4))
    trace_S1_2 <- sum(diag(Sigma1 %*% Sigma1))
    trace_S1 <- sum(diag(Sigma1))
    Q4_1 <- Q4_term1 * (trace_S1_2 + trace_S1^2)

    Q4_term2 <- nu2^2 / ((nu2 - 2)^2 * (nu2 - 4))
    trace_S2_2 <- sum(diag(Sigma2 %*% Sigma2))
    trace_S2 <- sum(diag(Sigma2))
    Q4_2 <- Q4_term2 * (trace_S2_2 + trace_S2^2)

    Q4 <- Q4_1 + Q4_2

    # Compute nu_star (Equation 6)
    nu_star <- numerator / Q4

    # Update Sigma_star (Section 4.3 step 3)
    Sigma_star_new <- ((nu_star - 2) / nu_star) * Sigma_init

    # Check convergence
    if (max(abs(Sigma_star_new - Sigma_star)) < tol) {
      Sigma_star <- Sigma_star_new
      converged <- TRUE
      break
    }
    Sigma_star <- Sigma_star_new
  }

  result <- list(
    mu_diff = mu1 - mu2,
    Sigma_star = Sigma_star,
    nu_star = nu_star,
    converged = converged,
    iterations = if(converged) iter else max_iter,
    method = "multivariate_general"
  )
  class(result) <- c("ws_tdiff_multivariate_general", "list")
  return(result)
}
