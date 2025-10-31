#' @export
print.ws_tdiff_univariate <- function(x, ...) {
  cat("Welch-Satterthwaite Approximation (Univariate)\n")
  cat("==============================================\n")
  cat(sprintf("Location (mu1 - mu2): %.4f\n", x$mu_diff))
  cat(sprintf("Effective scale (sigma*): %.4f\n", x$sigma_star))
  cat(sprintf("Effective df (nu*): %.4f\n", x$nu_star))
  invisible(x)
}

#' @export
print.ws_tdiff_multivariate_independent <- function(x, ...) {
  cat("Welch-Satterthwaite Approximation (Multivariate Independent)\n")
  cat("=============================================================\n")
  cat("Location difference:\n"); print(x$mu_diff)
  cat("\nEffective scale:\n"); print(x$sigma_star)
  cat("\nEffective df:\n"); print(x$nu_star)
  invisible(x)
}

#' @export
print.ws_tdiff_multivariate_general <- function(x, ...) {
  cat("Welch-Satterthwaite Approximation (General Multivariate)\n")
  cat("=========================================================\n")
  cat("Location difference:\n"); print(x$mu_diff)
  cat("\nEffective scale matrix:\n"); print(x$Sigma_star)
  cat(sprintf("\nEffective df: %.4f\n", x$nu_star))
  cat(sprintf("Converged: %s (%d iterations)\n", x$converged, x$iterations))
  invisible(x)
}
