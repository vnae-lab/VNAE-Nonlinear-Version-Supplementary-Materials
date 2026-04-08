# ============================================================================
# VNAE Fourth-Order Cancellation Test
# Based on: "Riemannian Manifolds of Asymmetric Equilibria" (Pereira, 2026)
# Demonstrates that the combination of second derivatives of the metric
# that contains fourth-order derivatives of V is numerically zero.
# ============================================================================

# Load required package for numerical derivatives
if (!require(numDeriv)) install.packages("numDeriv")
library(numDeriv)

# ----------------------------------------------------------------------------
# Core model functions (same as in the VNAE geometry)
# ----------------------------------------------------------------------------

# Interaction potential (cubic + bilinear)
V <- function(x, y) (x^3 / 3) + (y^3 / 3) + (x * y)

# Player-specific structural fields (cubic form)
phi_A <- function(x, theta_A) theta_A * (x^3 / 3)
phi_B <- function(y, theta_B) theta_B * (y^3 / 3)

# Hessian of (V + φ) (2x2 matrix)
hessian_Vphi <- function(x, y, theta_A, theta_B) {
  H11 <- 2*x + 2*theta_A*x   # ∂²(V+φ_A)/∂x²
  H22 <- 2*y + 2*theta_B*y   # ∂²(V+φ_B)/∂y²
  H12 <- 1                   # ∂²V/∂x∂y
  matrix(c(H11, H12, H12, H22), nrow = 2, byrow = TRUE)
}

# Riemannian metric g = ω_i δ_ij + β H_ij
metric <- function(x, y, theta_A, theta_B, beta, omega = c(1,1)) {
  H <- hessian_Vphi(x, y, theta_A, theta_B)
  g <- diag(omega) + beta * H
  # assign names for clarity
  colnames(g) <- rownames(g) <- c("x", "y")
  g
}

# ----------------------------------------------------------------------------
# Compute the residual that should be zero by Schwarz symmetry
# Residual = ∂²g_11/∂y² + ∂²g_22/∂x² - 2 * ∂²g_12/∂x∂y
# ----------------------------------------------------------------------------
cancellation_residual <- function(x, y, theta_A, theta_B, beta, omega = c(1,1)) {
 
  # Helper to extract a metric component at a point
  g_fun <- function(xx, yy, i, j) {
    metric(xx, yy, theta_A, theta_B, beta, omega)[i, j]
  }
 
  # Second partial derivative of g11 with respect to y (twice)
  d2g11_dy2 <- grad(function(p) {
    # first derivative with respect to y at p
    grad(function(q) g_fun(q[1], q[2], 1, 1), p)[2]
  }, c(x, y))[2]
 
  # Second partial derivative of g22 with respect to x (twice)
  d2g22_dx2 <- grad(function(p) {
    grad(function(q) g_fun(q[1], q[2], 2, 2), p)[1]
  }, c(x, y))[1]
 
  # Mixed second partial derivative of g12: ∂²g12/∂x∂y
  # We compute gradient of ∂g12/∂x w.r.t y
  d2g12_dxdy <- grad(function(p) {
    # ∂g12/∂x at p
    grad(function(q) g_fun(q[1], q[2], 1, 2), p)[1]
  }, c(x, y))[2]
 
  residual <- d2g11_dy2 + d2g22_dx2 - 2 * d2g12_dxdy
  return(residual)
}

# ----------------------------------------------------------------------------
# Example: positive curvature regime (θ_A=4, θ_B=6, β=0.1, point (0.5,0.4))
# ----------------------------------------------------------------------------
cat("\n================== VNAE CANCELLATION TEST ==================\n")
cat("Test point: (x, y) = (0.5, 0.4)\n")
cat("Parameters: θ_A = 4, θ_B = 6, β = 0.1\n")
cat("Weights ω = (1, 1)\n\n")

res <- cancellation_residual(x = 0.5, y = 0.4,
                             theta_A = 4.0, theta_B = 6.0,
                             beta = 0.1)

cat("Residual (∂²g11/∂y² + ∂²g22/∂x² - 2·∂²g12/∂x∂y):\n")
cat(sprintf("%.12e\n", res))

if (abs(res) < 1e-10) {
  cat("\nSTATUS: CANCELLATION CONFIRMED (residual within numerical precision).\n")
  cat("REASON: The fourth-order derivatives of V cancel analytically\n")
  cat("        due to the symmetry of mixed partials (Schwarz theorem).\n")
} else {
  cat("\nSTATUS: Residual not zero – check implementation or parameters.\n")
}
cat("===============================================================\n")
