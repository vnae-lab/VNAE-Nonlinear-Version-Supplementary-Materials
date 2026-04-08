# ============================================================================
# VNAE Fourth-Order Cancellation Test - 3 Players
# Based on: "Riemannian Manifolds of Asymmetric Equilibria" (Pereira, 2026)
# Demonstrates that the combination of second derivatives of the metric
# that contains fourth-order derivatives of V is numerically zero.
# ============================================================================

# Load required package for numerical derivatives
if (!require(numDeriv)) install.packages("numDeriv")
library(numDeriv)

# ----------------------------------------------------------------------------
# Core model functions for 3 players
# ----------------------------------------------------------------------------

# Interaction potential (cubic + all bilinear terms)
V <- function(x, y, z) {
  (x^3 / 3) + (y^3 / 3) + (z^3 / 3) + (x * y) + (x * z) + (y * z)
}

# Player-specific structural fields (cubic form)
phi_A <- function(x, theta_A) theta_A * (x^3 / 3)
phi_B <- function(y, theta_B) theta_B * (y^3 / 3)
phi_C <- function(z, theta_C) theta_C * (z^3 / 3)

# Hessian of (V + φ_A + φ_B + φ_C) (3x3 matrix)
hessian_Vphi <- function(x, y, z, theta_A, theta_B, theta_C) {
  # Diagonal elements
  H11 <- 2*x + 2*theta_A*x   # ∂²(V+φ_A)/∂x²
  H22 <- 2*y + 2*theta_B*y   # ∂²(V+φ_B)/∂y²
  H33 <- 2*z + 2*theta_C*z   # ∂²(V+φ_C)/∂z²
  # Off-diagonals (mixed partials of V only)
  H12 <- 1   # ∂²V/∂x∂y
  H13 <- 1   # ∂²V/∂x∂z
  H23 <- 1   # ∂²V/∂y∂z
  matrix(c(H11, H12, H13,
           H12, H22, H23,
           H13, H23, H33), nrow = 3, byrow = TRUE)
}

# Riemannian metric g = ω_i δ_ij + β H_ij
metric <- function(x, y, z, theta_A, theta_B, theta_C, beta, omega = c(1,1,1)) {
  H <- hessian_Vphi(x, y, z, theta_A, theta_B, theta_C)
  g <- diag(omega) + beta * H
  colnames(g) <- rownames(g) <- c("x", "y", "z")
  g
}

# ----------------------------------------------------------------------------
# Compute the residual for the (x,y) block:
# Residual = ∂²g_xx/∂y² + ∂²g_yy/∂x² - 2 * ∂²g_xy/∂x∂y
# This combination should be zero because it contains only fourth-order
# derivatives of V, which cancel due to Schwarz symmetry.
# ----------------------------------------------------------------------------
cancellation_residual_xy <- function(x, y, z, theta_A, theta_B, theta_C, beta, omega = c(1,1,1)) {
 
  # Helper to extract a metric component at a point
  g_fun <- function(xx, yy, zz, i, j) {
    metric(xx, yy, zz, theta_A, theta_B, theta_C, beta, omega)[i, j]
  }
 
  # Second partial derivative of g_xx with respect to y (twice)
  d2gxx_dy2 <- grad(function(p) {
    # first derivative with respect to y at p
    grad(function(q) g_fun(q[1], q[2], q[3], 1, 1), p)[2]
  }, c(x, y, z))[2]
 
  # Second partial derivative of g_yy with respect to x (twice)
  d2gyy_dx2 <- grad(function(p) {
    grad(function(q) g_fun(q[1], q[2], q[3], 2, 2), p)[1]
  }, c(x, y, z))[1]
 
  # Mixed second partial derivative of g_xy: ∂²g_xy/∂x∂y
  d2gxy_dxdy <- grad(function(p) {
    # ∂g_xy/∂x at p
    grad(function(q) g_fun(q[1], q[2], q[3], 1, 2), p)[1]
  }, c(x, y, z))[2]
 
  residual <- d2gxx_dy2 + d2gyy_dx2 - 2 * d2gxy_dxdy
  return(residual)
}

# ----------------------------------------------------------------------------
# Example: parameters from supplementary material (3-player positive curvature)
# ----------------------------------------------------------------------------
cat("\n================== VNAE CANCELLATION TEST (3 PLAYERS) ==================\n")
cat("Test point: (x, y, z) = (0.5, 0.5, 0.5)\n")
cat("Parameters: θ_A = 8, θ_B = 6, θ_C = 13, β = 0.1\n")
cat("Weights ω = (1, 1, 1)\n\n")

res <- cancellation_residual_xy(x = 0.5, y = 0.5, z = 0.5,
                                theta_A = 8.0, theta_B = 6.0, theta_C = 13.0,
                                beta = 0.1)

cat("Residual for (x,y) block (∂²g_xx/∂y² + ∂²g_yy/∂x² - 2·∂²g_xy/∂x∂y):\n")
cat(sprintf("%.12e\n", res))

if (abs(res) < 1e-10) {
  cat("\nSTATUS: CANCELLATION CONFIRMED (residual within numerical precision).\n")
  cat("REASON: The fourth-order derivatives of V cancel analytically\n")
  cat("        due to the symmetry of mixed partials (Schwarz theorem).\n")
} else {
  cat("\nSTATUS: Residual not zero – check implementation or parameters.\n")
}
cat("========================================================================\n")
