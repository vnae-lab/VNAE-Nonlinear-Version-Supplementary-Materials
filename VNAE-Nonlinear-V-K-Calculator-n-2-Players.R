# -------------------------------------------------------------------------
# GENERAL VNAE CALCULATOR (2-PLAYER NONLINEAR)
# Based on Pereira (2026)
# -------------------------------------------------------------------------

vnae_calc_2d <- function(thetas, omega, beta, s_point) {
 
  # 1. Elements of the Hessian Matrix (H)
  # Diagonals: Hi = 2 * s_i * (omega_i + theta_i)
  # Off-diagonals: Interaction coupling (standardized to 1.0)
  H11 <- 2 * s_point[1] * (omega[1] + thetas[1])
  H22 <- 2 * s_point[2] * (omega[2] + thetas[2])
  H12 <- 1.0  # Interaction term
  H21 <- 1.0
 
  H <- matrix(c(H11, H12, H21, H22), nrow = 2, byrow = TRUE)
 
  # 2. Jacobian Matrix (J) of the Flow
  # The Flow is s_dot = -F(s). The Jacobian is essentially -H
  J <- -H
 
  # 3. Riemannian Metric (g)
  # g = I + beta * H
  g <- diag(2) + beta * H
 
  # 4. Scalar Curvature (K) - Theorem 4.2
  # K = beta * |theta_A - theta_B| * (det(H) / det(g))
  weight <- abs(thetas[1] - thetas[2])
  K <- beta * weight * (det(H) / det(g))
 
  # -----------------------------------------------------------------------
  # PRINTING RESULTS
  # -----------------------------------------------------------------------
  cat("\n======================================================\n")
  cat("   GENERAL VNAE GEOMETRIC CALCULATOR (2-PLAYER)\n")
  cat("======================================================\n\n")
 
  cat("INPUTS:\n")
  cat("Asymmetries (theta):", thetas, "\n")
  cat("Weights (omega):    ", omega, "\n")
  cat("Rigidity (beta):    ", beta, "\n")
  cat("Current State (s):  ", s_point, "\n\n")
 
  cat("HESSIAN MATRIX (H):\n")
  print(round(H, 4))
 
  cat("\nJACOBIAN MATRIX (J):\n")
  print(round(J, 4))
 
  cat("\nRIEMANNIAN METRIC (g):\n")
  print(round(g, 4))
 
  cat("\n------------------------------------------------------\n")
  cat("SCALAR CURVATURE (K):", round(K, 6), "\n")
  cat("------------------------------------------------------\n")
 
  # Geometry Classification
  if (K > 0.0001) {
    cat("GEOMETRY:  SPHERICAL (K > 0)\n")
    cat("STABILITY: STABLE ATTRACTOR (VNAE Convergence)\n")
  } else if (K < -0.0001) {
    cat("GEOMETRY:  HYPERBOLIC (K < 0)\n")
    cat("STABILITY: UNSTABLE SADDLE (VNAE Divergence)\n")
  } else {
    cat("GEOMETRY:  FLAT (K = 0)\n")
    cat("STABILITY: NEUTRAL (Nash/Von Neumann Limit)\n")
  }
  cat("======================================================\n")
 
  return(invisible(list(K = K, H = H, g = g, J = J)))
}

# ---------------------------------------------------------
# EXAMPLES OF USE:
# ---------------------------------------------------------

# Case 1: K > 0 K = 0 or K < 0
#
vnae_calc_2d(thetas=c(4, 6), omega=c(1,1), beta=0.1, s_point=c(0.5, 0.4))
