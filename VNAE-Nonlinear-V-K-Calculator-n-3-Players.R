# -------------------------------------------------------------------------
# VICTORIA–NASH ASYMMETRIC EQUILIBRIUM (VNAE) CALCULATOR
# Based on Pereira (2026) Framework
# -------------------------------------------------------------------------

vnae_analysis <- function(thetas, beta, omega, points) {
 
  # Number of players
  n <- length(thetas)
 
  # -----------------------------------------------------
  # 1. HESSIAN MATRIX (H)
  # -----------------------------------------------------
  # Diagonals: 2 * s_i * (omega_i + theta_i)
  # Off-diagonals: 1.0 (interaction coupling)
 
  H <- matrix(1, nrow = n, ncol = n)
  for(i in 1:n) {
    H[i,i] <- 2 * points[i] * (omega[i] + thetas[i])
  }
 
  # -----------------------------------------------------
  # 2. RIEMANNIAN METRIC (g)
  # -----------------------------------------------------
  # g = I + beta * H
 
  g <- diag(1, n) + beta * H
 
  # -----------------------------------------------------
  # 3. SCALAR CURVATURE (K)
  # -----------------------------------------------------
  # K ≈ β Σ_{i<j} |θ_i - θ_j| det(H_ij) / det(g_ij)
 
  K_total <- 0
 
  if(n >= 2) {
    pairs <- combn(1:n, 2)
   
    for(p in 1:ncol(pairs)) {
      i <- pairs[1,p]
      j <- pairs[2,p]
     
      sub_H <- H[c(i,j), c(i,j)]
      sub_g <- g[c(i,j), c(i,j)]
     
      weight <- abs(thetas[i] - thetas[j])
     
      K_pair <- beta * weight * (det(sub_H) / det(sub_g))
      K_total <- K_total + K_pair
    }
  }
 
  # -----------------------------------------------------
  # OUTPUT
  # -----------------------------------------------------
 
  cat("====================================================\n")
  cat("   VNAE NONLINEAR GEOMETRIC CALCULATOR (n =", n, ")\n")
  cat("====================================================\n\n")
 
  cat("STRATEGIC STATE (s):", paste(points, collapse=", "), "\n")
  cat("ASYMMETRY PARAMETERS (theta):", paste(thetas, collapse=", "), "\n")
  cat("RIGIDITY COEFFICIENT (beta):", beta, "\n\n")
 
  cat("HESSIAN MATRIX (H):\n")
  print(round(H, 6))
 
  cat("\nRIEMANNIAN METRIC (g):\n")
  print(round(g, 6))
 
  cat("\n----------------------------------------------------\n")
  cat("FINAL SCALAR CURVATURE (K):", round(K_total, 6), "\n")
  cat("----------------------------------------------------\n")
 
  if(K_total > 0) {
    cat("GEOMETRY: SPHERICAL (K > 0)\n")
    cat("STABILITY: STABLE ATTRACTOR\n")
    cat("RESULT: Robust equilibrium (geometric contraction)\n")
  } else if(K_total < 0) {
    cat("GEOMETRY: HYPERBOLIC (K < 0)\n")
    cat("STABILITY: UNSTABLE REPELLOR\n")
    cat("RESULT: Divergent regime (geometric repulsion)\n")
  } else {
    cat("GEOMETRY: FLAT (K = 0)\n")
    cat("STABILITY: NEUTRAL (Nash/von Neumann limit)\n")
  }
 
  cat("====================================================\n")
 
  invisible(list(K = K_total, H = H, g = g))
}

# ---------------------------------------------------------
# TEST EXAMPLE
# ---------------------------------------------------------

vnae_analysis(
  thetas = c(8, 6, 13),
  beta   = 0.1,
  omega  = c(1, 1, 1),
  points = c(0.5, 0.5, 0.5)
)
