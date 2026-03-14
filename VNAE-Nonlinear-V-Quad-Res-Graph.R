# =============================================================================
# Residual and Curvature Proxy along Nullcline  (Nonlinear Cubic VNAE)
# =============================================================================
# INSTRUCTIONS: Change only the values in the "PARAMETERS" section below.
# =============================================================================

# --- PARAMETERS (change these) ---
omegaA <- 1.0
omegaB <- 1.0
thetaA <- 4.0      # asymmetry player A
thetaB <- 6.0      # asymmetry player B
beta   <- 0.1      # rigidity

# Domain for x along the nullcline
x_min <- -1.0
x_max <-  1.0
npoints <- 201

# --- DO NOT CHANGE BELOW THIS LINE ---

# Define the expectation field F (nonlinear)
F_x <- function(x, y) (omegaA + thetaA) * x^2 + omegaA * y
F_y <- function(x, y) (omegaB + thetaB) * y^2 + omegaB * x

# Nullcline of F_x: y = -((omegaA+thetaA)/omegaA) * x^2
slope_null <- - (omegaA + thetaA) / omegaA

# Generate points along the nullcline
xs <- seq(x_min, x_max, length.out = npoints)
ys_null <- slope_null * xs^2   # note: x^2, not x (because it's quadratic)

# Compute the residual |F| along this nullcline: we want F_y because F_x is zero by construction
resid <- sapply(seq_along(xs), function(i) {
  # F_x is zero on this curve (by definition), so residual is |F_y|
  abs(F_y(xs[i], ys_null[i]))
})

cat("Max residual |F_y| along the nullcline of F_x:", max(resid), "\n")

# Plot 1: Residual |F| along the nullcline
plot(xs, resid, type = 'l', col = "blue", lwd = 2,
     main = "Residual |F_y| Along the Nullcline (F_x = 0)",
     xlab = "x", ylab = "|F_y|")
grid()

# -----------------------------------------------------------------------------
# Approximate scalar curvature proxy (trace of Hessian of combined potential)
# -----------------------------------------------------------------------------
# Combined potential: U(x,y) = V + phiA + phiB
# V = x^3/3 + y^3/3 + x*y
# phiA = thetaA * x^3/3
# phiB = thetaB * y^3/3
# So U = (1+thetaA)*x^3/3 + (1+thetaB)*y^3/3 + x*y

U <- function(x, y) {
  (1 + thetaA) * x^3 / 3 + (1 + thetaB) * y^3 / 3 + x * y
}

# Finite-difference Hessian trace (Laplacian)
hessian_trace <- function(x, y, eps = 1e-4) {
  U_xx <- (U(x + eps, y) - 2 * U(x, y) + U(x - eps, y)) / eps^2
  U_yy <- (U(x, y + eps) - 2 * U(x, y) + U(x, y - eps)) / eps^2
  U_xx + U_yy
}

# Compute trace along the nullcline
trace_vals <- sapply(seq_along(xs), function(i) {
  hessian_trace(xs[i], ys_null[i])
})

# Plot 2: Hessian trace along the nullcline
plot(xs, trace_vals, type = 'l', col = "red", lwd = 2,
     main = "Hessian Trace Along the Nullcline (Curvature Proxy)",
     xlab = "x", ylab = "trace(Hessian)")
grid()

# Optionally, compute Jacobian at a specific point (e.g., origin)
jacobian_num <- function(x, y, eps = 1e-6) {
  F0 <- c(F_x(x, y), F_y(x, y))
  dFx_dx <- (F_x(x + eps, y) - F_x(x, y)) / eps
  dFx_dy <- (F_x(x, y + eps) - F_x(x, y)) / eps
  dFy_dx <- (F_y(x + eps, y) - F_y(x, y)) / eps
  dFy_dy <- (F_y(x, y + eps) - F_y(x, y)) / eps
  matrix(c(dFx_dx, dFy_dx, dFx_dy, dFy_dy), nrow = 2)
}

J0 <- jacobian_num(0, 0)
cat("\nJacobian at (0,0):\n")
print(J0)
cat("Eigenvalues:\n")
print(eigen(J0)$values)

# -----------------------------------------------------------------------------
# Notes:
# In this nonlinear cubic example, the nullcline of F_x is a parabola.
# The residual |F_y| along it indicates how far we are from the true equilibrium manifold.
# The Hessian trace is a rough proxy for curvature; the exact scalar curvature K
# is given by K = (det(H)/det(g)) * beta * |thetaA - thetaB|, which can be computed
# separately if desired (but we leave it as an exercise).
# -----------------------------------------------------------------------------
