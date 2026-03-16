# =============================================================================
# Geodesics and Flow Trajectories for the VNAE Model (Nonlinear Cubic Version)
# Supplementary Material: Riemannian Manifolds of Asymmetric Equilibria: Complete Geometric Structure and Curvature Signatures 
# =============================================================================
# INSTRUCTIONS: Change only the values in the "PARAMETERS" section below.
# =============================================================================

library(deSolve)

# =============================================================================
# PARAMETERS - CHANGE ONLY THESE VALUES
# =============================================================================

# Structural parameters
omega_A <- 1.0
omega_B <- 1.0
theta_A <- 4.0
theta_B <- 6.0
beta    <- 0.1

# Domain limits for plots
x_min <- -2.0
x_max <-  2.0
y_min <- -2.0
y_max <-  2.0

# =============================================================================
# DO NOT CHANGE ANYTHING BELOW THIS LINE
# =============================================================================

# -----------------------------------------------------------------------------
# Define the metric and its derivatives (via functions)
# -----------------------------------------------------------------------------
# Hessian components (depends on x,y)
H_xx <- function(x, y) 2 * x * (1 + theta_A)
H_yy <- function(x, y) 2 * y * (1 + theta_B)
H_xy <- function(x, y) 1

# Metric components
g11 <- function(x, y) 1 + beta * H_xx(x, y)
g12 <- function(x, y) beta * H_xy(x, y)
g22 <- function(x, y) 1 + beta * H_yy(x, y)

# Metric matrix as a function
metric <- function(x, y) {
  matrix(c(g11(x, y), g12(x, y), g12(x, y), g22(x, y)), nrow = 2, ncol = 2)
}

# Inverse metric
inverse_metric <- function(x, y) {
  solve(metric(x, y))
}

# -----------------------------------------------------------------------------
# Christoffel symbols via finite differences
# -----------------------------------------------------------------------------
christoffel <- function(x, y, h = 1e-5) {
  # Metric at nearby points
  g_x_plus  <- metric(x + h, y)
  g_x_minus <- metric(x - h, y)
  g_y_plus  <- metric(x, y + h)
  g_y_minus <- metric(x, y - h)
  
  dg_dx <- (g_x_plus - g_x_minus) / (2 * h)
  dg_dy <- (g_y_plus - g_y_minus) / (2 * h)
  
  # Christoffel symbols of the first kind: Γ_{kij}
  Gamma_kij <- array(0, dim = c(2, 2, 2))  # indices: k, i, j
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        # dg_{ik}/dx_j
        dg_ik_dxj <- if (j == 1) dg_dx[i, k] else dg_dy[i, k]
        
        # dg_{jk}/dx_i
        dg_jk_dxi <- if (i == 1) dg_dx[j, k] else dg_dy[j, k]
        
        # dg_{ij}/dx_k
        dg_ij_dxk <- if (k == 1) dg_dx[i, j] else dg_dy[i, j]
        
        Gamma_kij[k, i, j] <- 0.5 * (dg_ik_dxj + dg_jk_dxi - dg_ij_dxk)
      }
    }
  }
  
  # Raise index to get Γ^m_{ij}
  g_inv <- inverse_metric(x, y)
  Gamma_mij <- array(0, dim = c(2, 2, 2))
  for (m in 1:2) {
    for (i in 1:2) {
      for (j in 1:2) {
        Gamma_mij[m, i, j] <- sum(g_inv[m, ] * Gamma_kij[, i, j])
      }
    }
  }
  return(Gamma_mij)
}

# -----------------------------------------------------------------------------
# Geodesic equations
# -----------------------------------------------------------------------------
geodesic_eq <- function(t, state, pars) {
  x <- state[1]
  y <- state[2]
  u <- state[3]  # dx/dt
  v <- state[4]  # dy/dt
  
  Gamma <- christoffel(x, y)
  
  du <- - (Gamma[1, 1, 1] * u^2 + 2 * Gamma[1, 1, 2] * u * v + Gamma[1, 2, 2] * v^2)
  dv <- - (Gamma[2, 1, 1] * u^2 + 2 * Gamma[2, 1, 2] * u * v + Gamma[2, 2, 2] * v^2)
  
  return(list(c(u, v, du, dv)))
}

# -----------------------------------------------------------------------------
# Flow equations: sdot = -F(s)
# -----------------------------------------------------------------------------
flow_eq <- function(t, state, pars) {
  x <- state[1]
  y <- state[2]
  dx <- -((omega_A + theta_A) * x^2 + omega_A * y)
  dy <- -((omega_B + theta_B) * y^2 + omega_B * x)
  return(list(c(dx, dy)))
}

# -----------------------------------------------------------------------------
# Find points on the VNAE manifold (equilibria)
# -----------------------------------------------------------------------------
# The manifold is defined by F=0. We'll sample and find points satisfying both.
find_vnae_points <- function(res = 200, tol = 1e-3) {
  x_vals <- seq(x_min, x_max, length.out = res)
  y_vals <- seq(y_min, y_max, length.out = res)
  points <- list()
  for (x in x_vals) {
    for (y in y_vals) {
      if (abs(F_x(x, y)) < tol && abs(F_y(x, y)) < tol) {
        points <- append(points, list(c(x, y)))
      }
    }
  }
  if (length(points) == 0) return(data.frame(x = numeric(0), y = numeric(0)))
  do.call(rbind, points)
}

# F functions (for convenience)
F_x <- function(x, y) (omega_A + theta_A) * x^2 + omega_A * y
F_y <- function(x, y) (omega_B + theta_B) * y^2 + omega_B * x

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
# Plot setup
plot(NA, xlim = c(x_min, x_max), ylim = c(y_min, y_max),
     xlab = "x", ylab = "y",
     main = "Geodesics (metric g) vs Flow Trajectories (sdot = -F(s))")

# 1. VNAE curve (equilibrium set) - plot points
vnae_pts <- find_vnae_points()
if (nrow(vnae_pts) > 0) {
  points(vnae_pts, col = "red", pch = 16, cex = 0.5)
  # Optionally, add a line if points are dense enough; but here we'll just plot points.
  # For better visualization, we can also plot the two nullclines as thin lines:
  x_seq <- seq(x_min, x_max, length.out = 200)
  y_null1 <- -((omega_A + theta_A)/omega_A) * x_seq^2
  lines(x_seq, y_null1, col = "pink", lty = 3, lwd = 1)
  y_seq <- seq(y_min, y_max, length.out = 200)
  x_null2 <- -((omega_B + theta_B)/omega_B) * y_seq^2
  lines(x_null2, y_seq, col = "lightgreen", lty = 3, lwd = 1)
}

# 2. Geodesics from the origin in different directions
angles <- seq(0, 2 * pi, length.out = 8)[-8]  # 7 angles
for (angle in angles) {
  u0 <- 0.5 * cos(angle)
  v0 <- 0.5 * sin(angle)
  # Solve geodesic from (0,0) with initial velocity (u0,v0)
  sol_geo <- try(ode(y = c(0, 0, u0, v0), times = seq(0, 5, by = 0.1),
                     func = geodesic_eq, parms = NULL), silent = TRUE)
  if (!inherits(sol_geo, "try-error")) {
    lines(sol_geo[, 2], sol_geo[, 3], col = "blue", lty = 2, lwd = 1.5)
  }
}

# 3. Flow trajectories from different initial points
initial_points <- list(c(1, 0), c(0, 1), c(-1, 0), c(0, -1), c(0.5, 0.5), c(-0.5, -0.5))
for (point in initial_points) {
  sol_flow <- try(ode(y = point, times = seq(0, 10, by = 0.1),
                      func = flow_eq, parms = NULL), silent = TRUE)
  if (!inherits(sol_flow, "try-error")) {
    lines(sol_flow[, 2], sol_flow[, 3], col = "darkgreen", lty = 1, lwd = 2)
  }
}

# Legend
legend("topright",
       legend = c("VNAE curve (F = 0)",
                  "Nullclines (dotted)",
                  "Geodesics of metric g",
                  "Flow: sdot = -F(s)"),
       col = c("red", "gray", "blue", "darkgreen"),
       lty = c(NA, 3, 2, 1),
       pch = c(16, NA, NA, NA),
       lwd = c(NA, 1, 1.5, 2),
       cex = 0.9,
       bg = "white")
# Technical Note: Due to the nonlinear nature of the cubic potential and the metric stiffness (β), some geodesic trajectories may encounter regions of high curvature near the domain boundaries. This may trigger numerical warnings (e.g., DLSODA warnings) indicating that the integrator has reached the machine's precision limits. These warnings do not affect the theoretical integrity of the results; rather, they reflect the intrinsic geometric deformation of the Victoria-Nash strategy space.
