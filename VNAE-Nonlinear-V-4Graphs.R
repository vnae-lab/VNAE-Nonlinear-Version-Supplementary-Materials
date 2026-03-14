# =============================================================================
# Victoria-Nash Geometry Visualization (Nonlinear Cubic Formulation)
# =============================================================================
# INSTRUCTIONS: Change only the values in the "PARAMETERS" section below.
# =============================================================================

library(ggplot2)
library(reshape2)
library(gridExtra)
library(viridis)

# =============================================================================
# PARAMETERS - CHANGE ONLY THESE VALUES
# =============================================================================

# Structural parameters
omega_A <- 1.0    # Structural weight of player A
omega_B <- 1.0    # Structural weight of player B
theta_A <- 0.5    # Asymmetry of player A
theta_B <- 1.0    # Asymmetry of player B
beta    <- 0.1    # Rigidity coefficient

# Fixed point for curvature analysis
x_fixed <- 0.3    # x-coordinate of the fixed point
y_fixed <- -0.135 # y-coordinate of the fixed point

# Domain limits for plots
x_min <- -1.0     # Lower limit of x-axis
x_max <-  1.0     # Upper limit of x-axis
y_min <- -1.0     # Lower limit of y-axis
y_max <-  1.0     # Upper limit of y-axis

# =============================================================================
# DO NOT CHANGE ANYTHING BELOW THIS LINE
# =============================================================================

# Helper functions
F_x <- function(x, y) (omega_A + theta_A) * x^2 + omega_A * y
F_y <- function(x, y) (omega_B + theta_B) * y^2 + omega_B * x

H_xx <- function(x, y) 2 * x * (1 + theta_A)
H_yy <- function(x, y) 2 * y * (1 + theta_B)
H_xy <- function(x, y) 1

g_xx <- function(x, y) 1 + beta * H_xx(x, y)
g_yy <- function(x, y) 1 + beta * H_yy(x, y)
g_xy <- function(x, y) beta

curvature <- function(x, y) {
  H11 <- H_xx(x, y)
  H22 <- H_yy(x, y)
  H12 <- H_xy(x, y)
  detH <- H11 * H22 - H12^2
  g11 <- g_xx(x, y)
  g22 <- g_yy(x, y)
  g12 <- g_xy(x, y)
  detg <- g11 * g22 - g12^2
  (detH / detg) * beta * abs(theta_A - theta_B)
}

nullcline_x <- function(x) -((omega_A + theta_A) / omega_A) * x^2
nullcline_y <- function(y) -((omega_B + theta_B) / omega_B) * y^2

find_equilibria <- function(tol = 1e-4) {
  x_vals <- seq(x_min, x_max, length.out = 500)
  y_vals <- seq(y_min, y_max, length.out = 500)
  points <- list()
  
  for (x in x_vals) {
    y <- nullcline_x(x)
    if (y >= y_min && y <= y_max) {
      if (abs(F_y(x, y)) < tol) {
        points <- append(points, list(c(x, y)))
      }
    }
  }
  for (y in y_vals) {
    x <- nullcline_y(y)
    if (x >= x_min && x <= x_max) {
      if (abs(F_x(x, y)) < tol) {
        # Avoid duplicates
        dup <- any(sapply(points, function(p) abs(p[1]-x) < 1e-3 && abs(p[2]-y) < 1e-3))
        if (!dup) {
          points <- append(points, list(c(x, y)))
        }
      }
    }
  }
  if (length(points) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  do.call(rbind, points)
}

# -----------------------------------------------------------------------------
# Plot 1: Expectation field and nullclines
# -----------------------------------------------------------------------------
create_field_plot <- function() {
  x <- seq(x_min, x_max, length.out = 20)
  y <- seq(y_min, y_max, length.out = 20)
  grid <- expand.grid(x = x, y = y)
  grid$F_x <- F_x(grid$x, grid$y)
  grid$F_y <- F_y(grid$x, grid$y)
  
  x_curve <- seq(x_min, x_max, length.out = 200)
  y_curve1 <- nullcline_x(x_curve)
  y_curve2 <- seq(y_min, y_max, length.out = 200)
  x_curve2 <- nullcline_y(y_curve2)
  
  eq_df <- find_equilibria()
  
  p <- ggplot(grid, aes(x = x, y = y)) +
    geom_segment(aes(xend = x + F_x/8, yend = y + F_y/8),
                 arrow = arrow(length = unit(0.1, "cm")), alpha = 0.6) +
    geom_line(data = data.frame(x = x_curve, y = y_curve1),
              aes(x = x, y = y), color = "blue", size = 1, linetype = "dashed") +
    geom_line(data = data.frame(x = x_curve2, y = y_curve2),
              aes(x = x, y = y), color = "green", size = 1, linetype = "dashed") +
    labs(title = "VNAE: Expectation Field F(s;θ) and Nullclines",
         subtitle = paste0("θ_A = ", theta_A, ", θ_B = ", theta_B, ", β = ", beta),
         x = "Strategy x (Player A)", y = "Strategy y (Player B)") +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max))
  
  if (nrow(eq_df) > 0) {
    p <- p + geom_point(data = eq_df, aes(x = x, y = y), color = "red", size = 3)
  }
  p
}

# -----------------------------------------------------------------------------
# Plot 2: Curvature as function of θ_A (at fixed point)
# -----------------------------------------------------------------------------
create_curvature_plot <- function() {
  theta_A_range <- seq(0.5, 8, length.out = 100)
  K_values <- sapply(theta_A_range, function(thA) {
    # Local calculation for each theta_A
    H11_loc <- 2 * x_fixed * (1 + thA)
    H22_loc <- 2 * y_fixed * (1 + theta_B)
    H12_loc <- 1
    detH_loc <- H11_loc * H22_loc - H12_loc^2
    g11_loc <- 1 + beta * H11_loc
    g22_loc <- 1 + beta * H22_loc
    g12_loc <- beta
    detg_loc <- g11_loc * g22_loc - g12_loc^2
    (detH_loc / detg_loc) * beta * abs(thA - theta_B)
  })
  
  df <- data.frame(theta_A = theta_A_range, K = K_values)
  K_current <- curvature(x_fixed, y_fixed)  # with current theta_A
  
  ggplot(df, aes(x = theta_A, y = K)) +
    geom_line(color = "blue", size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(aes(x = theta_A, y = K_current), color = "red", size = 3) +
    labs(title = "VNAE Manifold Curvature Analysis",
         subtitle = paste0("K(θ_A = ", theta_A, ", θ_B = ", theta_B, ") = ", round(K_current, 4)),
         x = "θ_A (Structural Asymmetry)", y = "Gaussian Curvature K") +
    theme_minimal()
}

# -----------------------------------------------------------------------------
# Plot 3: Phase portrait (dynamics -F)
# -----------------------------------------------------------------------------
create_phase_portrait <- function() {
  x <- seq(x_min, x_max, length.out = 20)
  y <- seq(y_min, y_max, length.out = 20)
  grid <- expand.grid(x = x, y = y)
  grid$dx <- -F_x(grid$x, grid$y)
  grid$dy <- -F_y(grid$x, grid$y)
  grid$magnitude <- sqrt(grid$dx^2 + grid$dy^2)
  
  eq_df <- find_equilibria()
  
  p <- ggplot(grid, aes(x = x, y = y)) +
    geom_segment(aes(xend = x + dx/5, yend = y + dy/5, color = magnitude),
                 arrow = arrow(length = unit(0.08, "cm"))) +
    scale_color_viridis(name = "Field Strength") +
    labs(title = "VNAE Phase Portrait (ẋ = -F)",
         subtitle = paste0("θ_A = ", theta_A, ", θ_B = ", theta_B),
         x = "Strategy x", y = "Strategy y") +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max))
  
  if (nrow(eq_df) > 0) {
    p <- p + geom_point(data = eq_df, aes(x = x, y = y), color = "red", size = 3)
  }
  p
}

# -----------------------------------------------------------------------------
# Plot 4: Symmetric vs Asymmetric Comparison
# -----------------------------------------------------------------------------
create_comparison_plot <- function() {
  theta_sym <- 5.0
  # Local functions for symmetric case
  F_x_sym <- function(x, y) (omega_A + theta_sym) * x^2 + omega_A * y
  F_y_sym <- function(x, y) (omega_B + theta_sym) * y^2 + omega_B * x
  nullcline_x_sym <- function(x) -((omega_A + theta_sym)/omega_A) * x^2
  nullcline_y_sym <- function(y) -((omega_B + theta_sym)/omega_B) * y^2
  
  # Grid
  x <- seq(x_min, x_max, length.out = 15)
  y <- seq(y_min, y_max, length.out = 15)
  grid <- expand.grid(x = x, y = y)
  
  # Symmetric
  grid_sym <- grid
  grid_sym$F_x <- F_x_sym(grid_sym$x, grid_sym$y)
  grid_sym$F_y <- F_y_sym(grid_sym$x, grid_sym$y)
  
  # Asymmetric
  grid_asym <- grid
  grid_asym$F_x <- F_x(grid_asym$x, grid_asym$y)
  grid_asym$F_y <- F_y(grid_asym$x, grid_asym$y)
  
  # Curves
  x_curve <- seq(x_min, x_max, length.out = 200)
  y_curve1_sym <- nullcline_x_sym(x_curve)
  y_curve2_sym <- seq(y_min, y_max, length.out = 200)
  x_curve2_sym <- nullcline_y_sym(y_curve2_sym)
  
  y_curve1_asym <- nullcline_x(x_curve)
  x_curve2_asym <- nullcline_y(y_curve2_sym)
  
  # Equilibrium points (approximate for these cases)
  eq_sym <- data.frame(x = c(0, -0.408), y = c(0, -1))
  eq_asym <- data.frame(x = c(0, -0.5), y = c(0, -1.5))
  
  p1 <- ggplot(grid_sym, aes(x = x, y = y)) +
    geom_segment(aes(xend = x + F_x/8, yend = y + F_y/8),
                 arrow = arrow(length = unit(0.08, "cm")), alpha = 0.7) +
    geom_line(data = data.frame(x = x_curve, y = y_curve1_sym),
              aes(x = x, y = y), color = "blue", linetype = "dashed", size = 0.8) +
    geom_line(data = data.frame(x = x_curve2_sym, y = y_curve2_sym),
              aes(x = x, y = y), color = "green", linetype = "dashed", size = 0.8) +
    geom_point(data = eq_sym, aes(x = x, y = y), color = "red", size = 3) +
    labs(title = "Symmetric Case: θ_A = θ_B = 5",
         subtitle = "(K = 0 regime)",
         x = "Strategy x (Player A)", y = "Strategy y (Player B)") +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max))
  
  p2 <- ggplot(grid_asym, aes(x = x, y = y)) +
    geom_segment(aes(xend = x + F_x/8, yend = y + F_y/8),
                 arrow = arrow(length = unit(0.08, "cm")), alpha = 0.7) +
    geom_line(data = data.frame(x = x_curve, y = y_curve1_asym),
              aes(x = x, y = y), color = "blue", linetype = "dashed", size = 0.8) +
    geom_line(data = data.frame(x = x_curve2_sym, y = y_curve2_sym),
              aes(x = x, y = y), color = "green", linetype = "dashed", size = 0.8) +
    geom_point(data = eq_asym, aes(x = x, y = y), color = "red", size = 3) +
    labs(title = paste0("Asymmetric Case: θ_A = ", theta_A, ", θ_B = ", theta_B),
         subtitle = paste0("K = ", round(curvature(x_fixed, y_fixed), 4), " regime"),
         x = "Strategy x (Player A)", y = "Strategy y (Player B)") +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max))
  
  grid.arrange(p1, p2, ncol = 2)
}

# -----------------------------------------------------------------------------
# Execute - The plots will be generated with the parameters set above
# -----------------------------------------------------------------------------
print("Generating VNAE plots...")
print(paste("Current parameters: θ_A =", theta_A, "θ_B =", theta_B, "β =", beta))
print(paste("Fixed point for curvature: (", x_fixed, ",", y_fixed, ")"))

# Create and display the plots
plot1 <- create_field_plot()
plot2 <- create_curvature_plot()
plot3 <- create_phase_portrait()

# Display the plots
print(plot1)
print(plot2)
print(plot3)
create_comparison_plot()

# Optional: save the plots to files
# ggsave("VNAE_Field_Manifold.png", plot1, width = 8, height = 6, dpi = 300)
# ggsave("VNAE_Curvature_Analysis.png", plot2, width = 8, height = 6, dpi = 300)
# ggsave("VNAE_Phase_Portrait.png", plot3, width = 8, height = 6, dpi = 300)
# png("VNAE_Comparison.png", width = 12, height = 6, units = "in", res = 300)
# create_comparison_plot()
# dev.off()
