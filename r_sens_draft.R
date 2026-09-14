# ============================================================
# r x H SENSITIVITY SURFACE
#
# How does the 200-year population outcome depend jointly on:
#
#   1. population productivity (r)
#   2. sustained annual anthropogenic mortality (H)
#
# r <= 0 is not explored because, with H > 0, the model has
# no positive sustainable equilibrium under those conditions.
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(metR)


# ============================================================
# PARAMETERS
# ============================================================

N0 <- 4973
K  <- 6000

years <- 0:200


# ============================================================
# r SENSITIVITY RANGE
# ============================================================

r_values <- seq(
  0,
  0.12,
  by = 0.001
)


# ============================================================
# ANTHROPOGENIC MORTALITY RANGE
# ============================================================

H_grid <- seq(
  20,
  200,
  by = 2
)


# ============================================================
# MODEL FUNCTION
# ============================================================

run_r_sensitivity <- function(
    r_value,
    H,
    N0 = 4973,
    K = 6000,
    years = 0:200
) {
  
  N <- numeric(
    length(years)
  )
  
  N[1] <- N0
  
  
  for (t in 2:length(years)) {
    
    growth <-
      r_value *
      N[t - 1] *
      (
        1 -
          N[t - 1] / K
      )
    
    
    N[t] <-
      max(
        0,
        
        N[t - 1] +
          growth -
          H
      )
  }
  
  
  data.frame(
    
    r =
      r_value,
    
    H =
      H,
    
    final_population =
      last(N)
  )
}


# ============================================================
# RUN ALL r x H COMBINATIONS
# ============================================================

surface_results <- expand_grid(
  
  r =
    r_values,
  
  H =
    H_grid
  
) %>%
  
  rowwise() %>%
  
  do(
    
    run_r_sensitivity(
      
      r_value =
        .$r,
      
      H =
        .$H,
      
      N0 =
        N0,
      
      K =
        K,
      
      years =
        years
    )
  ) %>%
  
  ungroup()


# ============================================================
# POPULATION ISO-LINES
# ============================================================

iso_breaks <- c(
  500,
  1000,
  2000,
  3000,
  4000,
  4900
)


# ============================================================
# FIGURE
# ============================================================

p_r_surface <- ggplot(
  
  surface_results,
  
  aes(
    x = r,
    y = H,
    fill = final_population
  )
  
) +
  
  
  # ----------------------------------------------------------
# POPULATION SURFACE
# ----------------------------------------------------------

geom_tile() +
  
  
  # ----------------------------------------------------------
# POPULATION ISO-LINES
# ----------------------------------------------------------

geom_contour(
  
  aes(
    z =
      final_population
  ),
  
  breaks =
    iso_breaks,
  
  color =
    "grey20",
  
  linewidth =
    0.55
) +
  
  
  # ----------------------------------------------------------
# LABEL ALL POPULATION ISO-LINES
# ----------------------------------------------------------

metR::geom_text_contour(
  
  aes(
    z =
      final_population
  ),
  
  breaks =
    iso_breaks,
  
  stroke =
    0.2,
  
  size =
    3.6,
  
  color =
    "grey15",
  
  label.placer =
    metR::label_placer_fraction(
      frac =
        0.65
    )
) +
  
  
  # ----------------------------------------------------------
# FILL COLORS
#
# Same general muted red -> neutral -> blue logic
# as Figure 2.
# ----------------------------------------------------------

scale_fill_gradient2(
  
  low =
    "#B75C5C",
  
  mid =
    "#C6C1C7",
  
  high =
    "#5E83A3",
  
  midpoint =
    3000,
  
  limits =
    c(
      0,
      6000
    ),
  
  breaks =
    c(
      0,
      1000,
      2000,
      3000,
      4000,
      5000,
      6000
    ),
  
  labels =
    comma,
  
  name =
    "Population\nafter 200 years"
) +
  
  
  # ----------------------------------------------------------
# X AXIS
# ----------------------------------------------------------

scale_x_continuous(
  
  limits =
    c(
      0,
      0.12
    ),
  
  breaks =
    seq(
      0,
      0.12,
      by = 0.02
    ),
  
  labels =
    number_format(
      accuracy = 0.01
    ),
  
  expand =
    expansion(
      mult =
        c(
          0,
          0
        )
    )
) +
  
  
  # ----------------------------------------------------------
# Y AXIS
# ----------------------------------------------------------

scale_y_continuous(
  
  limits =
    c(
      20,
      200
    ),
  
  breaks =
    seq(
      20,
      200,
      by = 20
    ),
  
  expand =
    expansion(
      mult =
        c(
          0,
          0
        )
    )
) +
  
  
  # ----------------------------------------------------------
# LABELS
# ----------------------------------------------------------

labs(
  
  title =
    "Population persistence depends jointly on productivity and anthropogenic mortality",
  
  subtitle =
    paste0(
      "Population after 200 years; ",
      "N₀ = ",
      comma(N0),
      "; K = ",
      comma(K)
    ),
  
  x =
    "Annual population growth parameter (r)",
  
  y =
    "Annual anthropogenic mortality (whales/year)"
) +
  
  
  # ----------------------------------------------------------
# THEME
# ----------------------------------------------------------

theme_classic(
  base_size =
    13
) +
  
  theme(
    
    plot.title =
      element_text(
        face =
          "bold",
        size =
          16
      ),
    
    plot.subtitle =
      element_text(
        size =
          10.5
      ),
    
    legend.title =
      element_text(
        face =
          "bold"
      ),
    
    legend.position =
      "right"
  )


# ============================================================
# PRINT FIGURE
# ============================================================

p_r_surface

