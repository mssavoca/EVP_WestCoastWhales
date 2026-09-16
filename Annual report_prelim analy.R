# Humpback whale population models for EVP mid-year

library(ggplot2)
library(dplyr)
library(scales)

# ============================================================
# PARAMETERS
# ============================================================

N0 <- 4973
K <- 6000
r <- 0.118

years <- 0:200

# Constant annual anthropogenic mortality scenarios
takes <- seq(20, 200, by = 10)

# Replacement at starting population
replacement_N0 <- r * N0 * (1 - N0 / K)

# Maximum replacement under logistic growth
max_replacement <- r * K / 4
N_max_replacement <- K / 2


# ============================================================
# FUNCTION TO RUN ONE POPULATION TRAJECTORY
# ============================================================

run_model <- function(H, N0 = 4973, K = 6000, r = 0.119,
                      years = 0:200) {
  
  N <- numeric(length(years))
  N[1] <- N0
  
  for (t in 2:length(years)) {
    
    # Density-dependent population replacement
    growth <- r * N[t - 1] * (1 - N[t - 1] / K)
    
    # Population next year
    N[t] <- max(0, N[t - 1] + growth - H)
  }
  
  data.frame(
    year = years,
    population = N,
    take = H
  )
}


# ============================================================
# RUN ALL TAKE SCENARIOS
# ============================================================

results <- bind_rows(
  lapply(takes, run_model)
) %>%
  group_by(take) %>%
  mutate(
    final_population = last(population),
    trajectory = if_else(
      final_population >= N0,
      "Population ultimately grows",
      "Population ultimately declines"
    )
  ) %>%
  ungroup()


# Check final populations
results %>%
  group_by(take, trajectory) %>%
  summarise(
    final_population = last(population),
    .groups = "drop"
  )





# FIGURE 1; scenarios of constant anthropogenic take----

# Blue: trajectories ending above N0
# Red: trajectories ending below N0

blue_cols <- colorRampPalette(
  c("#08306B", "#9ECAE1")
)(length(unique(
  results$take[results$trajectory == "Population ultimately grows"]
)))

red_cols <- colorRampPalette(
  c("#FCAE91", "#99000D")
)(length(unique(
  results$take[results$trajectory == "Population ultimately declines"]
)))

# Give each take scenario its own color
take_colors <- c()

blue_takes <- takes[takes <= replacement_N0]
red_takes  <- takes[takes > replacement_N0]

take_colors <- c(
  setNames(blue_cols, blue_takes),
  setNames(red_cols, red_takes)
)


p1 <- ggplot(
  results,
  aes(
    x = year,
    y = population,
    group = take,
    color = factor(take)
  )
) +
  
  geom_line(linewidth = 0.8) +
  
  # Starting population
  geom_hline(
    yintercept = N0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey35"
  ) +
  
  # Carrying capacity
  geom_hline(
    yintercept = K,
    linetype = "dotted",
    linewidth = 0.5,
    color = "grey35"
  ) +
  
  
  scale_color_manual(
    values = take_colors,
    name = "Annual combined\ntake (whales/year)"
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 200, by = 20)
  ) +
  
  guides(
    color = guide_legend(
      ncol = 2,
      byrow = FALSE
    )
  ) +
  
  labs(
    title = "Population response to constant annual anthropogenic mortality",
    subtitle = paste0(
      "K = ", comma(K),
      "; r = ", r,
      "; N₀ = ", comma(N0),
      "; simulation = 200 years"
    ),
    x = "Years",
    y = "Humpback whale population"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )


p1



# ============================================================
# TIME TO QUASI-EXTINCTION
# ============================================================

quasi_extinction <- 500

# Take levels focused around and above the maximum
# replacement rate of 178.5 whales/year
take_levels <- seq(150, 250, by = 1)

time_to_qe <- data.frame(
  take = take_levels,
  years_to_quasi_extinction = NA
)

for (i in seq_along(take_levels)) {
  
  H <- take_levels[i]
  
  N <- N0
  
  for (t in 1:1000) {
    
    growth <- r * N * (1 - N / K)
    
    N <- max(0, N + growth - H)
    
    if (N <= quasi_extinction) {
      time_to_qe$years_to_quasi_extinction[i] <- t
      break
    }
  }
}






# TAKE x DURATION ANALYSIS ----


take_levels <- seq(100, 250, by = 5)
duration_levels <- seq(1, 150, by = 1)

duration_results <- expand.grid(
  take = take_levels,
  duration = duration_levels
)

duration_results$population <- NA

for (i in seq_len(nrow(duration_results))) {
  
  H <- duration_results$take[i]
  duration <- duration_results$duration[i]
  
  N <- N0
  
  for (t in 1:duration) {
    
    growth <- r * N * (1 - N / K)
    
    N <- max(0, N + growth - H)
  }
  
  duration_results$population[i] <- N
}





# FIGURE 2: population response to take and duration ----

library(metR)


contour_breaks <- c(500, 1000, 2000, 3000, 4000, 4900)

ggplot(
  duration_results,
  aes(
    x = duration,
    y = take,
    fill = population
  )
) +
  
  geom_tile() +
  
  # Maximum possible replacement
  geom_hline(
    yintercept = max_replacement,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey15"
  ) +
  
  # Population iso-lines
  geom_contour(
    aes(z = population),
    breaks = contour_breaks,
    color = "grey20",
    linewidth = 0.6
  ) +
  
  # ----------------------------------------------------------
# POPULATION ISO-LINE LABELS
# ----------------------------------------------------------

# 500
annotate(
  "label",
  x = 50,
  y = 248,
  label = "500",
  size = 3.5,
  label.size = 0.25,
  fill = alpha("white", 0.80)
) +
  
  # 1,000
  annotate(
    "label",
    x = 49,
    y = 239,
    label = "1,000",
    size = 3.5,
    label.size = 0.25,
    fill = alpha("white", 0.80)
  ) +
  
  # 2,000
  annotate(
    "label",
    x = 38,
    y = 237,
    label = "2,000",
    size = 3.5,
    label.size = 0.25,
    fill = alpha("white", 0.80)
  ) +
  
  # 3,000
  annotate(
    "label",
    x = 27,
    y = 230,
    label = "3,000",
    size = 3.5,
    label.size = 0.25,
    fill = alpha("white", 0.80)
  ) +
  
  # 4,000
  annotate(
    "label",
    x = 15,
    y = 215,
    label = "4,000",
    size = 3.5,
    label.size = 0.25,
    fill = alpha("white", 0.80)
  ) +
  
  # 4,900
  annotate(
    "label",
    x = 5,
    y = 155,
    label = "4,900",
    size = 3.5,
    label.size = 0.25,
    fill = alpha("white", 0.80)
  ) +
  
  # ----------------------------------------------------------
# PERSISTENCE / COLLAPSE ANNOTATIONS
# ----------------------------------------------------------

annotate(
  "text",
  x = 40,
  y = 135,
  label = "PERSISTENCE",
  fontface = "bold",
  size = 5,
  color = "grey15"
) +
  
  annotate(
    "text",
    x = 110,
    y = 225,
    label = "COLLAPSE /\nEXTINCTION",
    fontface = "bold",
    size = 5,
    color = "grey15"
  ) +
  
  # ----------------------------------------------------------
# COLOR SCALE
# ----------------------------------------------------------

scale_fill_gradient2(
  low = "#B45F5F",
  mid = "#B9A8AE",
  high = "#557A9B",
  midpoint = 2500,
  limits = c(0, N0),
  oob = scales::squish,
  name = "Population",
  labels = comma
) +
  
  scale_x_continuous(
    breaks = seq(0, 150, by = 25)
  ) +
  
  scale_y_continuous(
    breaks = seq(100, 250, by = 25)
  ) +
  
  labs(
    title = "Population response depends on both take level and duration",
    
    subtitle = paste0(
      "Starting population = ",
      comma(N0),
      "; K = ",
      comma(K),
      "; r = ",
      r
    ),
    
    x = "Years of elevated mortality",
    y = "Annual anthropogenic mortality (whales/year)"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    legend.title = element_text(
      face = "bold"
    )
  )





# FIGURE 3: effect of reduced reproductive rate (r) on population trajectories ----
# SENSITIVITY OF POPULATION OUTCOME TO PRODUCTIVITY (r)
# AND ANTHROPOGENIC MORTALITY (H)
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)


# ============================================================
# PARAMETERS
# ============================================================

N0 <- 4973
K  <- 6000

years <- 0:200


# ============================================================
# PRODUCTIVITY RANGE
#
# r <= 0 is not explored because, with positive H,
# the model has no positive sustainable equilibrium.
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
    
    # Density-dependent population growth
    growth <-
      r_value *
      N[t - 1] *
      (
        1 -
          N[t - 1] / K
      )
    
    
    # Population in following year
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

surface_results <-
  
  expand_grid(
    
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
# CLEAN SURFACE DATA
#
# Ensure exactly one population value for every r x H
# combination.
# ============================================================

surface_plot_data <-
  
  surface_results %>%
  
  group_by(
    r,
    H
  ) %>%
  
  summarise(
    
    final_population =
      mean(
        final_population,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) %>%
  
  arrange(
    H,
    r
  )


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
# ISO-LINE LABEL LOCATIONS
#
# These are deliberately positioned similarly to Figure 2.
# ============================================================

label_positions <- data.frame(
  
  population =
    iso_breaks,
  
  label_r = c(
    0.035,   # 500
    0.085,   # 1,000
    0.045,   # 2,000
    0.116,   # 3,000
    0.115,   # 4,000
    0.115    # 4,900
  ),
  
  label_H = c(
    80,     # 500
    160,     # 1,000
    95,      # 2,000
    165,     # 3,000
    135,     # 4,000
    80       # 4,900
  )
)


# ============================================================
# FIND APPROXIMATE POINT ON EACH ISO-LINE
#
# These points are used only as the starting locations for
# the short connector lines.
# ============================================================

connector_points <-
  
  bind_rows(
    
    lapply(
      
      seq_along(
        iso_breaks
      ),
      
      function(i) {
        
        this_population <-
          iso_breaks[i]
        
        this_r <-
          label_positions$label_r[i]
        
        
        surface_plot_data %>%
          
          mutate(
            
            # How close is this grid cell to the desired
            # population iso-line?
            population_distance =
              abs(
                final_population -
                  this_population
              ),
            
            
            # Prefer a point near the desired label r value
            r_distance =
              abs(
                r -
                  this_r
              )
          ) %>%
          
          arrange(
            
            population_distance,
            
            r_distance
          ) %>%
          
          slice(1) %>%
          
          mutate(
            
            population =
              this_population,
            
            label_r =
              label_positions$label_r[i],
            
            label_H =
              label_positions$label_H[i]
          )
      }
    )
  )


# ============================================================
# FIGURE 3
# ============================================================

p_r_surface <-
  
  ggplot(
    
    surface_plot_data,
    
    aes(
      
      x =
        r,
      
      y =
        H,
      
      fill =
        final_population
    )
    
  ) +
  
  
  # ==========================================================
# POPULATION SURFACE
#
# alpha = 1 ensures fully opaque colors.
# ==========================================================

geom_tile(
  alpha = 1
) +
  
  
  # ==========================================================
# POPULATION ISO-LINES
# ==========================================================

geom_contour(
  
  aes(
    
    x =
      r,
    
    y =
      H,
    
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
  
  
  # ==========================================================
# CONNECTOR LINES FROM LABELS TO ISO-LINES
# ==========================================================

geom_segment(
  
  data =
    connector_points,
  
  aes(
    
    x =
      r,
    
    y =
      H,
    
    xend =
      label_r,
    
    yend =
      label_H
  ),
  
  inherit.aes =
    FALSE,
  
  color =
    "grey20",
  
  linewidth =
    0.35
) +
  
  
  # ==========================================================
# POPULATION ISO-LINE LABELS
# ==========================================================

geom_label(
  
  data =
    label_positions,
  
  aes(
    
    x =
      label_r,
    
    y =
      label_H,
    
    label =
      comma(
        population
      )
  ),
  
  inherit.aes =
    FALSE,
  
  size =
    3.5,
  
  label.size =
    0.25,
  
  fill =
    alpha(
      "white",
      0.90
    ),
  
  color =
    "grey15",
  
  label.padding =
    unit(
      0.12,
      "lines"
    )
) +
  
  
  # ==========================================================
# COLLAPSE / EXTINCTION
#
# Upper-left red region
# ==========================================================

annotate(
  
  "text",
  
  x =
    0.012,
  
  y =
    175,
  
  label =
    "COLLAPSE /\nEXTINCTION",
  
  fontface =
    "bold",
  
  size =
    5,
  
  color =
    "grey15",
  
  hjust =
    0
) +
  
  
  # ==========================================================
# PERSISTENCE
#
# Lower-right blue region
# ==========================================================

annotate(
  
  "text",
  
  x =
    0.100,
  
  y =
    53,
  
  label =
    "PERSISTENCE",
  
  fontface =
    "bold",
  
  size =
    5,
  
  color =
    "grey15",
  
  hjust =
    0.5
) +
  
  
  # ==========================================================
# COLOR SCALE
#
# IMPORTANT:
# This now uses EXACTLY the same color logic as Figure 2.
#
# Values at or above N0 are squished to the darkest blue.
# ==========================================================

scale_fill_gradient2(
  
  low =
    "#B45F5F",
  
  mid =
    "#B9A8AE",
  
  high =
    "#557A9B",
  
  midpoint =
    2500,
  
  limits =
    c(
      0,
      N0
    ),
  
  oob =
    scales::squish,
  
  breaks =
    c(
      0,
      1000,
      2000,
      3000,
      4000,
      N0
    ),
  
  labels =
    c(
      "0",
      "1,000",
      "2,000",
      "3,000",
      "4,000",
      "4,973"
    ),
  
  name =
    "Population\nafter 200 years"
) +
  
  
  # ==========================================================
# X AXIS
# ==========================================================

scale_x_continuous(
  
  breaks =
    seq(
      0,
      0.12,
      by =
        0.02
    ),
  
  labels =
    number_format(
      accuracy =
        0.01
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
  
  
  # ==========================================================
# Y AXIS
# ==========================================================

scale_y_continuous(
  
  breaks =
    seq(
      20,
      200,
      by =
        20
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
  
  
  # ==========================================================
# PLOT WINDOW
# ==========================================================

coord_cartesian(
  
  xlim =
    c(
      0,
      0.12
    ),
  
  ylim =
    c(
      20,
      200
    ),
  
  expand =
    FALSE
) +
  
  
  # ==========================================================
# TITLES AND AXIS LABELS
# ==========================================================

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
  
  
  # ==========================================================
# THEME
# ==========================================================

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




# Fig.4 STOCHASTIC GOOD/BAD YEAR Simulation----
# Environmental productivity + anthropogenic mortality


library(ggplot2)
library(dplyr)
library(scales)


# PARAMETERS


N0 <- 4973
K  <- 6000

years <- 0:200
n_sims <- 10000


# GOOD YEAR
good_r_min <- 0.08
good_r_max <- 0.12

good_H_min <- 20
good_H_max <- 80



# BAD YEAR
bad_r_min <- 0.00
bad_r_max <- 0.05

bad_H_min <- 100
bad_H_max <- 200


# Frequencies of bad years to illustrate
bad_frequencies <- c(
  0.10,
  0.40,
  0.60
)

# Quasi-extinction threshold
quasi_extinction <- 500



# FUNCTION TO RUN STOCHASTIC MODEL

run_stochastic_model <- function(
    bad_frequency,
    N0 = 4973,
    K = 6000,
    years = 0:200,
    n_sims = 10000
) {
  
  results <- vector("list", n_sims)
  
  for (s in seq_len(n_sims)) {
    
    N <- numeric(length(years))
    N[1] <- N0
    
    # --------------------------------------------------------
    # Draw whether each year is good or bad
    # --------------------------------------------------------
    
    bad_year <- runif(
      length(years) - 1
    ) < bad_frequency
    
    
    # --------------------------------------------------------
    # Population dynamics
    # --------------------------------------------------------
    
    for (t in 2:length(years)) {
      
      if (bad_year[t - 1]) {
        
        # BAD YEAR:
        # low productivity + high anthropogenic mortality
        
        r <- runif(
          1,
          min = bad_r_min,
          max = bad_r_max
        )
        
        H <- runif(
          1,
          min = bad_H_min,
          max = bad_H_max
        )
        
      } else {
        
        # GOOD YEAR:
        # higher productivity + lower anthropogenic mortality
        
        r <- runif(
          1,
          min = good_r_min,
          max = good_r_max
        )
        
        H <- runif(
          1,
          min = good_H_min,
          max = good_H_max
        )
      }
      
      
      # Density-dependent population growth
      growth <- r *
        N[t - 1] *
        (1 - N[t - 1] / K)
      
      
      # Population update
      N[t] <- max(
        0,
        N[t - 1] + growth - H
      )
    }
    
    
    results[[s]] <- data.frame(
      year = years,
      population = N,
      simulation = s,
      bad_frequency = bad_frequency
    )
  }
  
  bind_rows(results)
}



# RUN STOCHASTIC SCENARIOS

set.seed(42)

stochastic_results <- bind_rows(
  
  lapply(
    bad_frequencies,
    
    function(p_bad) {
      
      run_stochastic_model(
        bad_frequency = p_bad,
        N0 = N0,
        K = K,
        years = years,
        n_sims = n_sims
      )
    }
  )
)



# LABEL SCENARIOS

stochastic_results <- stochastic_results %>%
  mutate(
    
    bad_frequency_label = factor(
      
      paste0(
        round(bad_frequency * 100),
        "% bad years"
      ),
      
      levels = c(
        "10% bad years",
        "40% bad years",
        "60% bad years"
      )
    )
  )



# SUMMARY STATISTICS THROUGH TIME

summary_results <- stochastic_results %>%
  group_by(
    bad_frequency_label,
    year
  ) %>%
  summarise(
    
    q10 = quantile(
      population,
      0.10
    ),
    
    median = median(
      population
    ),
    
    q90 = quantile(
      population,
      0.90
    ),
    
    .groups = "drop"
  )



# SUMMARY OF FINAL OUTCOMES

final_results <- stochastic_results %>%
  group_by(
    bad_frequency_label,
    simulation
  ) %>%
  summarise(
    
    final_population =
      last(population),
    
    minimum_population =
      min(population),
    
    quasi_extinction =
      any(
        population <= quasi_extinction
      ),
    
    .groups = "drop"
  )


scenario_summary <- final_results %>%
  group_by(
    bad_frequency_label
  ) %>%
  summarise(
    
    median_final_population =
      median(final_population),
    
    q10_final_population =
      quantile(
        final_population,
        0.10
      ),
    
    q90_final_population =
      quantile(
        final_population,
        0.90
      ),
    
    proportion_below_N0 =
      mean(
        final_population < N0
      ),
    
    probability_quasi_extinction =
      mean(
        quasi_extinction
      ),
    
    .groups = "drop"
  )

print(scenario_summary)



# FIGURE code----
# INDIVIDUAL STOCHASTIC TRAJECTORIES + MEDIAN

scenario_colors <- c(
  "10% bad years" = "#2F5D8A",
  "40% bad years" = "#8C6D31",
  "60% bad years" = "#A52A2A"
)


p_stochastic <- ggplot() +
  
# ----------------------------------------------------------
# Individual simulation trajectories

geom_line(
  data = stochastic_results,
  
  aes(
    x = year,
    y = population,
    group = interaction(
      bad_frequency_label,
      simulation
    )
  ),
  
  color = "grey40",
  alpha = 0.015,
  linewidth = 0.25
) +
  

# 10th–90th percentile envelope
geom_ribbon(
  data = summary_results,
  
  aes(
    x = year,
    ymin = q10,
    ymax = q90,
    group = bad_frequency_label
  ),
  
  fill = "grey50",
  alpha = 0.12
) +
  
  # ----------------------------------------------------------
# Median trajectory
# ----------------------------------------------------------

geom_line(
  data = summary_results,
  
  aes(
    x = year,
    y = median,
    color = bad_frequency_label
  ),
  
  linewidth = 1.3
) +
  
  # Starting population
  geom_hline(
    yintercept = N0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey35"
  ) +
  
  # Quasi-extinction threshold
  geom_hline(
    yintercept = quasi_extinction,
    linetype = "dotted",
    linewidth = 0.6,
    color = "grey45"
  ) +
  
  facet_wrap(
    ~ bad_frequency_label,
    ncol = 1,
    axes = "all_x",
    axis.labels = "all_x"
  ) +
  
  scale_color_manual(
    values = scenario_colors,
    guide = "none"
  ) +
  
  scale_x_continuous(
    limits = c(0, 200),
    breaks = seq(
      0,
      200,
      by = 20
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 6200),
    breaks = seq(
      0,
      6000,
      by = 2000
    ),
    labels = comma
  ) +
  
  labs(
    title = "Population response to environmental and anthropogenic stochasticity",
    
    subtitle = paste0(
      "Good year: r = 0.08–0.12, H = 20–80; ",
      "bad year: r = 0–0.05, H = 100–200"
    ),
    
    x = "Years",
    y = "Humpback whale population"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 10.5
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    )
  )


p_stochastic




# FIGURE 5: QUASI-EXTINCTION RISK CURVE----


# Colors match the 10%, 40%, and 60% scenarios
scenario_point_colors <- c(
  "0.1" = "#2F5D8A",   # 10% bad years - blue
  "0.4" = "#8C6D31",   # 40% bad years - ochre
  "0.6" = "#A52A2A"    # 60% bad years - red
)


# Pull out the three focal scenarios
risk_points <- risk_results %>%
  filter(
    bad_frequency %in% c(0.10, 0.40, 0.60)
  ) %>%
  mutate(
    scenario = factor(
      as.character(bad_frequency),
      levels = c("0.1", "0.4", "0.6")
    )
  )


# Plot
p_risk <- ggplot(
  risk_results,
  aes(
    x = bad_frequency,
    y = probability_quasi_extinction
  )
) +
  
  # Overall stochastic risk curve
  geom_line(
    linewidth = 1.1,
    color = "black"
  ) +
  
  # Highlight the three scenarios used in Fig. 4
  geom_point(
    data = risk_points,
    aes(
      color = scenario
    ),
    size = 3.5
  ) +
  
  scale_color_manual(
    values = scenario_point_colors,
    labels = c(
      "10% bad years",
      "40% bad years",
      "60% bad years"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    limits = c(0, 0.8),
    breaks = seq(
      0,
      0.8,
      by = 0.1
    ),
    labels = percent
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.1
    ),
    labels = percent
  ) +
  
  labs(
    title =
      "Quasi-extinction risk increases with the frequency of bad years",
    
    subtitle =
      "Quasi-extinction defined as N ≤ 500",
    
    x =
      "Frequency of bad years",
    
    y =
      "Probability of quasi-extinction within 200 years"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.text = element_text(
      size = 10
    )
  )


p_risk






# FIGURE 6 : EVENT-CENTERED RESPONSE TO A SINGLE OIL SPILL----
#
# Each simulation experiences:
#   Good years:
#       r = 0.08-0.12
#       H = 20-80
#
#   Bad years:
#       r = 0-0.05
#       H = 100-200
#
# Oil spill in relative year 0:
#       additional H = 500
#       r = 0 for one year
#
# Each simulation has 30 years of stochastic history BEFORE
# the spill and 100 years AFTER the spill.

library(ggplot2)
library(dplyr)
library(scales)



# PARAMETERS

N0 <- 4973
K  <- 6000


# Good years
good_r <- c(
  0.08,
  0.12
)

good_H <- c(
  20,
  80
)



# Bad years
bad_r <- c(
  0.00,
  0.05
)

bad_H <- c(
  100,
  200
)



# Oil spill
spill_H <- 1000
spill_r <- 0



# Event-centered time window
pre_years <- 30
post_years <- 5000

relative_years <- (
  -pre_years:
    post_years
)



# Number of simulations
n_sims <- 10000



# Bad-year scenarios
bad_frequencies <- c(
  0.10,
  0.40,
  0.60
)



# Quasi-extinction threshold
quasi_extinction <- 500


# ============================================================
# FUNCTION TO RUN EVENT-CENTERED STOCHASTIC MODEL
# ============================================================

run_event_centered_model <- function(
    bad_frequency,
    N0 = 4973,
    K = 6000,
    pre_years = 30,
    post_years = 100,
    n_sims = 10000
) {
  
  # ----------------------------------------------------------
  # Relative time scale
  # ----------------------------------------------------------
  
  relative_years <- c(
    -pre_years:post_years
  )
  
  n_years <- length(
    relative_years
  )
  
  
  # ----------------------------------------------------------
  # Storage matrices
  # ----------------------------------------------------------
  
  population_spill <- matrix(
    NA_real_,
    nrow = n_sims,
    ncol = n_years
  )
  
  population_no_spill <- matrix(
    NA_real_,
    nrow = n_sims,
    ncol = n_years
  )
  
  
  # ----------------------------------------------------------
  # All simulations begin 30 years before the spill at N0
  # ----------------------------------------------------------
  
  population_spill[, 1] <- N0
  population_no_spill[, 1] <- N0
  
  
  # ----------------------------------------------------------
  # Simulate through time
  # ----------------------------------------------------------
  
  for (t in 2:n_years) {
    
    # --------------------------------------------------------
    # Determine good vs. bad year
    # --------------------------------------------------------
    
    is_bad <- runif(
      n_sims
    ) < bad_frequency
    
    
    # --------------------------------------------------------
    # Draw r
    # --------------------------------------------------------
    
    r_year <- ifelse(
      
      is_bad,
      
      runif(
        n_sims,
        min = bad_r[1],
        max = bad_r[2]
      ),
      
      runif(
        n_sims,
        min = good_r[1],
        max = good_r[2]
      )
    )
    
    
    # --------------------------------------------------------
    # Draw H
    # --------------------------------------------------------
    
    H_year <- ifelse(
      
      is_bad,
      
      runif(
        n_sims,
        min = bad_H[1],
        max = bad_H[2]
      ),
      
      runif(
        n_sims,
        min = good_H[1],
        max = good_H[2]
      )
    )
    
    
    # ========================================================
    # NO-SPILL SCENARIO
    # ========================================================
    
    N_previous <-
      population_no_spill[, t - 1]
    
    growth_no_spill <-
      
      r_year *
      N_previous *
      (
        1 -
          N_previous / K
      )
    
    population_no_spill[, t] <-
      
      pmax(
        0,
        N_previous +
          growth_no_spill -
          H_year
      )
    
    
    # ========================================================
    # SPILL SCENARIO
    # ========================================================
    
    r_spill <- r_year
    H_spill <- H_year
    
    
    # Current relative year
    current_year <-
      relative_years[t]
    
    
    # Apply oil spill ONLY at relative year 0
    if (current_year == 0) {
      
      r_spill <- spill_r
      
      H_spill <-
        H_year +
        spill_H
    }
    
    
    N_previous <-
      population_spill[, t - 1]
    
    growth_spill <-
      
      r_spill *
      N_previous *
      (
        1 -
          N_previous / K
      )
    
    population_spill[, t] <-
      
      pmax(
        0,
        N_previous +
          growth_spill -
          H_spill
      )
  }
  
  
  # ==========================================================
  # CONVERT TO LONG FORMAT
  # ==========================================================
  
  # Give each column an explicit relative-year name
  colnames(population_spill) <-
    paste0(
      "yr_",
      relative_years
    )
  
  colnames(population_no_spill) <-
    paste0(
      "yr_",
      relative_years
    )
  
  
  # ----------------------------------------------------------
  # Spill trajectories
  # ----------------------------------------------------------
  
  spill_long <-
    
    as.data.frame(
      population_spill
    ) %>%
    
    mutate(
      simulation =
        row_number(),
      
      scenario =
        "Oil spill"
    ) %>%
    
    tidyr::pivot_longer(
      
      cols =
        starts_with("yr_"),
      
      names_to =
        "time",
      
      values_to =
        "population"
    ) %>%
    
    mutate(
      
      relative_year =
        as.numeric(
          sub(
            "yr_",
            "",
            time
          )
        )
    )
  
  
  # ----------------------------------------------------------
  # No-spill trajectories
  # ----------------------------------------------------------
  
  no_spill_long <-
    
    as.data.frame(
      population_no_spill
    ) %>%
    
    mutate(
      simulation =
        row_number(),
      
      scenario =
        "No spill"
    ) %>%
    
    tidyr::pivot_longer(
      
      cols =
        starts_with("yr_"),
      
      names_to =
        "time",
      
      values_to =
        "population"
    ) %>%
    
    mutate(
      
      relative_year =
        as.numeric(
          sub(
            "yr_",
            "",
            time
          )
        )
    )
  
  
  # ----------------------------------------------------------
  # Combine
  # ----------------------------------------------------------
  
  results <-
    bind_rows(
      spill_long,
      no_spill_long
    )
  
  
  # Add bad-year frequency
  results$bad_frequency <-
    bad_frequency
  
  
  results
}

# ============================================================
# RUN ALL THREE BAD-YEAR SCENARIOS
# ============================================================

set.seed(42)

event_results <-
  
  bind_rows(
    
    lapply(
      
      bad_frequencies,
      
      function(p) {
        
        run_event_centered_model(
          
          bad_frequency =
            p,
          
          N0 =
            N0,
          
          K =
            K,
          
          pre_years =
            pre_years,
          
          post_years =
            post_years,
          
          n_sims =
            n_sims
        )
      }
    )
  )


# ============================================================
# LABEL SCENARIOS
# ============================================================

event_results <-
  
  event_results %>%
  
  mutate(
    
    bad_frequency_label = factor(
      
      paste0(
        round(
          bad_frequency * 100
        ),
        "% bad years"
      ),
      
      levels = c(
        "10% bad years",
        "40% bad years",
        "60% bad years"
      )
    )
  )


# ============================================================
# SUMMARIZE POPULATION TRAJECTORIES
# ============================================================

event_summary <-
  
  event_results %>%
  
  group_by(
    
    bad_frequency_label,
    
    scenario,
    
    relative_year
    
  ) %>%
  
  summarise(
    
    q10 =
      quantile(
        population,
        0.10
      ),
    
    median =
      median(
        population
      ),
    
    q90 =
      quantile(
        population,
        0.90
      ),
    
    .groups = "drop"
  )


# ============================================================
# SPILL VS. NO-SPILL RETURN TIMES
# ============================================================

return_comparison <-
  
  event_results %>%
  
  filter(
    relative_year >= 0
  ) %>%
  
  group_by(
    bad_frequency_label,
    scenario,
    simulation
  ) %>%
  
  summarise(
    
    first_return =
      ifelse(
        any(
          population >= N0
        ),
        min(
          relative_year[
            population >= N0
          ]
        ),
        NA_real_
      ),
    
    .groups =
      "drop"
  ) %>%
  
  tidyr::pivot_wider(
    
    names_from =
      scenario,
    
    values_from =
      first_return
  ) %>%
  
  mutate(
    
    spill_delay =
      `Oil spill` -
      `No spill`
  )


return_comparison

# ============================================================
# COLORS FOR THE THREE BAD-YEAR SCENARIOS
# ============================================================

scenario_colors <- c(
  
  "10% bad years" =
    "#2F5D8A",
  
  "40% bad years" =
    "#8C6D31",
  
  "60% bad years" =
    "#A52A2A"
)


# ============================================================
# FIGURE
# ============================================================

p_event <- ggplot() +
  
  
  # ==========================================================
# INDIVIDUAL STOCHASTIC OIL-SPILL TRAJECTORIES
# ==========================================================

geom_line(
  
  data = event_results %>%
    
    filter(
      scenario ==
        "Oil spill"
    ),
  
  aes(
    
    x =
      relative_year,
    
    y =
      population,
    
    group =
      interaction(
        bad_frequency_label,
        simulation
      )
  ),
  
  color =
    "grey45",
  
  alpha =
    0.012,
  
  linewidth =
    0.25
) +
  
  
  # ==========================================================
# 10-90% ENVELOPE: OIL SPILL
# ==========================================================

geom_ribbon(
  
  data = event_summary %>%
    
    filter(
      scenario ==
        "Oil spill"
    ),
  
  aes(
    
    x =
      relative_year,
    
    ymin =
      q10,
    
    ymax =
      q90,
    
    group =
      bad_frequency_label
  ),
  
  fill =
    "grey55",
  
  alpha =
    0.15
) +
  
  
  # ==========================================================
# MEDIAN NO-SPILL COUNTERFACTUAL
# ==========================================================

geom_line(
  
  data = event_summary %>%
    
    filter(
      scenario ==
        "No spill"
    ),
  
  aes(
    
    x =
      relative_year,
    
    y =
      median,
    
    group =
      bad_frequency_label
  ),
  
  linetype =
    "dashed",
  
  linewidth =
    1.0,
  
  color =
    "grey35"
) +
  
  
  # ==========================================================
# MEDIAN OIL-SPILL TRAJECTORY
# ==========================================================

geom_line(
  
  data = event_summary %>%
    
    filter(
      scenario ==
        "Oil spill"
    ),
  
  aes(
    
    x =
      relative_year,
    
    y =
      median,
    
    color =
      bad_frequency_label,
    
    group =
      bad_frequency_label
  ),
  
  linewidth =
    1.4
) +
  
  
  # ==========================================================
# OIL-SPILL YEAR
# ==========================================================

geom_vline(
  
  xintercept =
    0,
  
  linetype =
    "dashed",
  
  linewidth =
    0.8,
  
  color =
    "black"
) +
  
  
  # ==========================================================
# STARTING POPULATION
# ==========================================================

geom_hline(
  
  yintercept =
    N0,
  
  linetype =
    "dashed",
  
  linewidth =
    0.45,
  
  color =
    "grey50"
) +
  
  
  # ==========================================================
# QUASI-EXTINCTION THRESHOLD
# ==========================================================

geom_hline(
  
  yintercept =
    quasi_extinction,
  
  linetype =
    "dotted",
  
  linewidth =
    0.6,
  
  color =
    "grey45"
) +
  
  
# ==========================================================
# FACETS
# ==========================================================

facet_wrap(
  
  ~ bad_frequency_label,
  
  ncol =
    1,
  
  axes =
    "all_x",
  
  axis.labels =
    "all"
) +
  
  
  # ==========================================================
# COLORS
# ==========================================================

scale_color_manual(
  
  values =
    scenario_colors,
  
  name =
    "Bad-year frequency"
) +
  
  
  # ==========================================================
# X AXIS
# ==========================================================

scale_x_continuous(
  
  limits =
    c(
      -30,
      100
    ),
  
  breaks =
    seq(
      -30,
      100,
      by = 10
    ),
  
  expand =
    expansion(
      mult =
        c(
          0,
          0.01
        )
    )
) +
  
  
  # ==========================================================
# Y AXIS
# ==========================================================

scale_y_continuous(
  
  limits =
    c(
      0,
      6200
    ),
  
  breaks =
    seq(
      0,
      6000,
      by = 2000
    ),
  
  labels =
    comma,
  
  expand =
    expansion(
      mult =
        c(
          0,
          0.02
        )
    )
) +
  
  
# ==========================================================
# LABELS
# ==========================================================

labs(
  
  title = "Population response to a single oil-spill event",
  x = "Year relative to oil spill",
  y = "Humpback whale population"
  ) +
  
  
# ==========================================================
# THEME
# ==========================================================

theme_classic(
  
  base_size =
    12
) +
  
  theme(
    
    plot.title =
      element_text(
        face =
          "bold",
        size =
          14
      ),

    strip.text =
      element_text(
        face =
          "bold",
        size =
          10
      ),
    
    legend.position =
      "none"
  )


# ============================================================
# PRINT FIGURE
# ============================================================

p_event



# ============================================================
# CORRECTED OIL-SPILL RECOVERY TABLE
#
# Two recovery metrics:
#
# 1. Return to N0:
#    First year the spill trajectory reaches N0 = 4,973 whales.
#
# 2. Recovery to matched no-spill trajectory:
#    First year the spill trajectory comes within 1% of the
#    matched no-spill trajectory AND remains within 1% thereafter,
#    provided the matched no-spill population remains above the
#    quasi-extinction threshold of 500 whales.
#
# Recovery times are summarized as:
#    median (25th-75th percentile)
#
# IMPORTANT:
# This requires event_results to have been rerun with a sufficiently
# long post-spill projection (e.g., post_years = 5000).
# ============================================================


library(dplyr)
library(tidyr)
library(scales)


# ============================================================
# PARAMETERS
# ============================================================

starting_population <- N0
quasi_extinction_threshold <- quasi_extinction


# ============================================================
# PUT SPILL AND NO-SPILL TRAJECTORIES ON SAME ROW
# ============================================================

paired_results <-
  
  event_results %>%
  
  filter(
    relative_year >= 0
  ) %>%
  
  select(
    bad_frequency_label,
    simulation,
    relative_year,
    scenario,
    population
  ) %>%
  
  pivot_wider(
    names_from =
      scenario,
    values_from =
      population
  ) %>%
  
  arrange(
    bad_frequency_label,
    simulation,
    relative_year
  )


# ============================================================
# CALCULATE DEFICIT RELATIVE TO MATCHED NO-SPILL TRAJECTORY
# ============================================================

paired_results <-
  
  paired_results %>%
  
  mutate(
    
    proportional_deficit = case_when(
      
      `No spill` > 0 ~
        pmax(
          0,
          (`No spill` - `Oil spill`) /
            `No spill`
        ),
      
      TRUE ~
        NA_real_
    ),
    
    # Spill has effectively recovered to its matched
    # no-spill trajectory.
    within_1pct =
      proportional_deficit <= 0.01 &
      `No spill` > quasi_extinction_threshold
  )


# ============================================================
# FUNCTION:
# FIRST YEAR OF SUSTAINED RECOVERY
#
# We require the population to be within 1% of the counterfactual
# and to remain there for the rest of the projection.
# ============================================================

first_sustained_recovery <- function(
    years,
    recovered
) {
  
  recovered[
    is.na(recovered)
  ] <- FALSE
  
  sustained_recovery <-
    
    rev(
      cummin(
        as.integer(
          rev(
            recovered
          )
        )
      )
    ) == 1
  
  possible_years <-
    
    years[
      sustained_recovery
    ]
  
  if (
    length(
      possible_years
    ) == 0
  ) {
    
    NA_real_
    
  } else {
    
    min(
      possible_years
    )
  }
}


# ============================================================
# CALCULATE RECOVERY TIME FOR EACH INDIVIDUAL SIMULATION
# ============================================================

recovery_by_sim <-
  
  paired_results %>%
  
  group_by(
    bad_frequency_label,
    simulation
  ) %>%
  
  summarise(
    
    # --------------------------------------------------------
    # A. RETURN TO STARTING POPULATION
    # --------------------------------------------------------
    
    return_to_N0 = {
      
      years_returned <-
        
        relative_year[
          `Oil spill` >=
            starting_population
        ]
      
      if (
        length(
          years_returned
        ) == 0
      ) {
        
        NA_real_
        
      } else {
        
        min(
          years_returned
        )
      }
    },
    
    
    # --------------------------------------------------------
    # B. RETURN TO MATCHED NO-SPILL TRAJECTORY
    # --------------------------------------------------------
    
    return_to_counterfactual =
      
      first_sustained_recovery(
        
        years =
          relative_year,
        
        recovered =
          within_1pct
      ),
    
    .groups =
      "drop"
  )


# ============================================================
# SUMMARIZE RECOVERY STATISTICS
# ============================================================

main_results_table <-
  
  recovery_by_sim %>%
  
  group_by(
    
    bad_frequency_label
    
  ) %>%
  
  summarise(
    
    # ========================================================
    # A. RETURN TO N0
    # ========================================================
    
    spill_return_rate =
      
      mean(
        !is.na(
          return_to_N0
        )
      ),
    
    spill_median_return =
      
      ifelse(
        
        any(
          !is.na(
            return_to_N0
          )
        ),
        
        median(
          return_to_N0,
          na.rm = TRUE
        ),
        
        NA_real_
      ),
    
    q25_spill_return =
      
      ifelse(
        
        any(
          !is.na(
            return_to_N0
          )
        ),
        
        quantile(
          return_to_N0,
          0.25,
          na.rm = TRUE
        ),
        
        NA_real_
      ),
    
    q75_spill_return =
      
      ifelse(
        
        any(
          !is.na(
            return_to_N0
          )
        ),
        
        quantile(
          return_to_N0,
          0.75,
          na.rm = TRUE
        ),
        
        NA_real_
      ),
    
    
    # ========================================================
    # B. RETURN TO MATCHED NO-SPILL TRAJECTORY
    # ========================================================
    
    counterfactual_return_rate =
      
      mean(
        !is.na(
          return_to_counterfactual
        )
      ),
    
    counterfactual_median_return =
      
      ifelse(
        
        any(
          !is.na(
            return_to_counterfactual
          )
        ),
        
        median(
          return_to_counterfactual,
          na.rm = TRUE
        ),
        
        NA_real_
      ),
    
    q25_counterfactual_return =
      
      ifelse(
        
        any(
          !is.na(
            return_to_counterfactual
          )
        ),
        
        quantile(
          return_to_counterfactual,
          0.25,
          na.rm = TRUE
        ),
        
        NA_real_
      ),
    
    q75_counterfactual_return =
      
      ifelse(
        
        any(
          !is.na(
            return_to_counterfactual
          )
        ),
        
        quantile(
          return_to_counterfactual,
          0.75,
          na.rm = TRUE
        ),
        
        NA_real_
      ),
    
    # --------------------------------------------------------
    # Number of simulations
    # --------------------------------------------------------
    
    n_simulations = n(),
    
    .groups =
      "drop"
  )


# ============================================================
# FORMAT FOR REPORT
# ============================================================

main_results_table_formatted <-
  
  main_results_table %>%
  
  mutate(
    
    `Bad-year frequency` =
      bad_frequency_label,
    
    
    # --------------------------------------------------------
    # SPILL RETURN TO N0
    # --------------------------------------------------------
    
    `Spill simulations returning to N₀` =
      percent(
        spill_return_rate,
        accuracy = 0.1
      ),
    
    
    `Median first return after spill (years)` =
      
      case_when(
        
        is.na(
          spill_median_return
        ) ~
          
          "Never",
        
        TRUE ~
          
          paste0(
            
            round(
              spill_median_return,
              0
            ),
            
            " (",
            
            round(
              q25_spill_return,
              0
            ),
            
            "–",
            
            round(
              q75_spill_return,
              0
            ),
            
            ")"
          )
      ),
    
    
    # --------------------------------------------------------
    # RECOVERY TO MATCHED NO-SPILL TRAJECTORY
    # --------------------------------------------------------
    
    `Simulations recovering to matched no-spill trajectory` =
      percent(
        counterfactual_return_rate,
        accuracy = 0.1
      ),
    
    
    `Median first return to matched no-spill trajectory (years)` =
      
      case_when(
        
        is.na(
          counterfactual_median_return
        ) ~
          
          "Never",
        
        TRUE ~
          
          paste0(
            
            round(
              counterfactual_median_return,
              0
            ),
            
            " (",
            
            round(
              q25_counterfactual_return,
              0
            ),
            
            "–",
            
            round(
              q75_counterfactual_return,
              0
            ),
            
            ")"
          )
      )
  ) %>%
  
  select(
    
    `Bad-year frequency`,
    
    `Spill simulations returning to N₀`,
    
    `Median first return after spill (years)`,
    
    `Simulations recovering to matched no-spill trajectory`,
    
    `Median first return to matched no-spill trajectory (years)`
    
  )


# ============================================================
# PRINT
# ============================================================

main_results_table_formatted


# ============================================================
# EXPORT TO WORD
# ============================================================

library(flextable)
library(officer)


ft <-
  
  flextable(
    main_results_table_formatted
  )


ft <-
  
  ft %>%
  
  # ----------------------------------------------------------
# HEADER
# ----------------------------------------------------------

bold(
  part = "header"
) %>%
  
  fontsize(
    size = 9,
    part = "header"
  ) %>%
  
  fontsize(
    size = 9,
    part = "body"
  ) %>%
  
  # ----------------------------------------------------------
# ALIGNMENT
# ----------------------------------------------------------

align(
  j = 2:5,
  align = "center",
  part = "all"
) %>%
  
  align(
    j = 1,
    align = "left",
    part = "all"
  ) %>%
  
  # ----------------------------------------------------------
# BORDERS
# ----------------------------------------------------------

border_outer(
  border = fp_border(
    color = "black",
    width = 1
  )
) %>%
  
  border_inner_h(
    border = fp_border(
      color = "grey70",
      width = 0.5
    )
  ) %>%
  
  border_inner_v(
    border = fp_border(
      color = "grey70",
      width = 0.5
    )
  ) %>%
  
  # ----------------------------------------------------------
# COLUMN WIDTHS / AUTOFIT
# ----------------------------------------------------------

autofit()


# ============================================================
# CREATE WORD DOCUMENT
# ============================================================

doc <-
  read_docx()


doc <-
  
  body_add_par(
    
    doc,
    
    "Table. Population recovery following a single oil-spill event under alternative frequencies of environmentally unfavorable years.",
    
    style =
      "Normal"
  )


doc <-
  
  body_add_flextable(
    
    doc,
    
    value =
      ft
  )


# ============================================================
# SAVE
# ============================================================

print(
  
  doc,
  
  target =
    "oil_spill_recovery_results_table.docx"
)