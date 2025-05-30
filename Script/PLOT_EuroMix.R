# --------------------------------------------------------------------------- #
# SCRIPT FOR COMPARING MEASURED AND PREDICTED DATA (EUROMIX)
# By: Jack Koster
# Date: 27-05-2025
# --------------------------------------------------------------------------- #

# rm(list=ls()) # to clear out the global environment #Input and Physio.c needed

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)
library(showtext)
library(purrr)
library(ggplot2)
library(stats)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

# Plot ideas
# Flow vs age (compare cmax, halflife) for each
# 
# mCP v sCP (measured observed conc in plasma)
# 
# 
# Correlation plot GFR and half-life
# 
# expAGE vs CP / expAGE vs Halflife

# Correlation plot duration of exposure and half-life (between same age/GFR)

# Data Extraction ####

# Load measured CP per subject
mPFOA_df <- read_csv(here("Input", "EM_PFOA_CP.csv"))

# Extract expCONC from Input
ext_expCONC <- function(expCONC_vec, idcode_vec) {
  tibble(
    Idcode = idcode_vec,
    expCONC = expCONC_vec
  )
}

expCONC_df <- ext_expCONC(Input$expCONC, Input$Idcode)


# Extract predicted CP at expSTOP
ext_pPFOA_CP <- function(result_list, expSTOP_vec) {
  map2_dfr(seq_along(result_list), expSTOP_vec, ~ {
    data_i <- result_list[[.x]]
    idx <- which.min(abs(data_i$time - .y))
    tibble(
      Idcode = .x,
      pPFOA_CP = data_i$CP[idx]
    )
  })
}

pPFOA_CP_flow_df <- ext_pPFOA_CP(RESULTS_flow$OUT_RAW_data, Input$expSTOP)
pPFOA_CP_age_df <- ext_pPFOA_CP(RESULTS_age$OUT_RAW_data, Input$expSTOP)

# OR

# # Extract predicted maximum CP (over predicts for old age Ids)
# ext_max_pPFOA_CP <- function(result_list) {
#   map_dfr(seq_along(result_list), ~ {
#     data_i <- result_list[[.x]]
#     tibble(
#       Idcode = .x,
#       pPFOA_CP = max(data_i$CP, na.rm = TRUE)
#     )
#   })
# }
# pPFOA_CP_flow_df <- ext_max_pPFOA_CP(RESULTS$OUT_RAW_data)


# Extract Halflives
ext_HalfLife <- function(analysed_list, idcode_vec) {
  map2_dfr(analysed_list, idcode_vec, ~ {
    half_life_c <- .x$HalfLife
    half_life_num <- as.numeric(sub("_years$", "", half_life_c))
    
    tibble(
      Idcode = .y,
      HalfLife = half_life_num
    )
  })
}

HalfLife_flow_df <- ext_HalfLife(RESULTS_flow$ANALYSED_data, Input$Idcode)
HalfLife_age_df <- ext_HalfLife(RESULTS_age$ANALYSED_data, Input$Idcode)

# Extract GFR for flow GFR
ext_GFR_flow <- function(calc_params_list, idcode_vec) {
  map_dfr(seq_along(calc_params_list), function(i) {
    tibble(
      Idcode = idcode_vec[[i]],
      GFR_flow = calc_params_list[[i]]$GFR
    )
  })
}

GFR_flow_df <- ext_GFR_flow(RESULTS_flow$CALC_Parameters, Input$Idcode)


# Extract GFR for age GFR
ext_GFR_age <- function(physio_df, expSTOP_vec, sex_vec, idcode_vec) {
  map_dfr(seq_along(expSTOP_vec), function(i) {
    stop_time <- expSTOP_vec[i]
    sex <- sex_vec[i]
    id <- idcode_vec[i]
    
    idx <- which.min(abs(physio_df$TIME - stop_time))
    
    gfr_value <- if (sex == "M") {
      physio_df$GFR_M[idx]
    } else if (sex == "F") {
      physio_df$GFR_F[idx]
    } 
    
    tibble(
      Idcode = id,
      GFR = gfr_value / 1.44    # (L/d -> ml/min)
    )
  })
}

GFR_age_df <- ext_GFR_age(Physio.c, Input$expSTOP, Input$sex, Input$Idcode)


# Merging dataframes and anonymizing ####  
merge_all <- function(
    mPFOA_df, 
    expCONC_df,
    pPFOA_CP_flow_df,
    pPFOA_CP_age_df,
    HalfLife_flow_df, 
    HalfLife_age_df,
    GFR_flow_df, 
    GFR_age_df
) {
  merged_df <- mPFOA_df %>%
    left_join(expCONC_df, by = "Idcode") %>%
    left_join(pPFOA_CP_flow_df, by = "Idcode") %>%
    left_join(pPFOA_CP_age_df, by = "Idcode") %>%
    left_join(HalfLife_flow_df, by = "Idcode") %>%
    left_join(HalfLife_age_df, by = "Idcode") %>%
    left_join(GFR_flow_df, by = "Idcode") %>%
    left_join(GFR_age_df, by = "Idcode")
  
  # Anonymize Idcodes
  set.seed(123)
  anon_labels <- tibble(
    Idcode = unique(merged_df$Idcode),
    AnonID = paste0("ID_", sample(seq_along(unique(merged_df$Idcode))))
  )
  
  merged_df <- merged_df %>%
    left_join(anon_labels, by = "Idcode")
  
  return(merged_df)
}

PFOA_OUT_df <- merge_all(
  mPFOA_df, 
  expCONC_df,
  pPFOA_CP_flow_df,
  pPFOA_CP_age_df,
  HalfLife_flow_df, 
  HalfLife_age_df,
  GFR_flow_df, 
  GFR_age_df
)


# Pivots for plotting ####

# Pivot for CP
plot_df <- PFOA_OUT_df %>%
  pivot_longer(
    cols = c(mPFOA_CP, pPFOA_CP.y),
    names_to = "Type",
    values_to = "Concentration"
  )


# Plots ####

# Plot using AnonID
ggplot(plot_df, aes(x = AnonID, y = Concentration, color = Type)) +
  geom_point(size = 3) +
  geom_smooth(
    aes(group = Type),
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    size = 1
  ) +
  scale_y_log10()+
  labs(
    y = "PFOA Concentration in Plasma (ng/mL)",
    color = "Measured vs Predicted"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank()
  )
  scale_color_manual(values = c("mPFOA_CP" = "darkorange", "pPFOA_CP" = "steelblue"))
  scale_y_log10()
  

# # Regression
# regression_mod <- lm(mPFOA_CP ~ pPFOA_CP, data = PFOA_CP_df)
# summary(regression_mod)
# 
# ggplot(PFOA_CP_df, aes(x = factor(Idcode))) +
#   geom_point(aes(y = mPFOA_df, color = "Measured"), alpha = 0.8) +
#   geom_point(aes(y = pPFOA_df, color = "Predicted"), alpha = 0.8) +
#   labs(
#     x = "Idcode",
#     y = "PFOA Concentration in Plasma (ng/mL)",
#     color = "Legend",
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

