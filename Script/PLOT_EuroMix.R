# --------------------------------------------------------------------------- #
# SCRIPT FOR COMPARING MEASURED AND PREDICTED DATA (EUROMIX)
# By: Jack Koster
# Date: 27-05-2025
# --------------------------------------------------------------------------- #

# rm(list=ls()) # to clear out the global environment

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


# Data Extraction ####

# Load RESULTS and rename

# RESULTS_flow <- RESULTS #rm(RESULTS)
# load(RESULTS.R)
# RESULTS_age <- RESULTS #rm(RESULTS)

# Load raw data
Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
Input <- read_csv(here("Input", "INPUT_EuroMix_20year.csv")) 
mPFOA_CP_df <- read_csv(here("Input", "EM_PFOA_CP.csv")) 

# Extract from Input

  ext_var <- function(var_name) {
    tibble(
      Idcode = Input$Idcode,
      !!var_name := Input[[var_name]]
    )
  }

sex_df <- ext_var("sex")
exptype_df <- ext_var("exposure_type")
expAGE_df <- ext_var("expAGE")
expCONC_df <- ext_var("expCONC")
  

# Extract predicted CP at expSTOP

  ext_pPFOA_CP <- function(result_list, expSTOP_vec) {
    map2_dfr(seq_along(result_list), expSTOP_vec, ~ {
      j <- result_list[[.x]]
      idx <- which.min(abs(j$time - .y))
      tibble(
        Idcode = .x,
        pPFOA_CP = j$CP[idx]
      )
    })
  }

pPFOA_CP_flow_df <- ext_pPFOA_CP(RESULTS_flow$OUT_RAW_data, Input$expSTOP)
pPFOA_CP_age_df <- ext_pPFOA_CP(RESULTS_age$OUT_RAW_data, Input$expSTOP)

# OR

# Extract predicted maximum CP (over predicts for old age Ids)

#   ext_max_pPFOA_CP <- function(result_list) {
#     map_dfr(seq_along(result_list), ~ {
#       j <- result_list[[.x]]
#       tibble(
#         Idcode = .x,
#         pPFOA_CP = max(j$CP, na.rm = TRUE)
#       )
#     })
#   }
# 
# pPFOA_CP_flow_df <- ext_max_pPFOA_CP(RESULTS_flow$OUT_RAW_data)
# pPFOA_CP_age_df <- ext_max_pPFOA_CP(RESULTS_age$OUT_RAW_data)


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
        GFR.x = calc_params_list[[i]]$GFR
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
        GFR.y = gfr_value / 1.44    # (L/d -> ml/min)
      )
    })
  }

GFR_age_df <- ext_GFR_age(Physio.c, Input$expSTOP, Input$sex, Input$Idcode)


# Merge dataframes ####  

  PFOA_OUT_df <- sex_df %>%
    left_join(exptype_df, by = "Idcode") %>%
    left_join(expAGE_df, by = "Idcode") %>%
    left_join(expCONC_df, by = "Idcode") %>%
    left_join(mPFOA_CP_df, by = "Idcode") %>%
    left_join(pPFOA_CP_flow_df, by = "Idcode") %>%
    left_join(pPFOA_CP_age_df, by = "Idcode") %>%
    left_join(HalfLife_flow_df, by = "Idcode") %>%
    left_join(HalfLife_age_df, by = "Idcode") %>%
    left_join(GFR_flow_df, by = "Idcode") %>%
    left_join(GFR_age_df, by = "Idcode")


## Sort IDs for plotting ####
ID_sort <- PFOA_OUT_df %>%
  distinct(Idcode, expAGE) %>%  # first by expAGE then Idcode
  arrange(expAGE, Idcode) %>%
  pull(Idcode)

PFOA_OUT_df <- PFOA_OUT_df %>%
  mutate(ageID = factor(Idcode, levels = ID_sort))


PFOA_Oral_df <- PFOA_OUT_df %>% filter(exposure_type == "Oral")
PFOA_OD_df <- PFOA_OUT_df %>% filter(exposure_type == "Oral_Dermal")

# write.csv(PFOA_OUT_df, here("Input", "PFOA_OUT_df.csv"), row.names = FALSE)


# Data Analysis ####

## Regression and Correlation ####

# Flow-based model
lm_flow <- lm(mPFOA_CP ~ pPFOA_CP.x, data = PFOA_OUT_df)
summary(lm_flow)

# Age-based model
lm_age <- lm(mPFOA_CP ~ pPFOA_CP.y, data = PFOA_OUT_df)
summary(lm_age)


# Halflife GFR
lm_hlgfr_flow <- lm(HalfLife.x ~ GFR.x, data = PFOA_OUT_df)


lm_hlgfr_age <- lm(HalfLife.y ~ GFR.y, data = PFOA_OUT_df)

summary(lm_hlgfr_flow)
summary(lm_hlgfr_age)

# HalfLife GFR Corr
corr_flow <- cor(PFOA_OUT_df$GFR.x, PFOA_OUT_df$HalfLife.x, use = "complete.obs", method = "pearson")
corr_age <- cor(PFOA_OUT_df$GFR.y, PFOA_OUT_df$HalfLife.y, use = "complete.obs", method = "pearson")




# Plots ####

## Fig 1a. Plot measured vs. predicted plasma concentrations ####

  # # Flow-based model
# lm_flow <- lm(mPFOA_CP ~ pPFOA_CP.x, data = PFOA_OUT_df)
# summary(lm_flow)


  
  CP_wide <- PFOA_OD_df %>%
    select(sex, mPFOA_CP, pPFOA_CP.x) %>%
    rename(
      Measured = "mPFOA_CP",
      Predicted = "pPFOA_CP.x") 

 CP_wide_F <- CP_wide %>% filter(sex == "F")
 CP_wide_M <- CP_wide %>% filter(sex == "M")
 
CRegression <- lm(Predicted~Measured, data=CP_wide_M)
summary(CRegression)

CRegression_plot <- CRegression %>%
  # Plot in log scale
  ggplot(aes(x = Measured, y = Predicted)) +
  
  # # Linear regression line
  # geom_smooth(method = 'lm', color = "black", se = TRUE) +
  # # reference lines
  # geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.5, color = "grey50") +  
  # geom_abline(intercept = log(1.1), slope = 1, linetype = "dashed", linewidth = 0.5, color = "grey50") +  # +10% line
  # geom_abline(intercept = log(0.9), slope = 1, linetype = "dashed", linewidth = 0.5, color = "grey50") +  # -10% line
  # geom_abline(intercept = log(2), slope = 1, linetype = "dotted", linewidth = 0.5, color = "grey50") +  # 2-fold upper
  # geom_abline(intercept = log(0.5), slope = 1, linetype = "dotted", linewidth = 0.5, color = "grey50") +  # 2-fold lower
  
  # add the actual points (these are the observed points)
  geom_point(color = "darkred", size = 1) +
  
  # theme etc
  theme_minimal()+
  theme(plot.title = element_text(size = 10, margin = margin(b = 20)),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

CRegression_plot


  p_PFOA_CPtest <- ggplot(CP_wide, aes(x = Measured, y = Predicted)) +
    geom_point(size = 3, color = "darkorchid3") +
    geom_smooth(method = "lm", linetype = "dashed", linewidth = 1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Measured plasma PFOA (ng/ml)",
      y = "Predicted plasma PFOA (ng/ml)",
      color = "Measured vs Predicted"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15))
    )

  ggsave(filename = here(OUTPUT, "Fig1a_CP.png"),
         dpi = 300,
         width = 24,
         height = 16,
         units = "cm")

p_PFOA_CPa



## Fig 1b. Plot age vs. measured and predicted plasma concentrations (by sex) ####

  CP_long <- PFOA_OUT_df %>%
    pivot_longer(
      cols = c(mPFOA_CP, pPFOA_CP.y), # .x for flow, .y for age
      names_to = "Type",
      values_to = "Concentration"
    ) %>%
    mutate(
      Type = recode(Type,
                    mPFOA_CP = "Measured",
                    pPFOA_CP.y = "Predicted"
      ),
      Sex = recode(sex, M = "Male", F = "Female")
    )
  
  CP_wide <- CP_long %>%
    pivot_wider(names_from = Type, values_from = Concentration)
  
  # # Correlation
  # corr_val <- cor(CP_wide$Measured, CP_wide$Predicted, use = "complete.obs", method = "pearson")
  
  # Correlation
  # corr_val <- cor(PFOA_OUT_df$GFR.x, PFOA_OUT_df$HalfLife.x, use = "complete.obs", method = "pearson")
  
  p_PFOA_CPb <- ggplot(CP_long, aes(x = ageID, y = Concentration, color = Type)) +
    geom_point(aes(shape = Sex), size = 3, stroke = 1, fill = "white") +
    geom_smooth(aes(group = Type), method = "lm", linetype = "dashed", size = 1) +
    scale_y_log10() +
    scale_shape_manual(values = c("Male" = 16, "Female" = 1)) +
    scale_color_manual(values = c("Measured" = "darkseagreen4", "Predicted" = "darkorchid3")) +
    labs(
      y = "PFOA Concentration in Plasma (ng/ml)",
      color = "Measured vs Predicted",
      shape = "Sex"
    ) +
    annotate("text", x = Inf, y = Inf, label = paste0("r = ", round(corr_val, 2)),
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
    theme_minimal(base_size = 20) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),

      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15))
    )
  
  ggsave(filename = here(OUTPUT, "Fig1b_CP_sex.png"),
         dpi = 300,
         width = 24,
         height = 16,
         units = "cm")
  
p_PFOA_CPb
  
  
## Fig 1c. Plot measured vs. predicted plasma concentrations (by exposure) ####
  
  CP_long <- PFOA_OUT_df %>%
    pivot_longer(
      cols = c(mPFOA_CP, pPFOA_CP.y),
      names_to = "Type",
      values_to = "Concentration"
    ) %>%
    mutate(
      Type = recode(Type,
                    mPFOA_CP = "Measured",
                    pPFOA_CP.y = "Predicted"),
      Exposure = recode(exposure_type,
                        Oral_Dermal = "Oral + Dermal",
                        Oral = "Oral")
    )
  
  CP_wide <- CP_long %>%
    pivot_wider(names_from = Type, values_from = Concentration)
  
  # Correlation
  corr_val <- cor(CP_wide$Measured, CP_wide$Predicted, use = "complete.obs", method = "pearson")
  
  p_PFOA_CPc <- ggplot(CP_long, aes(x = ageID, y = Concentration, color = Type)) +
    geom_point(aes(shape = Exposure), size = 3, stroke = 1, fill = "white") +
    geom_smooth(aes(group = Type), method = "lm", linetype = "dashed", size = 1) +
    scale_y_log10() +
    scale_shape_manual(values = c("Oral + Dermal" = 16, "Oral" = 1)) +
    scale_color_manual(values = c("Measured" = "darkseagreen4", "Predicted" = "darkorchid3")) +
    labs(
      y = "PFOA Concentration in Plasma (ng/ml)",
      color = "Measured vs Predicted",
      shape = "Exposure Type"
    ) +
    annotate("text", x = Inf, y = Inf, label = paste0("r = ", round(corr_val, 2)),
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
    theme_minimal(base_size = 20) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15))
    )
  
  ggsave(filename = here(OUTPUT, "Fig1c_CP_exp.png"),
         dpi = 300,
         width = 24,
         height = 16,
         units = "cm")
  
p_PFOA_CPc
  

## Fig 2. Plot half-lives vs. predicted plasma concentrations ####
  
  hl_pPFOA_long <- PFOA_OUT_df %>%
    select(ageID, 
           HalfLife_flow = HalfLife.x, 
           HalfLife_age = HalfLife.y,
           pPFOA_CP_flow = pPFOA_CP.x, 
           pPFOA_CP_age = pPFOA_CP.y) %>%
    pivot_longer(
      cols = c(HalfLife_flow, HalfLife_age, pPFOA_CP_flow, pPFOA_CP_age),
      names_to = c("Metric", "Type"),
      names_pattern = "(.+)_(.+)$",
      values_to = "Value"
    ) %>%
    pivot_wider(names_from = Metric, values_from = Value)

  p_hl_pPFOA <- ggplot(hl_pPFOA_long, aes(x = HalfLife, y = pPFOA_CP, color = Type)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_log10() +
    labs(
      x = "Half-Life (years)",
      y = "Predicted Plasma Concentration (ng/ml)",
      color = "GFR Type"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.title.x = element_text(margin = margin(t = 15)), 
      axis.title.y = element_text(margin = margin(r = 15))
    ) +
    scale_color_manual(values = c("flow" = "aquamarine3", "age" = "coral3"))
  
  ggsave(filename = here(OUTPUT, "Fig2_hl_CP.png"),
         dpi = 300,
         width = 24,
         height = 16,
         units = "cm")
  
p_hl_pPFOA
  
  
## Fig 3. Plot GFRs vs. predicted plasma concentrations ####  
  
  gfr_pPFOA_long <- PFOA_OUT_df %>%
    select(ageID, GFR_flow = GFR.x, GFR_age = GFR.y, pPFOA_flow = pPFOA_CP.x, pPFOA_age = pPFOA_CP.y) %>%
    pivot_longer(cols = -ageID, names_to = c("Metric", "Type"), names_sep = "_") %>%
    pivot_wider(names_from = Metric, values_from = value)
  
  p_GFR_pPFOA <- ggplot(gfr_pPFOA_long, aes(x = GFR, y = pPFOA, color = Type)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_y_log10() +
    labs(x = "GFR (ml/min)", y = "Predicted Plasma PFOA (ng/ml)", color = "GFR Type") +
    theme_minimal(base_size = 20) +
    scale_color_manual(values = c("flow" = "aquamarine3", "age" = "coral3"))
  
  ggsave(filename = here(OUTPUT, "Fig3_GFR_pPFOA.png"),
         dpi = 300,
         width = 24,
         height = 16,
         units = "cm")
  
p_GFR_pPFOA
  

## Fig 4. Plot GFRs vs. half-lives ####

  gfr_hl_long <- PFOA_OUT_df %>%
    select(ageID, GFR_flow = GFR.x, GFR_age = GFR.y,
           HalfLife_flow = HalfLife.x, HalfLife_age = HalfLife.y) %>%
    pivot_longer(cols = -ageID, names_to = c("Metric", "Type"), names_sep = "_") %>%
    pivot_wider(names_from = Metric, values_from = value)
  
  p_GFR_HL <- ggplot(gfr_hl_long, aes(x = HalfLife, y = GFR, color = Type)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(x = "Half-Life (years)", y = "GFR (ml/min)", color = "GFR Type") +
    theme_minimal(base_size = 20)+
    theme(
      axis.title.x = element_text(margin = margin(t = 5)), 
      axis.title.y = element_text(margin = margin(r = 5))
    ) +
    scale_color_manual(values = c("flow" = "aquamarine3", "age" = "coral3"))
  
  ggsave(filename = here(OUTPUT, "Fig4_GFR_hl.png"),
         dpi = 300,
         width = 12,
         height = 8,
         units = "cm")

p_GFR_HL


## Fig 5. Plot expCONC vs half-lives ####

  expconc_hl_long <- PFOA_OUT_df %>%
    select(ageID, expCONC, HalfLife_flow = HalfLife.x, HalfLife_age = HalfLife.y) %>%
    pivot_longer(cols = starts_with("HalfLife"), names_to = "Type", names_prefix = "HalfLife_", values_to = "HalfLife")

  p_expCONC_HL <- ggplot(expconc_hl_long, aes(x = expCONC, y = HalfLife, color = Type)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_x_log10() +
    labs( x = "Exposure Concentration (ug/kg/day)", y = "Half-Life (years)", color = "GFR Type") +
    theme_minimal(base_size = 20) +
    theme(
      axis.title.x = element_text(margin = margin(t = 15)), 
      axis.title.y = element_text(margin = margin(r = 15))
    ) +
    scale_color_manual(values = c("flow" = "aquamarine3", "age" = "coral3"))

  ggsave(filename = here(OUTPUT, "Fig5_expCONC_HL.png"),
         dpi = 300,
         width = 24,
         height = 16,
         units = "cm")

p_expCONC_HL
  
  
## Fig 6. Plot expCONC vs predicted PFOA concentration #### 

# Measured

  expconc_cp_long <- PFOA_OUT_df %>%
    select(ageID, expCONC, pPFOA_CP_flow = pPFOA_CP.x, pPFOA_CP_age = pPFOA_CP.y) %>%
    pivot_longer(
      cols = starts_with("pPFOA_CP"),
      names_to = "Type",
      names_prefix = "pPFOA_CP_",
      values_to = "Predicted_CP"
    )
  
  p_expCONC_CP <- ggplot(expconc_cp_long, aes(x = expCONC, y = Predicted_CP, color = Type)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Exposure Concentration (ug/kg/day)",
      y = "Predicted Plasma Concentration (ng/ml)",
      color = "GFR Type"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.title.x = element_text(margin = margin(t = 15)), 
      axis.title.y = element_text(margin = margin(r = 15))
    ) +
    scale_color_manual(values = c("flow" = "aquamarine3", "age" = "coral3"))
  
  ggsave(filename = here(OUTPUT, "Fig6_expCONC_CP.png"),
         dpi = 300,
         width = 24,
         height = 16,
         units = "cm")

p_expCONC_CP
  

# Data Analysis ####
  
## Regression and Correlation ####

# Flow-based model
lm_flow <- lm(mPFOA_CP ~ pPFOA_CP.x, data = PFOA_OUT_df)
summary(lm_flow)

# Age-based model
lm_age <- lm(mPFOA_CP ~ pPFOA_CP.y, data = PFOA_OUT_df)
summary(lm_age)


# Halflife GFR
lm_hlgfr_flow <- lm(HalfLife.x ~ GFR.x, data = PFOA_OUT_df)


lm_hlgfr_age <- lm(HalfLife.y ~ GFR.y, data = PFOA_OUT_df)

summary(lm_hlgfr_flow)
summary(lm_hlgfr_age)

# HalfLife GFR Corr
# Correlation
corr_flow <- cor(PFOA_OUT_df$GFR.x, PFOA_OUT_df$HalfLife.x, use = "complete.obs", method = "pearson")
corr_age <- cor(PFOA_OUT_df$GFR.y, PFOA_OUT_df$HalfLife.y, use = "complete.obs", method = "pearson")

# Plot Lms
GFR_hl_flow <- corr_flow %>%   ggplot(aes(x = GFR.x, y = HalfLife.x)) +   geom_smooth(method = 'lm', color = "black", se = TRUE) +   geom_point(color = "darkred", size = 1) +
  theme_minimal()+
  theme(plot.title = element_text(size = 10, margin = margin(b = 20)),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +   labs(title = "Correlation between GFR and Halflife", x = "GFR", y = "HalfLife") +
  annotate("text", x = Inf, y = Inf, label = paste0("r = ", round(corr_flow, 2)),
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
GFR_hl_flow


GFR_hl_age <- corr_age %>%   ggplot(aes(x = GFR.y, y = HalfLife.y)) +   geom_smooth(method = 'lm', color = "black", se = TRUE) +   geom_point(color = "darkred", size = 1) +
  theme_minimal()+
  theme(plot.title = element_text(size = 10, margin = margin(b = 20)),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +   labs(title = "Correlation between GFR and Halflife", x = "GFR", y = "HalfLife") +
  annotate("text", x = Inf, y = Inf, label = paste0("r = ", round(corr_age, 2)),
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
GFR_hl_age




## Fig 7. Plot measured and predicted and flow vs. age linear regression ####

  PFOA_lm_long <- PFOA_OUT_df %>%
    select(ageID, mPFOA_CP, pPFOA_flow = pPFOA_CP.x, pPFOA_age = pPFOA_CP.y) %>%
    pivot_longer(cols = starts_with("pPFOA"), names_to = "Model", values_to = "pPFOA_CP")
  
  p_PFOA_lm <- ggplot(PFOA_lm_long, aes(x = pPFOA_CP, y = mPFOA_CP, color = Model)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Predicted Plasma PFOA (ng/ml)",
      y = "Measured Plasma PFOA (ng/ml)",
      color = "GFR Type"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15))
    ) +
    scale_color_manual(values = c("pPFOA_flow" = "aquamarine3", "pPFOA_age" = "coral3"))

  ggsave(filename = here(OUTPUT, "Fig7_lm.png"),
         dpi = 300,
         width = 24,
         height = 16,
         units = "cm")
  
p_PFOA_lm


# Compare linear regression
broom::glance(lm_flow)
broom::glance(lm_age)



# # Anonymize Idcodes
# # Can be added in the merge_all before calling the function to randomize IDs
# set.seed(123)
# anons <- tibble(
#   Idcode = unique(merged_df$Idcode),
#   AnonID = paste0("ID_", sample(seq_along(unique(merged_df$Idcode))))
# )
# 
# merged_df <- merged_df %>%
#   left_join(anons, by = "Idcode")

