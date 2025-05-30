# --------------------------------------------------------------------------- #
# SCRIPT FOR COMPARING SIMULATED AGAINST OBSERVED DATA (EUROMIX)
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
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

# Evaluate against measured serum PFOA 
# create ggplot2, 2 scatters overlayed, one for PFOA_CP one for 
# max CP (from OUT_RAW_data) (x=id, y=CP ug/l)

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

# Measured CP per subject
mPFOA_df <- read_csv(here("Input", "EM_PFOA_CP.csv"))

ext_pPFOA_CP <- function(result_list, expSTOP_vec) {
  map2_dfr(result_list, expSTOP_vec, ~ {
    idx <- which.min(abs(.x$time - .y))
    tibble(
      Idcode = .y,   
    )
  }) %>%
    mutate(Idcode = row_number())  
}

pPFOA_CP_flow_df <- ext_pPFOA_CP(RESULTS$OUT_RAW_data, Input$expSTOP)


# Extract GFR at expSTOP
ext_GFR <- function(result_list) {
  map2_dfr(result_list, seq_along(result_list), ~ {
    exp_stop_time <- unique(.x$expSTOP)
    idx <- which.min(abs(.x$time - exp_stop_time))
    tibble(
      Idcode = .y,
      GFR_at_expSTOP = .x$GFR[idx],
      max_GFR = max(.x$GFR, na.rm = TRUE)
    )
  })
}

pPFOA_CP_flow_df <- ext_pPFOA_CP(RESULTS_flow$ANALYSED_data)


ext_time <- function(result_list) {
  map2_dfr(result_list, seq_along(result_list), ~ {
    Tstart <- unique(.x$Tstart)
    Tstop <- unique(.x$Tstop)
    tibble(
      Idcode = .y,
      duration = Tstop - Tstart
    )
  })
}

expCONC_df <- ext_expCONC_raw(RESULTS_flow$ANALYSED_data)


# Extract expCON (exposure concentration)
ext_expCONC <- function(data_list) {
  map2_dfr(data_list, seq_along(data_list), ~ {
    raw_val <- .x$expCONC
    if (is.null(raw_val)) raw_val <- NA_character_
    tibble(
      Idcode = .y,
      expCONC_raw = as.character(raw_val)
    )
  })
}

expCONC_df <- extract_expCONC_raw(RESULTS_flow$ANALYSED_data)


# Create flow DFs
# Create age DFs

# Create full DF
PFOA_OUT_df <- function(pPFOA_df, GFR_df, duration_df) {
  pPFOA_df %>%
    left_join(GFR_df, by = "Idcode") %>%
    left_join(duration_df, by = "Idcode")
}
# Example for RESULTS_flow
merged_flow <- merge_all(pPFOA_flow, GFR_flow, duration_flow)

# Example for RESULTS_flow_age
merged__age <- merge_all(pPFOA__age, GFR__age, duration_age)


view(PFOA_CP_full_df)

write.csv(PFOA_CP_full_df, "PFOA_CP_full_df.csv", row.names = FALSE)


# Check CPs only
merged_df <- mPFOA_df %>%
  left_join(pPFOA_CP_flow_df, by = "Idcode")




# Regression
regression_mod <- lm(mPFOA_CP ~ pPFOA_CP, data = PFOA_CP_df)
summary(regression_mod)

ggplot(PFOA_CP_df, aes(x = factor(Idcode))) +
  geom_point(aes(y = mPFOA_df, color = "Measured"), alpha = 0.8) +
  geom_point(aes(y = pPFOA_df, color = "Predicted"), alpha = 0.8) +
  labs(
    x = "Idcode",
    y = "PFOA Concentration in Plasma (ng/mL)",
    color = "Legend",
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot CP log on measured vs. predicted plot. 



# Weird plasma behavior

# Extract CALC_Parameters for each group
low_params  <- RESULTS[["CALC_Parameters"]][low_ids]
high_params <- RESULTS[["CALC_Parameters"]][high_ids]

# Collect all unique parameter names from each group
low_param_names  <- unique(unlist(lapply(low_params, names)))
high_param_names <- unique(unlist(lapply(high_params, names)))

# Parameters only in high CP group
only_in_high <- setdiff(high_param_names, low_param_names)

# Parameters only in low CP group
only_in_low <- setdiff(low_param_names, high_param_names)

# Result
list(
  only_in_high = only_in_high,
  only_in_low = only_in_low
)

low_exp <- lapply(low_ids, function(i) {
  conc <- RESULTS[["ANALYSED_data"]][[i]][["expCONC"]]
  data.frame(ID = i, Time = conc$time, expCONC = conc$conc)
})

high_exp <- lapply(high_ids, function(i) {
  conc <- RESULTS[["ANALYSED_data"]][[i]][["expCONC"]]
  data.frame(ID = i, Time = conc$time, expCONC = conc$conc)
})

# Combine into data frames
df_low <- do.call(rbind, low_exp)
df_high <- do.call(rbind, high_exp)

# Add group labels
df_low$Group <- "Low CP"
df_high$Group <- "High CP"