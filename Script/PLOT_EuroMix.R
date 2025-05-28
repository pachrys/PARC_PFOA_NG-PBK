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
PFOA_CP <- read_csv(here("Input", "EM_PFOA_CP.csv")) # (ng/ml)

# create ggplot2, 2 scatters overlayed, one for PFOA_CP one for max CP (from OUT_RAW_data) (x=id, y=CP ug/l)


## Regression on the plasma concentration starting point

# Import measured plasma concentrations
PFOA_CP <- read_csv(here("Input", "EM_PFOA_CP.csv")) # (ng/ml)

# Import RESULTS from the PBK simulation
PBK_OUT <- RESULTS$OUT_RAW_data


PredictedObserved.df <- PredictedObserved.df %>%
  mutate(
    CPatexpSTOP = map2_dbl( # Find the predicted concentration at the time of the stop of exposure to check if it's the same as the one measured
      PBK_OUT,
      expSTOP,
      ~ {
        idx <- which.min(abs(.x$time - .y))
        .x$CP[idx]
      }
    )
  ) %>%
  mutate(
    CP_final_predicted = map_dbl( # Find the predicted concentration at the end of the study
      PBK_OUT,
      ~ {
        last_row <- nrow(.x)
        .x$CP[last_row]
      }
    )
  )

CRegression <- lm(CP_final~CP_final_predicted, data=PredictedObserved.df)
summary(CRegression)



# V1


# Extract max predicted CP for each Idcode
predicted_df <- map_dfr(RESULTS$OUT_RAW_data, ~ {
  data_frame <- .x
  data_frame %>%
    summarise(
      Idcode = unique(Idcode),
      Predicted_CP = max(CP, na.rm = TRUE)
    )
})

measured_df <- read_csv(here("Input", "EM_PFOA_CP.csv")) %>%
  rename(Measured_CP = PFOA_CP)  # ng/ml = ug/L so same unit, no conversion needed

merged_df <- left_join(measured_df, predicted_df, by = "Idcode")

regression_model <- lm(Measured_CP ~ Predicted_CP, data = merged_df)
summary(regression_model)

ggplot(merged_df, aes(x = factor(Idcode))) +
  geom_point(aes(y = Measured_CP, color = "Measured"), alpha = 0.7) +
  geom_point(aes(y = Predicted_CP, color = "Predicted"), alpha = 0.7) +
  labs(
    x = "Idcode",
    y = "PFOA Concentration (ng/mL = ug/L)",
    color = "Legend",
    title = "Measured vs Predicted PFOA Plasma Concentration by Idcode"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# V2

