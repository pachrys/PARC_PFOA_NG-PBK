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
