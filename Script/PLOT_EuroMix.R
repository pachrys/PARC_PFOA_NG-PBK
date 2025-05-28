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


# Extract max predicted plasma concentrations
pPFOA_df <- map_dfr(RESULTS$OUT_RAW_data, ~ {
  data_frame <- .x
  data_frame %>%
    summarise(
      Idcode = unique(Idcode),
      pPFOA_CP = max(CP, na.rm = TRUE)
    )
})

mPFOA_df <- read_csv(here("Input", "EM_PFOA_CP.csv")) %>%

PFOA_CP_df <- left_join(mPFOA_df, pPFOA_df, by = "Idcode")

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

