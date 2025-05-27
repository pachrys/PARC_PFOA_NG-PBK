# EuroMix File Structure ####

# Packages
library(readxl)
library(dplyr)
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)
library(showtext)
library(readr)
library(scales)


exp_data <- read_xlsx(here("Input", "EM_PFOA_Exp.xlsx"))


mean_oral <- exp_data %>%
  select(mean_oral) %>%
  mutate(mean_oral = as.numeric(gsub("\"", "", mean_oral)),
         mean_oral = format(mean_oral, scientific = FALSE, trim = TRUE))

mean_dermal <- exp_data %>%
  select(mean_dermal) %>%
  mutate(mean_dermal = as.numeric(gsub("\"", "", mean_dermal)),
         mean_dermal = format(mean_dermal, scientific = FALSE, trim = TRUE))

write.csv(mean_oral, here("Input", "EM_PFOA_Exp_mean_oral.csv"), row.names = FALSE, quote = FALSE)
write.csv(mean_dermal, here("Input", "EM_PFOA_Exp_mean_dermal.csv"), row.names = FALSE, quote = FALSE)
