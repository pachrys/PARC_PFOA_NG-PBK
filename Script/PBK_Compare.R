# --------------------------------------------------------------------------- #
# SCRIPT FOR COMPARING PBK MODEL RESULTS BY PARAMETER
# By: Jack Koster
# Date: 08-05-2025
# --------------------------------------------------------------------------- #

# Packages
library(writexl)
library(ggplot2)
library(here)

# # Set output storage directory
# OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
# dir.create(OUTPUT, recursive = TRUE)

# Load results and rename

# Extract AUC and HalfLife #### 
extract_data <- function(dataset, dataset_name) {
  result_data <- data.frame(
    Individual = integer(),
    AUC = numeric(),
    HalfLife = numeric(),
    Condition = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(dataset$ANALYSED_data)) {
    indiv <- dataset$ANALYSED_data[[i]]
    
    if (!is.null(indiv) && all(c("AUC", "HalfLife") %in% names(indiv))) {
      # Extract and clean numeric part 
      auc_val <- as.numeric(sub("_.*", "", indiv[["AUC"]]))
      hl_val  <- as.numeric(sub("_.*", "", indiv[["HalfLife"]]))
      
      # Append to the result data frame
      result_data <- rbind(result_data, data.frame(
        Individual = i,
        AUC = auc_val,
        HalfLife = hl_val,
        Condition = dataset_name
      ))
    } else {
      warning("Missing or invalid data for individual ", i, " in ", dataset_name)
    }
  }
  
  return(result_data)
}

## Extract data for each data set ####
flow_data <- extract_data(RESULTS_Flow, "RESULTS_Flow")
age_data <- extract_data(RESULTS_Age, "RESULTS_Age")
flow_lt_data <- extract_data(RESULTS_Flow_LT, "RESULTS_Flow_LT")
age_lt_data <- extract_data(RESULTS_Age_LT, "RESULTS_Age_LT")

  final_results <- rbind(flow_data, age_data, flow_lt_data, age_lt_data)
  print(final_results)
  write_xlsx(final_results, path = here(OUTPUT, "final_results.xlsx"))

  
# Plots ####
# plot AUC 
p_AUC <- ggplot(final_results, aes(x = Individual, y = AUC, color = Condition)) +
  geom_point(size = 4) +  # Larger dots
  labs(title = "AUC by Condition",
       x = "Individual",
       y = "AUC") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    plot.title = element_text(size = 16),    
    axis.title = element_text(size = 14),   
    axis.text = element_text(size = 12)      
  )

# plot HalfLife
p_HalfLife <- ggplot(final_results, aes(x = Individual, y = HalfLife, color = Condition)) +
  geom_point(size = 4) +  # Larger dots
  labs(title = "HalfLife by Condition",
       x = "Individual",
       y = "HalfLife (days)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    plot.title = element_text(size = 16),    
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12)      
  )

ggsave(file.path(OUTPUT, "AUC_plot.png"), plot = p_AUC, width = 8, height = 6)
ggsave(file.path(OUTPUT, "HalfLife_plot.png"), plot = p_HalfLife, width = 8, height = 6)
