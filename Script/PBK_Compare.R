# PFOA PBK_Compare 8/5/2025 (JK) ####

# Packages
library(writexl)
library(ggplot2)
library(here)

# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

# Load results and rename

# Initialize the function to extract AUC and HalfLife from a given dataset
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
      # Extract and clean numeric part using regex
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

# Extract data for each data set
flow_data <- extract_data(RESULTS_Flow, "RESULTS_Flow")
age_data <- extract_data(RESULTS_Age, "RESULTS_Age")
flow_lt_data <- extract_data(RESULTS_Flow_LT, "RESULTS_Flow_LT")
age_lt_data <- extract_data(RESULTS_Age_LT, "RESULTS_Age_LT")

# Merge all results into one data frame
final_results <- rbind(flow_data, age_data, flow_lt_data, age_lt_data)
print(final_results)
write_xlsx(final_results, path = here(OUTPUT, "final_results.xlsx"))

# plot AUC
p_AUC <- ggplot(final_results, aes(x = Individual, y = AUC, color = Condition)) +
  geom_point(size = 4) +  # Larger dots
  labs(title = "AUC by Condition",
       x = "Individual",
       y = "AUC") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  # Larger legend title
    legend.text = element_text(size = 12),   # Larger legend text
    plot.title = element_text(size = 16),    # Larger plot title
    axis.title = element_text(size = 14),    # Larger axis titles
    axis.text = element_text(size = 12)      # Larger axis labels
  )

# plot HalfLife
p_HalfLife <- ggplot(final_results, aes(x = Individual, y = HalfLife, color = Condition)) +
  geom_point(size = 4) +  # Larger dots
  labs(title = "HalfLife by Condition",
       x = "Individual",
       y = "HalfLife (days)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  # Larger legend title
    legend.text = element_text(size = 12),   # Larger legend text
    plot.title = element_text(size = 16),    # Larger plot title
    axis.title = element_text(size = 14),    # Larger axis titles
    axis.text = element_text(size = 12)      # Larger axis labels
  )

ggsave(file.path(OUTPUT, "AUC_plot.png"), plot = p_AUC, width = 8, height = 6)
ggsave(file.path(OUTPUT, "HalfLife_plot.png"), plot = p_HalfLife, width = 8, height = 6)
