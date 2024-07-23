# Load necessary libraries
library(dplyr)

# Function to calculate standard curve parameters
calculate_standard_curve <- function(log_copies, cq_values) {
  fit <- lm(cq_values ~ log_copies)
  slope <- coef(fit)[2]
  intercept <- coef(fit)[1]
  efficiency <- 10^(-1/slope)
  return(list(slope = slope, intercept = intercept, efficiency = efficiency))
}

# Function to calculate copies from Cq values using standard curve
calculate_copies <- function(cq_values, slope, intercept) {
  log_copies <- (cq_values - intercept) / slope
  copies <- 10^log_copies
  return(copies)
}

# F9 standard curve data
f9_standard_data <- data.frame(
  log_copies = c(8.815644149, 7.815644149, 6.815644149, 5.815644149, 4.815644149, 3.815644149, 2.815644149, 1.815644149),
  cq_values = c(5.711062857, 10.0192069, 13.30810448, 16.56234853, 21.05404185, 24.0781887, 28.26, 31.92)
)

# Calculate F9 standard curve parameters
f9_standard_curve <- calculate_standard_curve(f9_standard_data$log_copies, f9_standard_data$cq_values)
f9_slope <- f9_standard_curve$slope
f9_intercept <- f9_standard_curve$intercept
f9_efficiency <- f9_standard_curve$efficiency

# Ppia standard curve data
ppia_standard_data <- data.frame(
  log_copies = c(6.485366466, 5.485366466, 4.485366466, 3.485366466, 2.485366466, 1.485366466),
  cq_values = c(16.0433524, 19.02719368, 22.57893792, 26.362526, 30.3438398, 33.652062)
)

# Calculate Ppia standard curve parameters
ppia_standard_curve <- calculate_standard_curve(ppia_standard_data$log_copies, ppia_standard_data$cq_values)
ppia_slope <- ppia_standard_curve$slope
ppia_intercept <- ppia_standard_curve$intercept
ppia_efficiency <- ppia_standard_curve$efficiency

# Print the standard curve parameters
cat("F9 Standard Curve Parameters:\n")
cat("Slope:", f9_slope, "\nIntercept:", f9_intercept, "\nEfficiency:", f9_efficiency, "\n\n")
cat("Ppia Standard Curve Parameters:\n")
cat("Slope:", ppia_slope, "\nIntercept:", ppia_intercept, "\nEfficiency:", ppia_efficiency, "\n\n")

# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# Load the RNA-seq data from the CSV file
file_path <- "C:/Users/harum/OneDrive/Documents/LAB/mhpQPCR/qpcrcsv/timepoint_copynumber.csv"
sample_data <- read_csv(file_path)

# Display the head of the data to verify it's loaded correctly
print(head(sample_data))

# Function to calculate copies from Cq values using standard curve
calculate_copies <- function(cq_values, slope, intercept) {
  log_copies <- (cq_values - intercept) / slope
  copies <- 10^log_copies
  return(copies)
}

# Standard curve parameters for F9
f9_slope <- -3.194
f9_intercept <- 36.990
f9_efficiency <- 10^(-1/f9_slope)

# Standard curve parameters for Ppia
ppia_slope <- -3.594
ppia_intercept <- 38.990
ppia_efficiency <- 10^(-1/ppia_slope)

# Calculate absolute copy numbers for F9 and Ppia for each time point
sample_data <- sample_data %>%
  pivot_longer(cols = starts_with("Day"), names_to = "Day", values_to = "Cq") %>%
  mutate(
    Copies = case_when(
      Target == "F9" ~ calculate_copies(Cq, f9_slope, f9_intercept),
      Target == "PPIA" ~ calculate_copies(Cq, ppia_slope, ppia_intercept)
    )
  )

# Print the sample data with calculated copy numbers
print(head(sample_data))

# Separate the data for F9 and Ppia
f9_data <- sample_data %>% filter(Target == "F9")
ppia_data <- sample_data %>% filter(Target == "PPIA")

# Join F9 and Ppia data by Content and Day
normalized_data <- f9_data %>%
  left_join(ppia_data, by = c("Content", "Day"), suffix = c("_F9", "_Ppia")) %>%
  mutate(Normalized_F9 = Copies_F9 / Copies_Ppia)

# Reshape data from wide to long format for replicates
new_data_long <- normalized_data %>%
  pivot_longer(cols = c(Cq_F9, Cq_Ppia, Copies_F9, Copies_Ppia, Normalized_F9), names_to = "Metric", values_to = "Value") %>%
  separate(Content, into = c("Sample", "Replicate"), sep = "-", remove = FALSE)

# Summarize the data for plotting
summary_data <- new_data_long %>%
  group_by(Sample, Day) %>%
  summarize(
    Mean = mean(Value[Metric == "Normalized_F9"], na.rm = TRUE),
    SEM = sd(Value[Metric == "Normalized_F9"], na.rm = TRUE) / sqrt(n())
  )

# Adjust the color scale based on the actual unique values in the Content column
color_values <- c("CoSCR" = "#162F51", "CKPL" = "#CD6D66", "MB" = "#15783D")

# Plot the data using ggplot2
ggplot(summary_data, aes(x = as.numeric(gsub("Day", "", Day)), y = Mean, color = Sample)) +
  geom_line() +  # Add lines
  geom_point(size = 3) +  # Add points
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 1) + # Add error bars
  labs(x = "Day", y = "Normalized F9 Copy Number", title = "Normalized F9 Copy Number Across Time Points") +
  scale_x_continuous(breaks = seq(3, 30, by = 6), labels = paste("Day", c(3, 9, 16, 21, 30))) +
  scale_color_manual(values = color_values) +
  theme_minimal()

######################################################################################################

#Without ppia_data

# Load necessary libraries
rm(list = ls())
cat("\014") 
library(dplyr)
library(ggplot2)
library(tidyr)  # Load tidyr for pivot_longer function
library(ggbreak)
library(readr)  # Load readr for reading CSV files

# Function to calculate standard curve parameters
calculate_standard_curve <- function(log_copies, cq_values) {
  fit <- lm(cq_values ~ log_copies)
  slope <- coef(fit)[2]
  intercept <- coef(fit)[1]
  efficiency <- 10^(-1/slope)
  return(list(slope = slope, intercept = intercept, efficiency = efficiency))
}

# Function to calculate copies from Cq values using standard curve
calculate_copies <- function(cq_values, slope, intercept) {
  log_copies <- (cq_values - intercept) / slope
  copies <- 10^log_copies
  return(copies)
}

# Example Cq value from a blank sample (average Cq of blank)
Cq_blank <- 32.3

# Calculate copy number for blank sample
standard_data <- data.frame(
  log_copies = c(8.815644149, 7.815644149, 6.815644149, 5.815644149, 4.815644149, 3.815644149, 2.815644149, 1.815644149),
  cq_values = c(5.711062857, 10.0192069, 13.30810448, 16.56234853, 21.05404185, 24.0781887, 28.26, 31.92)
)

standard_curve <- calculate_standard_curve(standard_data$log_copies, standard_data$cq_values)
slope <- standard_curve$slope
intercept <- standard_curve$intercept
efficiency <- standard_curve$efficiency

# Calculate copy number for blank sample
copy_number_blank <- calculate_copies(Cq_blank, slope, intercept)

# Define the new Cq values for the treated samples
# Define the path to the CSV file
file_path <- "C:/Users/harum/OneDrive/Documents/LAB/mhpQPCR/qpcrcsv/timepoint_copynumber.csv"  # Replace with the actual path

# Read the CSV file
treated_data <- read_csv(file_path)

# Extract data for target F9 and average the same Content
f9_data <- treated_data %>%
  filter(Target == "F9") %>%
  group_by(Content) %>%
  summarize(across(starts_with("Day"), mean))

# Reshape the data to long format
treated_data_long <- f9_data %>%
  pivot_longer(cols = starts_with("Day"), names_to = "Day", values_to = "Cq")

# Calculate copies per reaction
treated_data_long <- treated_data_long %>%
  mutate(copies_rxn = round(calculate_copies(Cq, slope, intercept)),
         copies_per_uL_sample = round(copies_rxn / 2), # Assuming 2 uL volume added per reaction
         cDNA_copies_per_ng = round(copies_per_uL_sample)) # Adjust according to specific requirements

# Normalize to blank
treated_data_long <- treated_data_long %>%
  mutate(CopyNumber_adjusted = cDNA_copies_per_ng - round(copy_number_blank / 2)) # Adjusting blank normalization for cDNA copies per ng

# Reshape data from wide to long format
new_data_long <- treated_data_long %>%
  separate(Content, into = c("Sample", "Replicate"), sep = "-")

# Calculate mean and standard error for each time point
new_data_summary <- new_data_long %>%
  group_by(Sample, Day) %>%
  summarize(Mean = mean(CopyNumber_adjusted, na.rm = TRUE),
            SEM = sd(CopyNumber_adjusted, na.rm = TRUE) / sqrt(n()))



custom_colors <- c("CKPL" = "#CC6D66", "CoSCR" = "#152F50", "MB" = "#12783D")

# Plot the data based on days
# Plot the data based on days with y-axis break
ggplot(new_data_summary, aes(x = as.numeric(gsub("Day", "", Day)), y = Mean, color = Sample)) +
  geom_line() +  # Add lines
  geom_point(size=3) +  # Add points
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 1) + # Add error bars
  labs(x = "Day", y = "Copy number per ng of mRNA", title = "Copy Number Across Time Points") +
  scale_x_continuous(breaks = seq(3, 30, by = 6), labels = paste("Day", c(3, 9, 16, 21, 30))) +
  scale_color_manual(values = custom_colors) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "top",
    #panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
    #panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "gray95")
  )
  #scale_y_break(c(4000, 7500), scales = c(0.6, 0.2))
# Summarize the data to get mean and standard deviation for each content and day
summary_data <- treated_data_long %>%
  group_by(Content, Day) %>%
  summarise(
    n = n(),
    mean_copy_number = mean(CopyNumber_adjusted, na.rm = TRUE),
    sd_copy_number = sd(CopyNumber_adjusted, na.rm = TRUE),
    se = sd_copy_number / sqrt(n)
  )

summary_data <- summary_data %>%
  mutate(
    lower = mean_copy_number - se,
    upper = mean_copy_number + se
  )

# Define your custom color palette
#custom_colors <- c("CKPL-1" = "#CC6D66", "CKPL-2" = "#CC6D10", "CoSCR-1" = "#152F50", "CoSCR-2" = "#152D90", "MB-1" = "#12789D", "MB-2" = "#12783D")
ggplot(summary_data, aes(x = as.numeric(gsub("Day", "", Day)), y = mean_copy_number, color = Content)) +
  geom_line() +  # Add lines
  geom_point(size=3) +  # Add points
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1) + # Add error bars
  labs(x = "Day", y = "Copy number per ng of mRNA", title = "Copy Number Across Time Points") +
  scale_x_continuous(breaks = seq(9, 30, by = 6), labels = paste("Day", c(9, 16, 21, 30))) +
  scale_color_manual(values = custom_colors) +
  theme_minimal()
  #scale_y_break(c(4000, 7500), scales = 0.5)


#############F9copy number with Ppia normalization#############
# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggbreak)

# Function to calculate standard curve parameters
calculate_standard_curve <- function(log_copies, cq_values) {
  fit <- lm(cq_values ~ log_copies)
  slope <- coef(fit)[2]
  intercept <- coef(fit)[1]
  efficiency <- 10^(-1/slope)
  return(list(slope = slope, intercept = intercept, efficiency = efficiency))
}

# Function to calculate copies from Cq values using standard curve
calculate_copies <- function(cq_values, slope, intercept) {
  log_copies <- (cq_values - intercept) / slope
  copies <- 10^log_copies
  return(copies)
}

# Standard curve data for F9
f9_standard_data <- data.frame(
  log_copies = c(8.815644149, 7.815644149, 6.815644149, 5.815644149, 4.815644149, 3.815644149, 2.815644149, 1.815644149),
  cq_values = c(5.711062857, 10.0192069, 13.30810448, 16.56234853, 21.05404185, 24.0781887, 28.26, 31.92)
)

# Calculate F9 standard curve parameters
f9_standard_curve <- calculate_standard_curve(f9_standard_data$log_copies, f9_standard_data$cq_values)
f9_slope <- f9_standard_curve$slope
f9_intercept <- f9_standard_curve$intercept
f9_efficiency <- f9_standard_curve$efficiency

# Load the data
file_path <- "C:/Users/harum/OneDrive/Documents/LAB/mhpQPCR/qpcrcsv/timepoint_copynumber.csv"
sample_data <- read_csv(file_path)

# Reshape the data to long format
sample_data_long <- sample_data %>%
  pivot_longer(cols = starts_with("Day"), names_to = "Day", values_to = "Cq")

# Calculate copies per reaction for F9 and relative expression levels for Ppia
sample_data_long <- sample_data_long %>%
  mutate(
    Copies = case_when(
      Target == "F9" ~ calculate_copies(Cq, f9_slope, f9_intercept),
      Target == "PPIA" ~ 2^(-Cq)
    )
  )

# Aggregate data to ensure unique Content and Day combinations
f9_data <- sample_data_long %>% filter(Target == "F9") %>% group_by(Content, Day) %>% summarize(Copies_F9 = mean(Copies, na.rm = TRUE), .groups = 'drop')
ppia_data <- sample_data_long %>% filter(Target == "PPIA") %>% group_by(Content, Day) %>% summarize(Relative_Ppia = mean(Copies, na.rm = TRUE), .groups = 'drop')

# Join F9 and Ppia data by Content and Day
normalized_data <- f9_data %>%
  left_join(ppia_data, by = c("Content", "Day")) %>%
  mutate(Normalized_F9 = Copies_F9 / Relative_Ppia)

# Summarize the data for plotting
summary_data <- normalized_data %>%
  group_by(Content, Day) %>%
  summarize(
    Mean = mean(Normalized_F9, na.rm = TRUE),
    SEM = sd(Normalized_F9, na.rm = TRUE) / sqrt(n())
  )

# Define custom color palette
custom_colors <- c("CKPL" = "#CD6D66", "MB" = "#15783D", "CoSCR" = "#162F51")

# Plot the data using ggplot2
ggplot(summary_data, aes(x = as.numeric(gsub("Day", "", Day)), y = Mean, color = Content)) +
  geom_line() +  # Add lines
  geom_point(size = 3) +  # Add points
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 1) + # Add error bars
  labs(x = "Day", y = "Normalized F9 Copy Number", title = "Normalized F9 Copy Number Across Time Points") +
  scale_x_continuous(breaks = seq(3, 30, by = 6), labels = paste("Day", c(3, 9, 16, 21, 30))) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  scale_y_break(c(2500, 7500), scales = 0.5)
