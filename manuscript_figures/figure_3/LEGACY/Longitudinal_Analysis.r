# Set up library paths and load libraries
myPaths <- .libPaths()
myPaths = c("/home/shivani/miniconda3/envs/cairo_env///lib/R/library", myPaths)
.libPaths(myPaths)
library(reticulate)
use_python("/home/shivani/miniconda3/envs/cairo_env//bin/python", required = TRUE)
use_condaenv("cairo_env", required = TRUE)
py_config()

dyn.load('/home/shivani/miniconda3/envs/cairo_env//lib/libhdf5_hl.so.310.0.3')
set.seed(665)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(qs)
library(BiocParallel)
library(data.table)
library(IRdisplay)
library(scales)
library(readxl)
library(tidyr)
library(stringr)
library(ggbeeswarm)

################# Initialization

result_dir <- "/data/work/Projects/MMRF/longitudinal/Analysis/StandardRisk/"

################# Reading data
cat("\nReading data ...\n")

counts <- readRDS("/data/shiny/shiny-server/apps/MMRF/data/mmrf_umi_counts.RDS")
metadata <- readRDS("/data/work/Projects/MMRF/data/scRNAseq/1110_per_cell_md.rds")
compartments <- readRDS("/data/work/Projects/MMRF/data/scRNAseq/1110_per_cell_md_refined_compartments.rds")
Clusters <- read.csv("/data/work/Projects/MMRF/longitudinal/Risk_DifferentialAbundance_DirichletReg_Cluster_Results.TSV - Risk_DifferentialAbundance_DirichletReg_Cluster_Results.TSV.csv")

metadata_info <- read_excel("/data/work/Projects/MMRF/longitudinal/616_metadata_by_sampleid.xlsx")

# Assuming we have two dataframes: 'compartment' and 'Clusters'
# 'subcluster_V03072023' is in 'compartment' and 'Cell type' is in 'Clusters'

# Rename the 'Cell type' column in Clusters to match the 'subcluster_V03072023' column in compartment
Clusters <- rename(Clusters, subcluster_V03072023 = `Cell.type`)

# Merge the dataframes
compartment_merged <- merge(compartments, Clusters[, c("subcluster_V03072023", "Description")], by = "subcluster_V03072023", all.x = TRUE)

# 'compartment_merged' will now have the 'Description' column from 'Clusters' added to it.
# Convert factors to characters (if they are factors)
compartment_merged$subcluster_V03072023 <- as.character(compartment_merged$subcluster_V03072023)
compartment_merged$Description <- as.character(compartment_merged$Description)

# Populate 'Description' with 'subcluster_V03072023' values where 'Description' is NA or empty
compartment_merged$Description <- ifelse(is.na(compartment_merged$Description) | compartment_merged$Description == "", compartment_merged$subcluster_V03072023, compartment_merged$Description)

MMRF_compartments <- compartment_merged[, c("aliquot_id", "sample_id","d_visit_specimen_id","public_id",'visit_type','VISIT','VJ_INTERVAL','VISITDY',"subcluster_V03072023","Description","RefinedCompartments", "collection_event","davies_based_risk","doublet_pred","predicted_cell_type","general_ct","compartment","seurat_clusters","siteXbatch",
                                      "d_tx_induction_cat")]

MMRF_compartments = MMRF_compartments[MMRF_compartments$davies_based_risk %in% c("high_risk", "standard_risk"), ]

# List of public_ids to change to high_risk
changed_risk <- c("MMRF_1371", "MMRF_1650", "MMRF_2126", "MMRF_2523")
change_to_standard_risk <- "MMRF_1424"
# Update the 'davies_based_risk' column for these IDs to 'high_risk'
MMRF_compartments$davies_based_risk[MMRF_compartments$public_id %in% changed_risk] <- 'high_risk'
MMRF_compartments$davies_based_risk[MMRF_compartments$public_id == change_to_standard_risk] <- 'standard_risk'

# Verify the changes (optional)
MMRF_compartments[MMRF_compartments$public_id %in% change_to_standard_risk, ]

length(unique(MMRF_compartments$sample_id))

result_dir <- "/data/work/Projects/MMRF/longitudinal/Analysis/StandardRisk/"
# Assuming MMRF_compartments is your data frame

# Define the result directory (ensure this directory exists or create it)
# result_dir <- "path/to/your/result_dir"  # Adjust this to your actual directory path

# Count unique samples per patient by public_id and davies_based_risk
sample_counts_per_patient_risk <- MMRF_compartments %>%
  group_by(public_id, davies_based_risk) %>%
  summarise(samples_count = n_distinct(sample_id), .groups = 'drop') %>%
  ungroup() %>%
  mutate(label_position = samples_count)  # Position for the label

# Plot histogram of the counts, divided by davies_based_risk
plot <- ggplot(sample_counts_per_patient_risk, aes(x = samples_count, fill = davies_based_risk)) +
  geom_histogram(aes(y = ..count..), binwidth = 1, position = "dodge", color = "black") +
  geom_text(aes(label = ..count.., y = ..count.. + 0.3), stat = "count", position = position_dodge(width = 1), vjust = 0, size=5.5) +
  theme_minimal() +
  labs(title = "Histogram of Samples Collected Per Patient \nby Davies Based Risk\n",
       x = "\nNumber of Samples Collected",
       y = "Number of Patients") +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face='bold'),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5),  # Center x-axis labels
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 18)) +
  scale_fill_manual(values = c("high_risk" = "coral", "standard_risk" = "cornflowerblue"), name = "Davies Based Risk")

# Display the plot
print(plot)

# Save the plot in the result_dir
plot_filename <- paste0(result_dir, "/samples_per_patient_by_risk_histogram_with_counts.png")  # Define the file name
ggsave(plot_filename, plot, width = 12, height = 8, dpi = 300)

plot_filename <- paste0(result_dir, "/samples_per_patient_by_risk_histogram_with_counts.svg")  # Define the file name
ggsave(plot_filename, plot, width = 12, height = 8, dpi = 300)

# Step 1: Identify and print patients affected by the change
affected_patients <- MMRF_compartments %>%
  group_by(public_id) %>%
  filter(all(c("Relapse/Progression", "Remission/Response") %in% collection_event)) %>%
  distinct(public_id)

print("Patients affected by the change (having both relapse and remission samples):")
print(affected_patients)

# Step 2: Prioritize Relapse Samples if both Relapse and Remission exist for a patient
MMRF_compartments <- MMRF_compartments %>%
  group_by(public_id) %>%
  mutate(collection_event = ifelse(
    all(c("Relapse/Progression", "Remission/Response") %in% collection_event) & collection_event == "Remission/Response",
    "Relapse/Progression",
    collection_event
  )) %>%
  ungroup()

# Step 3: Filter and Count Patients by Collection Event and Risk
patients_count <- MMRF_compartments %>%
  group_by(public_id, davies_based_risk) %>%
  # Check for Baseline and Relapse/Progression
  mutate(has_baseline_relapse = all(c("Baseline", "Relapse/Progression") %in% collection_event),
         # Check for Baseline and Remission/Response
         has_baseline_remission = all(c("Baseline", "Remission/Response") %in% collection_event)) %>%
  summarise(has_baseline_relapse = any(has_baseline_relapse),
            has_baseline_remission = any(has_baseline_remission),
            .groups = 'drop') %>%
  # Reshape for plotting
  pivot_longer(cols = starts_with("has"),
               names_to = "CollectionType",
               values_to = "Count") %>%
  filter(Count) %>%
  group_by(davies_based_risk, CollectionType) %>%
  summarise(PatientCount = n(), .groups = 'drop')

# Replace logical names with descriptive names
patients_count$CollectionType <- recode(patients_count$CollectionType,
                                        has_baseline_relapse = "Baseline & \nRelapse/\nProgression\n",
                                        has_baseline_remission = "Baseline & \nRemission/\nResponse\n")

# Step 4: Create the Plot
plot <- ggplot(patients_count, aes(x = CollectionType, y = PatientCount, fill = davies_based_risk)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +  # Adjust the dodge width
  geom_text(aes(label = PatientCount), position = position_dodge(width = 0.9), vjust = -0.250, size = 5.5) +  # Adjust the dodge width
  theme_minimal() +
  labs(title = "Number of Patients \nby Collection Event Combinations and Risk\n",
       x = "Collection Event Combination",
       y = "Number of Patients",
       fill = "Davies Based Risk") +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5),  # Center x-axis labels
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 18)) +
  scale_fill_manual(values = c("high_risk" = "coral", "standard_risk" = "cornflowerblue"), name = "Davies Based Risk")

# Display the plot
print(plot)

# If you need to save the plot
ggsave(filename = paste0(result_dir, "/patient_counts_by_event_risk_combinations.png"), plot = plot, width = 8, height = 8, dpi = 300)
ggsave(filename = paste0(result_dir, "/patient_counts_by_event_risk_combinations.svg"), plot = plot, width = 8, height = 8, dpi = 300)

# Filter data for patients with more than one sample
patients_sample_counts <- MMRF_compartments %>%
  group_by(public_id) %>%
  summarise(samples_count = n_distinct(sample_id)) %>%
  ungroup()

patients_with_multiple_samples <- patients_sample_counts %>%
  filter(samples_count > 1)

subset_MMRF_compartments <- MMRF_compartments %>%
  filter(public_id %in% patients_with_multiple_samples$public_id)

# Filter data to include only samples from public_ids with at least one "Relapse/Progression" event
public_ids_with_relapse <- MMRF_compartments %>%
  filter(collection_event == "Relapse/Progression") %>%
  distinct(public_id)

all_samples_from_relevant_patients <- subset_MMRF_compartments %>%
  filter(public_id %in% public_ids_with_relapse$public_id)

# Count the unique sample_ids from the filtered data
unique_samples_count <- all_samples_from_relevant_patients %>%
  summarise(number_of_samples = n_distinct(sample_id))

print(unique_samples_count)

risk_based_sample_counts <- all_samples_from_relevant_patients %>%
  group_by(davies_based_risk) %>%
  summarise(number_of_samples = n_distinct(sample_id))

print(risk_based_sample_counts)

cell_proportions <- subset_MMRF_compartments %>%
  group_by(sample_id, RefinedCompartments, Description) %>%
  summarise(count = n(), .groups = "keep") %>%
  ungroup() %>%
  group_by(sample_id) %>%
  mutate(total_count = sum(count)) %>%
  group_by(sample_id, RefinedCompartments) %>%
  mutate(proportion = ifelse(RefinedCompartments == "Fibroblast",
                             count / total_count,  # For Fibroblast, use total count across all compartments
                             count / sum(count))) %>%  # For other compartments, use total count within the compartment
  left_join(subset_MMRF_compartments, by = c("sample_id", "RefinedCompartments", "Description"))

print(head(cell_proportions))

cell_proportions[(cell_proportions$RefinedCompartments=='Fibroblast'),]

cell_proportions=cell_proportions[!cell_proportions$collection_event %in% c("Other"),]

timepoint_dir <- paste0(result_dir, "/TimepointAnalysis")

if(!dir.exists(timepoint_dir)) {
  dir.create(timepoint_dir, recursive = TRUE)
}

convert_to_months <- function(interval, days) {
  if (interval == "Baseline") {
    return(0)
  } else if (str_detect(interval, "^Month ")) {
    return(as.numeric(str_replace(interval, "Month ", "")))
  } else if (str_detect(interval, "^Year ")) {
    return(days / 30.44)
  } else {
    return(NA)  # Return NA for unrecognized formats
  }
}

high_risk_data <- cell_proportions %>%
  filter(davies_based_risk == "high_risk") %>%
  mutate(VJ_INTERVAL = mapply(convert_to_months, VJ_INTERVAL, VISITDY))

standard_risk_data <- cell_proportions %>%
  filter(davies_based_risk == "standard_risk") %>%
  mutate(VJ_INTERVAL = mapply(convert_to_months, VJ_INTERVAL, VISITDY))

combined_data <- bind_rows(high_risk_data %>% mutate(Risk_Group = "High Risk"),
                           standard_risk_data %>% mutate(Risk_Group = "Standard Risk"))

unique_combinations <- unique(high_risk_data[, c("RefinedCompartments", "Description")])

for(i in 1:nrow(unique_combinations)) {
  combination <- unique_combinations[i, ]
  combination_data <- high_risk_data %>%
    filter(RefinedCompartments == combination$RefinedCompartments,
           Description == combination$Description) %>%
    group_by(public_id) %>%
    mutate(num_visit_intervals = n_distinct(VJ_INTERVAL)) %>%
    filter(num_visit_intervals > 1) %>%
    ungroup()

  if(nrow(combination_data) > 0) {
    plot <- ggplot(combination_data, aes(x = VJ_INTERVAL, y = proportion, group = public_id, color = as.factor(public_id))) +
      geom_line() +
      theme_minimal() +
      labs(title = paste("Standard Risk - Cell Proportion Over Time for", combination$RefinedCompartments, "-", combination$Description),
           x = "Visit Interval (Months)",
           y = "Proportion of Cells",
           color = "Public ID") +
      theme(legend.position = "bottom")

    filename <- paste("StandardRisk", combination$RefinedCompartments, combination$Description, "timepoint_plot.png", sep = "_")
    filepath <- paste0(timepoint_dir, "/", filename)

    tryCatch({
      ggsave(filepath, plot, width = 10, height = 6)
    }, error = function(e) {
      message("Failed to save plot for ", combination$RefinedCompartments, " - ", combination$Description, ": ", e$message)
    })
  } else {
    message("Insufficient data to plot for High Risk - ", combination$RefinedCompartments, " - ", combination$Description)
  }
}

# Define colors for risk groups and collection events
risk_colors <- c("High Risk" = "coral", "Standard Risk" = "cornflowerblue")
event_colors <- c("Baseline" = "green", "Relapse/Progression" = "purple", "Remission/Response" = "midnightblue")

for(i in 1:nrow(unique_combinations)) {
  combination <- unique_combinations[i, ]
  combination_data <- combined_data %>%
    filter(RefinedCompartments == combination$RefinedCompartments,
           Description == combination$Description) %>%
    group_by(public_id) %>%
    mutate(num_visit_intervals = n_distinct(VJ_INTERVAL)) %>%
    filter(num_visit_intervals > 1) %>%
    ungroup()

  if(nrow(combination_data) > 0) {
    plot <- ggplot(combination_data, aes(x = VJ_INTERVAL, y = proportion, group = public_id)) +
      geom_line(aes(color = Risk_Group)) +
      geom_point(aes(color = collection_event, shape = collection_event)) +
      scale_color_manual(values = c(risk_colors, event_colors)) +
      theme_minimal() +
      labs(title = paste("Cell Proportion Over Time for", combination$RefinedCompartments, "-", combination$Description),
           x = "Visit Interval (Months)",
           y = "Proportion of Cells") +
      theme(legend.position = "bottom")

    filename <- paste(combination$RefinedCompartments, combination$Description, "timepoint_plot.svg", sep = "_")
    filepath <- paste0(timepoint_dir, "/", filename)

    tryCatch({
      ggsave(filepath, plot, width = 10, height = 6)
    }, error = function(e) {
      message("Failed to save plot for ", combination$RefinedCompartments, " - ", combination$Description, ": ", e$message)
    })
  } else {
    message("Insufficient data to plot for ", combination$RefinedCompartments, " - ", combination$Description)
  }
}

# Ensure that collection_event is correctly used to filter for Baseline and Relapse/Progression events
public_ids_with_both_events <- combined_data %>%
  filter(collection_event %in% c("Baseline", "Relapse/Progression")) %>%
  group_by(public_id, davies_based_risk) %>%
  filter(n_distinct(collection_event) == 2) %>%
  ungroup() %>%
  select(public_id, davies_based_risk) %>%
  distinct() %>%
  pull(public_id)

# Filter the original dataset to keep only relevant public_ids
filtered_data <- combined_data %>%
  filter(public_id %in% public_ids_with_both_events)

# Calculate fold change with inclusion of davies_based_risk to separate calculations for high and standard risk
fold_change_data <- filtered_data %>%
  group_by(public_id, davies_based_risk, Description, RefinedCompartments) %>%
  summarise(
    first_proportion = first(proportion[collection_event == "Baseline"]),
    last_proportion = last(proportion[collection_event == "Relapse/Progression"]),
    log2_fold_change = ifelse(!is.na(first_proportion) & !is.na(last_proportion), log2(last_proportion / first_proportion), NA_real_),
    .groups = 'drop'
  ) %>%
  filter(!is.na(log2_fold_change))

p_values <- fold_change_data %>%
  group_by(Description, RefinedCompartments, davies_based_risk) %>%
  summarise(
    first_proportions = list(first_proportion),
    last_proportions = list(last_proportion),
    .groups = 'drop'
  ) %>%
  rowwise() %>%
  mutate(
    p_value = ifelse(
      all(unlist(first_proportions) == unlist(last_proportions)),
      NA,  # Assign NA if all pairs are the same
      wilcox.test(unlist(first_proportions), unlist(last_proportions), paired = TRUE)$p.value
    )
  ) %>%
  ungroup() %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "" # not significant
  ))

# Now, p_values includes Description, RefinedCompartments, davies_based_risk, p_value, and significance

unique_refined_compartments <- unique(fold_change_data$RefinedCompartments)

for(ref_comp in unique_refined_compartments) {
  comp_data1 <- filter(fold_change_data, RefinedCompartments == ref_comp)
  
  # Calculate the median fold change per Description and adjust data
  medians <- comp_data1 %>%
    group_by(Description) %>%
    summarise(median_value = median(log2_fold_change), .groups = 'drop')
  
  comp_data1 <- left_join(comp_data1, medians, by = "Description")
  comp_data_with_p <- left_join(comp_data1, p_values, by = c("Description", "RefinedCompartments"))
  
  annotations <- comp_data_with_p %>%
    group_by(Description) %>%
    summarise(significance = first(significance), max_y = max(log2_fold_change, na.rm = TRUE), .groups = 'drop')
  
  extra_space <- max(annotations$max_y) * 0.1
  annotations$max_y <- annotations$max_y + extra_space
  
  # Plot with wrapped strip text and adjusted settings
  p <- ggplot(comp_data_with_p, aes(x = Description, y = log2_fold_change, fill = median_value)) +
      geom_boxplot(width=0.9) +
      geom_beeswarm() +
      geom_text(data = annotations, aes(x = Description, label = significance, y = max_y), vjust = -1, size = 3, inherit.aes = FALSE) +
      scale_fill_gradient2(low = "red", high = "blue", mid = "grey", midpoint = 0, name = "Median Fold Change") +
      theme_minimal() +
      labs(title = paste("Log2 Fold Change for Refined Compartments:", ref_comp),
           x = "Description", y = "Log2 Fold Change") +
      facet_wrap(~ Description, scales = "free_y", ncol=5, labeller = label_wrap_gen(width = 20)) + # Wrap labels to a specified width
      theme(strip.text.x = element_text(size = 12, face = "bold"),
            strip.background = element_rect(fill = "lightblue"),
            panel.spacing = unit(1, "lines"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      expand_limits(y = max(annotations$max_y)) +
      coord_cartesian(ylim = c(NA, max(annotations$max_y) + extra_space))

  file_path <- paste0(result_dir, "/Log2_FoldChange_", gsub("[^A-Za-z0-9]", "_", ref_comp), ".png")
  ggsave(file_path, p, width = 15, height = 10)
}

specific_descriptions <- c("Treg")
comp_data_specific <- filter(fold_change_data, Description %in% specific_descriptions)

# Calculate the median fold change per Description AND Risk Group
medians <- comp_data_specific %>%
  group_by(Description, davies_based_risk) %>%
  summarise(median_value = median(log2_fold_change), .groups = 'drop')

# Adjust data by joining with medians and p-values
comp_data_with_p <- left_join(comp_data_specific, medians, by = c("Description", "davies_based_risk")) %>%
  left_join(p_values, by = c("Description", "RefinedCompartments", "davies_based_risk"))

# Prepare annotations if needed
annotations <- comp_data_with_p %>%
  group_by(Description, davies_based_risk) %>%
  summarise(significance = first(significance), max_y = max(log2_fold_change, na.rm = TRUE), .groups = 'drop')

extra_space <- max(annotations$max_y, na.rm = TRUE) * 0.1
annotations$max_y <- annotations$max_y + extra_space

# Plot with High Risk and Standard Risk Box Plots Side by Side
p <- ggplot(comp_data_with_p, aes(x = interaction(Description, davies_based_risk), y = log2_fold_change, fill = median_value)) +
    geom_boxplot(width=0.5) +
    geom_beeswarm() +
    geom_text(data = annotations, aes(x = interaction(Description, davies_based_risk), label = significance, y = max_y), vjust = -1, size = 9, inherit.aes = FALSE) +
    scale_fill_gradient2(low = "coral", high = "cornflowerblue", mid = "grey", midpoint = 0, name = "Median Log2FC") +
    theme_minimal() +
    labs(x = "Description and Risk Group", y = "Log2 Fold Change") +
    theme(strip.text.x = element_text(size = 16, face = "bold"),
          strip.background = element_rect(fill = "lightblue"),
          panel.spacing = unit(1, "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=16),
          axis.text.y = element_text(size=16, face='bold'),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    expand_limits(y = max(annotations$max_y, na.rm = TRUE)) +
    coord_cartesian(ylim = c(NA, max(annotations$max_y, na.rm = TRUE) + extra_space))

# Adjust the file path as needed
file_path <- paste0(result_dir, "/TReg_Log2_FoldChange_by_Risk.png")
ggsave(file_path, p, width = 8, height = 10)
file_path <- paste0(result_dir, "/TReg_Log2_FoldChange_by_Risk.svg")
ggsave(file_path, p, width = 8, height = 10)

# Assuming 'convert_to_months' and 'arcsine_transform' functions are defined as provided
arcsine_transform <- function(p) {
  asin(sqrt(p))
}

# Prepare the data for both high_risk and standard_risk
all_risk_data <- cell_proportions %>%
  filter(davies_based_risk %in% c("high_risk", "standard_risk")) %>%
  mutate(
    months = mapply(convert_to_months, VJ_INTERVAL, VISITDY),
    arcsine_proportion = arcsine_transform(proportion),
    risk_group = davies_based_risk
  )

# Initialize a dataframe to store results
results_df <- data.frame(Description = character(), 
                         Risk_Group = character(),
                         Intercept = numeric(), 
                         Slope = numeric(), 
                         P_Value = numeric(), 
                         R_Squared = numeric(),
                         stringsAsFactors = FALSE)

# Fit linear models and store results
all_risk_data %>%
  group_by(Description, risk_group) %>%
  do({
    lm_fit <- lm(arcsine_proportion ~ months, data = .)
    tidy_lm <- tidy(lm_fit)
    summary_lm <- summary(lm_fit)
    
    data.frame(
      Description = unique(.$Description),
      Risk_Group = unique(.$risk_group),
      Intercept = tidy_lm$estimate[1],
      Slope = tidy_lm$estimate[2],
      P_Value = tidy_lm$p.value[2],
      R_Squared = summary_lm$r.squared
    )
  }) %>%
  ungroup() -> results_df

# Adjust p-values for multiple comparisons (FDR)
results_df$FDR <- p.adjust(results_df$P_Value, method = "BH")

unique_descriptions <- unique(all_risk_data$Description)

plot_save_directory <- result_dir  # Update this path

for (desc in unique_descriptions) {
  desc_data <- filter(all_risk_data, Description == desc)
  
  # Prepare statistics summary text for High Risk and Standard Risk
  high_risk_stats <- filter(results_df, Description == desc, Risk_Group == "high_risk")
  standard_risk_stats <- filter(results_df, Description == desc, Risk_Group == "standard_risk")
  
  # Example of how to format the statistics summary text correctly
  stats_summary <- paste(
      "High Risk - Slope:", format(high_risk_stats$Slope, scientific = TRUE, digits = 2),
      ", R²:", round(high_risk_stats$R_Squared, 3),
      ", FDR:", format(high_risk_stats$FDR, scientific = TRUE, digits = 2),
      "\nStandard Risk - Slope:", format(standard_risk_stats$Slope, scientific = TRUE, digits = 2),
      ", R²:", round(standard_risk_stats$R_Squared, 3),
      ", FDR:", format(standard_risk_stats$FDR, scientific = TRUE, digits = 2),
      sep = ""
    )

  # Base plot setup
  plot <- ggplot(desc_data, aes(x = months, y = arcsine_proportion)) +
    labs(title = paste("Trend for", desc), x = "Months", y = "Arcsine Transformed Proportion") +
    theme_minimal() +
    scale_color_manual(values = c("high_risk" = "coral", "standard_risk" = "cornflowerblue")) +
    geom_smooth(data = filter(desc_data, risk_group == "high_risk"),
                aes(x = months, y = arcsine_proportion, color = "high_risk", fill = "high_risk"),
                method = "lm", formula = 'y ~ x', se = TRUE, alpha = 0.2) +
    geom_smooth(data = filter(desc_data, risk_group == "standard_risk"),
                aes(x = months, y = arcsine_proportion, color = "standard_risk", fill = "standard_risk"),
                method = "lm", formula = 'y ~ x', se = TRUE, alpha = 0.2) +
    scale_fill_manual(values = c("high_risk" = "coral", "standard_risk" = "cornflowerblue")) +
    annotate("text", x = Inf, y = -Inf, label = stats_summary, hjust = 1.05, vjust = -0.5, size = 3.5, color = "black")

  # Save the plot
  file_name1 <- paste0(gsub("[[:space:]]+", "_", desc), "_trend_plot.png")
  file_name2 <- paste0(gsub("[[:space:]]+", "_", desc), "_trend_plot.svg")
  if(!dir.exists(plot_save_directory)) {
    dir.create(plot_save_directory, recursive = TRUE)
  }
  ggsave(filename = file.path(plot_save_directory, file_name1), plot = plot, width = 10, height = 6)
  ggsave(filename = file.path(plot_save_directory, file_name2), plot = plot, width = 10, height = 6)
}
