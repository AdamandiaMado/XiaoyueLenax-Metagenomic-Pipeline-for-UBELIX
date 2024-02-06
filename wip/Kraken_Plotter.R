#================================

# This R script is used for visualizing results from kraken.
# The actual pipeline has not implemented Kraken due to database storage.


# Install ggplot2 package if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

# Replace 'your_kraken_report.txt' with the actual path to your Kraken database report. placeholder. 
kraken_report_path <- 'Human1.report.txt'
# Read the Kraken report into a data frame
kraken_data <- read.table(kraken_report_path, header = FALSE, sep = '\t', stringsAsFactors = FALSE)


# Rename columns
colnames(kraken_data) <- c("Percentage", "Count1", "Count2", "Taxonomic Rank", "Code", "Taxon")
print(kraken_data["Taxon"])

# Create a pie chart - SUPER SLOW
#output <- ggplot(kraken_data, aes(x = "", y = Percentage, fill = Taxon)) +
#  geom_bar(stat = "identity", width = 1) +
#  coord_polar(theta = "y") +
#  labs(title = "Kraken Database Report") +
#  theme_void() +
#  theme(legend.position = "right")

# Save the plot as a PNG file - takes sooooo long
#ggsave("kraken_pie_chart_human1.png", output, width = 10, height = 6, units = "in")

# bar plot instead?

# Filter taxa with abundances greater than 10%
filtered_data <- kraken_data[as.numeric(gsub("%", "", kraken_data$Percentage)) > 10, ]

# Group small abundances under 'others'
filtered_data$Taxon <- ifelse(as.numeric(gsub("%", "", filtered_data$Percentage)) <= 10, 'Others', filtered_data$Taxon)

# Convert Percentage column to numeric (removing % sign)
filtered_data$Percentage <- as.numeric(gsub("%", "", filtered_data$Percentage))

# Create a bar plot
ggplot(filtered_data, aes(x = reorder(Taxon, Percentage), y = Percentage, fill = Taxon)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundance of Taxa", x = "Taxon", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# Now do the same for human2  =====================================================================
human2 <- 'Human2.report.txt'
k2<- read.table(human2, header = FALSE, sep = '\t', stringsAsFactors = FALSE)


# Rename columns
colnames(k2) <- c("Percentage", "Count1", "Count2", "Taxonomic Rank", "Code", "Taxon")
print(k2["Taxon"])
filtered_k2 <- k2[as.numeric(gsub("%", "", k2$Percentage)) > 10, ]

# Group small abundances under 'others'
filtered_k2$Taxon <- ifelse(as.numeric(gsub("%", "", filtered_k2$Percentage)) <= 10, 'Others', filtered_k2$Taxon)

# Convert Percentage column to numeric (removing % sign)
filtered_k2$Percentage <- as.numeric(gsub("%", "", filtered_k2$Percentage))

# Create a bar plot
ggplot(filtered_k2, aes(x = reorder(Taxon, Percentage), y = Percentage, fill = Taxon)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundance of Taxa", x = "Taxon", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# Now do the same for human3  =====================================================================
human3 <- 'Human3.report.txt'
k3<- read.table(human3, header = FALSE, sep = '\t', stringsAsFactors = FALSE)


# Rename columns
colnames(k3) <- c("Percentage", "Count1", "Count2", "Taxonomic Rank", "Code", "Taxon")
print(k3["Taxon"])
filtered_k3 <- k3[as.numeric(gsub("%", "", k3$Percentage)) > 10, ]

# Group small abundances under 'others'
filtered_k3$Taxon <- ifelse(as.numeric(gsub("%", "", filtered_k3$Percentage)) <= 10, 'Others', filtered_k3$Taxon)

# Convert Percentage column to numeric (removing % sign)
filtered_k3$Percentage <- as.numeric(gsub("%", "", filtered_k3$Percentage))

# Create a bar plot
ggplot(filtered_k3, aes(x = reorder(Taxon, Percentage), y = Percentage, fill = Taxon)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundance of Taxa", x = "Taxon", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))