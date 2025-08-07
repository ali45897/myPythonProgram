# 1. Create necessary subfolders
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")
data= read.csv(file.choose())

install.packages("dplyr")
install.packages("readr")
# 3. Load required libraries
library(readr)
library(dplyr)
# 5. Inspect the structure
cat("------ Dataset Structure ------\n")
str(data)
cat("------ Dataset Summary ------\n")
summary(data)
cat("------ First Few Rows ------\n")
print(head(data))
# Convert variables to appropriate types
data$gender <- as.factor(data$gender)
data$diagnosis <- as.factor(data$diagnosis)
data$smoker <- as.factor(data$smoker)
# Create new binary smoker variable
data$smoker_binary <- ifelse(data$smoker == "Yes", 1, 0)
data$smoker_binary <- as.factor(data$smoker_binary)
# Save the cleaned data
write.csv(data, "clean_data/patient_info_clean.csv", row.names = FALSE)
# Save all objects in the environment
save.image(file = "results/MaryamAhmedAli_class_Ib_Assignment.RData")

