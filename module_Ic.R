# ===================================================================
#               AI and Biotechnology / Bioinformatics
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Research Internship (2025)
# -------------------------------------------------------------------
#                Module I: Getting Started with R
# -------------------------------------------------------------------
# ===================================================================

# Topic: Syntax in R
# Includes:
#   1. Variables
#   2. Comments
#   3. Keywords

# --------------------------------------------------------------------------------------------------
# Syntax are the predefined structure and rules of any programming language
# It shows us how to right codes right way 
# Helps to write clear and error free codes &
# easy to debug and fix errors
# R have three syntax: Variables, Comments & Keywords 
# --------------------------------------------------------------------------------------------------

#### 1. VARIABLES ####
# Variables store values in R.
# Use the assignment operator `<-` to assign variable

# Example: store a gene name
# Use quotes (" ") or (' ') to write string data e.g "TP53", "Gene"
# variable name should not be in quotes

gene <- "TP53"
# Explanation:
# gene is the variable name
# <- is an assignment operator that stores the gene name "TP53" into variable
# "TP53" is now saved in variable gene

# Retrieve the value in the console
gene # type variable name and hit ctrl + entr 
print(gene) # or use print () function
# Store numeric values in one variable (vector)
expression_levels <- c(2.3, 4.6, 3.6, 7.2, 4.7)
# Import a CSV file as a variable
raw_data <- read.csv(file.choose())
raw_data
#### Rules for Variable Naming ####

# Must start with a letter
# 1gene <- 25  # Invalid, starts with a number
gene1 <- 25    # Valid, number at the end
# No spaces in names
sample id <- "s01"  # Invalid: contains space
# Use underscore (_) or dot (.) instead
sample_id <- "s01"
sample.id <- "s02"
# R is case-sensitive
Glucose_level <- 110
glucose_level <- 110  # Treated as a different variable
# Variables can be overwritten without warning
glucose_level <- c(110, 90, 120)
# Example of variable overwriting in a data
# Removing Columns

# First create a copy of raw_data to preserve the original dataset
data <- raw_data
# Remove the 'patient_id' column from the original raw_data
# This overwrites raw_data permanently, and the column cannot be restored
# unless a copy (like 'data') was created beforehand
raw_data$patient_id <- NULL
# Alternatively, create a new dataset without the first column
# This keeps the original 'data' intact and stores the modified version in clean_data
clean_data <- data[, -1]
# --------------------------------------------------------------------------------------------------

#### 2. COMMENTS ####

# Comments improve code readability.
# They are ignored during execution.
# This is only for your own understanding
# R doesn’t consider it as a code
# It starts with # and anything after that is ignored by R when the code runs
# If your code begin with # R treats the entire line as a comment and skips it completely

# Example:
# data_2 <- 23
data_2 <- 23
# Comments can be used as headings:
#### Heading 1 ####
#### Heading 2 ####

# Best practice: always add comment before a code or at the end for clarity”
# --------------------------------------------------------------------------------------------------

#### 3. KEYWORDS ####

# Reserved words have predefined functions and cannot be used as variable names.
# Examples: if, else, TRUE, FALSE, for

help("reserved")  # Check the list of reserved words
# You can check the documentation for any R function using help(function_name) or ?function_name
help(mean) # check other R function
?median
# --------------------------------------------------------------------------------------------------

#### Sorting Data with TRUE & FALSE keywords ####

# Sort 'age' from largest to smallest
sorted_age <- sort(raw_data$age, decreasing = TRUE)
sorted_age
# Sort 'age' from smallest to largest
sorted_age2 <- sort(raw_data$age, decreasing = FALSE)
sorted_age2
# TRUE and FALSE are logical constants in R
# They are written in uppercase without quotation marks

# --------------------------------------------------------------------------------------------------

#### if & else: Logical Conditions ####
gene_expression <- 50
# if statement:
# Executes the code inside the curly braces only if the condition is TRUE
# Here, the condition checks if gene_expression is greater than 50

# If TRUE, prints the message "Gene is highly expressed"
if (gene_expression > 30) {
  print("Gene is highly expressed")
}
# If FALSE, nothing happens because there is no else block
if (gene_expression > 50) {
  print("Gene is highly expressed")
}
# if-else statement:
# Executes the first block if the condition is TRUE
# If the condition is FALSE, executes the code in the else block

if (gene_expression > 50) {
  print("Gene is highly expressed")
} else {
  print("Gene expression is low")
}
# Invalid: cannot use 'if' as a variable name
# if <- 28
# --------------------------------------------------------------------------------------------------

#### for Loop: Automating Tasks ####


# Example: convert character columns to factors

# Check structure of raw_data
str(raw_data)
# Create a copy to avoid modifying the original data
clean_data <- raw_data
