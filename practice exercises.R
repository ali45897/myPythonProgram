# Practice Exercises 

# ----------------------------------------------------------------------------------------------------------------

# 1. Check Cholesterol level (using if) 
# Write an If statement to check cholesterol level is greater than 240, 
# if true, it will prints “High Cholesterol”
cholesterol <- 230
if(cholesterol>240){
  print("High Cholesterol")
}
# ----------------------------------------------------------------------------------------------------------------

# 2. Blood Pressure Status (using if...else)
# Write an if…else statement to check if blood pressure is normal.
# If it’s less than 120, print: “Blood Pressure is normal”
# If false then print: “Blood Pressure is high”
Systolic_bp <- 130
if(Systolic_bp < 120){
  print("Blood Pressure is normal")
}else{
  print("Blood Pressure is high")
}
# ----------------------------------------------------------------------------------------------------------------

# 3. Automating Data Type Conversion with for loop

# Use patient_info.csv data and metadata.csv
# Perform the following steps separately on each dataset (patient_info.csv data and metadata.csv)
# Create a copy of the dataset to work on.
# Identify all columns that should be converted to factor type.
# Store their names in a variable (factor_cols).
#### patient data ####
pdata<- read.csv(file.choose())
# check the str of data
str(pdata)
# Create a copy to avoid modifying the original data
clean_pdata <- pdata
str(clean_pdata)
factor_pcols<- c("gender", "diagnosis", "smoker")# Store their names in a variable (factor_cols).
for (col in factor_pcols) {
   clean_pdata[[col]] <- as.factor(clean_pdata[[col]])  # Replace 'data' with the name of your dataset
}
str(clean_pdata)
str(pdata)
#### meta Data ####
mdata <- read.csv(file.choose())
# check the str of data
str(mdata)
# Create a copy to avoid modifying the original data
clean_mdata <- mdata
str(clean_mdata)

factor_mcols<- c("height", "gender")# Store their names in a variable (factor_cols).
for (col in factor_mcols) {
  clean_mdata[[col]] <- as.factor(clean_mdata[[col]])  # Replace 'data' with the name of your dataset
}
str(clean_mdata)
str(mdata)
# ----------------------------------------------------------------------------------------------------------------

# 4. Converting Factors to Numeric Codes

# Choose one or more factor columns (e.g., smoking_status).
# Convert "Yes" to 1 and "No" to 0 using a for loop.

# Hint:
binary_cols <- c("gender")   # store column names in a vector
# use ifelse() condition inside the loop to replace Yes with 1 and No with 0
 for (col in binary_cols) {
   clean_mdata[[col]] <-  ifelse( clean_mdata[[col]] == "Male", 1,0)# insert your ifelse() code here
 }
str(clean_mdata)
binary_pcols <- c("smoker")
for (col in binary_pcols){
  clean_pdata[[col]]<- ifelse(clean_pdata[[col]] == "Yes", 1,0)
}

# ----------------------------------------------------------------------------------------------------------------

#  Verification:
#    Compare the original and modified datasets to confirm changes.
str(clean_pdata)
str(pdata)

# ----------------------------------------------------------------------------------------------------------------
