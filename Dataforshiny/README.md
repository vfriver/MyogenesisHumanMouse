# This folder contains the datasets required to run the visualization tool. Depending on which stage of myogenesis you wish to observe, you must follow these instructions:

# 🛠 Data Loading Instructions

To explore different time points of the myogenesis dataset, you can easily switch the data source directly within the R code. 

### **How to Switch Datasets**

1. **Locate the files**: In the repository, you will find several `.csv` files corresponding to different stages (e.g., log2TPM_Day0.csv, log2TPM_Day1.csv, log2TPM_Day6.csv).
2. **Open the script**: Open the `app.R` file in RStudio.
3. **Modify the code**: Look for the data loading line at the top of the script and update the filename within the `read.csv()` function.

#### **Examples:**

* **To observe Day 0:**
  ```R
  tableHM <- read.csv("log2TPM_Day0.csv")
