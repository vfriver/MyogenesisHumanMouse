# Load libraries
library(ggplot2)
library(ggpubr)

# df <- data.frame(log2TPMHuman_Mouse_data.csv) #upload the file average log2(TPM)

# Scatter plot with regression and correlation
ggplot(df, aes(x = human_val, y = mouse_val)) +
  
  # Each point represents one orthologous gene
  geom_point(alpha = 0.5, color = "steelblue") +
  
  # Linear regression line
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  
  # Pearson correlation coefficient
  stat_cor(method = "pearson") +
  
  # Identity line 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  
  # Clean theme
  theme_minimal() +
  
  # Axis labels and title
  labs(
    title = "Correlation of Gene Expression Between Human and Mouse",
    subtitle = "Each point represents a 1:1 orthologous gene",
    x = "Average log2(TPM + 1), Human",
    y = "Average log2(TPM + 1), Mouse"
  )
