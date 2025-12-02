# FWERvsN.r
# Code for plotting FWER vs Sample Size boxplot for UseCC=TRUE
library(tidyverse)

print(getwd())

# Read the data
data <- read_csv("internalData/FWERvsN.csv")

# Filter for UseCC=TRUE and create the boxplot
data %>%
  filter(UseCC == TRUE) %>%
  ggplot(aes(x = factor(SampleSize), y = FWER)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(
    title = "FWER vs Sample Size (UseCC = TRUE)",
    x = "Sample Size",
    y = "FWER"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Filter for UseCC=FALSE and create the boxplot
data %>%
  filter(UseCC == FALSE) %>%
  ggplot(aes(x = factor(SampleSize), y = FWER)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(
    title = "FWER vs Sample Size (UseCC = FALSE)",
    x = "Sample Size",
    y = "FWER"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Panel plot with both UseCC=TRUE and UseCC=FALSE (one below the other)
data %>%
  mutate(UseCC = ifelse(UseCC, "UseCC = TRUE", "UseCC = FALSE")) %>%
  ggplot(aes(x = factor(SampleSize), y = FWER)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ UseCC, ncol = 1) +
  labs(
    title = "FWER vs Sample Size by UseCC",
    x = "Sample Size",
    y = "FWER"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

# Panel plot by Selection variable (one below the other)
data %>%
  mutate(Selection = ifelse(Selection, "Selection = TRUE", "Selection = FALSE")) %>%
  ggplot(aes(x = factor(SampleSize), y = FWER)) +
  geom_boxplot(fill = "coral", alpha = 0.7) +
  facet_wrap(~ Selection, ncol = 1) +
  labs(
    title = "FWER vs Sample Size by Selection",
    x = "Sample Size",
    y = "FWER"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

# Panel plot by info_frac variable (one below the other)
data %>%
  ggplot(aes(x = factor(SampleSize), y = FWER)) +
  geom_boxplot(fill = "darkgreen", alpha = 0.7) +
  facet_wrap(~ info_frac, ncol = 1) +
  labs(
    title = "FWER vs Sample Size by Information Fraction",
    x = "Sample Size",
    y = "FWER"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )
