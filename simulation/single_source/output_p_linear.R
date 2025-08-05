options(bitmapType = "cairo")

library(ggplot2)
library(patchwork)

load("varying_p_single_source_linear.Rdata")

# Subset the data for each metric
df_beta_L1 <- subset(df_results, metric == "beta_L1")
df_beta_L2 <- subset(df_results, metric == "beta_L2")
df_eta_L1  <- subset(df_results, metric == "eta_L1")
df_eta_L2  <- subset(df_results, metric == "eta_L2")

# 1. Plot for Beta L1 Loss (x-axis is p)
p_beta_L1 <- ggplot(df_beta_L1, aes(x = p, y = value, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Beta L1 Loss", x = "p", y = "Loss") +
  theme_minimal()

# 2. Plot for Beta L2 Loss (x-axis is p)
p_beta_L2 <- ggplot(df_beta_L2, aes(x = p, y = value, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Beta L2 Loss", x = "p", y = "Loss") +
  theme_minimal()

# 3. Plot for Eta L1 Loss (x-axis is p)
p_eta_L1 <- ggplot(df_eta_L1, aes(x = p, y = value, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Eta L1 Loss", x = "p", y = "Loss") +
  theme_minimal()

# 4. Plot for Eta L2 Loss (x-axis is p)
p_eta_L2 <- ggplot(df_eta_L2, aes(x = p, y = value, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Eta L2 Loss", x = "p", y = "Loss") +
  theme_minimal()

# Arrange the four plots in a 2x2 grid using patchwork
combined_plot <- (p_beta_L1 | p_beta_L2) / (p_eta_L1 | p_eta_L2)
combined_plot


library(writexl)
load("varying_p_single_source_linear.Rdata")
df_results$value <- round(df_results$value, 4)
write_xlsx(df_results, "output/single_linear_p.xlsx")



