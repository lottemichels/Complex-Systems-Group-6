# ============================
# Load Libraries
# ============================
library(rEDM)
library(tidyverse)
library(scales)
library(Kendall)

# ============================
# Load and Prepare Data
# ============================
data <- D1_sensor_data[10000:11000, ] %>% drop_na(temperature, humidity, hive_power)

# Normalize columns
data$temperature_norm <- as.numeric(scale(data$temperature))
data$humidity_norm <- as.numeric(scale(data$humidity))
data$hive_power_norm <- as.numeric(scale(data$hive_power))

# ============================
# Define Helper Functions
# ============================

run_simplex <- function(df, var_name, E_range = 2:8) {
  df_numeric <- data.frame(time = 1:nrow(df), value = df[[var_name]])
  lib_range <- c(1, 500)
  pred_range <- c(501, 1000)
  stats <- lapply(E_range, function(E) {
    result <- Simplex(dataFrame = df_numeric, lib = lib_range, pred = pred_range,
                      E = E, columns = "value", target = "value")
    result <- result[!is.na(result$Predictions), ]
    rho <- cor(result$Observations, result$Predictions)
    mae <- mean(abs(result$Observations - result$Predictions))
    rmse <- sqrt(mean((result$Observations - result$Predictions)^2))
    data.frame(E = E, rho = rho, MAE = mae, RMSE = rmse)
  })
  stats_df <- bind_rows(stats)
  return(stats_df)
}

run_smap <- function(df, target, columns, E, title = NULL) {
  lib_range <- c(1, 500)
  pred_range <- c(501, 1000)
  result <- SMap(dataFrame = df, E = E, lib = lib_range, pred = pred_range,
                 columns = columns, target = target, embedded = TRUE)
  result$predictions <- result$predictions[!is.na(result$predictions$Predictions), ]
  result$stats <- list(
    rho = cor(result$predictions$Observations, result$predictions$Predictions),
    MAE = mean(abs(result$predictions$Observations - result$predictions$Predictions)),
    RMSE = sqrt(mean((result$predictions$Observations - result$predictions$Predictions)^2))
  )
  return(result)
}

run_ccm <- function(df, E, lib_column, target_column, lib_sizes) {
  result <- CCM(dataFrame = df, E = E, libSizes = lib_sizes, sample = 100, random = TRUE,
                columns = lib_column, target = target_column)
  return(as.data.frame(result))
}

# ============================
# Run Simplex for Each Variable
# ============================
simplex_temp <- run_simplex(data, "temperature_norm")
simplex_hum <- run_simplex(data, "humidity_norm")
simplex_hive <- run_simplex(data, "hive_power_norm")

# Optimal E values
optimal_E_temp <- simplex_temp$E[which.max(simplex_temp$rho)]
optimal_E_hum <- simplex_hum$E[which.max(simplex_hum$rho)]
optimal_E_hive <- simplex_hive$E[which.max(simplex_hive$rho)]

# ============================
# Multivariate SMap
# ============================
df_ht <- data.frame(time = 1:nrow(data),
                    humidity_norm = data$humidity_norm,
                    temperature_norm = data$temperature_norm)

df_hh <- data.frame(time = 1:nrow(data),
                    humidity_norm = data$humidity_norm,
                    hive_power_norm = data$hive_power_norm)

smap_ht <- run_smap(df_ht, "humidity_norm", c("humidity_norm", "temperature_norm"), optimal_E_hum)
smap_hh <- run_smap(df_hh, "humidity_norm", c("humidity_norm", "hive_power_norm"), optimal_E_hum)

# ============================
# CCM
# ============================
ccm_df <- data.frame(
  humidity_norm = data$humidity_norm,
  temperature_norm = data$temperature_norm
)

lib_sizes <- seq(10, 500, by = 50)

ccm_ht <- run_ccm(ccm_df, optimal_E_hum, "humidity_norm", "temperature_norm", lib_sizes)
ccm_th <- run_ccm(ccm_df, optimal_E_temp, "temperature_norm", "humidity_norm", lib_sizes)

# ============================
# Multiview
# ============================
mv_out <- Multiview(
  dataFrame = df_ht,
  lib = c(1, 500),
  pred = c(501, 1000),
  E = optimal_E_temp,
  target = "temperature_norm",
  columns = c("humidity_norm", "temperature_norm")
)

mv_out$predictions <- mv_out$predictions[!is.na(mv_out$predictions$Predictions), ]
mv_stats <- list(
  rho = cor(mv_out$predictions$Observations, mv_out$predictions$Predictions),
  MAE = mean(abs(mv_out$predictions$Observations - mv_out$predictions$Predictions)),
  RMSE = sqrt(mean((mv_out$predictions$Observations - mv_out$predictions$Predictions)^2))
)

# ============================
# Final Output
# ============================
print("Simplex Results (Temperature):")
print(simplex_temp)
print("Simplex Results (Humidity):")
print(simplex_hum)
print("Simplex Results (Hive Power):")
print(simplex_hive)

print("SMap - Humidity + Temperature:")
print(smap_ht$stats)
print("SMap - Humidity + Hive Power:")
print(smap_hh$stats)

print("CCM - Humidity → Temperature:")
print(ccm_ht)
print("CCM - Temperature → Humidity:")
print(ccm_th)

print("Multiview Stats:")
print(mv_stats)


# Humidity ➝ Temperature
ccm_ht <- CCM(
  dataFrame = ccm_df,
  E = optimal_E_hum,
  libSizes = lib_sizes,
  sample = 100,
  columns = "humidity_norm",
  target = "temperature_norm",
  random = TRUE
)




# Humidity ➝ Temperature
ccm_ht <- CCM(
  dataFrame = ccm_df,
  E = optimal_E_hum,
  libSizes = lib_sizes,
  sample = 100,
  columns = "humidity_norm",
  target = "temperature_norm",
  random = TRUE
)

# Temperature ➝ Humidity
ccm_th <- CCM(
  dataFrame = ccm_df,
  E = optimal_E_temp,
  libSizes = lib_sizes,
  sample = 100,
  columns = "temperature_norm",
  target = "humidity_norm",
  random = TRUE
)


df_ht <- as.data.frame(ccm_ht) %>%
  group_by(LibSize) %>%
  summarize(rho = mean(rho), .groups = 'drop')

df_th <- as.data.frame(ccm_th) %>%
  group_by(LibSize) %>%
  summarize(rho = mean(rho), .groups = 'drop')

df_ccm <- left_join(df_ht, df_th, by = "LibSize", suffix = c("_ht", "_th"))

ggplot(df_ccm, aes(x = LibSize)) +
  geom_line(aes(y = rho_ht, color = "Humidity ➝ Temperature")) +
  geom_line(aes(y = rho_th, color = "Temperature ➝ Humidity")) +
  labs(title = "Convergent Cross Mapping (CCM)",
       y = expression(rho),
       x = "Library Size") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()


mv_df <- data.frame(
  time = 1:nrow(data),
  temperature_norm = data$temperature_norm,
  humidity_norm = data$humidity_norm,
  hive_power_norm = data$hive_power_norm
)


# Build numeric dataframe for Multiview
numeric_mv_data <- data.frame(
  time = 1:nrow(data),
  humidity_norm = data$humidity_norm,
  temperature_norm = data$temperature_norm
)



mv_df <- df_ht
mv_out <- Multiview(dataFrame = mv_df, E = 3, lib = c(1, 500), pred = c(501, 1000),
                    target = "humidity_norm", columns = c("humidity_norm", "temperature_norm"))

mv_preds <- mv_out$Predictions %>% drop_na()
mv_stats <- list(
  rho = cor(mv_preds$Observations, mv_preds$Predictions),
  MAE = mean(abs(mv_preds$Observations - mv_preds$Predictions)),
  RMSE = sqrt(mean((mv_preds$Observations - mv_preds$Predictions)^2))
)

# ============================
# Output Summary
# ============================
print("Optimal E values:")
print(c(Temperature = optimal_E_temp, Humidity = optimal_E_hum, Hive = optimal_E_hive))

print("SMap Humidity ~ Temperature:")
print(smap_ht)

print("SMap Humidity ~ Hive Power:")
print(smap_hh)

print("Multiview Embedding:")
print(mv_stats)