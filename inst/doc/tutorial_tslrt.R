## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install_package, message = FALSE, warning = FALSE-------------------

#install.packages("devtools")
library(devtools)
#install_github("borealexander/tslrt")
library(tslrt)


## ----packages, message = FALSE, warning = FALSE--------------------------

library(data.table)
library(survival)
library(ggplot2)
library(dplyr)


## ----example_1-----------------------------------------------------------

# Set parameters

n_c <- 100
n_e <- 100
median_c_1 <- 10
median_c_2 <- 20
median_e <- 20
p_change <- 0.9
delay <- 3
end_event <- 150
rec_period <- 12
rec_power <- 1

# Model as input
model_delay_switch <- list(n_c = n_c,
                           n_e = n_e,
                           median_c_1 = median_c_1,
                           median_c_2 = median_c_2,
                           delay = delay,
                           median_e = median_e,
                           p_change = p_change,
                           end_event = end_event,
                           rec_period = rec_period,
                           rec_power = rec_power)

# Simulate data
ex_sim_delay <- sim_crossover_delay(model = model_delay_switch)

# print out the top of the simulated data
head(ex_sim_delay)


## ----plot_example_1, fig.width = 7, fig.height = 5-----------------------

# Plot data
fit_delay <- survfit( Surv(time, event) ~ group, data = ex_sim_delay)
plot(fit_delay, col = c("black", "blue"), mark.time = TRUE, pch = c(1, 2))
legend('topright', legend = c("Control", "Experimental"),
       col = c("black", "blue"), lty = c(1, 1), pch = c(1, 2))


## ----plot_HR_example_1, fig.width = 7, fig.height = 5--------------------

t <- seq(1,50)

# Hazard ratio
HR_switch <- HR_delay_switch(lambda_c_1 = log(2)/median_c_1,
                             lambda_c_2 = log(2)/median_c_2,
                             lambda_e = log(2)/median_e,
                             delay = delay,
                             p = p_change)

plot(t, HR_switch$HR(t), type = "l",
     xlab = "Time (months)", ylab = "HR")


## ----plot_survival_example_1, fig.width = 7, fig.height = 5--------------

plot(fit_delay, col = c("black", "blue"), mark.time = TRUE, pch = c(1, 2))
lines(t, exp(-t*log(2)/median_c_1), col = "red", lty = 2)
lines(t, HR_switch$S(t), col = "black", lty = 2)
lines(t, exp(-t*log(2)/median_e), col = "blue", lty = 2)
legend('topright', legend = c("Control", "Experimental",
                              "Control model - no switching",
                              "Control model - with switching",
                              "Experimental model"),
       col = c("black", "blue", "red", "black", "blue"),
       lty = c(1, 1, 2, 2, 2), 
       pch = c(1, 2, NA, NA, NA))


## ----risk_table_example_1------------------------------------------------

# Create risk table
delay_risk_table <- get_risk_table(ex_sim_delay)

head(delay_risk_table)


## ----weights_example_1, fig.width = 7, fig.height = 5--------------------

# create risk table
delay_risk_table <- get_risk_table(ex_sim_delay)

# calculate weights
logrank_weights <- calculate_weights(delay_risk_table, method = "logrank")
delay_weights <- calculate_weights(delay_risk_table, method = "theta", hr_fun = HR_switch$HR)
fh_weights <- calculate_weights(delay_risk_table, method = "fh", rho = 1, gamma = 0)


# data table for plotting the weights
weight_dt <- data.table(t  = c(logrank_weights$t, delay_weights$t, fh_weights$t),
                        w = c(logrank_weights$w, delay_weights$w, fh_weights$w),
                        Test = rep(c('logrank', 'new weights', 'fh(1,0)'),
                                   c(length(logrank_weights$w),length(delay_weights$w),length(fh_weights$w))))


ggplot(weight_dt, aes(x = t, y = w, color = Test)) +
  geom_line(size = 1.2) +
  theme_bw() +
  labs(x = "Time (months)", y = "w", title = "Comparison of weights for different tests") +
  theme(text = element_text(size = 18),
        panel.grid.minor = element_blank())


## ----calculate_Z_example_1-----------------------------------------------
# Calculate Z-scores

# Logrank
calculate_zs(logrank_weights)

# New weights
calculate_zs(delay_weights)

# Flemming-Harrington
calculate_zs(fh_weights)


