#' Monte Carlo simulation for treatment switching with delay time switching
#'
#' Runs a Monte Carlo simulation for the treatment switching model with switch after delay
#' @param n_c Number of patients in control
#' @param n_e Number of patients in experimental
#' @param median_c_1 Median in control before switching (and for patients that do not switch)
#' @param median_c_2 Median in control after switching (for the patients that switches)
#' @param median_e_1 Median in experimental
#' @param delay Delay time for when switching occurs
#' @param end_event Number of events in the study
#' @param rec_period Time for recruitment
#' @param rec_power Variable for recruitment that follows the power model
#' @param p Proportion(s) of patients that switches after the delay
#' @param alpha One-sided significance level
#' @param M Number of simulations (for each p)
#'
#' @export

MC_delay_crossover <- function(n_c = 100,
                               n_e = 100,
                               median_c_1 = 10,
                               median_c_2 = 10,
                               median_e = 15,
                               delay = 5,
                               end_event = 150,
                               rec_period = 12,
                               rec_power = 1,
                               p = c(0, 1),
                               alpha = 0.025,
                               M = 100)
{

  result <- data.table(p = p,
                       p_switched = rep(0, length(p)),
                       power_logrank = rep(0, length(p)),
                       power_weighted = rep(0, length(p)),
                       efficiency = rep(0, length(p)))

  for (i in 1:length(p)) {


    print(paste0("Simulation progress: ", i, " out of ", length(p)))

    # model parameters
    trt_switch_model <- list(n_c = n_c,
                             n_e = n_e,
                             median_c_1 = median_c_1,
                             median_c_2 = median_c_2,
                             median_e = median_e,
                             delay = delay,
                             p_change = p[i],
                             end_event = end_event,
                             rec_period = rec_period,
                             rec_power = rec_power)


    # HR function
    # change name ???
    # change to different parameters than in simulation ???
    HR_switch <- HR_delay_switch(lambda_c_1 = log(2)/median_c_1,
                                 lambda_c_2 = log(2)/median_c_2,
                                 lambda_e = log(2)/median_e,
                                 delay = delay,
                                 p = p[i])$HR


    # simulate data
    sim_data <- replicate(M, sim_crossover_delay(model = trt_switch_model), simplify = FALSE)

    # create risk table
    # change to own risk table ???
    risk_table <- lapply(sim_data, get_risk_table)

    # calculate standard logrank
    logrank_rt <- lapply(risk_table, calculate_weights, method = "logrank")
    logrank_Z <- lapply(logrank_rt, calculate_zs)

    # calculate with proposed weight from article
    weighted_logrank_rt <- lapply(risk_table, calculate_weights, method = "theta", hr_fun = HR_switch)
    weighted_logrank_Z <- lapply(weighted_logrank_rt, calculate_zs)

    # proportion of patients that have switched
    #switched <- lapply(sim_data, proportion_switched)
    switched <- data.table(t(as.data.table(lapply(sim_data, prop_switch))))
    result$p_switched[i] <- as.numeric(lapply(switched, mean))

    # calculate power
    result$power_logrank[i] <- mean(logrank_Z > qnorm(1-alpha))
    result$power_weighted[i] <- mean(weighted_logrank_Z > qnorm(1-alpha))

    result$efficiency[i] <- ((qnorm(1-alpha) + qnorm(result$power_weighted[i])) /
      (qnorm(1-alpha) + qnorm(result$power_logrank[i])))^2


  }

  # data table for plotting
  n <- length(result$p)
  plot_dt <- data.table(p = rep(result$p, 2),
                        p_switched = rep(result$p_switched, 2),
                        power = c(result$power_logrank, result$power_weighted),
                        test = rep(c("Logrank", "Weighted"), c(n, n)))


  p <- ggplot(aes(x = 100*p_switched, y = power, color = test), data = plot_dt) +
    geom_line(size = 1.5) +
    scale_x_continuous("% of switchers") +
    labs(y = "Power", color = "") +
    theme_bw() +
    theme(text = element_text(size = 12),
          legend.position = c(0.85, 0.85),
          legend.title = element_blank(),
          panel.grid = element_blank()) +
    scale_colour_manual(values = c("black", "grey80")) +
    ylim(0,1)


  p.eff <- ggplot(aes(x = 100*p_switched, y = efficiency), data = result) +
    geom_line(size = 1.5) +
    scale_x_continuous("% of switchers") +
    labs(y = "Efficiency", color = "") +
    theme_bw() +
    theme(text = element_text(size = 12),
          legend.position = c(0.85, 0.85),
          legend.title = element_blank(),
          panel.grid = element_blank()) +
    scale_colour_manual(values = c("black"))

  return(list(result = result, plot = p, plot.eff = p.eff))

}
