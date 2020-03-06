#' Monte Carlo simulation for treatment switching with switching after progression
#'
#' Runs a Monte Carlo simulation with treatment switching after progression
#' @param n_c Number of patients in control
#' @param n_e Number of patients in experimental
#' @param median_c Median in control arm
#' @param median_e Median in experimental arm
#' @param median_prog Median for progressions
#' @param end_event Number of events to stop the study
#' @param rec_period Time for recruitment
#' @param rec_power Variable for recruitment that follows the power model
#' @param p Proportion(s) of patients that switches after the delay
#' @param HR_fun Pre-specified HR function for calculating weights
#' @param alpha One-sided significance level
#' @param M Number of simulations (for each p)
#'
#' @export

MC_exp_prog_crossover_median <- function(n_c = 100,
                                         n_e = 100,
                                         median_c = c(5, 10),
                                         median_e = 15,
                                         median_prog = 5,
                                         end_event = 150,
                                         rec_period = 12,
                                         rec_power = 1,
                                         p = 1,
                                         HR_fun,
                                         alpha = 0.025,
                                         M = 100)
{

  result <- data.table(median_c = median_c,
                       p_switched = rep(0, length(median_c)),
                       Z_logrank = rep(0, length(median_c)),
                       Z_weighted = rep(0, length(median_c)),
                       Z_FH = rep(0, length(median_c)),
                       power_logrank = rep(0, length(median_c)),
                       power_weighted = rep(0, length(median_c)),
                       power_fh = rep(0, length(median_c)),
                       rel_eff_MWLRLR = rep(0, length(median_c)),
                       rel_eff_MWLRFH = rep(0, length(median_c)),
                       eff_MWLRLR = rep(0, length(median_c)),
                       eff_MWLRFH = rep(0, length(median_c)))

  for (i in 1:length(median_c)) {


    print(paste0("Simulation progress: ", i, " out of ", length(median_c)))

    # model parameters
    trt_switch_model <- list(n_c = n_c,
                             n_e = n_e,
                             median_c = median_c[i],
                             median_e = median_e,
                             median_progression = median_prog,
                             p_change = p,
                             end_event = end_event,
                             rec_period = rec_period,
                             rec_power = rec_power)
    # HR function
    HR_switch <- HR_fun

    # simulate data
    sim_data <- replicate(M, sim_exp_crossover(model = trt_switch_model), simplify = FALSE)

    # create risk table
    # change to own risk table ???
    risk_table <- lapply(sim_data, get_risk_table)

    # calculate standard logrank
    logrank_rt <- lapply(risk_table, calculate_weights, method = "logrank")
    logrank_Z <- lapply(logrank_rt, calculate_zs)
    result$Z_logrank[i] <- mean(unlist(logrank_Z))

    # calculate with proposed weight from article
    weighted_logrank_rt <- lapply(risk_table, calculate_weights, method = "theta", hr_fun = HR_switch)
    weighted_logrank_Z <- lapply(weighted_logrank_rt, calculate_zs)
    result$Z_weighted[i] <- mean(unlist(weighted_logrank_Z))

    #calculate with FH class of weights
    fh_logrank_rt <- lapply(risk_table, calculate_weights, method = "fh", rho = 1, gamma = 0)
    fh_logrank_Z <- lapply(fh_logrank_rt, calculate_zs)
    result$Z_FH[i] <- mean(unlist(fh_logrank_Z))

    # proportion of patients that have switched
    #switched <- lapply(sim_data, proportion_switched)
    switched <- data.table(t(as.data.table(lapply(sim_data, prop_switch))))
    result$p_switched[i] <- as.numeric(lapply(switched, mean))

    # calculate power
    result$power_logrank[i] <- mean(logrank_Z > qnorm(1-alpha))
    result$power_weighted[i] <- mean(weighted_logrank_Z > qnorm(1-alpha))
    result$power_fh[i] <- mean(fh_logrank_Z > qnorm(1-alpha))


  }

  #efficiency
  result$eff_MWLRLR <- (result$Z_weighted/result$Z_logrank)^2

  result$eff_MWLRFH <- (result$Z_weighted/result$Z_FH)^2

  # relative efficiency
  result$rel_eff_MWLRLR <- ((qnorm(1-alpha) + qnorm(result$power_weighted)) /
                                 (qnorm(1-alpha) + qnorm(result$power_logrank)))^2

  result$rel_eff_MWLRFH <- ((qnorm(1-alpha) + qnorm(result$power_weighted)) /
                                 (qnorm(1-alpha) + qnorm(result$power_fh)))^2

  # data table for plotting
  n <- length(result$median_c)
  plot_dt_power <- data.table(median_c = rep(result$median_c, 3),
                              power = c(result$power_logrank, result$power_weighted, result$power_fh),
                              test = rep(c("LR", "mWLR", "WLR-FH"), c(n, n, n)))


  p.power <- ggplot(aes(x = median_c, y = 100*power, color = test), data = plot_dt_power) +
    geom_line(size = 1) +
    scale_x_continuous("Median OS in control group (months)") +
    labs(y = "Power (%)", color = "") +
    theme_classic() +
    ylim(0, 100) +
    theme(text = element_text(size = 14),
          legend.position = c(0.75, 0.85),
          legend.title = element_blank(),
          panel.grid = element_blank()) +
    annotate("text", x = 9, y = 71, label = paste0("P(switching) = ", p))

  plot_dt_eff <- data.table(median_c = rep(result$median_c, 2),
                            efficiency = c(result$eff_MWLRLR, result$eff_MWLRFH),
                            test = rep(c("mWLR / LR", "mWLR / WLR-FH"), c(n, n)))

  p.eff <- ggplot(aes(x = median_c, y = 100*efficiency, color = test), data = plot_dt_eff) +
    geom_hline(yintercept=100, linetype="dashed", color = "black") +
    geom_line(size = 1) +
    scale_x_continuous("Median OS in control group (months)") +
    labs(y = "Efficiency (%)", color = "") +
    theme_classic() +
    theme(text = element_text(size = 14),
          legend.position = c(0.75, 0.85),
          legend.title = element_blank(),
          panel.grid = element_blank())

  return(list(result = result, plot.power = p.power, plot.eff = p.eff))

}
