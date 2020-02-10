#' Simulation with crossover effect
#'
#' \code{sim_exp_crossover} Simulating a trial with crossover effect.
#' Control and experimental arm are both exponentially distributed,
#' but control can have progression and might change to experimental arm
#' @param n_c Number of patients in control arm
#' @param n_e Number of patients in experimental arm
#' @param median_c Median in control arm
#' @param median_e Median in experimental arm
#' @param median_progression Median for getting a progression in control arm
#' @param p_change Probability that a patient in control that has a progression changes to experimental arm
#' @param end_event Number of events to stop the study
#' @param rec_period Recruitment period
#' @param rec_power Recruitment follows a power model, this variable how the power looks
#' @param model List containing all the variabels above
#' @return dt Data table containing time, event and group/arm
#'
#'@export

sim_exp_crossover <- function(n_c = 100,
                              n_e = 100,
                              median_c = 10,
                              median_e = 15,
                              median_progression = 5,
                              p_change = 0.5,
                              end_event = 150,
                              rec_period = 12,
                              rec_power = 1,
                              model = NULL)
{

  if( !is.null(model))
  {
    n_c <- model$n_c
    n_e <- model$n_e
    median_c <- model$median_c
    median_e <- model$median_e
    median_progression <- model$median_progression
    p_change <- model$p_change
    end_event <- model$end_event
    rec_period <- model$rec_period
    rec_power <- model$rec_power
  }

  # recruitment of patients (for control and experimental group)
  rec_t_c <- rec_period * runif(n_c) ^ (1 / rec_power)
  rec_t_e <- rec_period * runif(n_e) ^ (1 / rec_power)

  # events times for control with crossover after progression
  t_c_1 <- rexp(n = n_c, rate = log(2)/median_c)
  t_c_2 <- rexp(n = n_c, rate = log(2)/median_e)
  t_c_prog <- rexp(n = n_c, rate = log(2)/median_progression)
  prog <- rbinom(n_c, size = 1, p_change)

  switch_c <- ((t_c_1  > t_c_prog) & (prog == 1))

  t_c <- ifelse(switch_c, (t_c_2 + t_c_prog), t_c_1)

  # event times for experimental
  t_e <- rexp(n = n_e, rate = log(2)/median_e)

  # "calender" times for events
  cal_t_c = rec_t_c + t_c
  cal_t_e = rec_t_e + t_e

  # calender time for when our final event happens
  max_cal_t <- sort(c(cal_t_c, cal_t_e))[end_event]


  # events (if the event happens before the final event)
  event_c <- cal_t_c <= max_cal_t
  event_e <- cal_t_e <= max_cal_t

  # "observation" times for the study if event of censored
  obs_t_c <- ifelse(event_c, t_c, max_cal_t - rec_t_c)
  obs_t_e <- ifelse(event_e, t_e, max_cal_t - rec_t_e)

  dt <- data.table(time = c(obs_t_c, obs_t_e),
                   event = c(event_c, event_e),
                   group = factor(rep(1:2, c(n_c, n_e)),
                                  labels = c("control", "experimental")),
                   switch = c(switch_c, rep(FALSE, n_e)))

  dt

}

