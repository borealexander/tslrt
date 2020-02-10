
#' Calculates the hazard ratio for switching treatment at a delay time
#'
#' \code{HR_delay_switch} Hazard ratio function for switching treatment at a delay time
#' @param lambda_c_1 Hazard rate in control group before switching
#' @param lambda_c_2 Hazard rate in control group after switching
#' @param lambda_e Hazard rate in experimental group
#' @param delay Time where treatment changes
#' @param p Proportion of (alive) patients that switches treatment at time delay
#' @return Survival and hazard ratio function
#'
#' @export

HR_delay_switch <- function(lambda_c_1, lambda_c_2, lambda_e, delay, p){

  S <- function(t){

    S <- ifelse( t < delay,
                 exp(-lambda_c_1*t),
                 exp(-lambda_c_1*delay) * (p*exp(-lambda_c_2*(t-delay)) + (1-p)*exp(-lambda_c_1*(t-delay))) )

  }

  HR <- function(t){

    # hazard for control after treatment switching
    h_c <- (lambda_c_2*p*exp(-lambda_c_2*(t-delay)) + lambda_c_1*(1-p)*exp(-lambda_c_1*(t-delay)))/(p*exp(-lambda_c_2*(t-delay)) + (1-p)*exp(-lambda_c_1*(t-delay)))

    HR <- ifelse(t < delay, lambda_e / lambda_c_1, lambda_e/h_c)

  }

  return(list(S = S, HR = HR))

}



#################################################################################################


#' Creates numerically estimated survival and hazard ratio function for treatment switching with exponential model
#'
#' \code{numerical_prog_switch_exp} numerically estimates the survival and HR functions in a treatment
#' switching model with exponentially distributions
#' @param lambda_c rate for control group
#' @param lambda_e rate for experimental group
#' @param lambda_p rate for progression in control group
#' @param p proportion of patients that switches after progression
#' @param max_t maximum time in time interval
#' @param delta_t size of time intervals in numerical approximations
#' @return a list containing functions for the survival curve and HR and hazard for control and proportion of switchers
#'
#' @export

numerical_prog_switch_exp <- function(lambda_c, lambda_e, lambda_p, p, max_t = 100, delta_t = 0.01){

  # time points
  t <- seq(0, 100, by = delta_t)

  # Survival functions
  S_prog_switch <- rep(0, length(t))
  S_prog_no_switch <- rep(0, length(t))
  S_nonprog <- rep(0, length(t))
  S_nonprog[1] <-  1
  p_switch <- rep(0, length(t))


  for (i in 1:(length(t)-1)) {

    S_nonprog[i+1] <- S_nonprog[i] * exp(-lambda_p*(t[i+1] - t[i])) * exp(-lambda_c*(t[i+1] - t[i]))
    S_prog_switch[i+1] <- p*S_nonprog[i] * (1 - exp(-lambda_p*(t[i+1] - t[i]))) + S_prog_switch[i] * exp(-lambda_e*(t[i+1] - t[i]))
    S_prog_no_switch[i+1] <- (1-p)*S_nonprog[i] * (1 - exp(-lambda_p*(t[i+1] - t[i]))) + S_prog_no_switch[i] * exp(-lambda_c*(t[i+1] - t[i]))
    p_switch[i+1] <- p_switch[i] + p*S_nonprog[i] * (1 - exp(-lambda_p*(t[i+1] - t[i])))

  }

  # overall survival function
  S <- S_prog_no_switch + S_prog_switch + S_nonprog
  S_fun <- approxfun(x = t, y = S)

  # hazard function for control
  h_c <- rep(0, length(t))


  for (i in 1:(length(t)-1)) {

    h_c[i+1] <- -(S[i+1] - S[i]) / (S[i+1] *(t[i+1] - t[i]))

  }

  h_c[1] <- h_c[2]
  h_c_fun <- approxfun(x = t, y = h_c)

  h_e <- rep(lambda_e, length(t))

  HR_fun <- approxfun(x = t, y = h_e/h_c)

  p_fun <- approxfun(x = t, y = p_switch)

  return(list(S = S_fun, h = h_c_fun, HR = HR_fun, p_switch = p_fun))


}
