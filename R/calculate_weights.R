#' Function for calculating weights for logrank test
#'
#' Calculates weights for treatment switching logrank test
#' @param risk_table Risk table from the function \code{}
#' @param method Which method for calculating the weights
#' \code{logrank} Standard logrank test (w=1 for all)
#' \code{theta} w = -log(hr(t))
#' \code{hr_weight} Weight from article
#' \code{fh} Flemming-Harrington weights
#' @param hr_fun Hazard ratio function for method \code{theta} and \code{hr_weight}
#' @param rho Variable for method \code{fh}
#' @param gamma Varaiable for method \code{fh}
#' @return Risk table updated with weights
#'
#' @export

calculate_weights <- function(risk_table, method = "logrank", hr_fun = NULL, rho = 0, gamma = 0)
{

  if(method == "logrank")
  {
    risk_table$w <- 1
  }

  else if(method == "theta")
  {
    risk_table$w <- -log(hr_fun(risk_table$t))
  }

  else if(method == "hr_weight")
  {
    #risk_table$w <- (hr_fun(risk_table$t)^(-1) -hr_fun(max(risk_table$t))^(-1)) / (hr_fun(0)^(-1) -hr_fun(max(risk_table$t))^(-1))
    risk_table$w <- (hr_fun(risk_table$t)^(-1) - 1) / (hr_fun(0)^(-1) - 1)

  }

  else if(method == "fh"){

    # Kaplan-Meier estimate of survival in the lumped data:
    risk_table$s <- exp(cumsum(log(1 - risk_table$d / risk_table$n)))
    # Fleming-Harrington (rho, gamma) weights:
    risk_table$w <- risk_table$s ^ rho * (1 - risk_table$s) ^ gamma

  }



  risk_table
}
