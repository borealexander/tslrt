#' Get (standardized) score and weighted log-rank statistics.
#'
#' \code{get_zs} Calculate (standardized) score and weighted log-rank statistics from risk table (with weights).
#' @param risk_table A risk table with weights produced by the functions \code{get_risk_table} and \code{calculate_weights}.
#' @return The standardized weighted log-rank statistic.
#'
#'@export


# Dominics code from modestly

calculate_zs <- function(risk_table){

  #if (class(risk_table) == "list") risk_table = risk_table$risk_table

  n_e <- max(risk_table$n_e)
  n_c <- max(risk_table$n_c)


  # formulas for U and V[U] in Leton and Zuluaga pg. 596.

  u <- with(risk_table, sum(w * (d_c - d * n_c / n)))
  v_u <- with(risk_table, sum(w^2 * n_c * n_e * d * (n - d) / n / n / (n - 1), na.rm = TRUE))

  z_u <- u / sqrt(v_u)

  z_u


}
