#' Calculates the propoertion of patients that have switched treatment in 1 simulation
#'
#' @param dt A simulation with treatment switching
#' @return Proportion of patients that have switched treatment
#' @export

prop_switch <- function(dt)
{

  mean(dt[dt$group == "control",]$switch)

  # why does this not work???
  # mean(dt[group == "control"]$switch)


}
