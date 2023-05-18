#' Update params
#' Updata parameters such that they balance. Changing popoulation size will
#' lead to the chanages in the initials values of susceptibles, infecteds, etc.
#' Similarly, changing the proportion of immune in the beginning of the population
#' will lead to chantges in initial values of the state variables
#' Change in the date of introduction change the simulation period
#' Updating occurs by way of the name of the parameter
#'
#' @param pars # parameter values to be varied
#' @param pars_baseline # This represents a complete set of parameters needed to run the model including initial values for the state variables. Some of the parameter values may be updated based on the pars
#' @param model # model (e.g., seirw)
#' @param output_time #
#' @param output_state
#'
#' @return
#' @export
#'
#' @examples
#'
update_params <- function(pars=NULL,
                          pars_baseline=NULL) {

  if (is.null(pars) | is.null(pars_baseline)) {
    stop("Both pars and pars_baseline must not be NULL")
  }

  params = pars_baseline

  nms = names(pars)
  for (nm in nms) {
    params[[nm]] = pars[[nm]]
  }

  # simulation times changes by Day 1 (introduction of the virus)
  # index starts from 1 and therefore one must be added to get the correct
  params[["ndays"]] <- round(params[["day1"]]) +
    (params[["obs_length"]] * params[["output_days"]]) + 1

  pop = params[["population"]]
  prop_eff = params[["prop_eff_pop"]]
  prop_r = params[["prop_immune"]]
  i0 = 1
    # initial values for the state variables
  params$Susceptible <- pop*prop_eff*(1-prop_r) - i0
  params$Exposed <- 0
  params$Infectious<- i0
  params$Asymptomatic <- 0
  params$Recovered <- pop*prop_eff*prop_r
  params$Water <- 0
  params$CumulExposed <- 0
  params$CumulInfectious <- 0

  return (params)
}
