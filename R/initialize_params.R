initialize_params <- function(agetwo=FALSE, erlang=FALSE, ...){

  params <- list() # input parameters for the SEIR model

  # parameters to be estimates based on the outbreak data
  params$day1 <- 14 #
  params$prop_immune <- 0.0 #
  params$R0 <- 3.0 #
  params$prop_report <- 0.2 #
  # baseline population size is critical for the outbreak size
  # and the currently is the the population size of the admin in which the
  # outbreak was reported.  Actual population is likely to be a fraction
  params$prop_eff_pop <- 0.1 #

  # time window over which the number of cases is tracked (i.e., weekly reported)
  params$output_days <- 7
  # obs_length refers the duration during which data are available in days
  # this is updated for each data set
  params$obs_length <- 20 # 20 weeks
  # this refers to the total simulation days and has to be larger than or equal
  # to the obs_length
  params$ndays <- 365 # total number of days for output

  params$model <- seiarw
  params$erlang <- erlang
  params$tau <- 0.01 # time step size for numerical integration


  params$epsilon <- 1/1.4 # mean latent period = 1/epsilon
  params$gamma <- 1/2 # mean infectious period = 1/gamma
  params$fA <- 0.5 # fraction of asymptomatic state
  params$bA <- 0.05 # relative infectiousness of asymptomatic state
  params$kappa <- 575 # excretion rate cells per person per day

  params$prop_children <- 0.193 #

  params$xi <- 1/21 # mean decay rate of Vibrio cholerae
  params$K <- 10000 # half-infective bacteria dose (10,000 cells/ml)

  params$R0W <- 0.8 #
  params$sigma <- 1/(4*365) # 1/sigma = mean duration of natural immunity
  # 1/sigma_v1 = mean duration of OCV-induced immunity (1st dose)
  params$sigma_v1 <- 1/(2*365)
  # 1/sigma_v2 = mean duration of OCV-induced immunity (2nd dose)
  params$sigma_v2 <- 1/(4*365)
  params$vacc_cov <- c(0.6, 0)
  params$campaign_dur <- 5 # 14 days of vaccination campaign
  params$delay_until_2nd_campaign <- 14 # 30 days of delay between the 1st and the 2nd campaign
  params$vacc_eff_v1 <- c(0.3, 0.5)
  params$vacc_eff_v2 <- c(0.5, 0.7)

  params$alpha <- 0.0; # proportional reduction in R0
  params$case_track_window <- 7;
  # case threshold over which intervention will be implemented
  #
  params$case_threshold <- 10

  i0=10
  pop=1e5
  prop_r = 0.4
  # initial values for the state variables
  params$population <- pop
  params$susceptible <- params$population*(1-prop_r) - i0
  params$exposed <- 0
  params$infectious<- i0
  params$asymptomatic <- 0
  params$recovered <- params$population*prop_r
  params$water <- 0
  params$cumul_exposed <- 0
  params$cumul_infectious <- 0

  if (agetwo) {
    params$model <- seiarw_2ag_erlang
    # params$init <- NULL
    # params$init$S1 <- params$prop_children*(pop - I0)
    # params$init$S2 <- (1 - params$prop_children)*(pop - I0)
    # params$init$E1 <- 0
    # params$init$E2 <- 0
    # params$init$I1 <- params$prop_children*I0
    # params$init$I2 <- (1 - params$prop_children)*I0
    # params$init$A1 <- 0
    # params$init$A2 <- 0
    # params$init$R1 <- 0
    # params$init$R2 <- 0
    # params$init$CE1 <- 0
    # params$init$CI1 <- 0
    # params$init$CE2 <- 0
    # params$init$CI2 <- 0
    # params$init$W <- 0
  }

  params = update_params(list(...), params)

  return(params)
}
