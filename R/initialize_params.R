initialize_params <- function(I0=10, pop=1e5,
                        agetwo=FALSE,
                        erlang=FALSE){
  # 25970000 North Korea population size
  # since SEPIAR model is a stochastic model stochastic die-out possible
  y0 <- c(S=pop-I0, E=0, A=0, I=I0, R=0, W=0, CE=0, CI=0)

  params <- list() # input parameters for the SEIR model
  params$measure_var <- "CI" # daily diff of CI gives daily symptomatic case
  params$model <- seiarw
  params$erlang <- erlang
 # (Leung PLoS NTD 2022)
  params$prop_children <- 0.193 # mean latent period = 1/epsilon
  params$epsilon <- 1/2 # mean latent period = 1/epsilon
  params$gamma <- 1/4 # mean infectious period = 1/gamma
  params$kappa <- 575 # excretion rate cells per person per day

  params$xi <- 1/21 # mean decay rate of Vibrio cholerae
  params$K <- 10000 # half-infective bacteria dose (10,000 cells/ml)
  params$R0 <- 3.0
  params$R0W <- 0.8 #
  params$sigma <- 1/(4*365) # 1/sigma = mean duration of natural immunity
  # 1/sigma_v1 = mean duration of OCV-induced immunity (1st dose)
  params$sigma_v1 <- 1/(2*365)
  # 1/sigma_v2 = mean duration of OCV-induced immunity (2nd dose)
  params$sigma_v2 <- 1/(4*365)
  params$vacc_cov <- c(0.6, 0)
  params$campaign_dur <- 14 # 14 days of vaccination campaign
  params$delay_until_2nd_campaign <- 30 # 30 days of delay between the 1st and the 2nd campaign
  params$vacc_eff_v1 <- c(0.3, 0.5)
  params$vacc_eff_v2 <- c(0.5, 0.7)

  params$day_intervention <- 100.0

  params$fA <- 0.9 # fraction of asymptomatic state
  params$bA <- 1 # relative infectiousness of asymptomatic state
  params$tau <- 0.1 # time step size
  params$ndays <- 365 # number of days for output
  # North Korea data have >80 observations, ignoring some zero incidence
  params$obslength <- 80

  params$alpha <- 0.0; # proportional reduction in R0
  # time window over which the number of cases is tracked
  params$case_track_window <- 7
  # case threshold over which intervention will be implemented
  params$case_threshold <- 10

  params$init$S <- y0[["S"]]
  params$init$E <- y0[["E"]]
  params$init$I <- y0[["I"]]
  params$init$A <- y0[["A"]]
  params$init$R <- y0[["R"]]
  params$init$W <- y0[["W"]]
  params$init$CE <- y0[["CE"]]
  params$init$CI <- y0[["CI"]]

  if (agetwo) {
    params$model <- seiar_erlang
    params$init <- NULL
    params$init$S1 <- params$prop_children*(pop - I0)
    params$init$S2 <- (1 - params$prop_children)*(pop - I0)
    params$init$E1 <- 0
    params$init$E2 <- 0
    params$init$I1 <- params$prop_children*I0
    params$init$I2 <- (1 - params$prop_children)*I0
    params$init$A1 <- 0
    params$init$A2 <- 0
    params$init$R1 <- 0
    params$init$R2 <- 0
    params$init$CE1 <- 0
    params$init$CI1 <- 0
    params$init$CE2 <- 0
    params$init$CI2 <- 0
    params$init$W <- 0
  }

  return(params)
}
