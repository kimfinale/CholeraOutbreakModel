#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List seiarw(List params) {
  double tau = params["tau"]; // time step size
  double ndays = params["ndays"]; // number of days for output
  int nsteps = ceil(ndays / tau) + 1;
  // vectors for state variables
  NumericVector S(nsteps);
  NumericVector E(nsteps);
  NumericVector I(nsteps);
  NumericVector A(nsteps); //asymptomatic
  NumericVector R(nsteps); // recovered
  NumericVector W(nsteps); // recovered
  NumericVector CE(nsteps); // cumulative infection
  NumericVector CI(nsteps); // cumulative symptomatic
  // vaccination assumes that vaccination does not change the status of E, I, R
  NumericVector V1(nsteps); // susceptibles who received the first dose of OCV
  NumericVector V2(nsteps); // susceptibles who received the second dose of OCV
  NumericVector CV1(nsteps); // cumulative number of the first dose of OCV
  NumericVector CV2(nsteps); // cumulative number of the second dose of OCV

  NumericVector time(nsteps);// fill time w/zeros

  // initial values
  S(0) = params["susceptible"];
  E(0) = params["exposed"];
  I(0) = params["infectious"];
  A(0) = params["asymptomatic"];
  R(0) = params["recovered"];
  W(0) = params["water"];
  CE(0) = params["cumul_exposed"];
  CI(0) = params["cumul_infectious"];
  V1(0) = 0;
  V2(0) = 0;
  CV1(0) = 0;
  CV2(0) = 0;

  double prop_report = params["prop_report"]; // fraction
  double epsilon = params["epsilon"]; // 1 / latent period
  double kappa = params["kappa"]; // excretion rate
  double xi = params["xi"]; // decay rate
  double K = params["K"]; // half-infective dose (eg, 10,000 doses/ml)
  double gamma = params["gamma"]; // 1 / recovery period
  double sigma = params["sigma"]; // 1 / duration of natural immunity
  double sigma_v1 = params["sigma_v1"]; // 1 / duration of OCV-induced immunity (1st dose)
  double sigma_v2 = params["sigma_v2"]; // 1 / duration of OCV-induced immunity (2nd dose)
  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double fA = params["fA"]; // fraction asymptomatic
  double bA = params["bA"]; // relative infectivity of A to I
  double R0 = params["R0"];
  double R0W = params["R0W"]; // transmission rate arising from water
  double alpha = params["alpha"]; //proportional reduction in R0
  // time window over which the number of cases is tracked
  double case_track_window = params["case_track_window"];
  // case threshold over which intervention will be implemented
  double case_threshold = params["case_threshold"];

  double campaign_dur = params["campaign_dur"]; // vaccine duration
  // double delay_until_2nd = params["delay_until_2nd_campaign"];
  NumericVector vacc_cov_ = params["vacc_cov"]; // vaccine coverage
  NumericVector VE_v1_ = params["vacc_eff_v1"]; // vaccine efficacy (1st dose) for 5+ yo
  NumericVector VE_v2_ = params["vacc_eff_v2"]; // vaccine efficacy (2nd dose) for 5+ yo

  double vacc_cov = vacc_cov_[0]; // vaccine coverage
  double VE_v1 = VE_v1_[0]; // vaccine efficacy (1st dose) for 5+ yo
  double VE_v2 = VE_v2_[0]; // vaccine efficacy (2nd dose) for 5+ yo

  // double vacc_eff_u5 = params["vacc_eff_u5"]; // vaccine efficacy for 5+ yo
  // daily vaccination rate during a campaign_dur is given by the following relationship
  //  -log(1-vacc_cov)/campaign_dur, where vacc_cov is between 0 and 1
  // double day_intervention = params["day_intervention"];

  // double beta = R0 / (bA*fA/gamma + (1-fA)/gamma);
  double beta = R0 * gamma;
  // double betaW = beta;
  double pop_init = S[0] + E[0] + I[0] + A[0] + R[0] + V1[0] + V2[0];
  double betaW = R0W * gamma * (xi / kappa) * K / pop_init; // needs
  double vacc_rate = 0;
  // double vacc_rate = -log(1-vacc_cov)/campaign_dur;

  // Rprintf("the value of nsteps : %.1d \n", nsteps);
  // Rprintf("the value of R0 : %.3f \n", R0);
  // Rprintf("the value of R0W : %.3f \n", R0W);
  // Rprintf("the value of beta : %.3f \n", beta);
  // Rprintf("the value of betaW : %.3f \n", betaW);
  // Rprintf("the value of kappa : %.3f \n", kappa);
  // Rprintf("the value of xi : %.3f \n", xi);
  // Rprintf("the value of epsilon : %.3f \n", epsilon);
  // Rprintf("the value of gamma : %.3f \n", gamma);

  // Calculate the number of events for each step, update state vectors
  double case_tracked = 0;
  bool intervention_in_place = false;
  // arbitrary small number (at least abs value is arger than the simulation days)
  double campaign_start = -1e5;

  for (int step = 0; step < nsteps - 1; step++) {
    // --------------- Intervention ----------------------------
    if (case_tracked >= case_threshold) {
      beta = (1 - alpha) * R0 / (bA*fA/gamma + (1-fA)/gamma);
      intervention_in_place = true;
      if(campaign_start < 0){ // initial values assigned to be zero
        campaign_start = step*tau;
      }
    }
    if((step % (int) ceil(case_track_window/tau)) == 0 && !intervention_in_place) {
      case_tracked = 0;
    }

    // vaccination related
    double time_since_campaign = (step+1)*tau - campaign_start;
    if (0 < time_since_campaign && time_since_campaign <= campaign_dur) {
      vacc_rate = -log(1-vacc_cov)/campaign_dur;
    }
    else {
      vacc_rate = 0;
    }
    // -----------------------------------------------------------------
    // test vacc_rate
    // if (step % 5 == 0) {
    //   Rprintf("step: %d, time: %.1f, vacc_rate: %.3f \n",
    //           step, (step+1)*tau, vacc_rate);
    // }
    double St = S[step];
    double Et = E[step];
    double It = I[step];
    double At = A[step];
    double Rt = R[step];
    double Wt = W[step];

    double CEt = CE[step];
    double CIt = CI[step];

    double V1t = V1[step];
    double V2t = V2[step];
    double CV1t = CV1[step];
    double CV2t = CV2[step];

    // State Equations
    double N = St + Et + It + At + Rt + V1t + V2t;
    double foi = beta * (bA * At + It) / N + betaW * Wt / (K + Wt);
    // Rprintf("the value of N : %.1f \n", N);
    // double new_infection = (iS + (1-VE_v1)*iV1 + (1-VE_v2)*iV2) * foi * tau;
    double StoE = St * foi * tau;
    double RtoS = Rt * sigma * tau;
    double V2toS = V2t * sigma_v2 * tau; // immunity wane (2nd dose)
    double V1toS = V1t * sigma_v1 * tau;
    // new vacc doses NEED TO BE ADJUSTED
    double new_vacc_first_dose = (St+Et+At+Rt) * vacc_rate * tau;
    double new_vacc_second_dose = (St+Et+At+Rt+V1t) * vacc_rate * tau;
    // Rprintf("the value of new infection : %f \n", new_infection);

    // Rprintf("the value of new_infection : %.1f \n", new_infection);
    // Rprintf("the value of EtoP : %.1f \n", EtoP);
    double EtoI = Et * epsilon * (1-fA) * tau;//
    double EtoA = Et * epsilon * fA * tau;//

    // Rprintf("the value of PtoI : %.1f \n", PtoI);
    // Rprintf("the value of PtoA : %.1f \n", PtoA);
    double ItoR = It * gamma * tau;
    double AtoR = At * gamma * tau;//

    double ItoW = At * kappa * tau;//
    double AtoW = At * bA * kappa * tau;//
    double fromW = Wt * xi * tau;//

    double StoV1 = St * vacc_rate * tau;//
    double V1toV2 = V1t * vacc_rate * tau;//
    double V1toE = V1t * (1-VE_v1) * foi * tau;
    double V2toE = V2t * (1-VE_v2) * foi * tau;
    // Rprintf("the value of ItoR : %.1f \n", ItoR);
    // Rprintf("the value of AtoR : %.1f \n", AtoR);

    // Calculate the change in each state variable
    double dS = - StoE - StoV1 + RtoS + V1toS + V2toS;
    double dE = StoE + V1toE + V2toE - EtoA - EtoI;
    double dI = EtoI - ItoR;
    double dA = EtoA - AtoR;
    double dR = ItoR + AtoR - RtoS;
    double dW = ItoW + AtoW - fromW;
    double dV1 = StoV1 - V1toS - V1toV2 - V1toE;
    double dV2 = V1toV2 - V2toS - V2toE;

    // Update next timestep
    S[step + 1] = St + dS;
    E[step + 1] = Et + dE;
    I[step + 1] = It + dI;
    A[step + 1] = At + dA;
    R[step + 1] = Rt + dR;
    W[step + 1] = Wt + dW;

    V1[step + 1] = V1t + dV1;
    V2[step + 1] = V2t + dV2;

    CE[step + 1] = CEt + StoE + V1toE + V2toE;// cumulative infection
    CI[step + 1] = CIt + prop_report * EtoI;// OBSERVED cumulative symptomatic
    // this will later need to be adjusted to account for the delay of observation
    case_tracked = case_tracked + EtoI;

    CV1[step + 1] = CV1t + new_vacc_first_dose;
    CV2[step + 1] = CV2t + new_vacc_second_dose;
    time[step + 1] = (step + 1) * tau;// time in fractional years
  }
// Return results as data.frame
  DataFrame result = DataFrame::create(
    Named("time") = time,
    Named("S") = S,
    Named("E") = E,
    Named("I") = I,
    Named("A") = A,
    Named("R") = R,
    Named("W") = W,
    Named("V1") = V1,
    Named("V2") = V2,
    Named("CE") = CE,
    Named("CI") = CI,
    Named("CV1") = CV1,
    Named("CV2") = CV2);

  return result;
}

// [[Rcpp::export]]
List seiarw_2ag_erlang(List params) {
  double tau = params["tau"]; // time step size
  double ndays = params["ndays"]; // number of days for output
  int nsteps = ceil(ndays / tau) + 1;
  // vectors for state variables
  // waiting times are modeled as Erlang with shape=2
  // two age groups implemented as two columns of the matrix
  NumericMatrix S(nsteps,2); //
  NumericMatrix E1(nsteps,2); //
  NumericMatrix E2(nsteps,2); //
  NumericMatrix I1(nsteps,2); //
  NumericMatrix I2(nsteps,2); //
  NumericMatrix A1(nsteps,2); //
  NumericMatrix A2(nsteps,2); //
  NumericMatrix R1(nsteps,2); //
  NumericMatrix R2(nsteps,2); //

  NumericMatrix CE(nsteps,2); // cumulative infection
  NumericMatrix CI(nsteps,2); // cumulative symptomatic
  // vaccination assumes that vaccination does not change the status of E, I, R
  NumericMatrix O1(nsteps,2); // susceptibles who received the one dose of OCV
  NumericMatrix O2(nsteps,2); // susceptibles who received the dose of OCV
  NumericMatrix T1(nsteps,2); // susceptibles who received the second dose of OCV
  NumericMatrix T2(nsteps,2); // susceptibles who received the second dose of OCV

  NumericVector W(nsteps); // recovered
  NumericVector time(nsteps);// fill time w/zeros


  double s0 = params["Susceptible"];
  double i0 = params["Infectious"];
  double r0 = params["Recovered"];
  double prop_children = params["prop_children"];

  S(0,0) = s0 * prop_children;
  S(0,1) = s0 * (1-prop_children);
  I1(0,0) = i0 * prop_children;
  I1(0,1) = i0 * (1-prop_children);
  R1(0,0) = r0 * prop_children;
  R1(0,1) = r0 * (1-prop_children);

  // initial values
  double prop_report = params["prop_report"]; // 1 / latent period
  double epsilon = params["epsilon"]; // 1 / latent period
  double kappa = params["kappa"]; // excretion rate
  double xi = params["xi"]; // decay rate
  double K = params["K"]; // half-infective dose (eg, 10,000 doses/ml)
  double gamma = params["gamma"]; // 1 / recovery period
  double sigma = params["sigma"]; // 1 / duration of natural immunity
  double sigma_v1 = params["sigma_v1"]; // 1 / duration of OCV-induced immunity (1st dose)
  double sigma_v2 = params["sigma_v2"]; // 1 / duration of OCV-induced immunity (2nd dose)
  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double fA = params["fA"]; // fraction asymptomatic
  double bA = params["bA"]; // relative infectivity of A to I
  double R0 = params["R0"];
  double R0W = params["R0W"]; // transmission rate arising from water
  double alpha = params["alpha"]; //proportional reduction in R0
  // time window over which the number of cases is tracked
  double case_track_window = params["case_track_window"];
  // case threshold over which intervention will be implemented
  double case_threshold = params["case_threshold"];

  NumericVector vacc_cov = params["vacc_cov"]; // vaccine coverage
  // Rprintf("vacc coverage : %.2f \n", vacc_cov[0]);
  // Rprintf("vacc coverage : %.2f \n", vacc_cov[1]);
  // Rprintf("the value of R0 : %.3f \n", R0);
  double campaign_dur = params["campaign_dur"]; // duration of a vaccination campaign
  double delay_until_2nd = params["delay_until_2nd_campaign"];

  NumericVector VE_v1 = params["vacc_eff_v1"]; // vaccine efficacy (1st dose) for 5+ yo
  NumericVector VE_v2 = params["vacc_eff_v2"]; // vaccine efficacy (2nd dose) for 5+ yo
  // double vacc_eff_u5 = params["vacc_eff_u5"]; // vaccine efficacy for 5+ yo
  // daily vaccination rate during a campaign_dur is given by the following relationship
  // -log(1-vacc_cov)/campaign_dur, where vacc_cov is between 0 and 1
  // double day_intervention = params["day_intervention"];

  // double beta = R0 / (bA*fA/gamma + (1-fA)/gamma);
  double beta = R0 * gamma;
  // double betaW = beta;
  double pop_init = sum(S(0,_)) + sum(E1(0,_)) + sum(E2(0,_)) + sum(I1(0,_))
    + sum(I2(0,_)) + sum(A1(0,_)) + sum(A2(0,_)) + sum(R1(0,_)) + sum(R2(0,_))
    + sum(O1(0,_)) + sum(O2(0,_)) + sum(T1(0,_)) + sum(T2(0,_));
  double betaW = R0W * gamma * (xi / kappa) * K / pop_init; // needs

  NumericVector vacc_rate(2); // vaccination rate
  // double vacc_rate_per_tau = 0;

  double case_tracked = 0;
  bool intervention_in_place = false;
  // arbitrary small number (at least abs value is larger than the simulation days)
  double campaign_1_start = -1e5;
  double campaign_2_start = -1e5;
  double campaign_1_end = -1e5;
  double campaign_2_end = -1e5;

  for (int step = 0; step < nsteps - 1; step++) {
    // --------------- Intervention ----------------------------
    double current_t = step*tau;
    if (case_tracked >= case_threshold) {
      beta = (1 - alpha) * R0 / (bA*fA/gamma + (1-fA)/gamma);
      intervention_in_place = true;
      if(campaign_1_start < 0){ // initial values assigned to be zero
        campaign_1_start = current_t;
        campaign_2_start = campaign_1_start + delay_until_2nd;
        campaign_1_end = campaign_1_start + campaign_dur;
        campaign_2_end = campaign_2_start + campaign_dur;
        // Rprintf("the first campaign : %.1f\n", campaign_1_start);
        // Rprintf("the second campaign : %.1f\n", campaign_2_start);
      }
      // Rprintf("case tracked: %.1f \n", case_tracked);
    }
    // reset to case_tracked to zero every case_track_window
    if((step % (int) ceil(case_track_window/tau)) == 0 && !intervention_in_place) {
      case_tracked = 0;
    }

    // vaccination related
    if (campaign_1_start < current_t && current_t <= campaign_1_end) {
      vacc_rate[0] = -log(1-vacc_cov[0])/campaign_dur;
      // vacc_rate_per_tau = 1 - exp(log(1-vacc_cov)/(campaign_dur/tau));
    }
    else if (campaign_2_start < current_t && current_t <= campaign_2_end) {
      vacc_rate[1] = -log(1-vacc_cov[1])/campaign_dur;
      // vacc_rate_per_tau = 1 - exp(log(1-vacc_cov)/(campaign_dur/tau));
    }
    else {
      vacc_rate[0] = 0;
      vacc_rate[1] = 0;
    }
    // ----------------------------------------------------------
    // current state variables
    NumericVector St = S(step, _);
    NumericVector E1t = E1(step, _);
    NumericVector I1t = I1(step, _);
    NumericVector A1t = A1(step, _);
    NumericVector R1t = R1(step, _);
    NumericVector E2t = E2(step, _);
    NumericVector I2t = I2(step, _);
    NumericVector A2t = A2(step, _);
    NumericVector R2t = R2(step, _);

    double iW = W[step];
    NumericVector CEt = CE(step, _);
    NumericVector CIt = CI(step, _);

    NumericVector O1t = O1(step, _);
    NumericVector O2t = O2(step, _);
    NumericVector T1t = T1(step, _);
    NumericVector T2t = T2(step, _);

    double N = sum(St) + sum(E1t) + sum(I1t) + sum(A1t)
      + sum(R1t) + sum(E2t) + sum(I2t) + sum(A2t)
      + sum(R2t) + sum(O1t) + sum(O2t) + sum(T1t) + sum(T2t);
    double isum = sum(I1t) + sum(I2t);
    double asum = sum(A1t) + sum(A2t);
    double foi = beta * (isum + bA * asum) / N + betaW * iW / (K + iW);
    //
    NumericVector StoE1 = St * foi * tau;
    NumericVector E1toE2 = E1t * 2 * epsilon * tau;
    NumericVector E2toI1 = E2t * (1-fA) * 2 * epsilon * tau;
    NumericVector E2toA1 = E2t * fA * 2 * epsilon * tau;

    NumericVector I1toI2 = I1t * 2 * gamma * tau;
    NumericVector A1toA2 = A1t * 2 * gamma * tau;
    NumericVector I2toR1 = I2t * 2 * gamma * tau;
    NumericVector A2toR1 = A2t * 2 * gamma * tau;
    NumericVector R1toR2 = R1t * 2 * sigma * tau;
    NumericVector R2toS = R2t * 2 * sigma * tau;

    NumericVector O1toE1 = O1t * (1 - VE_v1) * foi * tau;
    NumericVector O2toE1 = O2t * (1 - VE_v1) * foi * tau;
    NumericVector T1toE1 = T1t * (1 - VE_v2) * foi * tau;
    NumericVector T2toE1 = T2t * (1 - VE_v2) * foi * tau;
    NumericVector O1toO2 = O1t * 2 * sigma_v1 * tau;
    NumericVector O2toS = O2t * 2 * sigma_v1 * tau;
    NumericVector T1toT2 = T1t * 2 * sigma_v2 * tau;
    NumericVector T2toS = T2t * 2 * sigma_v2 * tau;

    NumericVector StoO1 = St * (vacc_rate[0] + vacc_rate[1]) * tau;
    NumericVector O1toT1 = O1t * vacc_rate[1] * tau;
    NumericVector O2toT1 = O2t * vacc_rate[1] * tau;

    // Calculate the change in each state variable
    // Update next timestep
    S(step + 1, _) = St - StoE1 - StoO1 + R2toS + O2toS + T2toS;
    E1(step + 1, _) = E1t + StoE1 + O1toE1 + O2toE1 + T1toE1 + T2toE1 - E1toE2;
    E2(step + 1, _) = E2t + E1toE2 - E2toI1 - E2toA1;
    I1(step + 1, _) = I1t + E2toI1 - I1toI2;
    I2(step + 1, _) = I2t + I1toI2 - I2toR1;
    A1(step + 1, _) = A1t + E2toA1 - A1toA2;
    A2(step + 1, _) = A2t + A1toA2 - A2toR1;
    R1(step + 1, _) = R1t + I2toR1 + A2toR1 - R1toR2;
    R2(step + 1, _) = R2t + R1toR2 - R2toS;

    O1(step + 1, _) = O1t + StoO1 - O1toE1 - O1toO2 - O1toT1;
    O2(step + 1, _) = O2t + O1toO2 - O2toE1 - O2toS - O2toT1;

    T1(step + 1, _) = T1t + O1toT1 + O2toT1 - T1toE1 - T1toT2;
    T2(step + 1, _) = T2t + T1toT2 - T2toS - T2toE1;

    CE(step + 1, _) = CEt + StoE1;
    CI(step + 1, _) = CIt + prop_report * E2toI1;

    case_tracked = case_tracked + sum(E2toI1);
    // this will later need to be adjusted to account for the delay of observation
    time[step + 1] = (step + 1) * tau;// time in fractional years
  }

  // Return results as data.frame
  DataFrame result = DataFrame::create(
    Named("time") = time,
    Named("S") = S,
    Named("E1") = E1,
    Named("E2") = E2,
    Named("I1") = I1,
    Named("I2") = I2,
    Named("A1") = A1,
    Named("A2") = A2,
    Named("R1") = R1,
    Named("R2") = R2,
    Named("O1") = O1,
    Named("O2") = O2,
    Named("T1") = T1,
    Named("T2") = T2,
    Named("CE") = CE,
    Named("CI") = CI);

  return result;
}


List seiarw_2ag_erlang2(List params) {
  double tau = params["tau"]; // time step size
  double ndays = params["ndays"]; // number of days for output
  int nsteps = ceil(ndays / tau) + 1;
  // vectors for state variables
  // waiting times are modeled as Erlang with shape=2
  // two age groups implemented as two columns of the matrix
  // infectious stages are modeled as two stage to account for
  NumericMatrix S(nsteps,2); //
  NumericMatrix E1(nsteps,2); // exposed, the first stage
  NumericMatrix E2(nsteps,2); // exposed, the second stage
  NumericMatrix I1(nsteps,2); // infectious (also means symptomatic)
  NumericMatrix I2(nsteps,2); //
  NumericMatrix A1(nsteps,2); // Asymptomatic
  NumericMatrix A2(nsteps,2); //
  NumericMatrix R1(nsteps,2); // Recovered, immune, and so removed
  NumericMatrix R2(nsteps,2); //

  NumericMatrix CE(nsteps,2); // cumulative infection
  NumericMatrix CI(nsteps,2); // cumulative symptomatic
  // vaccination assumes that vaccination does not change the status of E, I, R
  NumericMatrix O1(nsteps,2); // susceptibles who received the one dose of OCV (first stage)
  NumericMatrix O2(nsteps,2); // susceptibles who received the one dose of OCV
  NumericMatrix T1(nsteps,2); // susceptibles who received the second dose of OCV
  NumericMatrix T2(nsteps,2); // susceptibles who received the second dose of OCV

  NumericVector W(nsteps); // recovered
  NumericVector time(nsteps);// fill time w/zeros


  double s0 = params["Susceptible"];
  double i0 = params["Infectious"];
  double r0 = params["Recovered"];
  double prop_children = params["prop_children"];

  S(0,0) = s0 * prop_children;
  S(0,1) = s0 * (1-prop_children);
  I1(0,0) = i0 * prop_children;
  I1(0,1) = i0 * (1-prop_children);
  R1(0,0) = r0 * prop_children;
  R1(0,1) = r0 * (1-prop_children);

  // initial values
  double epsilon = params["epsilon"]; // 1 / latent period
  double kappa = params["kappa"]; // excretion rate
  double xi = params["xi"]; // decay rate
  double K = params["K"]; // half-infective dose (eg, 10,000 doses/ml)
  double gamma = params["gamma"]; // 1 / recovery period
  double sigma = params["sigma"]; // 1 / duration of natural immunity
  double sigma_v1 = params["sigma_v1"]; // 1 / duration of OCV-induced immunity (1st dose)
  double sigma_v2 = params["sigma_v2"]; // 1 / duration of OCV-induced immunity (2nd dose)
  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double fA = params["fA"]; // fraction asymptomatic
  double bA = params["bA"]; // relative infectivity of A to I
  double R0 = params["R0"];
  double R0W = params["R0W"]; // transmission rate arising from water
  double alpha = params["alpha"]; //proportional reduction in R0
  // time window over which the number of cases is tracked
  double case_track_window = params["case_track_window"];
  // case threshold over which intervention will be implemented
  double case_threshold = params["case_threshold"];

  NumericVector vacc_cov = params["vacc_cov"]; // vaccine coverage
  // Rprintf("vacc coverage : %.2f \n", vacc_cov[0]);
  // Rprintf("vacc coverage : %.2f \n", vacc_cov[1]);
  // Rprintf("the value of R0 : %.3f \n", R0);
  double campaign_dur = params["campaign_dur"]; // duration of a vaccination campaign
  double delay_until_2nd = params["delay_until_2nd_campaign"];

  NumericVector VE_v1 = params["vacc_eff_v1"]; // vaccine efficacy (1st dose) for 5+ yo
  NumericVector VE_v2 = params["vacc_eff_v2"]; // vaccine efficacy (2nd dose) for 5+ yo
  // double vacc_eff_u5 = params["vacc_eff_u5"]; // vaccine efficacy for 5+ yo
  // daily vaccination rate during a campaign_dur is given by the following relationship
  // -log(1-vacc_cov)/campaign_dur, where vacc_cov is between 0 and 1
  // double day_intervention = params["day_intervention"];

  // double beta = R0 / (bA*fA/gamma + (1-fA)/gamma);
  double beta = R0 * gamma;
  // double betaW = beta;
  double pop_init = sum(S(0,_)) + sum(E1(0,_)) + sum(E2(0,_)) + sum(I1(0,_))
    + sum(I2(0,_)) + sum(A1(0,_)) + sum(A2(0,_)) + sum(R1(0,_)) + sum(R2(0,_))
    + sum(O1(0,_)) + sum(O2(0,_)) + sum(T1(0,_)) + sum(T2(0,_));
    double betaW = R0W * gamma * (xi / kappa) * K / pop_init; // needs

    NumericVector vacc_rate(2); // vaccination rate
    // double vacc_rate_per_tau = 0;

    double case_tracked = 0;
    bool intervention_in_place = false;
    // arbitrary small number (at least abs value is larger than the simulation days)
    double campaign_1_start = -1e5;
    double campaign_2_start = -1e5;
    double campaign_1_end = -1e5;
    double campaign_2_end = -1e5;

    for (int step = 0; step < nsteps - 1; step++) {
      // --------------- Intervention ----------------------------
      double current_t = step*tau;
      if (case_tracked >= case_threshold) {
        beta = (1 - alpha) * R0 / (bA*fA/gamma + (1-fA)/gamma);
        intervention_in_place = true;
        if(campaign_1_start < 0){ // initial values assigned to be zero
          campaign_1_start = current_t;
          campaign_2_start = campaign_1_start + delay_until_2nd;
          campaign_1_end = campaign_1_start + campaign_dur;
          campaign_2_end = campaign_2_start + campaign_dur;
          // Rprintf("the first campaign : %.1f\n", campaign_1_start);
          // Rprintf("the second campaign : %.1f\n", campaign_2_start);
        }
        // Rprintf("case tracked: %.1f \n", case_tracked);
      }
      // reset to case_tracked to zero every case_track_window
      if((step % (int) ceil(case_track_window/tau)) == 0 && !intervention_in_place) {
        case_tracked = 0;
      }

      // vaccination related
      if (campaign_1_start < current_t && current_t <= campaign_1_end) {
        vacc_rate[0] = -log(1-vacc_cov[0])/campaign_dur;
        // vacc_rate_per_tau = 1 - exp(log(1-vacc_cov)/(campaign_dur/tau));
      }
      else if (campaign_2_start < current_t && current_t <= campaign_2_end) {
        vacc_rate[1] = -log(1-vacc_cov[1])/campaign_dur;
        // vacc_rate_per_tau = 1 - exp(log(1-vacc_cov)/(campaign_dur/tau));
      }
      else {
        vacc_rate[0] = 0;
        vacc_rate[1] = 0;
      }
      // ----------------------------------------------------------
      // current state variables
      NumericVector St = S(step, _);
      NumericVector E1t = E1(step, _);
      NumericVector I1t = I1(step, _);
      NumericVector A1t = A1(step, _);
      NumericVector R1t = R1(step, _);
      NumericVector E2t = E2(step, _);
      NumericVector I2t = I2(step, _);
      NumericVector A2t = A2(step, _);
      NumericVector R2t = R2(step, _);

      double iW = W[step];
      NumericVector CEt = CE(step, _);
      NumericVector CIt = CI(step, _);

      NumericVector O1t = O1(step, _);
      NumericVector O2t = O2(step, _);
      NumericVector T1t = T1(step, _);
      NumericVector T2t = T2(step, _);

      double N = sum(St) + sum(E1t) + sum(I1t) + sum(A1t)
        + sum(R1t) + sum(E2t) + sum(I2t) + sum(A2t)
        + sum(R2t) + sum(O1t) + sum(O2t) + sum(T1t) + sum(T2t);
        double isum = sum(I1t) + sum(I2t);
        double asum = sum(A1t) + sum(A2t);
        double foi = beta * (isum + bA * asum) / N + betaW * iW / (K + iW);
        //
        NumericVector StoE1 = St * foi * tau;
        NumericVector E1toE2 = E1t * 2 * epsilon * tau;
        NumericVector E2toI1 = E2t * (1-fA) * 2 * epsilon * tau;
        NumericVector E2toA1 = E2t * fA * 2 * epsilon * tau;

        NumericVector I1toI2 = I1t * 2 * gamma * tau;
        NumericVector A1toA2 = A1t * 2 * gamma * tau;
        NumericVector I2toR1 = I2t * 2 * gamma * tau;
        NumericVector A2toR1 = A2t * 2 * gamma * tau;
        NumericVector R1toR2 = R1t * 2 * sigma * tau;
        NumericVector R2toS = R2t * 2 * sigma * tau;

        NumericVector O1toE1 = O1t * (1 - VE_v1) * foi * tau;
        NumericVector O2toE1 = O2t * (1 - VE_v1) * foi * tau;
        NumericVector T1toE1 = T1t * (1 - VE_v2) * foi * tau;
        NumericVector T2toE1 = T2t * (1 - VE_v2) * foi * tau;
        NumericVector O1toO2 = O1t * 2 * sigma_v1 * tau;
        NumericVector O2toS = O2t * 2 * sigma_v1 * tau;
        NumericVector T1toT2 = T1t * 2 * sigma_v2 * tau;
        NumericVector T2toS = T2t * 2 * sigma_v2 * tau;

        NumericVector StoO1 = St * (vacc_rate[0] + vacc_rate[1]) * tau;
        NumericVector O1toT1 = O1t * vacc_rate[1] * tau;
        NumericVector O2toT1 = O2t * vacc_rate[1] * tau;

        // Calculate the change in each state variable
        // Update next timestep
        S(step + 1, _) = St - StoE1 - StoO1 + R2toS + O2toS + T2toS;
        E1(step + 1, _) = E1t + StoE1 + O1toE1 + O2toE1 + T1toE1 + T2toE1 - E1toE2;
        E2(step + 1, _) = E2t + E1toE2 - E2toI1 - E2toA1;
        I1(step + 1, _) = I1t + E2toI1 - I1toI2;
        I2(step + 1, _) = I2t + I1toI2 - I2toR1;
        A1(step + 1, _) = A1t + E2toA1 - A1toA2;
        A2(step + 1, _) = A2t + A1toA2 - A2toR1;
        R1(step + 1, _) = R1t + I2toR1 + A2toR1 - R1toR2;
        R2(step + 1, _) = R2t + R1toR2 - R2toS;

        O1(step + 1, _) = O1t + StoO1 - O1toE1 - O1toO2 - O1toT1;
        O2(step + 1, _) = O2t + O1toO2 - O2toE1 - O2toS - O2toT1;

        T1(step + 1, _) = T1t + O1toT1 + O2toT1 - T1toE1 - T1toT2;
        T2(step + 1, _) = T2t + T1toT2 - T2toS - T2toE1;

        CE(step + 1, _) = CEt + StoE1;
        CI(step + 1, _) = CIt + E2toI1;

        case_tracked = case_tracked + sum(E2toI1);
        // this will later need to be adjusted to account for the delay of observation
        time[step + 1] = (step + 1) * tau;// time in fractional years
    }

    // Return results as data.frame
    DataFrame result = DataFrame::create(
      Named("time") = time,
      Named("S") = S,
      Named("E1") = E1,
      Named("E2") = E2,
      Named("I1") = I1,
      Named("I2") = I2,
      Named("A1") = A1,
      Named("A2") = A2,
      Named("R1") = R1,
      Named("R2") = R2,
      Named("O1") = O1,
      Named("O2") = O2,
      Named("T1") = T1,
      Named("T2") = T2,
      Named("CE") = CE,
      Named("CI") = CI);

    return result;
}
