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
  List init = params["init"];
  S(0) = init["S"];
  E(0) = init["E"];
  I(0) = init["I"];
  A(0) = init["A"];
  R(0) = init["R"];
  W(0) = init["W"];
  CE(0) = init["CE"];
  CI(0) = init["CI"];
  V1(0) = 0;
  V2(0) = 0;
  CV1(0) = 0;
  CV2(0) = 0;

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
  double vacc_cov = params["vacc_cov"]; // vaccine coverage
  double vacc_dur = params["vacc_dur"]; // vaccine duration
  double VE_v1 = params["vacc_eff_v1"]; // vaccine efficacy (1st dose) for 5+ yo
  double VE_v2 = params["vacc_eff_v2"]; // vaccine efficacy (2nd dose) for 5+ yo
  // double vacc_eff_u5 = params["vacc_eff_u5"]; // vaccine efficacy for 5+ yo
  // daily vaccination rate during a vacc_dur is given by the following relationship
  //  -log(1-vacc_cov)/vacc_dur, where vacc_cov is between 0 and 1
  // double day_intervention = params["day_intervention"];

  // double beta = R0 / (bA*fA/gamma + (1-fA)/gamma);
  double beta = R0 * gamma;
  // double betaW = beta;
  double pop_init = S[0] + E[0] + I[0] + A[0] + R[0] + V1[0] + V2[0];
  double betaW = R0W * gamma * (xi / kappa) * K / pop_init; // needs
  double vacc_rate = 0;
  // double vacc_rate = -log(1-vacc_cov)/vacc_dur;

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

  for (int istep = 0; istep < nsteps - 1; istep++) {
    // --------------- Intervention ----------------------------
    if (case_tracked >= case_threshold) {
      beta = (1 - alpha) * R0 / (bA*fA/gamma + (1-fA)/gamma);
      intervention_in_place = true;
      if(campaign_start < 0){ // initial values assigned to be zero
        campaign_start = istep*tau;
      }
    }
    if((istep % (int) ceil(case_track_window/tau)) == 0 && !intervention_in_place) {
      case_tracked = 0;
    }

    // vaccination related
    double time_since_campaign = (istep+1)*tau - campaign_start;
    if (0 < time_since_campaign && time_since_campaign <= vacc_dur) {
      vacc_rate = -log(1-vacc_cov)/vacc_dur;
    }
    else {
      vacc_rate = 0;
    }
    // -----------------------------------------------------------------
    // test vacc_rate
    // if (istep % 5 == 0) {
    //   Rprintf("step: %d, time: %.1f, vacc_rate: %.3f \n",
    //           istep, (istep+1)*tau, vacc_rate);
    // }
    double iS = S[istep];
    double iE = E[istep];
    double iI = I[istep];
    double iA = A[istep];
    double iR = R[istep];
    double iW = W[istep];

    double iCE = CE[istep];
    double iCI = CI[istep];

    double iV1 = V1[istep];
    double iV2 = V2[istep];
    double iCV1 = CV1[istep];
    double iCV2 = CV2[istep];

    // State Equations
    double N = iS + iE + iI + iA + iR + iV1 + iV2;
    double foi = beta * (bA * iA + iI) / N + betaW * iW / (K + iW);
    // Rprintf("the value of N : %.1f \n", N);
    // double new_infection = (iS + (1-VE_v1)*iV1 + (1-VE_v2)*iV2) * foi * tau;
    double StoE = iS * foi * tau;
    double RtoS = iR * sigma * tau;
    double V2toS = iV2 * sigma_v2 * tau; // immunity wane (2nd dose)
    double V1toS = iV1 * sigma_v1 * tau;
    // new vacc doses NEED TO BE ADJUSTED
    double new_vacc_first_dose = (iS+iE+iA+iR) * vacc_rate * tau;
    double new_vacc_second_dose = (iS+iE+iA+iR+iV1) * vacc_rate * tau;
    // Rprintf("the value of new infection : %f \n", new_infection);

    // Rprintf("the value of new_infection : %.1f \n", new_infection);
    // Rprintf("the value of EtoP : %.1f \n", EtoP);
    double EtoI = iE * epsilon * (1-fA) * tau;//
    double EtoA = iE * epsilon * fA * tau;//

    // Rprintf("the value of PtoI : %.1f \n", PtoI);
    // Rprintf("the value of PtoA : %.1f \n", PtoA);
    double ItoR = iI * gamma * tau;
    double AtoR = iA * gamma * tau;//

    double ItoW = iA * kappa * tau;//
    double AtoW = iA * bA * kappa * tau;//
    double fromW = iW * xi * tau;//

    double StoV1 = iS * vacc_rate * tau;//
    double V1toV2 = iV1 * vacc_rate * tau;//
    double V1toE = iV1 * (1-VE_v1) * foi * tau;
    double V2toE = iV2 * (1-VE_v2) * foi * tau;
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
    S[istep + 1] = iS + dS;
    E[istep + 1] = iE + dE;
    I[istep + 1] = iI + dI;
    A[istep + 1] = iA + dA;
    R[istep + 1] = iR + dR;
    W[istep + 1] = iW + dW;

    V1[istep + 1] = iV1 + dV1;
    V2[istep + 1] = iV2 + dV2;

    CE[istep + 1] = iCE + StoE + V1toE + V2toE;// cumulative infection
    CI[istep + 1] = iCI + EtoI;// cumulative symptomatic
    // this will later need to be adjusted to account for the delay of observation
    case_tracked = case_tracked + EtoI;

    CV1[istep + 1] = iCV1 + new_vacc_first_dose;
    CV2[istep + 1] = iCV2 + new_vacc_second_dose;
    time[istep + 1] = (istep + 1) * tau;// time in fractional years
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
List seiarw_2ag(List params) {
  double tau = params["tau"]; // time step size
  double ndays = params["ndays"]; // number of days for output
  int nsteps = ceil(ndays / tau) + 1;
  // vectors for state variables
  NumericVector S1(nsteps);
  NumericVector E1(nsteps);
  NumericVector I1(nsteps);
  NumericVector A1(nsteps);//asymptomatic
  NumericVector R1(nsteps);//recovered

  NumericVector S2(nsteps);
  NumericVector E2(nsteps);
  NumericVector I2(nsteps);
  NumericVector A2(nsteps);//asymptomatic
  NumericVector R2(nsteps);//recovered

  NumericVector W(nsteps); //recovered

  NumericVector CE1(nsteps); //cumulative infection
  NumericVector CI1(nsteps); //cumulative symptomatic

  NumericVector CE2(nsteps); //cumulative infection
  NumericVector CI2(nsteps); //cumulative infection

  NumericVector time(nsteps);//fill time w/zeros

  // initial values
  List init = params["init"];
  S1(0) = init["S1"];
  E1(0) = init["E1"];
  I1(0) = init["I1"];
  A1(0) = init["A1"]; //asymptomatic
  R1(0) = init["R1"]; // recovered

  S2(0) = init["S2"];
  E2(0) = init["E2"];
  I2(0) = init["I2"];
  A2(0) = init["A2"]; //asymptomatic
  R2(0) = init["R2"]; // recovered

  W(0) = init["W"];
  CE1(0) = init["CE1"]; // cumulative infection
  CI1(0) = init["CI1"]; // cumulative symptomatic
  CE2(0) = init["CE2"]; // cumulative infection
  CI2(0) = init["CI2"]; // cumulative symptomatic

  double epsilon = params["epsilon"]; // 1 / latent period
  double kappa = params["kappa"]; // excretion rate
  double xi = params["xi"]; // decay rate
  double K = params["K"]; // half-infective dose
  double gamma = params["gamma"]; // 1 / recovery period
  double sigma = params["sigma"]; // 1 / duration of natural immunity

  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double fA = params["fA"]; // fraction asymptomatic infectious
  double bA = params["bA"]; // relative infectivity of A to I
  double R0 = params["R0"];
  double R0W = params["R0W"]; // fraction of R0 when intervention is in place
  // double day_intervention = params["day_intervention"];

  double beta = R0 / (bA*fA/gamma + (1-fA)/gamma);
  // double betaW = R0W * (xi / kappa) * gamma; // needs
  double pop_init = S1[0]+E1[0]+I1[0]+A1[0]+R1[0]+S2[0]+E2[0]+I2[0]+A2[0]+R2[0];
  double betaW = R0W * gamma * (xi / kappa) * K / pop_init;
  // Rprintf("the value of delta : %.2f \n", delta);
  // Rprintf("the value of nsteps : %.1d \n", nsteps);
  // Rprintf("the value of kappa : %.2f \n", kappa);
  // Rprintf("the value of epsilon : %.2f \n", epsilon);
  // Rprintf("the value of gamma : %.2f \n", gamma);
  // Rprintf("the value of rate_P_I : %.2f \n", rate_P_I);
  // Rprintf("the value of P_rates[0] : %.2f \n", P_rates[0]);
  // Rprintf("the value of P_rates[1] : %.2f \n", P_rates[1]);

  // Calculate the number of events for each step, update state vectors
  for (int istep = 0; istep < nsteps - 1; istep++) {
    // if (istep*tau >= day_intervention) {
    //   beta = R0 / (bA*fA/gamma + (1-fA)/gamma);
    // }
    double iS1 = S1[istep];
    double iE1 = E1[istep];
    double iI1 = I1[istep];
    double iA1 = A1[istep];
    double iR1 = R1[istep];

    double iS2 = S2[istep];
    double iE2 = E2[istep];
    double iI2 = I2[istep];
    double iA2 = A2[istep];
    double iR2 = R2[istep];

    double iW = W[istep];

    double iCE1 = CE1[istep];
    double iCI1 = CI1[istep];
    double iCE2 = CE2[istep];
    double iCI2 = CI2[istep];

    // State Equations
    double N = iS1 + iE1 + iI1 + iA1 + iR1 + iS2 + iE2 + iI2 + iA2 + iR2;
    double foi = beta * (bA * iA1 + iI1 + bA * iA2 + iI2) / N + betaW * iW;
    // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
    double new_infection1 = iS1 * foi * tau;
    double new_infection2 = iS2 * foi * tau;

    // Rprintf("the value of new_infection : %.1f \n", new_infection);
    // Rprintf("the value of EtoP : %.1f \n", EtoP);
    double E1toI1 = iE1 * epsilon * (1-fA) * tau;//
    double E1toA1 = iE1 * epsilon * fA * tau;//
    double E2toI2 = iE2 * epsilon * (1-fA) * tau;//
    double E2toA2 = iE2 * epsilon * fA * tau;//

        // Rprintf("the value of PtoI : %.1f \n", PtoI);
    // Rprintf("the value of PtoA : %.1f \n", PtoA);
    double I1toR1 = iI1 * gamma * tau;
    double A1toR1 = iA1 * gamma * tau;//
    double I2toR2 = iI2 * gamma * tau;
    double A2toR2 = iA2 * gamma * tau;//

    double I1toW = iA1 * kappa * tau;//
    double A1toW = iA1 * bA * kappa * tau;//
    double I2toW = iA2 * kappa * tau;//
    double A2toW = iA2 * bA * kappa * tau;//


    double fromW = iW * xi * tau;//
    // Rprintf("the value of ItoR : %.1f \n", ItoR);
    // Rprintf("the value of AtoR : %.1f \n", AtoR);

    // Calculate the change in each state variable
    double dS1 = - new_infection1;
    double dE1 = new_infection1 - E1toA1 - E1toI1;
    double dI1 = E1toI1 - I1toR1;
    double dA1 = E1toA1 - A1toR1;
    double dR1 = I1toR1 + A1toR1;
    double dS2 = - new_infection2;
    double dE2 = new_infection2 - E2toA2 - E2toI2;
    double dI2 = E2toI2 - I2toR2;
    double dA2 = E2toA2 - A2toR2;
    double dR2 = I2toR2 + A2toR2;

    double dW = I1toW + A1toW + I2toW + A2toW - fromW;
    // Update next timestep
    S1[istep + 1] = iS1 + dS1;
    E1[istep + 1] = iE1 + dE1;
    I1[istep + 1] = iI1 + dI1;
    A1[istep + 1] = iA1 + dA1;
    R1[istep + 1] = iR1 + dR1;
    S2[istep + 1] = iS2 + dS2;
    E2[istep + 1] = iE2 + dE2;
    I2[istep + 1] = iI2 + dI2;
    A2[istep + 1] = iA2 + dA2;
    R2[istep + 1] = iR2 + dR2;
    W[istep + 1] = iW + dW;

    CE1[istep + 1] = iCE1 + new_infection1;// cumulative infection
    CI1[istep + 1] = iCI1 + E1toI1;// cumulative symptomatic
    CE2[istep + 1] = iCE2 + new_infection2;// cumulative infection
    CI2[istep + 1] = iCI2 + E2toI2;// cumulative symptomatic
    time[istep + 1] = (istep + 1) * tau;// time in fractional years
  }
  // Return results as data.frame
  DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S1") = S1,
    Named("E1") = E1,
    Named("I1") = I1,
    Named("A1") = A1,
    Named("R1") = R1,
    Named("CE1") = CE1,
    Named("CI1") = CI1,
    Named("S2") = S2,
    Named("E2") = E2,
    Named("I2") = I2,
    Named("A2") = A2,
    Named("R2") = R2,
    Named("CE2") = CE2,
    Named("CI2") = CI2,
    Named("W") = W);

  return sim;
}

