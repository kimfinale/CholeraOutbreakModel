#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector reulermultinom(double size, NumericVector rate, double dt) {
  int ncol = rate.size();
  NumericVector trans(ncol); // transition events
  double p = sum(rate); //total event rate
  double tmpp = p;
  double tmpsize = R::rbinom(size, (1 - exp(-tmpp*dt))); // total number of events
  for (int k = 0; k < (ncol-1); k++) {
    double tr = R::rbinom(tmpsize, rate(k)/p);
    trans(k) = tr;
    tmpsize -= tr;
    tmpp -= rate(k);
  }
  trans(ncol-1) = tmpsize;
  return(trans);
}


// [[Rcpp::export]]
List seiarw_euler(List params) {
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

  NumericVector time(nsteps);// fill time w/zeros

  // initial values
  List init = params["init"];
  S(0) = init["S"];
  E(0) = init["E"];
  I(0) = init["I"];
  A(0) = init["A"]; //asymptomatic
  R(0) = init["R"]; // recovered
  W(0) = init["W"];
  CE(0) = init["CE"]; // cumulative infection
  CI(0) = init["CI"]; // cumulative symptomatic

  double epsilon = params["epsilon"]; // 1 / latent period
  double kappa = params["kappa"]; // excretion rate
  double xi = params["xi"]; // decay rate
  double gamma = params["gamma"]; // 1 / recovery period

  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double fA = params["fA"]; // fraction asymptomatic
  double bA = params["bA"]; // relative infectivity of A to I
  double R0 = params["R0"];
  double R0W = params["R0W"]; // fraction of R0 when intervention is in place
  // double day_intervention = params["day_intervention"];

  double beta = R0 / (bA*fA/gamma + (1-fA)/gamma);
  double betaW = R0W * (xi / kappa) * gamma; // needs

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
    double iS = S[istep];
    double iE = E[istep];
    double iI = I[istep];
    double iA = A[istep];
    double iR = R[istep];
    double iW = W[istep];

    double iCE = CE[istep];
    double iCI = CI[istep];

    // State Equations
    double N = iS + iE + iI + iA + iR;
    double foi = beta * (bA * iA + iI) / N + betaW * iW;
    // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
    double new_infection = iS * foi * tau;

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
    // Rprintf("the value of ItoR : %.1f \n", ItoR);
    // Rprintf("the value of AtoR : %.1f \n", AtoR);

    // Calculate the change in each state variable
    double dS = - new_infection;
    double dE = new_infection - EtoA - EtoI;
    double dI = EtoI - ItoR;
    double dA = EtoA - AtoR;
    double dR = ItoR + AtoR;
    double dW = ItoW + AtoW - fromW;
    // Update next timestep
    S[istep + 1] = iS + dS;
    E[istep + 1] = iE + dE;
    I[istep + 1] = iI + dI;
    A[istep + 1] = iA + dA;
    R[istep + 1] = iR + dR;
    W[istep + 1] = iW + dW;

    CE[istep + 1] = iCE + new_infection;// cumulative infection
    CI[istep + 1] = iCI + EtoI;// cumulative symptomatic
    time[istep + 1] = (istep + 1) * tau;// time in fractional years
  }
// Return results as data.frame
  DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S") = S,
    Named("E") = E,
    Named("I") = I,
    Named("A") = A,
    Named("R") = R,
    Named("W") = W,
    Named("CE") = CE,
    Named("CI") = CI);

  return sim;
}

