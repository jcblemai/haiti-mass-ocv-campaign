double foi, foi_stoc; // force of infection and its stochastic version
double dw;            // extra-demographic stochasticity on foi
double dB;            // deterministic forward time difference of bacteria in the environment
double k1, k2, k3, k4;  // coefficients of  the Runge-Kutta method
double rate[39];      // vector of all rates in model
double dN[94];        // vector of transitions between classes during integration timestep

double thetaA = thetaI * XthetaA;
double rhoI = rhoA * XrhoI;

int previous_vacc_campaign = TRUE ; /* flag that indicate if we are on the first or second campain */
double r_v_wdn = 0.0;       // rate of vaccination: 0 if out of time window, r_v if not
double p1d = 0;
double mobility = 0;

int scenario =  cases_ext;

double pdd = 0;
double t_eff =  0;
double t_eff_alt = 0;

