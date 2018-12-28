double foi, foi_stoc; // force of infection and its stochastic version
double dw;            // extra-demographic stochasticity on foi
double dB;            // deterministic forward time difference of bacteria in the environment
double k1, k2, k3, k4;  // coefficients of  the Runge-Kutta method
double r_v_wdn;       // rate of vaccination: 0 if out of time window, r_v if not
double rate[19];      // vector of all rates in model
double dN[19];        // vector of transitions between classes during integration timestep

double thetaA = thetaI * XthetaA;
double rhoI = rhoA * XrhoI;

// force of infection
foi = betaB * (B / (1 + B)) + foi_add;

if(std_W > 0.0) {
  // white noise (extra-demographic stochasticity)
  dw = rgammawn(std_W, dt);
  // apply stochasticity
  foi_stoc = foi * dw/dt;
} else {
  foi_stoc = foi;
}

// vaccination window
r_v_wdn = 0.0;

// define transition rates for each type of event
// S compartment
rate[0] = sigma * foi_stoc;   // infections
rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
// I compartment
rate[2] = mu;           // natural deaths
rate[3] = alpha;        // cholera-induced deaths
rate[4] = gammaI;       // recovery from infection
// A compartment (not in order because was added after initial model formulation)
rate[5] = mu;           // natural death
rate[6] = gammaA;       // symptoms development
// RI1,2,3 compartment
rate[7] = 3*rhoI;        // loss of natural immunity
rate[8] = mu;            // natural death
// RI2 compartment
rate[9] = 3*rhoI;        // loss of natural immunity
rate[10] = mu;
// RI3 compartment
rate[11] = 3*rhoI;        // loss of natural immunity
rate[12] = mu;
// RA1,2,3 compartment
rate[13] = 3*rhoA;        // loss of natural immunity
rate[14] = mu;            // natural death
// RA2 compartment
rate[15] = 3*rhoA;        // loss of natural immunity
rate[16] = mu;
// RA3 compartment
rate[17] = 3*rhoA;        // loss of natural immunity
rate[18] = mu;


// simulate all transitions
reulermultinom(2, S, &rate[0], dt, &dN[0]);
reulermultinom(3, I, &rate[2], dt, &dN[2]);
reulermultinom(2, A, &rate[5], dt, &dN[5]);
reulermultinom(2, RI1, &rate[7], dt, &dN[7]);
reulermultinom(2, RI2, &rate[9], dt, &dN[9]);
reulermultinom(2, RI3, &rate[11], dt, &dN[11]);
reulermultinom(2, RA1, &rate[13], dt, &dN[13]);
reulermultinom(2, RA2, &rate[15], dt, &dN[15]);
reulermultinom(2, RA3, &rate[17], dt, &dN[17]);


// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k2 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k3 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k4 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;

// update state variables

I   += dN[0] - dN[2] - dN[3] - dN[4];
A   += dN[1] - dN[5] - dN[6];
RI1 += dN[4] - dN[7] - dN[8];
RI2 += dN[7] - dN[9] - dN[10];
RI3 += dN[9] - dN[11] - dN[12];
RA1 += dN[6] - dN[13] - dN[14];
RA2 += dN[13] - dN[15] - dN[16];
RA3 += dN[15] - dN[17] - dN[18];
C   +=  dN[0];
W   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
B += (((dB) < -B) ? (-B + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S = nearbyint(H - I - A - RI1 - RI2 - RI3 - RA1 - RA2 - RA3);