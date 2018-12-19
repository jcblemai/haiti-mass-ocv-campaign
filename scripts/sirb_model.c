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

if(std_W > 0.0) 
{
    dw = rgammawn(std_W, dt);  // white noise (extra-demographic stochasticity)
    foi_stoc = foi * dw/dt;      // apply stochasticity
} else 
{
    foi_stoc = foi;
}

// vaccination window 
/*
if (t >= t_vacc_start && t <= (t_vacc_end + dt)) 
{
    r_v_wdn = (r_v / (S + A + RI1 + RI2 + RI3 + RA1 + RA2 + RA3));
}
else 
{*/
    r_v_wdn = 0.0;
//}


// define transition rates for each type of event (i.e what multplies the thing)
// S compartment
rate[0] = sigma * foi_stoc;   // infections
rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
rate[2] = r_v_wdn
// I compartment
rate[3] = mu;           // natural deaths
rate[4] = alpha;        // cholera-induced deaths
rate[5] = gammaI;       // recovery from infection
// A compartment
rate[6] = mu;           // natural death
rate[7] = gammaA;       // symptoms development
rate[8] = r_v_wdn
// RI1,2,3 compartment
rate[9] = 3*rhoI;        // loss of natural immunity
rate[10] = mu;            // natural death
rate[11] = r_v_wdn
// RA1,2,3 compartment
rate[12] = 3*rhoA;        // loss of natural immunity
rate[13] = mu;            // natural death
rate[15] = r_v_wdn



// simulate all transitions
reulermultinom(3, S,   &rate[0], dt, &dN[0]);
reulermultinom(3, I,   &rate[3], dt, &dN[3]);
reulermultinom(3, A,   &rate[6], dt, &dN[6]);
reulermultinom(3, RI1, &rate[9], dt, &dN[9]);
reulermultinom(3, RI2, &rate[9], dt, &dN[9]);
reulermultinom(3, RI3, &rate[9], dt, &dN[9]);
reulermultinom(3, RA1, &rate[12], dt, &dN[12]);
reulermultinom(3, RA2, &rate[12], dt, &dN[12]);
reulermultinom(3, RA3, &rate[12], dt, &dN[12]);



// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambda, lambdaR, rain, r, D);
k2 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambda, lambdaR, rain, r, D);
k3 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambda, lambdaR, rain, r, D);
k4 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambda, lambdaR, rain, r, D);
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