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

if (t >= t_vacc_start && t <= (t_vacc_end + dt))
{
    r_v_wdn = (r_v / (S + A + RI1 + RI2 + RI3 + RA1 + RA2 + RA3));
}
else
{
    r_v_wdn = 0.0;
}

// time in the vacc_eff referential. We assume different timing for 1d and 2d
t_eff = t - (t_vacc_start + (t_vacc_end - t_vacc_start)/2)

// define transition rates for each type of event (i.e what multplies the thing)
// S compartment
rate[0] = sigma * foi_stoc;   // infections
rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
rate[2] = p1d * r_v_wdn
rate[3] = pdd * r_v_wdn
// I compartment
rate[4] = mu;           // natural deaths
rate[5] = alpha;        // cholera-induced deaths
rate[6] = gammaI;       // recovery from infection
// A compartment
rate[7] = mu;           // natural death
rate[8] = gammaA;       // symptoms development
rate[9] = p1d * r_v_wdn
rate[10] = pdd * r_v_wdn
// RI1,2,3 compartment
rate[11] = 3*rhoI;        // loss of natural immunity
rate[12] = mu;            // natural death
rate[13] = p1d * r_v_wdn
rate[14] = pdd * r_v_wdn
// RA1,2,3 compartment
rate[15] = 3*rhoA;        // loss of natural immunity
rate[16] = mu;            // natural death
rate[17] = p1d * r_v_wdn
rate[18] = pdd * r_v_wdn
// V1d_S compartments
rate[19] = sigma       * (1 - eff_v_1d(t_eff)) * foi_stoc; // symptomatic infections
rate[20] = (1 - sigma) * (1 - eff_v_1d(t_eff)) * foi_stoc; // asymptomatic infections
rate[21] = mu;          // natural death
// V1d_RI1,2,3 compartment and V2d_RI1,2,3 compartment
rate[22] = mu;          // natural death
rate[23] = 3*rhoI
// V1d_RA1,2,3 compartment and V2d_RA1,2,3 compartment
rate[24] = mu;          // natural death
rate[25] = 3*rhoA
// V2d_S compartments
rate[26] = sigma * (1 - eff_v_2d(t_eff)) * foi_stoc; // symptomatic infections
rate[27] = (1 - sigma) * (1 - eff_v_2d(t_eff)) * foi_stoc; // asymptomatic infections
rate[28] = mu;          // natural death

/* For previous vacc campagain */
// V1d_S compartments
rate[29] = sigma       * (1 - eff_v_1d(t_eff_alt)) * foi_stoc; // symptomatic infections
rate[30] = (1 - sigma) * (1 - eff_v_1d(t_eff_alt)) * foi_stoc; // asymptomatic infections
rate[31] = mu;          // natural death
// V1d_RI1,2,3 compartment and V2d_RI1,2,3 compartment
rate[32] = mu;          // natural death
rate[33] = 3*rhoI
// V1d_RA1,2,3 compartment and V2d_RA1,2,3 compartment
rate[34] = mu;          // natural death
rate[35] = 3*rhoA
// V2d_S compartments
rate[36] = sigma *       (1 - eff_v_2d(t_eff_alt)) * foi_stoc; // symptomatic infections
rate[37] = (1 - sigma) * (1 - eff_v_2d(t_eff_alt)) * foi_stoc; // asymptomatic infections
rate[38] = mu;          // natural death



// simulate all transitions
/* Probably we can reuse the rates (because const in C function)
but the dN should be different */
reulermultinom(4, S,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I,     &rate[4],  dt, &dN[4]);
reulermultinom(4, A,     &rate[7],  dt, &dN[7]);
reulermultinom(4, RI1,   &rate[11], dt, &dN[11]);
reulermultinom(4, RI2,   &rate[11], dt, &dN[15]);
reulermultinom(4, RI3,   &rate[11], dt, &dN[19]);
reulermultinom(4, RA1,   &rate[15], dt, &dN[23]);
reulermultinom(4, RA2,   &rate[15], dt, &dN[27]);
reulermultinom(4, RA3,   &rate[15], dt, &dN[31]);
/* Vaccinated 1 dose */
reulermultinom(3, VSd,   &rate[19], dt, &dN[35]);
reulermultinom(2, VRI1d, &rate[22], dt, &dN[38]);
reulermultinom(2, VRI2d, &rate[22], dt, &dN[40]);
reulermultinom(2, VRI3d, &rate[22], dt, &dN[42]);
reulermultinom(2, VRA1d, &rate[24], dt, &dN[44]);
reulermultinom(2, VRA2d, &rate[24], dt, &dN[46]);
reulermultinom(2, VRA3d, &rate[24], dt, &dN[48]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd,  &rate[26], dt, &dN[50]);
reulermultinom(2, VRI1dd,&rate[22], dt, &dN[53]);
reulermultinom(2, VRI2dd,&rate[22], dt, &dN[55]);
reulermultinom(2, VRI3dd,&rate[22], dt, &dN[57]);
reulermultinom(2, VRA1dd,&rate[24], dt, &dN[59]);
reulermultinom(2, VRA2dd,&rate[24], dt, &dN[61]);
reulermultinom(2, VRA3dd,&rate[24], dt, &dN[63]);

/* For the previous vaccination campain */
/* Vaccinated 1 dose */
reulermultinom(3, VSd,   &rate[29], dt, &dN[65]);
reulermultinom(2, VRI1d, &rate[32], dt, &dN[68]);
reulermultinom(2, VRI2d, &rate[32], dt, &dN[70]);
reulermultinom(2, VRI3d, &rate[32], dt, &dN[72]);
reulermultinom(2, VRA1d, &rate[34], dt, &dN[74]);
reulermultinom(2, VRA2d, &rate[34], dt, &dN[76]);
reulermultinom(2, VRA3d, &rate[34], dt, &dN[78]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd,  &rate[36], dt, &dN[80]);
reulermultinom(2, VRI1dd,&rate[32], dt, &dN[83]);
reulermultinom(2, VRI2dd,&rate[32], dt, &dN[85]);
reulermultinom(2, VRI3dd,&rate[32], dt, &dN[87]);
reulermultinom(2, VRA1dd,&rate[34], dt, &dN[89]);
reulermultinom(2, VRA2dd,&rate[34], dt, &dN[91]);
reulermultinom(2, VRA3dd,&rate[34], dt, &dN[93]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k2 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k3 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k4 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;


// update state variables
S = -dN[0] - dN[1] /* FOI */
- dN[2] - dN[3] /* VACC */
+ dN[4] + dN[7] + dN[12] + dN[16] + dN[20] + dN[24] + dN[28] + dN[32] /* Mortality of nn vacc: rebirth*/
+ dN[19] + dN[31]  /* Recovery */
+ dN[37] + dN[38] + dN[40] + dN[42] + dN[44] + dN[46] + dN[48]        /* Mortality of 1d vacc: rebirth*/
+ dN[52] + dN[53] + dN[55] + dN[57] + dN[59] + dN[61] + dN[63]		   /* Mortality of dd vacc: rebirth*/
+ dN[37+20] + dN[38+20] + dN[40+20] + dN[42+20] + dN[44+20] + dN[46+20] + dN[48+20]        /* Prev campain*/
+ dN[52+20] + dN[53+20] + dN[55+20] + dN[57+20] + dN[59+20] + dN[61+20] + dN[63+20];		   

I   += dN[0] + dN[35] + dN[50] + dN[35+20] + dN[50+20] - dN[4] - dN[5] - dN[6];
A   += dN[1] + dN[36] + dN[51] + dN[36+20] + dN[51+20] - dN[7] - dN[8] - dN[9] - dN[10];
RI1 += dN[6] -  dN[11] - dN[12] - dN[13] - dN[14] ;
RI2 += dN[11] - dN[15] - dN[16] - dN[17] - dN[18];
RI3 += dN[15] - dN[19] - dN[20] - dN[21] - dN[22];
RA1 += dN[8]  - dN[23] - dN[24] - dN[25] - dN[26];
RA2 += dN[23] - dN[27] - dN[28] - dN[29] - dN[30];
RA3 += dN[27] - dN[31] - dN[32] - dN[33] - dN[34];

if (previous_vacc_campaign){
	VSd_alt    += dN[2];
	VRI1d_alt  += dN[13];
	VRI2d_alt  += dN[17];
	VRI3d_alt  += dN[21];
	VRA1d_alt  += dN[9] + dN[25] ;
	VRA2d_alt  += dn[29];
	VRA3d_alt  += dN[33];
	
	VSdd_alt    += dN[3];
	VRI1dd_alt  += dN[14];
	VRI2dd_alt  += dN[18];
	VRI3dd_alt  += dN[22];
	VRA1dd_alt  += dN[10] + dN[26];
	VRA2dd_alt  += dN[30];
	VRA3dd_alt  += dN[34];
} else {
	VSd    += dN[2];
	VRI1d  += dN[13];
	VRI2d  += dN[17];
	VRI3d  += dN[21];
	VRA1d  += dN[9] + dN[25] ;
	VRA2d  += dn[29];
	VRA3d  += dN[33];
	
	VSdd    += dN[3];
	VRI1dd  += dN[14];
	VRI2dd  += dN[18];
	VRI3dd  += dN[22];
	VRA1dd  += dN[10] + dN[26];
	VRA2dd  += dN[30];
	VRA3dd  += dN[34];

}

VSd    += dN[43] + dN[49] - dN[35] - dN[36] - dN[37];
VRI1d  += - dN[38] - dN[39];
VRI2d  += dN[39] - dN[40] - dN[41];
VRI3d  += dN[41] - dN[42] - dN[43];
VRA1d  += - dN[44] - dN[45];
VRA2d  += dN[45] - dN[46] - dN[47];
VRA3d  += dN[47] - dN[48] - dN[49];

VSdd    +=  dN[58] + dN[64] - dN[50] - dN[51] - dN[52];
VRI1dd  += - dN[53] - dN[54];
VRI2dd  +=  dN[54] - dN[55] - dN[56] ;
VRI3dd  +=  dN[56] - dN[57] - dN[58];
VRA1dd  += - dN[59] - dN[60];
VRA2dd  +=  dN[60] - dN[61] - dN[62];
VRA3dd  +=  dN[62] - dN[63] - dN[64];

/* *previous* vacccination campain */

VSd_alt    += dN[43+20] + dN[49+20] - dN[35+20] - dN[36+20] - dN[37+20];
VRI1d_alt  += - dN[38+20] - dN[39+20];
VRI2d_alt  += dN[39+20] - dN[40+20] - dN[41+20];
VRI3d_alt  += dN[41+20] - dN[42+20] - dN[43+20];
VRA1d_alt  += - dN[44+20] - dN[45+20];
VRA2d_alt  += dN[45+20] - dN[46+20] - dN[47+20];
VRA3d_alt  += dN[47+20] - dN[48+20] - dN[49+20];

VSdd_alt    +=  dN[58+20] + dN[64+20] - dN[50+20] - dN[51+20] - dN[52+20];
VRI1dd_alt  += - dN[53+20] - dN[54+20];
VRI2dd_alt  +=  dN[54+20] - dN[55+20] - dN[56+20] ;
VRI3dd_alt  +=  dN[56+20] - dN[57+20] - dN[58+20];
VRA1dd_alt  += - dN[59+20] - dN[60+20];
VRA2dd_alt  +=  dN[60+20] - dN[61+20] - dN[62+20];
VRA3dd_alt  +=  dN[62+20] - dN[63+20] - dN[64+20];






C   +=  dN[0] + dN[35] + dN[50] + dN[35+20] + dN[50+20];
W   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
B += (((dB) < -B) ? (-B + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
//S = nearbyint(H - I - A - RI1 - RI2 - RI3 - RA1 - RA2 - RA3);
