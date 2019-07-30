previous_vacc_campaign = TRUE ; 
r_v_wdn = 0.0;    
p1d = 0;
mobility = 0;
pdd = 0;
t_eff =  0;
t_eff_alt = 0;
dw = 0;

mobility = BArtibonite / (1 + BArtibonite) +
           BCentre / (1 + BCentre) +
           BGrande_Anse / (1 + BGrande_Anse) +
           BNippes / (1 + BNippes) +
           BNord / (1 + BNord) +
           BNord_Est / (1 + BNord_Est) +
           BNord_Ouest / (1 + BNord_Ouest) +
           BOuest / (1 + BOuest) +
           BSud / (1 + BSud) +
           BSud_Est / (1 + BSud_Est) -
           B%s / (1 + B%s);

foi = betaB%s * ( (1 - m%s) * B%s / (1 + B%s) + m%s * mobility);


if(std_W > 0.0)
{
    dw = rgammawn(std_W, dt);   // white noise (extra-demographic stochasticity)
    foi_stoc = foi * dw/dt;      // apply stochasticity
} else
{
    foi_stoc = foi;
}

if (t <= (t_vacc_end_alt%s + dt)){
	previous_vacc_campaign = TRUE;
	if (t >= t_vacc_start_alt%s && t <= (t_vacc_end_alt%s + dt)) {
    	r_v_wdn = (r_v_year_alt%s / (S%s + A%s + R1%s + R2%s + R3%s));
	}
	p1d = p1d_reg_alt%s;
} else {
	previous_vacc_campaign = FALSE;
	if (t >= t_vacc_start%s && t <= (t_vacc_end%s + dt)) {
    	r_v_wdn = (r_v_year%s / (S%s + A%s + R1%s + R2%s + R3%s));
	}
	p1d = p1d_reg%s;
}
pdd = 1 - p1d;



// time in the vacc_eff referential. We assume different timing for 1d and 2d
t_eff =     t - (t_vacc_start%s + (t_vacc_end%s - t_vacc_start%s)/2);
t_eff_alt = t - (t_vacc_start_alt%s + (t_vacc_end_alt%s - t_vacc_start_alt%s)/2);


// define transition rates for each type of event (i.e what multplies the thing)
// S compartment
rate[0] = sigma * foi_stoc;         // infections
rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
rate[2] = p1d * r_v_wdn;
rate[3] = pdd * r_v_wdn;
// I compartment
rate[4] = mu;                       // natural deaths
rate[5] = alpha;                    // cholera-induced deaths
rate[6] = gamma;                    // recovery from infection
// A compartment
rate[7] = mu;                       // natural death
rate[8] = gamma;                    // symptoms development
rate[9] = p1d * r_v_wdn;
rate[10] = pdd * r_v_wdn;
// R1,2,3 compartment
rate[11] = 3*rho;                  // loss of natural immunity
rate[12] = mu;                      // natural death
rate[13] = p1d * r_v_wdn;
rate[14] = pdd * r_v_wdn;
// V1d_S compartments
rate[15] = sigma       * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc; // symptomatic infections
rate[16] = (1 - sigma) * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc; // asymptomatic infections
rate[17] = mu;                                                       // natural death
// V1d_R1,2,3 + V2d_R1,2,3 + alt(V1d_R1,2,3  and V2d_R1,2,3 compartment)
rate[18] = mu;                                                       // natural death
rate[19] = 3*rho;
// V2d_S compartments
rate[20] = sigma * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc;       // symptomatic infections
rate[21] = (1 - sigma) * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc; // asymptomatic infections
rate[22] = mu;          // natural death
/* For previous vacc campagain */
// V1d_S compartments
rate[23] = sigma       * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // symptomatic infections
rate[24] = (1 - sigma) * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // asymptomatic infections
rate[25] = mu;                                                           // natural death
// V2d_S compartments
rate[26] = sigma *       (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // symptomatic infections
rate[27] = (1 - sigma) * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // asymptomatic infections
rate[28] = mu;          // natural death


// simulate all transitions
/* Probably we can reuse the rates (because const in C function)
but the dN should be different */
reulermultinom(4, S%s,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I%s,     &rate[4],  dt, &dN[4]);
reulermultinom(4, A%s,     &rate[7],  dt, &dN[7]);
reulermultinom(4, R1%s,   &rate[11], dt, &dN[11]);
reulermultinom(4, R2%s,   &rate[11], dt, &dN[15]);
reulermultinom(4, R3%s,   &rate[11], dt, &dN[19]);
/* Vaccinated 1 dose */
reulermultinom(3, VSd%s,  &rate[15], dt, &dN[23]);
reulermultinom(2, VR1d%s, &rate[18], dt, &dN[26]);
reulermultinom(2, VR2d%s, &rate[18], dt, &dN[28]);
reulermultinom(2, VR3d%s, &rate[18], dt, &dN[30]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd%s, &rate[20], dt, &dN[32]);
reulermultinom(2, VR1dd%s,&rate[18], dt, &dN[35]);
reulermultinom(2, VR2dd%s,&rate[18], dt, &dN[37]);
reulermultinom(2, VR3dd%s,&rate[18], dt, &dN[39]);
/* For the previous vaccination campain */
/* Vaccinated 1 dose */
reulermultinom(3, VSd_alt%s,  &rate[23], dt, &dN[41]);
reulermultinom(2, VR1d_alt%s, &rate[18], dt, &dN[44]);
reulermultinom(2, VR2d_alt%s, &rate[18], dt, &dN[46]);
reulermultinom(2, VR3d_alt%s, &rate[18], dt, &dN[48]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd_alt%s,  &rate[36], dt, &dN[50]);
reulermultinom(2, VR1dd_alt%s, &rate[18], dt, &dN[53]);
reulermultinom(2, VR2dd_alt%s, &rate[18], dt, &dN[55]);
reulermultinom(2, VR3dd_alt%s, &rate[18], dt, &dN[57]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k2 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k3 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k4 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;

I%s   += dN[0] + dN[23] + dN[32] + dN[23+18] + dN[32+18] - dN[4] - dN[5] - dN[6];
A%s   += dN[1] + dN[24] + dN[33] + dN[24+18] + dN[33+18] - dN[7] - dN[8] - dN[9] - dN[10];
R1%s += dN[6]  + dN[8]  -  dN[11] - dN[12] - dN[13] - dN[14] ;
R2%s += dN[11] - dN[15] - dN[16] - dN[17] - dN[18];
R3%s += dN[15] - dN[19] - dN[20] - dN[21] - dN[22];

if (previous_vacc_campaign){
	VSd_alt%s   += dN[2];
	VR1d_alt%s  += dN[13]+ dN[9];
	VR2d_alt%s  += dN[17];
	VR3d_alt%s  += dN[21];

	VSdd_alt%s   += dN[3];
	VR1dd_alt%s  += dN[14]+ dN[10];
	VR2dd_alt%s  += dN[18];
	VR3dd_alt%s  += dN[22];
} else {
	VSd_alt%s   += dN[2];
	VR1d_alt%s  += dN[13]+ dN[9];
	VR2d_alt%s  += dN[17];
	VR3d_alt%s  += dN[21];

	VSdd_alt%s   += dN[3];
	VR1dd_alt%s  += dN[14]+ dN[10];
	VR2dd_alt%s  += dN[18];
	VR3dd_alt%s  += dN[22];
}

VSd%s   += dN[31] - dN[23] - dN[24] - dN[25];
VR1d%s  += - dN[26] - dN[27];
VR2d%s  += dN[27] - dN[28] - dN[29];
VR3d%s  += dN[29] - dN[30] - dN[31];

VSdd%s   +=  dN[40]  - dN[32] - dN[33] - dN[34];
VR1dd%s  += - dN[35] - dN[36];
VR2dd%s  +=  dN[36] - dN[37] - dN[38] ;
VR3dd%s  +=  dN[38] - dN[39] - dN[40];

/* *previous* vacccination campain */

VSd%s   += dN[31+18] - dN[23+18] - dN[24+18] - dN[25+18];
VR1d%s  += - dN[26+18] - dN[27+18];
VR2d%s  += dN[27+18] - dN[28+18] - dN[29+18];
VR3d%s  += dN[29+18] - dN[30+18] - dN[31+18];

VSdd%s   +=  dN[40+18]  - dN[32+18] - dN[33+18] - dN[34+18];
VR1dd%s  += - dN[35+18] - dN[36+18];
VR2dd%s  +=  dN[36+18] - dN[37+18] - dN[38+18] ;
VR3dd%s  +=  dN[38+18] - dN[39+18] - dN[40+18];


C%s   +=  dN[0] + dN[15] + dN[20] + dN[15+18] + dN[20+18];
W%s   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
B%s += (((dB) < -B%s) ? (-B%s + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S%s = nearbyint(H%s - I%s - A%s - R1%s - R2%s - R3%s -
	VSd%s - VR1d%s - VR2d%s - VR3d%s  -
	VSdd%s- VR1dd%s -VR2dd%s -VR3dd%s  -
	VSd_alt%s - VR1d_alt%s - VR2d_alt%s - VR3d_alt%s -
	VSdd_alt%s - VR1dd_alt%s - VR2dd_alt%s - VR3dd_alt%s);



IncidenceAll +=  dN[0] + dN[15] + dN[20] + dN[15+18] + dN[20+18] + dN[1] + dN[16] + dN[21] + dN[16+18] + dN[21+18];
if (!previous_vacc_campaign)
{
	DosesAll  += dN[2] + dN[13] + dN[17] + dN[21] + dN[9] + 2*(dN[3] + dN[14] + dN[18] + dN[22] + dN[10]);
}
CasesAll  +=  dN[0] + dN[15] + dN[20] + dN[15+18] + dN[20+18];

