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
    	r_v_wdn = (r_v_year_alt%s / (S%s + A%s + RI1%s + RI2%s + RI3%s + RA1%s + RA2%s + RA3%s));
	}
	p1d = p1d_reg_alt%s;
} else {
	previous_vacc_campaign = FALSE;
	if (t >= t_vacc_start%s && t <= (t_vacc_end%s + dt)) {
    	r_v_wdn = (r_v_year%s / (S%s + A%s + RI1%s + RI2%s + RI3%s + RA1%s + RA2%s + RA3%s));
	}
	p1d = p1d_reg%s;
}
pdd = 1 - p1d;



// time in the vacc_eff referential. We assume different timing for 1d and 2d
t_eff =     t - (t_vacc_start%s + (t_vacc_end%s - t_vacc_start%s)/2);
t_eff_alt = t - (t_vacc_start_alt%s + (t_vacc_end_alt%s - t_vacc_start_alt%s)/2);


// define transition rates for each type of event (i.e what multplies the thing)
// S compartment
rate[0] = sigma * foi_stoc;   // infections
rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
rate[2] = p1d * r_v_wdn;
rate[3] = pdd * r_v_wdn;
// I compartment
rate[4] = mu;           // natural deaths
rate[5] = alpha;        // cholera-induced deaths
rate[6] = gammaI;       // recovery from infection
// A compartment
rate[7] = mu;           // natural death
rate[8] = gammaA;       // symptoms development
rate[9] = p1d * r_v_wdn;
rate[10] = pdd * r_v_wdn;
// RI1,2,3 compartment
rate[11] = 3*rhoI;        // loss of natural immunity
rate[12] = mu;            // natural death
rate[13] = p1d * r_v_wdn;
rate[14] = pdd * r_v_wdn;
// RA1,2,3 compartment
rate[15] = 3*rhoA;        // loss of natural immunity
rate[16] = mu;            // natural death
rate[17] = p1d * r_v_wdn;
rate[18] = pdd * r_v_wdn;
// V1d_S compartments
rate[19] = sigma       * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc; // symptomatic infections
rate[20] = (1 - sigma) * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc; // asymptomatic infections
rate[21] = mu;          // natural death
// V1d_RI1,2,3 compartment and V2d_RI1,2,3 compartment
rate[22] = mu;          // natural death
rate[23] = 3*rhoI;
// V1d_RA1,2,3 compartment and V2d_RA1,2,3 compartment
rate[24] = mu;          // natural death
rate[25] = 3*rhoA;
// V2d_S compartments
rate[26] = sigma * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc; // symptomatic infections
rate[27] = (1 - sigma) * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc; // asymptomatic infections
rate[28] = mu;          // natural death

/* For previous vacc campagain */
// V1d_S compartments
rate[29] = sigma       * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // symptomatic infections
rate[30] = (1 - sigma) * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // asymptomatic infections
rate[31] = mu;          // natural death
// V1d_RI1,2,3 compartment and V2d_RI1,2,3 compartment
rate[32] = mu;          // natural death
rate[33] = 3*rhoI;
// V1d_RA1,2,3 compartment and V2d_RA1,2,3 compartment
rate[34] = mu;          // natural death
rate[35] = 3*rhoA;
// V2d_S compartments
rate[36] = sigma *       (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // symptomatic infections
rate[37] = (1 - sigma) * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // asymptomatic infections
rate[38] = mu;          // natural death


// simulate all transitions
/* Probably we can reuse the rates (because const in C function)
but the dN should be different */
reulermultinom(4, S%s,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I%s,     &rate[4],  dt, &dN[4]);
reulermultinom(4, A%s,     &rate[7],  dt, &dN[7]);
reulermultinom(4, RI1%s,   &rate[11], dt, &dN[11]);
reulermultinom(4, RI2%s,   &rate[11], dt, &dN[15]);
reulermultinom(4, RI3%s,   &rate[11], dt, &dN[19]);
reulermultinom(4, RA1%s,   &rate[15], dt, &dN[23]);
reulermultinom(4, RA2%s,   &rate[15], dt, &dN[27]);
reulermultinom(4, RA3%s,   &rate[15], dt, &dN[31]);
/* Vaccinated 1 dose */
reulermultinom(3, VSd%s,   &rate[19], dt, &dN[35]);
reulermultinom(2, VRI1d%s, &rate[22], dt, &dN[38]);
reulermultinom(2, VRI2d%s, &rate[22], dt, &dN[40]);
reulermultinom(2, VRI3d%s, &rate[22], dt, &dN[42]);
reulermultinom(2, VRA1d%s, &rate[24], dt, &dN[44]);
reulermultinom(2, VRA2d%s, &rate[24], dt, &dN[46]);
reulermultinom(2, VRA3d%s, &rate[24], dt, &dN[48]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd%s,  &rate[26], dt, &dN[50]);
reulermultinom(2, VRI1dd%s,&rate[22], dt, &dN[53]);
reulermultinom(2, VRI2dd%s,&rate[22], dt, &dN[55]);
reulermultinom(2, VRI3dd%s,&rate[22], dt, &dN[57]);
reulermultinom(2, VRA1dd%s,&rate[24], dt, &dN[59]);
reulermultinom(2, VRA2dd%s,&rate[24], dt, &dN[61]);
reulermultinom(2, VRA3dd%s,&rate[24], dt, &dN[63]);
/* For the previous vaccination campain */
/* Vaccinated 1 dose */
reulermultinom(3, VSd_alt%s,   &rate[29], dt, &dN[65]);
reulermultinom(2, VRI1d_alt%s, &rate[32], dt, &dN[68]);
reulermultinom(2, VRI2d_alt%s, &rate[32], dt, &dN[70]);
reulermultinom(2, VRI3d_alt%s, &rate[32], dt, &dN[72]);
reulermultinom(2, VRA1d_alt%s, &rate[34], dt, &dN[74]);
reulermultinom(2, VRA2d_alt%s, &rate[34], dt, &dN[76]);
reulermultinom(2, VRA3d_alt%s, &rate[34], dt, &dN[78]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd_alt%s,  &rate[36], dt, &dN[80]);
reulermultinom(2, VRI1dd_alt%s,&rate[32], dt, &dN[83]);
reulermultinom(2, VRI2dd_alt%s,&rate[32], dt, &dN[85]);
reulermultinom(2, VRI3dd_alt%s,&rate[32], dt, &dN[87]);
reulermultinom(2, VRA1dd_alt%s,&rate[34], dt, &dN[89]);
reulermultinom(2, VRA2dd_alt%s,&rate[34], dt, &dN[91]);
reulermultinom(2, VRA3dd_alt%s,&rate[34], dt, &dN[93]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k2 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k3 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k4 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;


I%s   += dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30] - dN[4] - dN[5] - dN[6];
A%s   += dN[1] + dN[36] + dN[51] + dN[36+30] + dN[51+30] - dN[7] - dN[8] - dN[9] - dN[10];
RI1%s += dN[6] -  dN[11] - dN[12] - dN[13] - dN[14] ;
RI2%s += dN[11] - dN[15] - dN[16] - dN[17] - dN[18];
RI3%s += dN[15] - dN[19] - dN[20] - dN[21] - dN[22];
RA1%s += dN[8]  - dN[23] - dN[24] - dN[25] - dN[26];
RA2%s += dN[23] - dN[27] - dN[28] - dN[29] - dN[30];
RA3%s += dN[27] - dN[31] - dN[32] - dN[33] - dN[34];

if (previous_vacc_campaign){
	VSd_alt%s    += dN[2];
	VRI1d_alt%s  += dN[13];
	VRI2d_alt%s  += dN[17];
	VRI3d_alt%s  += dN[21];
	VRA1d_alt%s  += dN[9] + dN[25] ;
	VRA2d_alt%s  += dN[29];
	VRA3d_alt%s  += dN[33];

	VSdd_alt%s    += dN[3];
	VRI1dd_alt%s  += dN[14];
	VRI2dd_alt%s  += dN[18];
	VRI3dd_alt%s  += dN[22];
	VRA1dd_alt%s  += dN[10] + dN[26];
	VRA2dd_alt%s  += dN[30];
	VRA3dd_alt%s  += dN[34];
} else {
	VSd%s    += dN[2];
	VRI1d%s  += dN[13];
	VRI2d%s  += dN[17];
	VRI3d%s  += dN[21];
	VRA1d%s  += dN[9] + dN[25] ;
	VRA2d%s  += dN[29];
	VRA3d%s  += dN[33];

	VSdd%s    += dN[3];
	VRI1dd%s  += dN[14];
	VRI2dd%s  += dN[18];
	VRI3dd%s  += dN[22];
	VRA1dd%s  += dN[10] + dN[26];
	VRA2dd%s  += dN[30];
	VRA3dd%s  += dN[34];

}

VSd%s    += dN[43] + dN[49] - dN[35] - dN[36] - dN[37];
VRI1d%s  += - dN[38] - dN[39];
VRI2d%s  += dN[39] - dN[40] - dN[41];
VRI3d%s  += dN[41] - dN[42] - dN[43];
VRA1d%s  += - dN[44] - dN[45];
VRA2d%s  += dN[45] - dN[46] - dN[47];
VRA3d%s  += dN[47] - dN[48] - dN[49];

VSdd%s    +=  dN[58] + dN[64] - dN[50] - dN[51] - dN[52];
VRI1dd%s  += - dN[53] - dN[54];
VRI2dd%s  +=  dN[54] - dN[55] - dN[56] ;
VRI3dd%s  +=  dN[56] - dN[57] - dN[58];
VRA1dd%s  += - dN[59] - dN[60];
VRA2dd%s  +=  dN[60] - dN[61] - dN[62];
VRA3dd%s  +=  dN[62] - dN[63] - dN[64];

/* *previous* vacccination campain */

VSd_alt%s    += dN[43+30] + dN[49+30] - dN[35+30] - dN[36+30] - dN[37+30];
VRI1d_alt%s  += - dN[38+30] - dN[39+30];
VRI2d_alt%s  += dN[39+30] - dN[40+30] - dN[41+30];
VRI3d_alt%s  += dN[41+30] - dN[42+30] - dN[43+30];
VRA1d_alt%s  += - dN[44+30] - dN[45+30];
VRA2d_alt%s  += dN[45+30] - dN[46+30] - dN[47+30];
VRA3d_alt%s  += dN[47+30] - dN[48+30] - dN[49+30];

VSdd_alt%s    +=  dN[58+30] + dN[64+30] - dN[50+30] - dN[51+30] - dN[52+30];
VRI1dd_alt%s  += - dN[53+30] - dN[54+30];
VRI2dd_alt%s  +=  dN[54+30] - dN[55+30] - dN[56+30] ;
VRI3dd_alt%s  +=  dN[56+30] - dN[57+30] - dN[58+30];
VRA1dd_alt%s  += - dN[59+30] - dN[60+30];
VRA2dd_alt%s  +=  dN[60+30] - dN[61+30] - dN[62+30];
VRA3dd_alt%s  +=  dN[62+30] - dN[63+30] - dN[64+30];


C%s   +=  dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30];
W%s   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
B%s += (((dB) < -B%s) ? (-B%s + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S%s = nearbyint(H%s - I%s - A%s - RI1%s - RI2%s - RI3%s - RA1%s - RA2%s - RA3%s -
	VSd%s - VRI1d%s - VRI2d%s - VRI3d%s - VRA1d%s - VRA2d%s -VRA3d%s -
	VSdd%s- VRI1dd%s -VRI2dd%s -VRI3dd%s -VRA1dd%s-VRA2dd%s-VRA3dd%s -
	VSd_alt%s - VRI1d_alt%s - VRI2d_alt%s - VRI3d_alt%s - VRA1d_alt%s - VRA2d_alt%s - VRA3d_alt%s -
	VSdd_alt%s - VRI1dd_alt%s - VRI2dd_alt%s - VRI3dd_alt%s - VRA1dd_alt%s- VRA2dd_alt%s - VRA3dd_alt%s);



IncidenceAll +=  dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30] + dN[1] + dN[36] + dN[51] + dN[36+30] + dN[51+30];
if (!previous_vacc_campaign)
{
	DosesAll  += dN[2] + dN[13] + dN[17] + dN[21] + dN[9] + dN[25] + dN[29] + dN[33] + 2*(dN[3] + dN[14] + dN[18] + dN[22] + dN[10] + dN[26] + dN[30] + dN[34]);
}
CasesAll  +=  dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30];