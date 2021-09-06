/***************************************************************************************
 Standard Human Cellular Functions
 
***************************************************************************************/

#include "Global_variables_CELL.h"
#include <math.h>

void Calculate_Reset (Cell_param *Cell_ptr, double t, double tTemp){
	if (t>=waitTime && tTemp <=(0.001)) {				//0.001
		Cell_ptr->I_Na_peak = 0;
		Cell_ptr->t_90_old = Cell_ptr->t_90;
		Cell_ptr->flag2 = 0; 
		Cell_ptr->flag_EAD=0;
		Cell_ptr->flag_EAD2=0;
		Cell_ptr->peak_slope = 0;
		Cell_ptr->V_min = -88.654973;
		Cell_ptr->t_min = 0;
		Cell_ptr->V_thr = -88.654973;
		Cell_ptr->t_thr = 0;
		Cell_ptr->V_max = -88.654973;
		Cell_ptr->t_max = t;
		Cell_ptr->t_EAD= t;
		Cell_ptr->t_EAD2=t;
		Cell_ptr->t_90=t;
		Cell_ptr->V_EAD= -88.654973;
		Cell_ptr->V_EAD2=-88.654973;
		Cell_ptr->V_90 = -88.654973;
		Cell_ptr->dV_old = 0;
		
		Cell_ptr->flagI_total = 0;
	}
	
}

void Calculate_Update (Cell_param *Cell_ptr, double t, double tTemp){
	Cell_ptr->dV_old = Cell_ptr->dV;	
	Cell_ptr->dV = -1*(Cell_ptr->I_total)*dt;
	Cell_ptr->V = Cell_ptr->V + Cell_ptr->dV;
}

void Calculate_I_total(Cell_param *Cell_ptr, double t, double tTemp){ 
	Cell_ptr->I_total_old = Cell_ptr->I_total;
	Cell_ptr->dI_total_old = Cell_ptr->dI_total;
	
	if (t>waitTime && tTemp <=I_duration) {Cell_ptr->I_stim=stimulus;}
	else {Cell_ptr->I_stim=0;}		

	Cell_ptr->I_Na_ion_total = (1-mutant)*(Cell_ptr->I_Na) + (mutant)*(Cell_ptr->mI_Na) + Cell_ptr->I_Na_b + 3*Cell_ptr->I_Na_K + 3*Cell_ptr->I_Na_Ca;
	Cell_ptr->I_K_ion_total = Cell_ptr->I_Kr + Cell_ptr->I_Ks + Cell_ptr->I_K1 + Cell_ptr->I_Kp - 2*Cell_ptr->I_Na_K + Cell_ptr->I_to;
	Cell_ptr->I_Ca_ion_total = Cell_ptr->I_Ca_L + Cell_ptr->I_p_Ca + Cell_ptr->I_Ca_b - 2*Cell_ptr->I_Na_Ca;
	Cell_ptr->I_total = Cell_ptr->I_stim + Cell_ptr->I_Na_ion_total + Cell_ptr->I_K_ion_total + Cell_ptr->I_Ca_ion_total;

	Cell_ptr->dI_total = (Cell_ptr->I_total - Cell_ptr->I_total_old)/dt;	
}

void Calculate_I_Na(Cell_param *Cell_ptr, double t, double tTemp, double *Input_param){
	return( WT_SCN5A_function (Cell_ptr, t, tTemp, Input_param) );
}

void Calculate_mI_Na(Cell_param *Cell_ptr, double t, double tTemp){
	
#ifdef KPQ
	return( KPQ_SCN5A_function (Cell_ptr, t, tTemp) );
#endif
	
#ifdef D1790G
	return( DG_SCN5A_function (Cell_ptr, t, tTemp) );
#endif

#ifdef Y1795C
	return( YC_SCN5A_function (Cell_ptr, t, tTemp) );
#endif	
}

void Calculate_I_Na_Ca (Cell_param *Cell_ptr, double t, double tTemp){
    double k_Na_Ca = 1000;			// pA/pF
	double Km_Na = 87.5;			// Na_in half-saturation concentration of NaCa exhanger (mM)
	double Km_Ca = 1.38;			// Ca_in half-saturation concentration of NaCa exhanger (mM)
	double k_Na_Ca_sat = 0.1;		// Saturation factor 
	double gamma = 0.35;			// Position of energy barrier controlling voltage dependance of inaca
	double alpha = 2.5;				// Factor enhancing outward nature of I_Na_Ca
	
	Cell_ptr->I_Na_Ca =  k_Na_Ca*(1./(Km_Na*Km_Na*Km_Na + Na_out* Na_out* Na_out))*(1./(Km_Ca+ Ca_out))*
	(1./(1+k_Na_Ca_sat*exp((gamma-1)* Cell_ptr->V *F/(R*T))))*
	(exp(gamma* Cell_ptr->V *F/(R*T))* Cell_ptr->Na_in* Cell_ptr->Na_in* Cell_ptr->Na_in* Ca_out -
	 exp((gamma-1)* Cell_ptr->V *F/(R*T))* Na_out* Na_out* Na_out* Cell_ptr->Ca_in *alpha);      
	
}

void Calculate_I_Ca_L(Cell_param *Cell_ptr, double t, double tTemp ){
	double d_inf, a_d, b_d, c_d, tau_d;				
	double f_inf, a_f, b_f, c_f, tau_f;
	double f2_inf, a_f2, b_f2, c_f2, tau_f2;
	double f_Ca_inf, tau_f_Ca;
	const double G_Ca_L = 3.980E-5;		// cm/(ms*uF)
	
	d_inf = 1./(1.+exp((-8-Cell_ptr->V)/7.5));
	a_d=1.4/(1.+exp((-35-Cell_ptr->V)/13))+0.25;
	b_d=1.4/(1.+exp((Cell_ptr->V+5)/5));
	c_d=1./(1.+exp((50-Cell_ptr->V)/20));
	
	tau_d = a_d*b_d + c_d; 
	Cell_ptr->d = d_inf - (d_inf - Cell_ptr->d)*exp(-dt/tau_d);
	
	f_inf = 1./(1.+exp((Cell_ptr->V+20)/7));
	a_f=1102.5*exp(-(Cell_ptr->V+27)*(Cell_ptr->V+27)/225);
	b_f=200./(1+exp((13-Cell_ptr->V)/10.));
	c_f=(180./(1+exp((Cell_ptr->V+30)/10)))+20;
	
	tau_f = a_f + b_f + c_f;
	Cell_ptr->f = f_inf - (f_inf - Cell_ptr->f)*exp(-dt/tau_f);
	
	f2_inf = 0.67/(1.+exp((Cell_ptr->V+35)/7))+0.33;
	a_f2= 600*exp(-(Cell_ptr->V+25)*(Cell_ptr->V+25)/170);
	b_f2= 31/(1.+exp((25-Cell_ptr->V)/10));
	c_f2= 16/(1.+exp((Cell_ptr->V+30)/10));
	
	tau_f2 = a_f2 + b_f2 + c_f2;
	Cell_ptr->f2 = f2_inf - (f2_inf - Cell_ptr->f2)*exp(-dt/tau_f2);
	
	f_Ca_inf = 0.6/(1+(Cell_ptr->Ca_ss/0.05)*(Cell_ptr->Ca_ss/0.05))+0.4;
	tau_f_Ca = 80./(1+(Cell_ptr->Ca_ss/0.05)*(Cell_ptr->Ca_ss/0.05))+2.;
	
	Cell_ptr->f_Ca = f_Ca_inf - (f_Ca_inf - Cell_ptr->f_Ca)*exp(-dt/tau_f_Ca);
	
	Cell_ptr->I_Ca_L = G_Ca_L*Cell_ptr->d* Cell_ptr->f* Cell_ptr->f2* Cell_ptr->f_Ca*4*(Cell_ptr->V-15)*(F*F/(R*T))*
	(0.25*exp(2*(Cell_ptr->V-15)*F/(R*T))* Cell_ptr->Ca_ss- Ca_out)/(exp(2*(Cell_ptr->V-15)*F/(R*T))-1.);
	
}	

void Calculate_I_Kr(Cell_param *Cell_ptr, double t, double tTemp) {
	double a_xr1, b_xr1, xr1_ss, tau_xr1;				
	double a_xr2, b_xr2, xr2_ss, tau_xr2;	
	const double G_Kr = 0.153;		// nS/pF
	double E_Kr = ((R*T)/(z_K*F))*log(K_out/Cell_ptr->K_in);
	
	a_xr1=450./(1.+exp((-45.-Cell_ptr->V)/10.));
    b_xr1=6./(1.+exp((Cell_ptr->V-(-30.))/11.5));
	xr1_ss = 1./(1.+exp((-26.-Cell_ptr->V)/7.));
	tau_xr1 = a_xr1*b_xr1;
	Cell_ptr->xr1 = xr1_ss - (xr1_ss - Cell_ptr->xr1)*exp(-dt/tau_xr1);
	
	a_xr2=3./(1.+exp((-60.-Cell_ptr->V)/20.));
    b_xr2=1.12/(1.+exp((Cell_ptr->V-60.)/20.));
	xr2_ss = 1./(1.+exp((Cell_ptr->V-(-88.))/24.));
	tau_xr2 = a_xr2*b_xr2;
	Cell_ptr->xr2 = xr2_ss - (xr2_ss - Cell_ptr->xr2)*exp(-dt/tau_xr2);
	
	Cell_ptr->I_Kr = G_Kr*sqrt(K_out/5.4)*Cell_ptr->xr1*Cell_ptr->xr2*(Cell_ptr->V-E_Kr);
}

void Calculate_I_Ks(Cell_param *Cell_ptr, double t, double tTemp ){
	double xs_ss, ax_s, bx_s, tau_xs, G_Ks;
	const double PR_NaK = 0.03;
	
	if (Cell_ptr->Cell_type == 1)			// Cell type = (1) for endo; (2) for M,  (3) epi
	{G_Ks = 0.392;}			// nS/pF
	else if (Cell_ptr->Cell_type == 2)
	{G_Ks = 0.098;}			// nS/pF
	else if (Cell_ptr->Cell_type == 3)
	{G_Ks = 0.392;}			// nS/pF		//0.392
	
	double E_Ks = ((R*T)/(z_K*F))*(log((K_out+PR_NaK*Na_out)/(Cell_ptr->K_in + PR_NaK*Cell_ptr->Na_in)));
    
	xs_ss = 1./(1.+exp((-5.-Cell_ptr->V)/14.));
	ax_s =(1400./(sqrt(1.+exp((5.-Cell_ptr->V)/6))));
    bx_s =(1./(1.+exp((Cell_ptr->V-35.)/15.)));
    tau_xs =(ax_s*bx_s) + 80;
	
	Cell_ptr->xs = xs_ss - (xs_ss - Cell_ptr->xs)*exp(-dt/tau_xs);
	
	Cell_ptr->I_Ks = G_Ks*Cell_ptr->xs*Cell_ptr->xs*(Cell_ptr->V-E_Ks);
}

void Calculate_I_K1(Cell_param *Cell_ptr, double t, double tTemp ){
	double a_K1, b_K1, K1_s;
	const double G_K1 = 5.405;			// nS/pF
	double E_K1 = ((R*T)/(z_K*F))*log(K_out/Cell_ptr->K_in);
	
	a_K1 =0.1/(1.+exp(0.06*(Cell_ptr->V-E_K1-200)));
	b_K1 = (3.*exp(0.0002*(Cell_ptr->V-E_K1+100)) + exp(0.1*(Cell_ptr->V-E_K1-10)))/(1.+exp(-0.5*(Cell_ptr->V-E_K1)));
	
	K1_s = a_K1/(a_K1 + b_K1); 
	
	Cell_ptr->I_K1=G_K1*K1_s*(Cell_ptr->V-E_K1);
}

void Calculate_I_to(Cell_param *Cell_ptr, double t, double tTemp ){
	double s_inf, tau_s, r_inf, tau_r, G_to;	
	if (Cell_ptr->Cell_type == 1)			// Cell type = (1) for endocardial
	{G_to = 0.073;			// nS/pF
		s_inf = 1./(1.+exp((Cell_ptr->V+28)/5.));
		tau_s = 1000.*exp(-(Cell_ptr->V+67)*(Cell_ptr->V+67)/1000.)+8.;
	}		
	
	else if (Cell_ptr->Cell_type == 2)	// Cell type = (2) for M Cell
	{G_to = 0.294;			// nS/pF
		s_inf = 1./(1.+exp((Cell_ptr->V+20)/5.));
		tau_s = 85.*exp(-(Cell_ptr->V+45.)*(Cell_ptr->V+45.)/320.)+5./(1.+exp((Cell_ptr->V-20.)/5.))+3.;
	}	
	
	else if (Cell_ptr->Cell_type == 3)	// Cell type = (3) for epicardial
	{G_to = 0.294;			// nS/pF
		s_inf=1./(1.+exp((Cell_ptr->V+20)/5.));
		tau_s = 85.*exp(-(Cell_ptr->V+45.)*(Cell_ptr->V+45.)/320.)+5./(1.+exp((Cell_ptr->V-20.)/5.))+3.;
	}			
	
	double E_to = ((R*T)/(z_K*F))*log(K_out/Cell_ptr->K_in); 
	
	r_inf = 1./(1.+exp((20-Cell_ptr->V)/6.));	
	tau_r = 9.5*exp(-(Cell_ptr->V+40.)*(Cell_ptr->V+40.)/1800.)+0.8;
	
	Cell_ptr->r = r_inf-(r_inf-Cell_ptr->r)*exp(-dt/tau_r);
	Cell_ptr->s = s_inf-(s_inf-Cell_ptr->s)*exp(-dt/tau_s);
	
	Cell_ptr->I_to = G_to*Cell_ptr->r*Cell_ptr->s*(Cell_ptr->V-E_to); 
} 


void Calculate_I_p_Ca(Cell_param *Cell_ptr, double t, double tTemp){
	const double G_p_Ca = 0.1238;		// nS/pF
	const double Km_p_Ca = 0.0005;		// Half saturation constant (mM)
	
	Cell_ptr->I_p_Ca=G_p_Ca*(Cell_ptr->Ca_in/(Km_p_Ca + Cell_ptr->Ca_in));
}

void Calculate_I_Kp(Cell_param *Cell_ptr, double t, double tTemp ){
	const double G_Kp = 0.0146;			// nS/pF
	
	double E_Kp=(R*T/F)*log(K_out/Cell_ptr->K_in);
	
	Cell_ptr->I_Kp = G_Kp*(1./(1.+exp((25-Cell_ptr->V)/5.98)))*(Cell_ptr->V-E_Kp);	
}

void Calculate_I_Ca_b(Cell_param *Cell_ptr, double t, double tTemp ){
	const double G_Ca_b =  0.000592;
	
	double E_Ca_b = ((R*T)/(z_Ca*F))*log(Ca_out/Cell_ptr->Ca_in);
	Cell_ptr->I_Ca_b = G_Ca_b*(Cell_ptr->V-E_Ca_b);
}

void Calculate_I_Na_b(Cell_param *Cell_ptr, double t, double tTemp ){
	const double G_Na_b = 0.000290;		// nS/pF
	
	double E_Na_b=((R*T)/(z_Na*F))*log(Na_out/Cell_ptr->Na_in);			
	Cell_ptr->I_Na_b = G_Na_b*(Cell_ptr->V-E_Na_b);
}

void Calculate_I_Na_K(Cell_param *Cell_ptr, double t, double tTemp){
	const double I_Na_K_bar = 2.724;	// Maximal I_Na_K (pA/pF)
	const double Km_Na_in = 40;			// Na_in half saturation constant
	const double Km_K_out = 1;			// K_out half saturation constant
	
	Cell_ptr->I_Na_K= I_Na_K_bar*(K_out*Cell_ptr->Na_in/((K_out+Km_K_out)*(Cell_ptr->Na_in+Km_Na_in)))*
	(1./(1.+0.1245*exp(-0.1*Cell_ptr->V*F/(R*T))+0.0365*exp(-Cell_ptr->V*F/(R*T))));
}


void Calculate_I_rel (Cell_param *Cell_ptr, double t, double tTemp){
	double k_Ca_sr, k1_rel, k2_rel, dR_bar;
	const double G_rel = 0.102;				// Max. rate constant of Ca release from JSR due to overload (mM/ms)
	const double k1_rel_ = 0.15;			// R to O and RI to I I_rel transition rate (mM^2/ms)
	const double k2_rel_ = 0.045;			// O to I and R to RI I_rel transition rate (mM^2/ms)
	const double k3_rel = 0.060;			// O to R and I to RI I_rel transition rate (ms^-1)
	const double k4_rel = 0.005;			// I to O and RI to I I_rel transition rate (ms^-1)
	const double EC_sr = 1.5;				// Ca_sr half-saturation constant of k_Ca_sr
	const double max_sr = 2.5;				// Max value of k_Ca_sr (dimensionless)
	const double min_sr = 1;				// Min value of k_Ca_sr (dimensionless)
	
	k_Ca_sr = max_sr-((max_sr-min_sr)/(1+(EC_sr/Cell_ptr->Ca_sr)*(EC_sr/Cell_ptr->Ca_sr)));
	
	k1_rel = k1_rel_ / k_Ca_sr;
	k2_rel = k2_rel_ * k_Ca_sr;
	
	dR_bar = (-k2_rel*Cell_ptr->Ca_ss* Cell_ptr->R_bar + k4_rel*(1- Cell_ptr->R_bar))*dt;
	Cell_ptr->R_bar = Cell_ptr->R_bar + dR_bar;
	
	Cell_ptr->OO = (k1_rel*Cell_ptr->Ca_ss* Cell_ptr->Ca_ss* Cell_ptr->R_bar)/(k3_rel + k1_rel* Cell_ptr->Ca_ss* Cell_ptr->Ca_ss);
	
	Cell_ptr->I_rel = G_rel*Cell_ptr->OO*(Cell_ptr->Ca_sr - Cell_ptr->Ca_ss);
}

void Calculate_I_up(Cell_param *Cell_ptr, double t, double tTemp ){
	const double G_up = 0.006375;		// mM/ms
	const double Km_up = 0.00025;			//mM
	
	Cell_ptr->I_up = G_up*((Cell_ptr->Ca_in* Cell_ptr->Ca_in)/((Cell_ptr->Ca_in* Cell_ptr->Ca_in)+(Km_up*Km_up)));
}

void Calculate_I_leak(Cell_param *Cell_ptr, double t, double tTemp ){
	const double G_leak = 0.00036;		// mM/ms
	
	Cell_ptr->I_leak = G_leak*(Cell_ptr->Ca_sr- Cell_ptr->Ca_in); 
}

void Calculate_I_tr( Cell_param *Cell_ptr, double t, double tTemp){
	const double G_tr = 0.0038;		// mM/ms
	
	Cell_ptr->I_tr = G_tr*(Cell_ptr->Ca_ss- Cell_ptr->Ca_in);
}

/************************************** Dynamic Ion concentrations ************************************/

//New K concentration in myoplasm
void Calculate_K_in(Cell_param *Cell_ptr, double t, double tTemp ){
	double dK_in;
	
	dK_in =  -dt*((Cell_ptr->I_K_ion_total + Cell_ptr->I_stim)*CAP)/(V_cyto*z_K*F);
	Cell_ptr->K_in = Cell_ptr->K_in + dK_in;
}

//New Na concentration in myoplasm
void Calculate_Na_in(Cell_param *Cell_ptr, double t, double tTemp ){
	double dNa_in;
	dNa_in = -1*dt*(Cell_ptr->I_Na_ion_total*CAP)/(V_cyto*z_Na*F);
	Cell_ptr->Na_in = Cell_ptr->Na_in + dNa_in;
}

//New Ca concentration in cytoplasm
void Calculate_Ca_in(Cell_param *Cell_ptr, double t, double tTemp ){
	double b_Ca, c_Ca, dCa_in;
	const double buffer = 0.2;				//Total cytoplasmic buffer concentration (mM)
	const double K_buffer = 0.001;			// Ca_in half-saturation constant for cytoplasmic buffer (mM)
	
	Cell_ptr->Ca_in_buffer = (Cell_ptr->Ca_in* buffer)/(Cell_ptr->Ca_in + K_buffer);
	dCa_in = dt*((-(Cell_ptr->I_Ca_b+ Cell_ptr->I_p_Ca - 2*Cell_ptr->I_Na_Ca)*CAP/(V_cyto*z_Ca*F))+(V_sr/V_cyto)*(Cell_ptr->I_leak- Cell_ptr->I_up)+ Cell_ptr->I_tr);
	
	b_Ca = buffer - Cell_ptr->Ca_in_buffer - dCa_in - Cell_ptr->Ca_in + K_buffer;
	c_Ca = K_buffer*(Cell_ptr->Ca_in_buffer + dCa_in + Cell_ptr->Ca_in);
	
	Cell_ptr->Ca_in = (sqrt(b_Ca*b_Ca + 4*c_Ca) - b_Ca)/2;
}

//New Ca concentration in SR
void Calculate_Ca_sr(Cell_param *Cell_ptr, double t, double tTemp){
	double b_jsr, c_jsr, dCa_sr;
	const double buffer_sr = 10;				//Total SR buffer concentration (mM)
	const double K_buffer_sr = 0.3;				// Ca_sr half-saturation constant for SR buffer (mM)
	
	Cell_ptr->Ca_sr_buffer = (Cell_ptr->Ca_sr* buffer_sr)/(Cell_ptr->Ca_sr + K_buffer_sr);
	dCa_sr = dt*(Cell_ptr->I_up - Cell_ptr->I_rel - Cell_ptr->I_leak);
	
	b_jsr = buffer_sr - Cell_ptr->Ca_sr_buffer - dCa_sr - Cell_ptr->Ca_sr + K_buffer_sr;
	c_jsr = K_buffer_sr*(Cell_ptr->Ca_sr_buffer + dCa_sr + Cell_ptr->Ca_sr);
	
	Cell_ptr->Ca_sr = (sqrt(b_jsr*b_jsr + 4*c_jsr) - b_jsr)/2;
}

//New Ca concentration in SS
void Calculate_Ca_ss(Cell_param *Cell_ptr, double t, double tTemp){
	double b_Ca_ss, c_Ca_ss, dCa_ss;
	const double buffer_ss = 0.4;				//Total SS buffer concentration (mM)
	const double K_buffer_ss = 0.00025;			// Ca_ss half-saturation constant for SR buffer (mM)
	
	Cell_ptr->Ca_ss_buffer = (Cell_ptr->Ca_ss*buffer_ss)/(Cell_ptr->Ca_ss + K_buffer_ss);
	dCa_ss = dt*((-Cell_ptr->I_Ca_L*CAP)/(V_ss*z_Ca*F)+(V_sr/V_ss)* Cell_ptr->I_rel-(V_cyto/V_ss)* Cell_ptr->I_tr);
	
	b_Ca_ss = buffer_ss - Cell_ptr->Ca_ss_buffer - dCa_ss - Cell_ptr->Ca_ss + K_buffer_ss;
    c_Ca_ss = K_buffer_ss*(Cell_ptr->Ca_ss_buffer + dCa_ss + Cell_ptr->Ca_ss);
	
	Cell_ptr->Ca_ss = (sqrt(b_Ca_ss*b_Ca_ss + 4*c_Ca_ss) - b_Ca_ss)/2;
}

void Calculate_Points(Cell_param *Cell_ptr, double t, double tTemp){	
	if (fabs((1-mutant)*Cell_ptr->I_Na + mutant*Cell_ptr->mI_Na) > fabs(Cell_ptr->I_Na_peak) ){
		Cell_ptr->I_Na_peak = (1-mutant)*Cell_ptr->I_Na + mutant*Cell_ptr->mI_Na;
		Cell_ptr->t_I_Na_peak = t;}
	
	
	//Minimum voltage catcher
	if (t>=waitTime && tTemp<0.01){	
		Cell_ptr->t_min = t;
		Cell_ptr->V_min = Cell_ptr->V;
	}
	
	//Threshold voltage point catcher needed for t_min
	if (t>=waitTime && (Cell_ptr->dV/dt) > Cell_ptr->peak_slope) {
		Cell_ptr->peak_slope = Cell_ptr->dV/dt;
		Cell_ptr->V_thr = Cell_ptr->V;
		Cell_ptr->t_thr = t;
	}
	
	//Maximum voltage point catcher	
	if (t>=waitTime && Cell_ptr->dV_old >0 && Cell_ptr->dV < 0 && Cell_ptr->V >10) {
		Cell_ptr->V_max = Cell_ptr->V;
		Cell_ptr->t_max = t;	
	}

	//I_total minimum point
	if (t>=waitTime && Cell_ptr->dI_total_old < 0 && Cell_ptr->dI_total> 0 && tTemp >200 && Cell_ptr->flagI_total ==0){	
		Cell_ptr->flagI_total = 1;
		Cell_ptr->t_I_total = t;
		Cell_ptr->I_total_pt = Cell_ptr->I_total;
		Cell_ptr->V_I_total = Cell_ptr->V;
	}
	
	//Local minima	
	if (t>=waitTime && Cell_ptr->dV_old < 0 && Cell_ptr->dV > 0 && Cell_ptr->V > -40 && tTemp > 200 && Cell_ptr->flag_EAD==0) {
		Cell_ptr->flag_EAD=1;
		Cell_ptr->V_EAD = Cell_ptr->V;
		Cell_ptr->t_EAD = t;
	}
	
	//Local maxima	
	if (t > Cell_ptr->t_EAD && Cell_ptr->dV_old > 0 && Cell_ptr->dV < 0 && tTemp > 200 && Cell_ptr->flag_EAD2==0) {
		Cell_ptr->flag_EAD2=1;
		Cell_ptr->V_EAD2 = Cell_ptr->V;
		Cell_ptr->t_EAD2 = t;
	}
	
	//V_90 voltage point catcher, cycle number and APD_90 calculation
	Cell_ptr->V_90 = Cell_ptr->V_max - .90*(Cell_ptr->V_max - Cell_ptr->V_min);
	
	//Diastolic Interval	
	Cell_ptr->DI = Cell_ptr->t_thr - Cell_ptr->t_90_old;			
	
	//Height (or length) of the EAD	
	Cell_ptr->L_EAD=Cell_ptr->V_EAD2 - Cell_ptr->V_EAD;
	if (t>=waitTime && fabs(Cell_ptr->V - Cell_ptr->V_90) < (0.015) && Cell_ptr->dV < 0 && tTemp > 30 && t > Cell_ptr->t_max && Cell_ptr->flag2 ==0) {	
		Cell_ptr->flag2 = 1;
		Cell_ptr->t_90 = t;
		Cell_ptr->APD_90 = Cell_ptr->t_90 - Cell_ptr->t_thr;
		//Cell_ptr->peak_slope = 0;
	}
	
}