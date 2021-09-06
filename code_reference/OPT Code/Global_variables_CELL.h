/************************************************************************************************************
											Global Variables for Human Cable Simulations
By: Jonathan D. Moreno
Colleen Clancy Laboratory

August 2007
************************************************************************************************************/


#include <Carbon/Carbon.h>

#ifndef GLOBAL_H
#define GLOBAL_H


/*********************************************************************************************************
Universal Constants
*********************************************************************************************************/

//Ion Valences and Universal Constants	
	const double R = 8314.472;				// J/mol*K
	const double T = 310;					// K
	const double F = 96485.3415;			 
	const double Cm = 2.0;					// uF/ cm^2
	const double CAP = 0.185;				//Cellular capacitance
	const double rho = 162;					// ohm*cm
	const double z_Na = 1;						
	const double z_Ca = 2;
	const double z_K = 1;	
	const double dx= 0.01;					//.01 //Can go from 0.01 - 0.02 cm

//Cell Geometry
	const double pi = 3.141592;						
	const double S_cg = 0.2;				//Surface to volume ratio (um^-1)
	const double mutant = 0.0;

//Intracellular volumes
	const double V_cyto=0.016404;			//16404 uL
	const double V_sr=0.001094;
	const double V_ss=0.00005468;

//Extraceullular concentrations
	const double Na_out = 140;
	const double Ca_out = 2;
	const double K_out = 5.4;
	
//Beating parameters
	const double I_duration = 1.0; //0.5
	const double stimulus = -80;
	const double dt = 0.01; //0.01
	const double waitTime = 1000;


//Drug Conditions and Mutant
#define FLECAINIDE
//#define LIDOCAINE

//#define KPQ
//#define D1790G
//#define Y1795C

	const double drug=2.0*(1e-6);	//****** 10uM for Flecainide, 300uM for Lidocaine



/*********************************************************************************************************
Structures
*********************************************************************************************************/		
typedef struct cell_param{
	double V, dV, V_new;
	double Na_in, K_in, Ca_in, Ca_sr, Ca_ss, Ca_in_buffer, Ca_ss_buffer, Ca_sr_buffer;
	double I_Na, mI_Na, I_Ca_L, I_Kr, I_Ks, I_K1, I_Kp, I_to, I_Na_Ca, I_Na_K, I_p_Ca, I_Ca_b, I_Na_b, I_stim;
	double I_Na_ion_total, I_Ca_ion_total, I_K_ion_total, I_total, I_axial;
	double I_tr, I_leak, I_up, I_rel;
	int Cell_type;
	double m, h, j, d, f, f2, f_Ca, r, s, xr1, xr2, xs, OO, R_bar;
	
	double IC3, IC2, IF, IM1, IM2, C3, C2, O, OS, C1;
	double DIC3, DIC2, DIF, DIM1, DIM2, DC3, DC2, DC1, DO, DOS;
	double D_IC3, D_IC2, D_IF, D_IM1, D_IM2, D_C3, D_C2, D_C1, D_O, D_OS;

	double mIC3, mIC2, mIF, mIM1, mIM2, mC3, mC2, mO, mOS, mC1, mBC3, mBC2, mBC1, mBO;
	double DmIC3, DmIC2, DmIF, DmIM1, DmIM2, DmC3, DmC2, DmC1, DmO, DmOS, DmBC3, DmBC2, DmBC1, DmBO;
	double D_mIC3, D_mIC2, D_mIF, D_mIM1, D_mIM2, D_mC3, D_mC2, D_mC1, D_mO, D_mOS, D_mBC3, D_mBC2, D_mBC1, D_mBO;
	double mBIF, DmBIF, D_mBIF;
	
	double peak_slope, t_min, V_min, t_thr, V_thr, t_max, V_max, t_EAD, V_EAD, I_Na_peak, t_I_Na_peak;
	double t_EAD2, V_EAD2, t_90, V_90, dV_old, flag2;
	double APD_90, DI, L_EAD;
	double t_90_old, flag_EAD, flag_EAD2;
	double CV, CV_EAD;
	
	double I_total_old, dI_total, dI_total_old, t_I_total, I_total_pt, V_I_total, flagI_total;
	
    double Drug;
    
    } Cell_param;


/*********************************************************************************************************
Current Function Prototypes
*********************************************************************************************************/			
	void Calculate_I_Na(Cell_param *Cell_ptr, double t, double tTemp, double *Input_param);	//Fast sodium current
	void Calculate_mI_Na(Cell_param *Cell_ptr, double t, double tTemp);
	void Calculate_I_K1(Cell_param *Cell_ptr, double t, double tTemp );						//Time independent potassium current
	void Calculate_I_Kp(Cell_param *Cell_ptr, double t, double tTemp );						//Potassium plateau current, time independent, K_out insensitive	
	void Calculate_I_Na_Ca(Cell_param *Cell_ptr, double t, double tTemp );					//Sodium calcium exchanger
	void Calculate_I_Na_K(Cell_param *Cell_ptr, double t, double tTemp );					//Sodium potassium pump
	void Calculate_I_Ca_L(Cell_param *Cell_ptr, double t, double tTemp);					//L type calcium channel Ca contribution
	void Calculate_I_Kr(Cell_param *Cell_ptr, double t, double tTemp );						//Rapid rectifier
	void Calculate_I_Ks(Cell_param *Cell_ptr, double t, double tTemp);						//Slow rectifier
	void Calculate_I_to(Cell_param *Cell_ptr, double t, double tTemp );						//Transient outward current
	void Calculate_I_p_Ca(Cell_param *Cell_ptr, double t, double tTemp );					//Sarcolemmal calcium pump
	void Calculate_I_Ca_b(Cell_param *Cell_ptr, double t, double tTemp );					//Calcium background current
	void Calculate_I_Na_b(Cell_param *Cell_ptr, double t, double tTemp );					//Sodium background current
	void Calculate_I_total(Cell_param *Cell_ptr, double t, double tTemp);

	void Calculate_I_tr(Cell_param *Cell_ptr, double t, double tTemp);
	void Calculate_I_leak(Cell_param *Cell_ptr, double t, double tTemp);					//Leak from SR
	void Calculate_I_up(Cell_param *Cell_ptr, double t, double tTemp);						//Uptake to SR
	void Calculate_I_rel(Cell_param *Cell_ptr, double t, double tTemp);						//Release from SR

//Dynamic Concentrations			
	void Calculate_Na_in(Cell_param *Cell_ptr, double t, double tTemp );					//Dynamic [Na] in Myoplasm
	void Calculate_K_in(Cell_param *Cell_ptr, double t, double tTemp );						//Dynamic [K] in Myoplasm
	void Calculate_Ca_in(Cell_param *Cell_ptr, double t, double tTemp );					//Dynamic [Ca] in Myoplasm
	void Calculate_Ca_sr(Cell_param *Cell_ptr, double t, double tTemp );					//Dynamic [Ca] in JSR
	void Calculate_Ca_ss(Cell_param *Cell_ptr, double t, double tTemp );					//Dynamic [Ca] in NSR
	
	void Calculate_Reset (Cell_param *Cell_ptr, double t, double tTemp);
	void Calculate_Update (Cell_param *Cell_ptr, double t, double tTemp);
	void Calculate_Points (Cell_param *Cell_ptr, double t, double tTemp );

#endif