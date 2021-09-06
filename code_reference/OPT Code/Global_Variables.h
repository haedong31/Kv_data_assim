/************************************************************************************************************
 Global Variables for Human Cable Simulations
 By: Jonathan D. Moreno
 Colleen Clancy Laboratory
 
 August 2007
 ************************************************************************************************************/
#include <iostream>
#include <fstream>
#include <math.h>
#include <unistd.h>


#ifndef GLOBAL_H
#define GLOBAL_H


/*********************************************************************************************************
 Universal Constants
 *********************************************************************************************************/
const double pi = 3.141592;	
const double R = 8314.472;				// J/mol*K
const double F = 96485.3415;			 
const double Na_out = 138;
const double Na_in = 10;	

//#define PHYSIOLOGIC
#ifdef PHYSIOLOGIC
	const double V_duration = 400;
	const double BCL_1 = 1000;
	const double T = 310;
#else




	const double V_duration = 25;				//25
	const double BCL_1 = 100;	//For concentration dependent UDB, BCL_1 = 200ms; change this for frequency dependent UDB, 1Hz = 1000, 2Hz = 500, 5Hz = 200, 10Hz = 100
	const double T = 295;
#endif

//Beating parameters
const double V_o=-100;	
const double V_step= -10;
const double waitTime = 10000;  //***************You MUST hold for 30,000ms to get block of late current!!!!

const double total_cycles = 300;
int cycle, cycle2;
double interval;
const double BCL_3 = 40;
const int which = 2; // (1) KPQ, (2) WT, (3) DG, (4) YC

//Inactivation Parameters
double V_initial;
const double t_hold = 5000;

//Activation parameters
double V_temp;

//Testing Parameters
double t_peak, I_Na_peak, norm;
double I_Na_peak_0, I_Na_peak_end, late_avg, late_avg_0, late_avg_end;
double I_Na_peak_tau, t_tau, flag_tau, tau;

#define FLECAINIDE
//#define LIDOCAINE

/*********************************************************************************************************
 Structures
 *********************************************************************************************************/		
typedef struct cell_param{
	double V, I_Na, sum;
	
	double IC3, IC2, IF, IM1, IM2, C3, C2, O, OS, C1;
	double DIC3, DIC2, DIF, DIM1, DIM2, DC3, DC2, DO, DOS, DC1;
	double D_IC3, D_IC2, D_IF, D_IM1, D_IM2, D_C3, D_C2, D_O, D_OS, D_C1;
	
	double mIC3, mIC2, mIF, mIM1, mIM2, mC3, mC2, mO, mOS, mC1, mBC3, mBC2, mBC1, mBO;
	double DmIC3, DmIC2, DmIF, DmIM1, DmIM2, DmC3, DmC2, DmO, DmOS, DmC1, DmBC3, DmBC2, DmBC1, DmBO;
	double D_mIC3, D_mIC2, D_mIF, D_mIM1, D_mIM2, D_mC3, D_mC2, D_mO, D_mOS, D_mC1, D_mBC3, D_mBC2, D_mBC1, D_mBO;
	double mBIF, DmBIF, D_mBIF;
	
	double Drug;
	
} Cell_param;


/*********************************************************************************************************
 Current Function Prototypes
 *********************************************************************************************************/			
void Calculate_I_Na(Cell_param *Cell_ptr, double t, double dt, double *Input_param);	//Fast sodium current
void Calculate_peak (Cell_param *Cell_ptr, double t, double dt);

#endif