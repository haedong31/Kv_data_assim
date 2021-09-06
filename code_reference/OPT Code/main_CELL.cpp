/************************************************************************************************************
 Updated Single Cell Formulation SINGLE CELL SINGLE CELL SINGLE CELL SINGLE CELL
 Human Cable Model from ten Tusscher 2006 based off of the LR Model formulation.
 
 Na channel functions as well as additional functions separated into seperate header files
 
 By: Jonathan D. Moreno 
 Colleen Clancy Laboratory
 January 27, 2009
 *************************************************************************************************************/
#include "Global_variables_CELL.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include "mex.h"

#ifdef FLECAINIDE
	//#include "WT_Flec_Implicit.h"
    #include "WT_Flec_Implicit_CELL.h"
	/*
	#ifdef KPQ 
		#include "KPQ_Flecainide.h"
	#endif
	
	#ifdef D1790G
		#include "DG_Flecainide.h"
	#endif

	#ifdef Y1795C 
		#include "YC_Flecainide.h"
	#endif
	 */
#endif


#ifdef LIDOCAINE
	#include "WT_Lidocaine.h"
	/*
	#ifdef KPQ 
		#include "KPQ_Lidocaine.h"
	#endif

	#ifdef D1790G
		#include "DG_Lidocaine.h"
	#endif

	#ifdef Y1795C 
		#include "YC_Flecainide.h"
	#endif
	 */
#endif


#include "Functions_CELL.h"

using namespace std;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	double *Input_param;
	int i;
	Input_param = mxGetPr(prhs[0]);
	
	
	Cell_param Cell;
	int S1_cycle;
	double t = 0, tTemp = 0, interval;
	double Upstroke_Velocity_0, Upstroke_Velocity_400, Upstroke_Velocity_500;
	//const char *stateFileName = "initial_100_beats.j";
	string str = "Cell";
	int counter = 0;											
			
	ofstream result_file((str+".txt").c_str()   );			// This file contains all of the voltages for each cell (the master-esque file)
	ofstream result_file0((str+"_param.txt").c_str()   );		// These files hold data for individual cells to compare V_min, V_max, V_90 etc. to compare between cells
	ofstream result_file2((str+"_state.txt").c_str()   );
	ofstream result_file3((str+"_UV.txt").c_str()   );
	//ofstream result_fileS19((str+"_S19.txt").c_str()   );
	
	Cell.Cell_type = 3;		// Cell type = (1) for endo; (2) for M,  (3) epi
	{
	//Loading in initial conditions for each cell in the array    	
			Cell.V=-86.2;
			Cell.V_new=-86.2;
			Cell.dV=0;
			Cell.Na_in = 7.67;	
			Cell.K_in = 138.3;
			Cell.Ca_in = 0.00007;
			Cell.Ca_sr = 1.3;
			Cell.Ca_ss = 0.00007;
			Cell.Ca_in_buffer = 0;
			Cell.Ca_ss_buffer = 0;
			Cell.Ca_sr_buffer = 0;
			Cell.m = 0 ;
			Cell.h = 0.75;
			Cell.j = 0.75;
			Cell.d = 0 ;
			Cell.f = 1;
			Cell.f2 = 1;
			Cell.f_Ca = 1;
			Cell.r = 0;
			Cell.s = 1;
			Cell.xr1 = 0;
			Cell.xr2 = 1;
			Cell.xs = 0;
			Cell.OO = 0;
			Cell.R_bar = 1;
			
			Cell.I_Na = 0;
			Cell.mI_Na = 0;
			Cell.I_Ca_L = 0;
			Cell.I_Kr = 0;
			Cell.I_Ks = 0;
			Cell.I_K1 = 0;
			Cell.I_Kp = 0;
			Cell.I_to = 0;
			Cell.I_Na_Ca = 0;
			Cell.I_Na_K = 0;
			Cell.I_p_Ca = 0;
			Cell.I_Ca_b = 0;
			Cell.I_Na_b = 0;
			Cell.I_stim = 0;
			Cell.I_tr = 0;
			Cell.I_leak = 0;
			Cell.I_up = 0;
			Cell.I_rel = 0;
			
			Cell.I_Na_ion_total = 0;
			Cell.I_K_ion_total = 0;
			Cell.I_Ca_ion_total = 0;
			Cell.I_total = 0.0;
			Cell.I_axial=0.0;									
			
			Cell.peak_slope = 0;
			Cell.V_min = -88.654973;
			Cell.t_min = 0;
			Cell.V_thr = -88.654973;
			Cell.t_thr = 0;
			Cell.V_max = -88.654973;
			Cell.t_max = 0;
			Cell.t_EAD= 0;
			Cell.V_EAD= -88.654973;
			Cell.t_EAD2= 0;
			Cell.V_EAD2= -88.654973;
			Cell.V_90 = -88.654973;
			Cell.t_90 = 0;
			Cell.t_90_old = 0;
			Cell.dV_old = 0;
			Cell.flag2 = 0;
			Cell.flag_EAD=0;
			Cell.flag_EAD2=0;

//WT Channel Markov Initial Conditions		
		//Drug Free States of WT Channel
			Cell.IC3 = 0;
			Cell.IC2 = 0;
			Cell.IF = 0;
			Cell.IM1 = 0;
			Cell.IM2 = 0;
			Cell.C3 = 1;
			Cell.C2 = 0;
			Cell.O = 0;
			Cell.OS =0;
			//Charged Drug States of WT Channel
			Cell.DIC3 = 0;
			Cell.DIC2 =0;
			Cell.DIF = 0;
			Cell.DIM1 = 0;
			Cell.DIM2 = 0;
			Cell.DC3 = 0;
			Cell.DC2 =0;
			Cell.DC1 = 0;
			Cell.DO = 0;
			Cell.DOS = 0;
			//Neutral Drug States of WT Channel
			Cell.D_IC3 = 0;
			Cell.D_IC2 =0;
			Cell.D_IF = 0;
			Cell.D_IM1 = 0;
			Cell.D_IM2 = 0;
			Cell.D_C3 = 0;
			Cell.D_C2 =0;
			Cell.D_O = 0;
			Cell.D_OS = 0;
			Cell.D_C1 = 0;
			
			Cell.C1 = 1- (Cell.O + Cell.OS + Cell.C3 + Cell.C2 + Cell.IC3 + Cell.IC2 + Cell.IF + Cell.IM1 + Cell.IM2 +
				   Cell.DO + Cell.DOS + Cell.DC1 + Cell.DC2 + Cell.DC3 + Cell.DIC3 + Cell.DIC2 + Cell.DIF + Cell.DIM1 + Cell.DIM2 +
				   Cell.D_O + Cell.D_OS + Cell.D_C1 + Cell.D_C2 + Cell.D_C3 + Cell.D_IC3 + Cell.D_IC2 + Cell.D_IF + Cell.D_IM1 + Cell.D_IM2);

//Mutant Channel Markov Initial Conditions	
			//Drug Free States of Mutant Channel
			Cell.mIC3 = 0;
			Cell.mIC2 = 0;
			Cell.mIF = 0;
			Cell.mIM1 = 0;
			Cell.mIM2 = 0;
			Cell.mC3 = 1;
			Cell.mC2 = 0;
			Cell.mO = 0;
			Cell.mOS =0;
			Cell.mBC3 =0;
			Cell.mBC2 = 0;
			Cell.mBC1 = 0;
			Cell.mBO = 0;
			//Charged Drug States of Mutant Channel
			Cell.DmIC3 = 0;
			Cell.DmIC2 =0;
			Cell.DmIF = 0;
			Cell.DmIM1 = 0;
			Cell.DmIM2 = 0;
			Cell.DmC3 = 0;
			Cell.DmC2 =0;
			Cell.DmC1 = 0;
			Cell.DmO = 0;
			Cell.DmOS = 0;
			Cell.DmBC3 =0;
			Cell.DmBC2 = 0;
			Cell.DmBC1 = 0; 
			Cell.DmBO = 0;
			//Neutral Drug States of Mutant Channel
			Cell.D_mIC3 = 0;
			Cell.D_mIC2 =0;
			Cell.D_mIF = 0;
			Cell.D_mIM1 = 0;
			Cell.D_mIM2 = 0;
			Cell.D_mC3 = 0;
			Cell.D_mC2 =0;
			Cell.D_mO = 0;
			Cell.D_mOS = 0;
			Cell.D_mC1 = 0;
			Cell.D_mBC3 =0;
			Cell.D_mBC2 = 0;
			Cell.D_mBC1 = 0; 
			Cell.D_mBO = 0;		
			
			//Cell.mBIF = 0;
			//Cell.DmBIF = 0;
			//Cell.D_mBIF = 0;
		
			Cell.mC1 = 1- (Cell.mO + Cell.mOS + Cell.mC3 + Cell.mC2 + Cell.mIC3 + Cell.mIC2 + Cell.mIF + Cell.mIM1 + Cell.mIM2 + Cell.mBC3 + Cell.mBC2 + Cell.mBC1 + Cell.mBO +
							  Cell.DmO + Cell.DmOS + Cell.DmC1 + Cell.DmC2 + Cell.DmC3 + Cell.DmIC3 + Cell.DmIC2 + Cell.DmIF + Cell.DmIM1 + Cell.DmIM2 + Cell.DmBC3 + Cell.DmBC2 + Cell.DmBC1 + Cell.DmBO +
							  Cell.D_mO + Cell.D_mOS + Cell.D_mC1 + Cell.D_mC2 + Cell.D_mC3 + Cell.D_mIC3 + Cell.D_mIC2 + Cell.D_mIF + Cell.D_mIM1 + Cell.D_mIM2 + Cell.D_mBC3 + Cell.D_mBC2 + Cell.D_mBC1 + Cell.D_mBO + Cell.mBIF + Cell.DmBIF + Cell.D_mBIF);
	
	} //Initial Conditions
	
	/*********************************************************************************************************************************
	 Begin Time Loop Here
	 *********************************************************************************************************************************/		
	
	for (S1_cycle = 0; S1_cycle <=250; S1_cycle = S1_cycle + 1){		
		
        if (S1_cycle ==0){interval=waitTime;}
		/*
		 else if (S1_cycle <=20){interval=800;
		 test= Cell.APD_90;
		 cout <<"test = " << test << endl;
		 result_file_SS << S1_cycle <<", " <<  Cell.APD_90 << endl;}
		 
		 else if (S1_cycle == 21){interval = (test+50);
		 cout << "interval is now" << interval << endl;}
		 
		 else if (S1_cycle == 22){interval = 3500;}
		 //	test = Cell.APD_90;
		 
		 else if (S1_cycle == 23){interval = (3500);
		 cout << "interval is now" << interval << endl;}
		 */
		else {interval=500;}
		
		tTemp = fmod((t - waitTime) , (interval) );
		for (tTemp = 0; tTemp <= interval; tTemp=tTemp+dt){ 
			
				//Resets cell specific parameters every beat such as APD, DI, V90 etc.					
				Calculate_Reset (&Cell, t, tTemp);
				
				//Updating current calculations for each cell	
				Calculate_I_Na(&Cell, t, tTemp, Input_param );			
				//Calculate_mI_Na(&Cell, t, tTemp );
				Calculate_I_Ca_L(&Cell, t, tTemp );
				Calculate_I_Kr(&Cell, t, tTemp );
				Calculate_I_Ks(&Cell, t, tTemp );
				Calculate_I_K1(&Cell, t, tTemp );
				Calculate_I_Kp(&Cell, t, tTemp );
				Calculate_I_to(&Cell, t, tTemp );
				Calculate_I_Na_Ca(&Cell, t, tTemp );		
				Calculate_I_Na_K(&Cell, t, tTemp );		
				Calculate_I_p_Ca(&Cell, t, tTemp );
				Calculate_I_Ca_b(&Cell, t, tTemp );
				Calculate_I_Na_b(&Cell, t, tTemp );
				Calculate_I_total(&Cell, t, tTemp);
				
				//Updating ionic concentrations	
				Calculate_Na_in(&Cell, t, tTemp );
				Calculate_K_in(&Cell, t, tTemp ); 
				Calculate_I_tr(&Cell, t, tTemp);
				Calculate_I_leak(&Cell, t, tTemp);
				Calculate_I_up(&Cell, t, tTemp);
				Calculate_I_rel(&Cell, t, tTemp);
				Calculate_Ca_sr(&Cell, t, tTemp );
				Calculate_Ca_ss(&Cell, t, tTemp );
				Calculate_Ca_in(&Cell, t, tTemp );
				
				//Calculates t_min, V_min, t_max, V_max etc. etc.; Update calculates axial currents and updates the voltage
				Calculate_Points (&Cell, t, tTemp); 	
				Calculate_Update (&Cell, t, tTemp);
				
				t = t+dt;
			
						
			/*******************************************************************************************************************************************************
						Beginning of result files			 
			 *******************************************************************************************************************************************************/	
			/*
			if (t>=59999.999) {
				
				cout << "time = " << t << endl;
				
				cout << "IC3 = " << Cell.IC3 << endl;
				cout << "IC2 = " << Cell.IC2 << endl;
				cout << "IF = " << Cell.IF << endl;
				cout << "IM1 = " << Cell.IM1 << endl;
				cout << "IM2 = " << Cell.IM2 << endl;
				cout << "C3 = " << Cell.C3 << endl;
				cout << "C2 = " << Cell.C2 << endl;
				cout << "C1 = " << Cell.C1 <<endl;
				cout << "O = " << Cell.O << endl;
				cout << "OS = " << Cell.OS << endl;;
				//Charged Drug States
				cout << "DIC3 = " << Cell.DIC3 << endl;
				cout << "DIC2 = " << Cell.DIC2 << endl;
				cout << "DIF = " << Cell.DIF << endl;
				cout << "DIM1 = " << Cell.DIM1 << endl;
				cout << "DIM2 = " << Cell.DIM2 << endl;
				cout << "DC3 = " << Cell.DC3 << endl;
				cout << "DC2 = " << Cell.DC2 << endl;
				cout << "DC1 = " << Cell.DC1 << endl;
				cout << "DO = " << Cell.DO << endl;
				cout << "DOS = " << Cell.DOS << endl;
				//Neutral Drug States
				cout << "D_IC3 = " << Cell.D_IC3 << endl;
				cout << "D_IC2 = " << Cell.D_IC2 << endl;
				cout << "D_IF = " << Cell.D_IF << endl;
				cout << "D_IM1 = " << Cell.D_IM1 << endl;
				cout << "D_IM2 = " << Cell.D_IM2 << endl;
				cout << "D_C3 = " << Cell.D_C3 << endl;
				cout << "D_C2 = " << Cell.D_C2 << endl;
				cout << "D_O = " << Cell.D_O << endl;
				cout << "D_OS = " << Cell.D_OS << endl;
				cout << "D_C1 = " << Cell.D_C1 << endl;	
				
				//Drug Free States of Mutant Channel
				cout << "mIC3 = " << Cell.mIC3 << endl;
				cout << "mIC2 = " << Cell.mIC2 << endl;
				cout << "mIF = " << Cell.mIF << endl;
				cout << "mIM1 = " << Cell.mIM1 << endl;
				cout << "mIM2 = " << Cell.mIM2 << endl;
				cout << "mC3 = " << Cell.mC3 << endl;
				cout << "mC2 = " << Cell.mC2 << endl;
				cout << "mC1 = " << Cell.mC1 <<endl;
				cout << "mO = " << Cell.mO << endl;
				cout << "mOS = " << Cell.mOS << endl;;
				cout << "mBC3 = " << Cell.mBC3 << endl;;
				cout << "mBC2 = " << Cell.mBC2 << endl;
				cout << "mBC1 = " << Cell.mBC1 << endl;
				cout << "mBO = " << Cell.mBO << endl;
				//Charged Drug States
				cout << "DmIC3 = " << Cell.DmIC3 << endl;
				cout << "DmIC2 = " << Cell.DmIC2 << endl;
				cout << "DmIF = " << Cell.DmIF << endl;
				cout << "DmIM1 = " << Cell.DmIM1 << endl;
				cout << "DmIM2 = " << Cell.DmIM2 << endl;
				cout << "DmC3 = " << Cell.DmC3 << endl;
				cout << "DmC2 = " << Cell.DmC2 << endl;
				cout << "DmC1 = " << Cell.DmC1 << endl;
				cout << "DmO = " << Cell.DmO << endl;
				cout << "DmOS = " << Cell.DmOS << endl;
				cout << "DmBC3 = " << Cell.DmBC3 << endl;
				cout << "DmBC2 = " << Cell.DmBC2 << endl;
				cout << "DmBC1 = " << Cell.DmBC1 << endl; 
				cout << "DmBO = " << Cell.DmBO << endl;
				//Neutral Drug States
				cout << "D_mIC3 = " << Cell.D_mIC3 << endl;
				cout << "D_mIC2 = " << Cell.D_mIC2 << endl;
				cout << "D_mIF = " << Cell.D_mIF << endl;
				cout << "D_mIM1 = " << Cell.D_mIM1 << endl;
				cout << "D_mIM2 = " << Cell.D_mIM2 << endl;
				cout << "D_mC3 = " << Cell.D_mC3 << endl;
				cout << "D_mC2 = " << Cell.D_mC2 << endl;
				cout << "D_mO = " << Cell.D_mO << endl;
				cout << "D_mOS = " << Cell.D_mOS << endl;
				cout << "D_mC1 = " << Cell.D_mC1 << endl;
				cout << "D_mBC3 = " << Cell.D_mBC3 << endl;
				cout << "D_mBC2 = " << Cell.D_mBC2 << endl;
				cout << "D_mBC1 = " << Cell.D_mBC1 << endl; 
				cout << "D_mBO = " << Cell.D_mBO << endl;		
				
				
				cout << "f = " << Cell.f << endl;
				cout << "f2 = " << Cell.f2 << endl;
				cout << "f_Ca = " << Cell.f_Ca << endl;
				cout << "xs = " << Cell.xs << endl;
				cout << "xr1 = " << Cell.xr1 << endl;
				cout << "xr2 = " << Cell.xr2 << endl;
				cout << "r = " << Cell.r << endl;
				cout << "d = " << Cell.d <<endl;
				cout << "s = " << Cell.s << endl;
				cout << "OO = " << Cell.OO << endl;;
				cout << "R_bar = " << Cell.R_bar << endl;;
				cout << "Na_in = " << Cell.Na_in << endl;
				cout << "K_in = " << Cell.K_in << endl;
				cout << "Ca_in = " << Cell.Ca_in << endl;
				cout << "Ca_sr = " << Cell.Ca_sr << endl;
				cout << "Ca_ss = " << Cell.Ca_ss << endl;
				cout << "V = " << Cell.V << endl;
				cout << "V_new = " << Cell.V_new << endl;
				
				
				exit(0);
			}
			*/
			 //Gate Outputs
			
			
			//This file is for individual cells that calculates all of the min, max, V_90 parameters. 	
			/*
			if ( t>=(waitTime-500) && counter%10==0){result_file << t << ", " << Cell.V << ", " << Cell.I_Na << ", " << Cell.I_Ca_L<< ", " << Cell.I_Na_b << ", " << Cell.I_Ca_b << ", " << Cell.I_Kr << ", "<< Cell.I_Ks << ", " << Cell.I_K1 << ", " << Cell.I_Kp << ", " 
			<< Cell.I_p_Ca << ", " << Cell.I_Na_Ca << ", " << Cell.I_Na_K << ", " << Cell.I_to << ", "<< Cell.Na_in << ", " << Cell.K_in << ", " << Cell.Ca_in << ", " << Cell.I_total << ", " << Cell.Ca_ss << ", "
			<< Cell.d << ", " << Cell.f << ", "<< Cell.f2 << ", " << Cell.f_Ca << ", " << Cell.s << ", " << Cell.xr1 << ", " << Cell.xr2 << ", " << Cell.xs << ", " << Cell.R_bar << ", " << Cell.OO << endl;}
			*/
			
			if (counter%100==0 && S1_cycle >=0){result_file << t << ", " << Cell.V << ", " << (      (1-mutant)*(Cell.I_Na) + (mutant)*(Cell.mI_Na)   ) << ", " << Cell.I_Ca_L << ", " << Cell.I_Na_Ca << ", "<< Cell.Ca_in << ", "<< Cell.Na_in << ", "  << Cell.I_Ks << endl;}
			if (counter%1000==0 && S1_cycle >=0){result_file2 << t << ", " << Cell.C3 + Cell.C2 + Cell.C1 << ", " << Cell.IC3 + Cell.IC2 << ", " << Cell.IF << ", " << Cell.OS << ", "<< Cell.IM1 + Cell.IM2 << ", "<< Cell.O << ", "  
															<< Cell.DC3 + Cell.DC2 + Cell.DC1 << ", " << Cell.DIC3 + Cell.DIC2 << ", " << Cell.DIF << ", " << Cell.DOS << ", "<< Cell.DIM1 + Cell.DIM2 << ", "<< Cell.DO << ", "
															<< Cell.D_C3 + Cell.D_C2 + Cell.D_C1 << ", " << Cell.D_IC3 + Cell.D_IC2 << ", " << Cell.D_IF << ", " << Cell.D_OS << ", "<< Cell.D_IM1 + Cell.D_IM2 << ", "<< Cell.D_O << ", "<< endl;}
			
			if (Cell.flag2==1) {			
				if (S1_cycle ==1){Upstroke_Velocity_0 = Cell.peak_slope;}
				if (S1_cycle ==200){Upstroke_Velocity_400 = Cell.peak_slope;}
				if (S1_cycle ==250){Upstroke_Velocity_500 = Cell.peak_slope;}
				
				
				result_file0 << Cell.t_min << ", " << Cell.V_min << ", " << Cell.t_thr<< ", " << Cell.V_thr << ", " << Cell.t_max << ", " << Cell.V_max << ", "<< Cell.t_90 << ", " << Cell.V_90 << ", "<< Cell.t_EAD << ", " << Cell.V_EAD << ", "
				<< Cell.t_EAD2 << ", "<<Cell.V_EAD2<< ", " << Cell.L_EAD << ", " << Cell.APD_90 << ", " << Cell.DI << ", " << S1_cycle<< ", "<< Cell.peak_slope<< ", " << Cell.t_I_Na_peak << ", " << Cell.I_Na_peak << endl;
				//if (S1_cycle == 19){result_fileS19 << Cell.APD_90 << ", " << Cell.DI << ", " << Cell.L_EAD << ", " << Cell.t_EAD << ", " << Cell.V_EAD <<  ", " << Cell.CV_EAD << ", " << Cell.CV<< endl;} cout <<"S19"<< endl;								
				
				if (S1_cycle == 250){result_file3 << S1_cycle << ", " << Upstroke_Velocity_0  << ", " << Upstroke_Velocity_400 << ", " << Upstroke_Velocity_500 << ", " << (Upstroke_Velocity_500/Upstroke_Velocity_0)*100 << ", " << fabs(Upstroke_Velocity_500 - Upstroke_Velocity_400) << endl;}
				
				Cell.flag2 = 2;
			} // End of individual result files
			
			
			

			/*******************************************************************************************************************************************************
						End of result files
			*******************************************************************************************************************************************************/
			
			counter++;
			
			//if (fmod(t,5000)<0.001) {cout << "S1_cycle=" <<S1_cycle<< ", " << "t=" << t << ", " "interval =" << interval <<endl;}
			
		}//End of For tTemp ==
		tTemp = 0;
		
	}//End of for S1
	S1_cycle = 1;	//Resets S1 cycle number
	
	//}	//End of for interval 
	
    return;    
	
} /******************************************************************************************************************************************************
 End of Mex
 ********************************************************************************************************************************************************/

