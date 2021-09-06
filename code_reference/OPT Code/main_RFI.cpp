/**********************************************************************************
 Single Channel Formulation with Flecainide Drug Gating
 By: Jonathan D. Moreno
 Colleen Clancy Laboratory 
 **********************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include "mex.h"

#include "Global_Variables.h"

#ifdef FLECAINIDE
#include "WT_Flec_Implicit.h"
#endif

#ifdef LIDOCAINE 
#include "WT_Lido_Implicit.h"
#endif

#ifdef RANOLAZINE
#include "Mutant_Implicit.h"
#endif

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	double *Input_param;
	int i;
	Input_param = mxGetPr(prhs[0]);
	
	
	Cell_param Cell;
	int counter = 0;
	double t = 0, dt = 0.01, tTemp = 0;
	double V_cycle, V_rec;
	
	//char word[20];
	//cout << "To begin, enter file names:" << endl<< endl;
	//scanf ("%s", word);
	
	string str = "test";
	//cout << "Thanks!"<< endl;
	//chdir ("results");	
	
	ofstream result_file((str+".txt").c_str()   );
	ofstream result_file2((str+"RFI120_param.txt").c_str()   );
	ofstream result_file3((str+"RFI100_param.txt").c_str()   );
	ofstream result_file4((str+"RFI80_param.txt").c_str()   );
	
	ofstream result_file5((str+"RFI_L_param.txt").c_str()   );
	

	int Protocol = 1;
	for (Protocol = 1; Protocol <=1; Protocol = Protocol + 1){
		t = 0, dt = 0.01, tTemp = 0; counter = 0;
	
	
	
	//RFI Protocol*****************************************************************************************************************************
	if (Protocol ==1){
		double t_hold;
		t = 0;
		Cell.Drug = 0.0*(1e-6);
		for (V_cycle = 1; V_cycle <=1; V_cycle = V_cycle + 1){
			
			if (V_cycle == 0){V_rec = -120; t_hold = 100;}
			else if (V_cycle == 1){V_rec = -100; t_hold = 100;}
			else if (V_cycle ==2) {V_rec = -80; t_hold = 100;}			//Oginosawa is hold for 500ms, our protocol is hold for 100ms
			else {V_rec = -100;}
			
			
			for (cycle = 0; cycle <=1; cycle = cycle + 1){
				
				if (cycle ==0){interval=waitTime;}
				else {interval=10501;}
				tTemp=fmod((t-waitTime), (interval));
				
				for (tTemp = 0; tTemp <= (interval); tTemp = tTemp + dt){
					
					if (t<= waitTime) {dt = 0.005, Cell.V = V_rec;}
					else if (tTemp <= t_hold ) {dt = 0.005, Cell.V = -20;}
					else {dt = 0.005, Cell.V= V_rec;}				//This is the recovery Voltage
					
					Calculate_I_Na(&Cell, t, dt, Input_param);
					
					//if (t>=0 && counter%100==0)	{result_file << t << ", " << Cell.V << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << ", " << Cell.IF << ", " << Cell.IM1 <<", " << Cell.IM2 << ", " << Cell.OS << ", " << (Cell.IC3 + Cell.IC2) << ", " << Cell.I_Na << Cell.sum << endl;}
					
					if (V_cycle ==0){
						//if (cycle ==1 && counter%10==0 && (tTemp - t_hold) >=0.1)	{result_file2 << (tTemp-t_hold) << ", " << (Cell.mC1 + Cell.mC2 + Cell.mC3) << ", " << (Cell.DIF + Cell.DIC2 + Cell.DIC3) << ", "  <<Cell.DIM1 << ", " << Cell.DIM2 << ", " << Cell.sum << endl;}
						
						if (cycle ==1 && (tTemp-t_hold) >=0) {
							if (  (tTemp >=(t_hold + 1) && tTemp < (t_hold + 1)+dt) || (tTemp >=(t_hold + 5) && tTemp < (t_hold + 5)+dt) || (tTemp >=(t_hold + 10) && tTemp < (t_hold + 10)+dt) || (tTemp >=(t_hold + 25) && tTemp < (t_hold + 25)+dt) || (tTemp >=(t_hold + 50) && tTemp < (t_hold + 50)+dt) || (tTemp >=(t_hold + 100) && tTemp < (t_hold + 100)+dt) || 
								(tTemp >=(t_hold + 150) && tTemp < (t_hold + 150)+dt) || (tTemp >=(t_hold + 250) && tTemp < (t_hold + 250)+dt) || (tTemp >=(t_hold + 500) && tTemp < (t_hold + 500)+dt) || (tTemp >=(t_hold + 1000) && tTemp < (t_hold + 1000)+dt) || (tTemp >=(t_hold + 2000) && tTemp < (t_hold + 2000)+dt) || (tTemp >=(t_hold + 3000) && tTemp < (t_hold + 3000)+dt) || 
								(tTemp >=(t_hold + 5000) && tTemp < (t_hold + 5000)+dt) || (tTemp >=(t_hold + 7500) && tTemp < (t_hold + 7500) +dt) || (tTemp >=(t_hold + 10000) && tTemp < (t_hold + 10000)+dt) )
								
							{result_file2 << (tTemp-t_hold) << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << endl;}
							
						}
						
					}
					
					else if (V_cycle ==1){
						//if (cycle ==1 && counter%10==0 && (tTemp - t_hold) >=0.1)	{result_file3 << (tTemp-t_hold) << ", " << (Cell.mC1 + Cell.mC2 + Cell.mC3) << ", " << (Cell.DIF + Cell.DIC2 + Cell.DIC3) << ", "  <<Cell.DIM1 << ", " << Cell.DIM2 << ", " << Cell.sum << endl;}
						
						if (cycle ==1 && (tTemp-t_hold) >=0) {
							if (  (tTemp >=(t_hold + 1) && tTemp < (t_hold + 1)+dt) || (tTemp >=(t_hold + 5) && tTemp < (t_hold + 5)+dt) || (tTemp >=(t_hold + 10) && tTemp < (t_hold + 10)+dt) || (tTemp >=(t_hold + 25) && tTemp < (t_hold + 25)+dt) || (tTemp >=(t_hold + 50) && tTemp < (t_hold + 50)+dt) || (tTemp >=(t_hold + 100) && tTemp < (t_hold + 100)+dt) || 
								(tTemp >=(t_hold + 150) && tTemp < (t_hold + 150)+dt) || (tTemp >=(t_hold + 250) && tTemp < (t_hold + 250)+dt) || (tTemp >=(t_hold + 500) && tTemp < (t_hold + 500)+dt) || (tTemp >=(t_hold + 1000) && tTemp < (t_hold + 1000)+dt) || (tTemp >=(t_hold + 2000) && tTemp < (t_hold + 2000)+dt) || (tTemp >=(t_hold + 3000) && tTemp < (t_hold + 3000)+dt) || 
								(tTemp >=(t_hold + 5000) && tTemp < (t_hold + 5000)+dt) || (tTemp >=(t_hold + 7500) && tTemp < (t_hold + 7500) +dt) || (tTemp >=(t_hold + 10000) && tTemp < (t_hold + 10000)+dt) )
								
							{result_file3 << (tTemp-t_hold) << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << endl;}
							
						}
						
					}
		
					else if (V_cycle ==2){
						//if (cycle ==1 && counter%10==0 && (tTemp - t_hold) >=0.1)	{result_file4 << (tTemp-t_hold) << ", " << (Cell.mC1 + Cell.mC2 + Cell.mC3) << ", " << (Cell.DIF + Cell.DIC2 + Cell.DIC3) << ", "  <<Cell.DIM1 << ", " << Cell.DIM2 << ", " << Cell.sum << endl;}
						
						if (cycle ==1 && (tTemp-t_hold) >=0) {
							if (  (tTemp >=(t_hold + 1) && tTemp < (t_hold + 1)+dt) || (tTemp >=(t_hold + 5) && tTemp < (t_hold + 5)+dt) || (tTemp >=(t_hold + 10) && tTemp < (t_hold + 10)+dt) || (tTemp >=(t_hold + 25) && tTemp < (t_hold + 25)+dt) || (tTemp >=(t_hold + 50) && tTemp < (t_hold + 50)+dt) || (tTemp >=(t_hold + 100) && tTemp < (t_hold + 100)+dt) || 
								(tTemp >=(t_hold + 150) && tTemp < (t_hold + 150)+dt) || (tTemp >=(t_hold + 250) && tTemp < (t_hold + 250)+dt) || (tTemp >=(t_hold + 500) && tTemp < (t_hold + 500)+dt) || (tTemp >=(t_hold + 1000) && tTemp < (t_hold + 1000)+dt) || (tTemp >=(t_hold + 2000) && tTemp < (t_hold + 2000)+dt) || (tTemp >=(t_hold + 3000) && tTemp < (t_hold + 3000)+dt) || 
								(tTemp >=(t_hold + 5000) && tTemp < (t_hold + 5000)+dt) || (tTemp >=(t_hold + 7500) && tTemp < (t_hold + 7500) +dt) || (tTemp >=(t_hold + 10000) && tTemp < (t_hold + 10000)+dt) )
								
							{result_file4 << (tTemp-t_hold) << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << endl;}
							
						}
						
					}
					
					
					
					//Counter for t so you know where the program is
					//if (fmod(t, 5000) < dt) {cout << t  << ", " << cycle << ", " << tTemp << endl;}
					
					t = t+ dt;
					counter++; 
				} // End for tTemp
				
			}//End Cycle
			
			t = 0;
			
		} //End for V_cycle
		
		
		
	} //End of protocol 1 ****************************************************************************************************************
		
	//Recovery from Inactivation For Late Current****************************************************************************************************************
	else if (Protocol ==2){
		double t_hold;
		double Control_peak, Test_peak, t_control, t_test, t_rec;
		
		t = 0;
		Cell.Drug = 0*(1e-6);
		
		V_rec = -140; t_hold = 50;
		
		for (cycle = 0; cycle <=15; cycle = cycle + 1){
			
			if (cycle ==0){interval=waitTime;}
			else if (cycle ==1) {t_rec = 1, interval= t_hold + t_rec + t_hold + 5000;}		//Note, 5000 is the interpulse interval
			else if (cycle ==2) {t_rec = 5, interval= t_hold + t_rec + t_hold+ 5000;}
			else if (cycle ==3) {t_rec = 10, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==4) {t_rec = 25, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==5) {t_rec = 50, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==6) {t_rec = 100, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==7) {t_rec = 150, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==8) {t_rec = 250, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==9) {t_rec = 500, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==10) {t_rec = 1000, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==11) {t_rec = 2000, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==12) {t_rec = 3000, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==13) {t_rec = 5000, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==14) {t_rec = 7500, interval= t_hold + t_rec + t_hold + 5000;}
			else if (cycle ==15) {t_rec = 10000, interval= t_hold + t_rec + t_hold + 5000;}
			
			
			tTemp=fmod((t-waitTime), (interval));
			
			for (tTemp = 0; tTemp <= (interval); tTemp = tTemp + dt){
				
				if (t<= waitTime) {dt = 0.01, Cell.V = -140;}													//Wait Time
				else if (tTemp <= t_hold ) {dt = 0.01, Cell.V = -20;}											//First control pulse
				else if (tTemp > t_hold && tTemp < t_hold + t_rec){dt = 0.01, Cell.V = V_rec;}					//Recovery Interval
				else if (tTemp >=(t_hold + t_rec) && tTemp <=(t_hold+t_rec+t_hold) ){dt = 0.01, Cell.V = -20;}	//Second Pulse
				else {dt = 0.01, Cell.V= V_rec;}																//End
				
				Calculate_I_Na(&Cell, t, dt, Input_param);
				Calculate_peak(&Cell, t, dt);
				
				if (t>=0 && counter%100==0)	{result_file << t << ", " << Cell.V << ", " << Cell.I_Na << ", "  << Cell.sum << endl;}		 
				
				if ( t>waitTime && tTemp >=48 && tTemp <48+dt) {Control_peak = Cell.I_Na, t_control = t;}
				if ( t>waitTime && tTemp >= (t_hold+t_rec+48.0) && tTemp < (t_hold + t_rec + 48.0 + dt)  ) {
					Test_peak = Cell.I_Na, t_test = t;
					//result_file3 << t_rec << ", " << (Test_peak/Control_peak) << t_control << ", " << Control_peak << ", " << t_test << ", " << Test_peak  << endl;
				}
				
				
				
				t = t+ dt;
				counter++; 
			} // End for tTemp
			
			//{result_file7 << t_rec << ", " << t_control << ", " << Control_peak << ", " << t_test << ", " << Test_peak << ", " << (Test_peak/Control_peak)*100 << endl;}
			
		}//End Cycle
		
		
		
		//t = 0;
		
		
		
	}//end if protocol ==2//**************************************************************************************************************************************
		
		
		
		
} //End of For protocol loop		
	
	result_file2.close();	
	
	return;
} //End of Mex Function

/***************************************************************************************************************************************************************
 Function calls begin here
 ***************************************************************************************************************************************************************/ 

void Calculate_I_Na(Cell_param *Cell_ptr, double t, double dt, double *Input_param){	
	
	if(which == 2) 
		return(WT_SCN5A_function(Cell_ptr, t, dt, Input_param)  );
	
	else {cout << "Invalid which" << endl; exit(1);}	
	
} // end calculate I_Na



void Calculate_peak (Cell_param *Cell_ptr, double t, double dt){
	
	if (fabs(Cell_ptr->I_Na) > fabs(I_Na_peak) ){
		I_Na_peak = Cell_ptr->I_Na;
		I_Na_peak_tau = I_Na_peak*0.5;		//0.3679
		flag_tau = 0;
	t_peak = t;}
	
	
}


