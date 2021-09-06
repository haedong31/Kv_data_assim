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
	ofstream result_file2((str+"SSA_param.txt").c_str()   );
	

	int Protocol = 1;
	for (Protocol = 1; Protocol <=1; Protocol = Protocol + 1){
		t = 0, dt = 0.01, tTemp = 0; counter = 0;
	
	
	
	//Steady State Inactivation Protocol*****************************************************************************************************************************
	if (Protocol ==1){
		t = 0;
		Cell.Drug = 0.0*(1e-6);
		for (cycle = 0; cycle <= 20; cycle = cycle + 1){		//10
			
			if (cycle ==0){interval=waitTime;}
			else {interval=5025;}
			
			V_initial = -130 + cycle*(5.);						//-135 , 10
			tTemp=fmod((t-waitTime), (interval)); 
			
			for (tTemp = 0; tTemp <= (interval); tTemp = tTemp + dt){
				
				if (t<= waitTime) {dt = 0.005, Cell.V = -100;}  //0.05
				else if (tTemp <= t_hold) {dt = 0.005, Cell.V = V_initial;}
				else if (tTemp> t_hold && tTemp <= t_hold + 25) {dt = 0.005, Cell.V = -10;}
				else {dt = 0.005, Cell.V= -100;}
				
				Calculate_I_Na(&Cell, t, dt, Input_param);
				
				if (tTemp>=(t_hold-dt) && tTemp < t_hold){I_Na_peak = 0;} //Resets I_Na_peak just before the pulse to -10. Just in case you get some current at your holding step potential!
				Calculate_peak(&Cell, t, dt);
			
				//if (t>=0 && counter%1000==0)	{result_file << t << ", " << Cell.V << ", " << Cell.I_Na << ", " << Cell.C3 + Cell.C2 + Cell.C1 << ", "<< Cell.D_IC3 + Cell.D_IC2 << ", " << Cell.D_IF << ", " << Cell.D_IM1 + Cell.D_IM2 + Cell.D_OS << ", " << Cell.D_O << ", "<< Cell.DC3 + Cell.DC2 + Cell.DC1<< ", " << Cell.DIC3 + Cell.DIC2 << ", " << Cell.DIF << ", " << Cell.DIM1 + Cell.DIM2 + Cell.DOS << endl;}
				
				if (tTemp >= (interval-dt) ){
					if (cycle == 1) {I_Na_peak_0 = I_Na_peak;}
					norm = (I_Na_peak/I_Na_peak_0);
					if (cycle >=1){result_file2 << t_peak << ", " << I_Na_peak << ", " << V_initial << ", " << norm << endl;}		//norm is the normalized current. Plot V_initial vs. norm to get the SSA curve
					
				}
				//Counter for t so you know where the program is
				//if (fmod(t, 5000) < dt) {cout << t  << ", " << cycle << ", " << tTemp << endl;}
				
				t = t+ dt;
				counter++; 
			} // End for tTemp
			I_Na_peak = 0;
			
		}// End for cycle
		
	}//*****************************************************************************************************************************************
	
			
	} //End of for protocol loop
	
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


