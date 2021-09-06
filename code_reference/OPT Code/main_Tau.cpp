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
	ofstream result_file2((str+"Tau_param.txt").c_str()   );
	

	int Protocol = 1;
	for (Protocol = 1; Protocol <=1; Protocol = Protocol + 1){
		t = 0, dt = 0.01, tTemp = 0; counter = 0;
	
	
	
	//Tau to 50% Protocol*****************************************************************************************************************************
	if (Protocol ==1){
		t = 0;
		Cell.Drug = 0.0*(1e-6);
		for (cycle = 0; cycle <= 9; cycle = cycle + 1){					//6
			
			if (cycle ==0){interval=waitTime;}
			else {interval=6000;}
			
			V_temp = 25 - cycle*(5);									//30, 10
			tTemp=fmod((t-waitTime), (interval)); 
			
			for (tTemp = 0; tTemp <= (interval); tTemp = tTemp + dt){
				
				if (t<= waitTime) {dt = 0.005, Cell.V = -90;}							//-90
				else if (tTemp <= t_hold) {dt = 0.005, Cell.V = -90;}					//-90
				else if (tTemp> t_hold && tTemp <= t_hold + 40) {dt = 0.005, Cell.V = V_temp;}			//40
				else {dt = 0.005, Cell.V= -90;}											//-90
				
				
				Calculate_I_Na(&Cell, t, dt, Input_param);
				if (tTemp>=(t_hold-dt) && tTemp < t_hold){I_Na_peak = 0; flag_tau = 1;} //Resets I_Na_peak just before the pulse to -10. Just in case you get some current at your holding step potential!
				if (t>waitTime && fabs(Cell.I_Na) - fabs(I_Na_peak_tau) < 0.001 && tTemp > t_hold && flag_tau == 0 ){t_tau = t; flag_tau = 1; tau = t_tau - t_peak;}
				Calculate_peak(&Cell, t, dt);
				
				//Result files begin here	
				if (t>=0 && counter%10==0)	{result_file << t << ", " << Cell.V << ", " << Cell.I_Na << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << ", " << Cell.IF << ", " << Cell.IM1 <<", " << Cell.IM2 << ", " << Cell.OS << ", " << (Cell.IC3 + Cell.IC2) << ", " << Cell.sum <<endl;}
				
				if (tTemp >= (interval-dt) ){
					if (cycle == 1) {I_Na_peak_0 = I_Na_peak;}
					
					if (cycle >=1){result_file2 << cycle << ", " << V_temp << ", " << t_peak << ", " << I_Na_peak << ", " << t_tau << ", " << I_Na_peak_tau << ", " << tau << endl;}		//norm is the normalized current. Plot V_initial vs. norm to get the SSA curve
					flag_tau = 0;
				}
				
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
		I_Na_peak_tau = I_Na_peak*0.50;		//0.3679
		flag_tau = 0;
	t_peak = t;}
	
	
}


