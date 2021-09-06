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
#include "WT_Lido_MATLAB.h"
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
	ofstream result_file12((str+"FDUDB_NUFL_param.txt").c_str()  );
	ofstream result_file14((str+"Onset_param.txt").c_str()  );
	
	

	int Protocol = 1;
	for (Protocol = 1; Protocol <=1; Protocol = Protocol + 1){
		t = 0, dt = 0.01, tTemp = 0; counter = 0;
	
//Full FD UDB*******************************************************************************************************************************************
	if (Protocol ==1){
		t = 0;
		int cycle_BCL, Hz, BCL_test;
		//Cell.Drug=10.0*(1e-6);
		
		for(cycle_BCL = 1; cycle_BCL <=2; cycle_BCL = cycle_BCL + 1){
			
			if (cycle_BCL == 1){BCL_test = 100; Hz = 10; Cell.Drug=10.0*(1e-6); }
			else if (cycle_BCL == 2){BCL_test = 100; Hz = 10; Cell.Drug=100.0*(1e-6);}
			else {BCL_test = 200;}
			
			for (cycle = 0; cycle <=total_cycles; cycle = cycle + 1){		
				
				if (cycle ==0){interval=waitTime;}
				else {interval=BCL_test;}
				
				tTemp = fmod((t - waitTime) , (interval) );
				for (tTemp = 0; tTemp <= interval; tTemp=tTemp+dt){ 
					
					if (t<= waitTime) {dt = 0.005, Cell.V = V_o;}					//0.01
					else if (tTemp <= V_duration) {dt = 0.005, Cell.V = V_step;}	//0.005
					else {dt = 0.005, Cell.V = V_o;}								//0.01
					
					Calculate_I_Na(&Cell, t, dt, Input_param);
					Calculate_peak(&Cell, t, dt);
					
					late_avg = Cell.I_Na;
					
					
					if (cycle == 1 && tTemp >=23.5 && tTemp < (23.5+dt)){ I_Na_peak_0 = I_Na_peak, late_avg_0 = late_avg; /*cout << "peak_0 is " << I_Na_peak_0 << ", " << "late_0 is " << late_avg_0 << endl;*/}
					if (cycle == total_cycles && tTemp >=23.5 && tTemp < (23.5+dt)){ I_Na_peak_end = I_Na_peak, late_avg_end = late_avg; /*cout << "peak_end is " << I_Na_peak_end << ", " << "late_end is " << late_avg_end << endl;*/}
					
					if ( (cycle_BCL ==1 && cycle>=1 && fmod(cycle, 50)==0 )  || (cycle_BCL ==1 && cycle ==1)    ){  
						
						if (tTemp >=23.5 && tTemp <(23.5 + dt)  ) {result_file14 << cycle << ", " << (I_Na_peak / I_Na_peak_0)*100 << endl; }
					}
					
					//Result files begin here
					//if (t>=0 && counter%1000==0){result_file << t << ", " << Cell.V << ", " << Cell.I_Na << ", " << Cell.C3 + Cell.C2 + Cell.C1 << ", "<< Cell.DIC3 + Cell.DIC2 << ", " << Cell.DIF << ", " << Cell.DIM1 + Cell.DIM2 + Cell.DOS << ", " << Cell.DO << ", "<< Cell.DC3 + Cell.DC2 + Cell.DC1<< endl;}
					
					//Counter for t so you know where the program is
					//if (fmod(t, 5000) < dt) {cout << t << ", " << "Sum = " << Cell.sum << ", " << "Cycle = " << cycle << ", " << "Interval = " << interval << ", " << "Drug = " << Cell.Drug << endl;}
					
					t = t+ dt;
					counter++;
				}//End for tTemp	
				
				I_Na_peak = 0;
				late_avg = 0;
				
			}//End for cycle	
			//{cout << "Fractional UDB of Peak = " << (1-(I_Na_peak_end/I_Na_peak_0))*100 << ",   " << "Fractional UDB of Late = " << (1-(late_avg_end/late_avg_0))*100 << endl;}
			
			{result_file12 << Hz << ", " <<  (1-(I_Na_peak_end/I_Na_peak_0))*100 << ", " << (1-(late_avg_end/late_avg_0))*100 << ", " << (I_Na_peak_end/I_Na_peak_0)*100 << ", " <<  (late_avg_end/late_avg_0)*100 << endl;}			//FDUDB would be 1, 2
			
			t = 0;
			
		}//End for cycle_BCL
		
	}//**************************************************************************************************************************************************************	
		
			
	} //End of for protocol loop
	
	result_file12.close(), result_file14.close();
	
	
	
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


