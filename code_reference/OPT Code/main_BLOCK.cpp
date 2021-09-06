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
//#include "WT_Flec_New.h"
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
	ofstream result_file13((str+"Block_param.txt").c_str()  );
	
	
	

	int Protocol = 1;
	for (Protocol = 1; Protocol <=1; Protocol = Protocol + 1){
		t = 0, dt = 0.01, tTemp = 0; counter = 0;
	
//Full UDB protocol and TB protocol****************************************************************************************************************************** 
	if (Protocol ==1){
		t = 0;
		int cycle_drug;
		double Peak_0_TB;
		
		for(cycle_drug = 1; cycle_drug <=4; cycle_drug = cycle_drug + 1){
			
			if (cycle_drug == 1){Cell.Drug = 0.0*(1e-6);}
			else if (cycle_drug == 2){Cell.Drug = 10*(1e-6);}
			else if (cycle_drug == 3){Cell.Drug = 100*(1e-6);}
			else if (cycle_drug == 4){Cell.Drug = 1000*(1e-6);}
			
            else if (cycle_drug == 5){Cell.Drug = 0.1*(1e-6);}
            else if (cycle_drug == 6){Cell.Drug = 1*(1e-6);}
			else if (cycle_drug == 7){Cell.Drug = 3*(1e-6);}
			else if (cycle_drug == 8){Cell.Drug = 30*(1e-6);}
			else if (cycle_drug == 9){Cell.Drug = 60*(1e-6);}
			else if (cycle_drug == 10){Cell.Drug = 100*(1e-6);}
			else if (cycle_drug == 11){Cell.Drug = 200*(1e-6);}
			else if (cycle_drug == 12){Cell.Drug = 300*(1e-6);}
			else if (cycle_drug == 13){Cell.Drug = 400*(1e-6);}
			else if (cycle_drug == 14){Cell.Drug = 500*(1e-6);}
			else if (cycle_drug == 15){Cell.Drug = 600*(1e-6);}
			else if (cycle_drug == 16){Cell.Drug = 700*(1e-6);}
			else if (cycle_drug == 17){Cell.Drug = 800*(1e-6);}
			else if (cycle_drug == 18){Cell.Drug = 900*(1e-6);}
			else if (cycle_drug == 19){Cell.Drug = 1000*(1e-6);}
            
			else if (cycle_drug == 20){Cell.Drug = 1000*(1e-6);}
			else {Cell.Drug = 0;}
			
			for (cycle = 0; cycle <=total_cycles; cycle = cycle + 1){		
				
				if (cycle ==0){interval=waitTime;}
				else {interval=200;}
				
				tTemp = fmod((t - waitTime) , (interval) );
				for (tTemp = 0; tTemp <= interval; tTemp=tTemp+dt){ 
					
					if (t<= waitTime) {dt = 0.005, Cell.V = V_o;}					//0.005
					else if (tTemp <= V_duration) {dt = 0.005, Cell.V = V_step;}	//0.005
					else {dt = 0.005, Cell.V = V_o;}								//0.005
					
					Calculate_I_Na(&Cell, t, dt, Input_param);
					Calculate_peak(&Cell, t, dt);
					
					late_avg = Cell.I_Na;
					
					
					if (cycle == 1 && tTemp >=23.5 && tTemp < (23.5+dt)){ 
						I_Na_peak_0 = I_Na_peak, late_avg_0 = late_avg; //cout << "peak_0 is " << I_Na_peak_0 << ", " << "late_0 is " << late_avg_0 << endl;
						if (cycle_drug ==1){Peak_0_TB = I_Na_peak_0;}
						
					}
					
					
					if (cycle == total_cycles && tTemp >=23.5 && tTemp < (23.5+dt)){ I_Na_peak_end = I_Na_peak, late_avg_end = late_avg; /*cout << "peak_end is " << I_Na_peak_end << ", " << "late_end is " << late_avg_end << endl;*/}
					
					//Result files begin here
					//if (t>=0 && counter%1000==0)	{result_file << t << ", " << Cell.V << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << ", " << Cell.DIF << ", " << Cell.DIM1 <<", " << Cell.DIM2 << ", " << Cell.DOS << ", " << Cell.DIC3 + Cell.DIC2 << ", " << Cell.DC1 + Cell.DC2 + Cell.DC3 << Cell.sum << endl;}
					
					//Counter for t so you know where the program is
					//if (fmod(t, 5000) < dt) {cout << t << ", " << "Sum = " << Cell.sum << ", " << "Cycle = " << cycle << ", " << "Interval = " << interval << ", " << "Drug = " << Cell.Drug<< endl;}
					
					t = t+ dt;
					counter++;
				}//End for tTemp	
				
				I_Na_peak = 0;
				late_avg = 0;
				
			}//End for cycle	
			//{cout << "Fractional UDB of Peak = " << (1-(I_Na_peak_end/I_Na_peak_0))*100 << ",   " << "Fractional UDB of Late = " << (1-(late_avg_end/late_avg_0))*100 << endl;}
			
			{result_file13 << Cell.Drug*(1e6) << ", " <<    /*(I_Na_peak_end/Peak_0_TB)*100   */       (I_Na_peak_end/I_Na_peak_0)*100    << ", " <<  (I_Na_peak_0/Peak_0_TB)*100 << ", " <<(late_avg_end/late_avg_0)*100 << endl;}			//UDB would be 1, 2, TB would be 1, 3
			
			t =0;
			
		}//End for cycle_drug
		
		
	}//**************************************************************************************************************************************************************
		
		
		
	} //End of for protocol loop
	
	result_file13.close();
	
	
	
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


