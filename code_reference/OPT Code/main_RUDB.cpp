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
	//ofstream result_file2((str+"RUDB80_param.txt").c_str()   );
	//ofstream result_file3((str+"RUDB100_param.txt").c_str()   );
	ofstream result_file4((str+"RUDB120_param.txt").c_str()   );
	
	ofstream result_file2((str+"RUDB100_param.txt").c_str()   );
	ofstream result_file3((str+"RUDB100_NUFL_param.txt").c_str()   );
	
	int Protocol = 1;
	for (Protocol = 1; Protocol <=1; Protocol = Protocol + 1){
		t = 0, dt = 0.01, tTemp = 0; counter = 0;
	
//Recovery after UDB ********************************************************************************************************************************************
	if (Protocol ==1){
		t = 0;
		Cell.Drug = 0.0*(1e-6);
		for (V_cycle = 0; V_cycle <1; V_cycle = V_cycle + 1){
			
		//	if (V_cycle == 0){V_rec = -80;}
		//	else if (V_cycle == 1){V_rec = -100;}
		//	else if (V_cycle ==2) {V_rec = -120;}
         //   else {V_rec = -100;}
			
			if (V_cycle == 0){V_rec = -100; Cell.Drug = 0.0*(1e-6);}
			else if (V_cycle == 1){V_rec = -100; Cell.Drug = 100.0*(1e-6);}
			else if (V_cycle ==2) {V_rec = -120;}
			else {V_rec = -100;}
			
			
		
		for (cycle = 0; cycle <= 101; cycle = cycle + 1){
			//Input BCL of 40 (25Hz)
			if (cycle ==0){interval=waitTime;}
			else if (cycle >=1 && cycle<=100){interval = BCL_3;}
			else {interval=10026;}
			tTemp=fmod((t-waitTime), (interval));
			
			for (tTemp = 0; tTemp <= (interval); tTemp = tTemp + dt){
				
				if (t<= waitTime) {dt = 0.005, Cell.V = -100;}  //0.005
				else if (tTemp <=V_duration ) {dt = 0.005, Cell.V = -10;} //0.005
				else {dt = 0.005, Cell.V= V_rec;}				//This is the recovery Voltage, 0.005
				
				Calculate_I_Na(&Cell, t, dt, Input_param);
				
				//Result files begin here	
				//if (t>=0 && counter%100==0)	{result_file << t << ", " << Cell.V << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << ", " << Cell.DIF << ", " << Cell.DIM1 <<", " << Cell.DIM2 << ", " << Cell.DOS << ", " << Cell.DIC3 + Cell.DIC2 << ", " << Cell.DC1 + Cell.DC2 + Cell.DC3 << Cell.sum << endl;}
				
				if (V_cycle ==0){
					//if (cycle ==101 && counter%10==0 && (tTemp - V_duration) >=0)	{result_file2 << (tTemp-V_duration) << ", " << (Cell.mC1 + Cell.mC2 + Cell.mC3) << ", " << (Cell.DIF + Cell.DIC2 + Cell.DIC3) << ", "  <<Cell.DIM1 << ", " << Cell.DIM2 << ", " << Cell.sum << endl;}
					
					if (cycle ==101 && (tTemp-V_duration) >=0) {
						if (  (tTemp >=26 && tTemp < 26+dt) || (tTemp >=30 && tTemp < 30+dt) || (tTemp >=35 && tTemp < 35+dt) || (tTemp >=50 && tTemp < 50+dt) || (tTemp >=75 && tTemp < 75 +dt) || (tTemp >=125 && tTemp < 125+dt) || (tTemp >=175 && tTemp < 175+dt) || (tTemp >=275 && tTemp < 275+dt) || 
							(tTemp >=525 && tTemp < 525+dt) || (tTemp >=1025 && tTemp < 1025+dt) || (tTemp >=2025 && tTemp < 2025+dt) || (tTemp >=3025 && tTemp < 3025+dt) || (tTemp >=5025 && tTemp < 5025 +dt) || (tTemp >=7525 && tTemp < 7525+dt) || (tTemp >=10025 && tTemp < 10025+dt) )
							
						{result_file2 << (tTemp-V_duration) << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << endl;}
						
					}
					
				}
				
				
				
				if (V_cycle ==1){
					//if (cycle ==101 && counter%10==0 && (tTemp - V_duration) >=0)	{result_file3 << (tTemp-V_duration) << ", " << (Cell.mC1 + Cell.mC2 + Cell.mC3) << ", " << (Cell.DIF + Cell.DIC2 + Cell.DIC3) << ", "  <<Cell.DIM1 << ", " << Cell.DIM2 << ", " << Cell.sum << endl;}
					
					if (cycle ==101 && (tTemp-V_duration) >=0) {
						if (  (tTemp >=26 && tTemp < 26+dt) || (tTemp >=30 && tTemp < 30+dt) || (tTemp >=35 && tTemp < 35+dt) || (tTemp >=50 && tTemp < 50+dt) || (tTemp >=75 && tTemp < 75 +dt) || (tTemp >=125 && tTemp < 125+dt) || (tTemp >=175 && tTemp < 175+dt) || (tTemp >=275 && tTemp < 275+dt) || 
							(tTemp >=525 && tTemp < 525+dt) || (tTemp >=1025 && tTemp < 1025+dt) || (tTemp >=2025 && tTemp < 2025+dt) || (tTemp >=3025 && tTemp < 3025+dt) || (tTemp >=5025 && tTemp < 5025 +dt) || (tTemp >=7525 && tTemp < 7525+dt) || (tTemp >=10025 && tTemp < 10025+dt) )
								
						{result_file3 << (tTemp-V_duration) << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << endl;}
					
					}
					
				}
				
				
				
				if (V_cycle ==2){
				//	if (cycle ==101 && counter%10==0 && (tTemp - V_duration) >=0)	{result_file4 << (tTemp-V_duration) << ", " << (Cell.mC1 + Cell.mC2 + Cell.mC3) << ", " << (Cell.DIF + Cell.DIC2 + Cell.DIC3) << ", "  <<Cell.DIM1 << ", " << Cell.DIM2 << ", " << Cell.sum << endl;}
					
					if (cycle ==101 && (tTemp-V_duration) >=0) {
						if (  (tTemp >=26 && tTemp < 26+dt) || (tTemp >=30 && tTemp < 30+dt) || (tTemp >=35 && tTemp < 35+dt) || (tTemp >=50 && tTemp < 50+dt) || (tTemp >=75 && tTemp < 75 +dt) || (tTemp >=125 && tTemp < 125+dt) || (tTemp >=175 && tTemp < 175+dt) || (tTemp >=275 && tTemp < 275+dt) || 
							(tTemp >=525 && tTemp < 525+dt) || (tTemp >=1025 && tTemp < 1025+dt) || (tTemp >=2025 && tTemp < 2025+dt) || (tTemp >=3025 && tTemp < 3025+dt) || (tTemp >=5025 && tTemp < 5025 +dt) || (tTemp >=7525 && tTemp < 7525+dt) || (tTemp >=10025 && tTemp < 10025+dt) )
								
						{result_file4 << (tTemp-V_duration) << ", " << (Cell.C1 + Cell.C2 + Cell.C3) << endl;}
						
					}
					
				}
				
				//Counter for t so you know where the program is
				//if (fmod(t, 5000) < dt) {cout << t  << ", " << cycle << ", " << tTemp << endl;}
				
				t = t+ dt;
				counter++; 
			} // End for tTemp
		
		}// End for cycle
		
			t = 0;
			
		}//End for V_cycle
		
	
	}//******************************************************************************************************************************************************
		
	} //End of for protocol loop
	
	result_file2.close(), result_file3.close(), result_file4.close();
	
	
	
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


