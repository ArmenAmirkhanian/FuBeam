/*
The MIT License(MIT)

Copyright(c) 2015 Armen Amirkhanian

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// FuBeam (*Fu*damental Beam)
// Author: Armen Amirkhanian

// Calculates the deflection of a beam on an elastic foundation with possible separation.
// Output is half of the beam length because symmetry is assumed

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <exception>
#include <cmath>

double C1, C2, C3, C4, C11, C13;

// Reset all array elements to zero
void reset_zero(double *in_Arry, int Num_elements){
	for (int i = 0; i < Num_elements; i++){
		in_Arry[i] = 0;
	}
}

// Deflection due to catilever condition at the separation point
double cant_defl(double dT, double beta, double gamma, double UW, double LCant, double E, double I, double localX, double t, double slope_in, int return_type){
	double Ra, Rb, Ma, Mb, Mt, thetaA, thetaB, yA, yB;

	Rb = UW*LCant;
	Mb = (-UW*pow(LCant, 2)) / 2;
	thetaA = (UW*pow(LCant, 3)) / (6 * E*I);
	yA = (-UW*pow(LCant, 3) * 3 * LCant) / (24 * E*I);
	Ma = 0.0;

	// Simply supported and guided
	//yA = (-UW * 5 * LCant*LCant*LCant*LCant) / (24 * E*I);
	//Ma = (UW*LCant*LCant) / 2;
	//thetaA = 0;


	thetaB = 0.0;
	yB = 0.0;
	Ra = 0.0;
	// Moment from temperature differential
	Mt = (-gamma*E*I*dT)/t;

	double y, slope, moment, shear;

	switch (return_type)
	{
	case 1:
		// Deflection due only to cantilever
		//y = yA + thetaA + thetaA*localX + (Ma*pow(localX, 2)) / (2 * E*I) + (Ra*pow(localX, 3)) / (6 * E*I) - (UW*pow(localX, 4)) / (24 * E*I);
		y = (-UW*(LCant - localX)*(LCant - localX)*(LCant - localX)*(LCant - localX)) / (24 * E*I) + (UW*LCant*(LCant - localX)*(LCant - localX)*(LCant - localX)) / (6 * E*I) + (-dT*gamma*(C1*C2 + C3*C4 - C2)*(LCant - localX)) / (beta*t*C11);
		// Add in deflection due to temperature differential
		//y = y + (gamma*dT*(LCant-localX)*(LCant-localX))/(2*t);
		return y;
		break;
	case 2:
		slope = (-UW*(LCant - localX)*(LCant - localX)*(LCant - localX)) / (6 * E*I) + (UW*LCant*(LCant - localX)*(LCant - localX)) / (2 * E*I) + (-dT*gamma*(C1*C2 + C3*C4 - C2)) / (beta*t*C11);
		return slope;
		break;
	case 3:
		moment = UW*(LCant*(LCant - localX) - ((LCant - localX)*(LCant - localX))/2);
		return moment;
		break;
	case 4:
		shear = UW*(-(LCant - localX)+LCant);
		return shear;
		break;
	default: return 9999999999.9;
	}
}

// Equations for free/fixed end cantilever solely due to temperature differential
double free_temp(double gamma, double E, double I, double t, double dT, double LCant, double localX, double slope_in, int return_type){
	double Ra, Ma, thetaA, yA, Rb, Mb, thetaB, yB, tempMo;

	Ra = 0.0;
	Ma = 0.0;
	Rb = 0.0;
	Mb = 0.0;
	thetaB = 0.0;
	yB = 0.0;

	thetaA = (-gamma*dT*LCant)/t + slope_in;
	yA = (gamma*dT*LCant*LCant) / (2 * t);

	double y, slope, moment, shear;

	switch (return_type)
	{
	case 1:
		y = yA + thetaA*localX + (Ma*localX*localX) / (2 * E*I) + (Ra*pow(localX, 3)) / (6 * E*I) + (gamma*dT*localX*localX) / (2 * t) - (yA + thetaA*LCant + (Ma*LCant*LCant) / (2 * E*I) + (Ra*pow(LCant, 3)) / (6 * E*I) + (gamma*dT*LCant*LCant) / (2 * t));
		return y;
		break;
	case 2:
		slope = thetaA + (Ma*localX) / (E*I) + (Ra*localX*localX) / (2 * E*I) + (gamma*dT*localX) / t;
		return slope;
		break;
	case 3:
		moment = Ma + Ra*localX;
		return moment;
		break;
	case 4:
		shear = Ra;
		return shear;
		break;
	default: return 9999999999.9;
	}
}

// Deflection due to uniform load and Winkler foundation
// Return value is dictated by return flag
double uni_defl(double beta, double UW, double C1, double C2, double C3, double C4, double C11, double E, double I, double x, int return_type){

	double Ra, Ma, thetaA, F1, F2, F3, F4, F5, yA;

	//Boundary conditions for both ends free
	Ra = 0;
	Ma = 0;
	thetaA = 0;
	yA = (UW*(C4*C2 - 2 * pow(C3, 2))) / (4 * E*I*pow(beta, 4)*C11);

	//Boundary condition for middle of beam free and separation point fixed
	//thetaA = 0.0;
	//yA = 0.0;
	//Ra = (UW*(C1*C2 + C4*C3)) / (beta*(2 + C11));
	//Ma = (UW*(2 * C1*C3 - C2*C2)) / (2 * beta*beta*(2 + C11));

	F1 = cosh(beta*x)*cos(beta*x);
	F2 = cosh(beta*x)*sin(beta*x) + sinh(beta*x)*cos(beta*x);
	F3 = sinh(beta*x)*sin(beta*x);
	F4 = cosh(beta*x)*sin(beta*x) - sinh(beta*x)*cos(beta*x);
	F5 = 1 - cosh(beta*x)*cos(beta*x);

	
	
	double y, slope, moment, shear;

	switch (return_type)
	{
	case 1:
		y = yA*F1 + (thetaA*F2) / (2 * beta) + (Ma*F3) / (2 * E*I*pow(beta, 2)) + (Ra*F4) / (4 * E*I*pow(beta, 3)) - (UW*F5) / (4 * E*I*pow(beta, 4));
		return y;
		break;
	case 2:
		slope = thetaA*F1 + (Ma*F2)/(2*E*I*beta) + (Ra*F3)/(2*E*I*beta*beta) - yA*beta*F4 - (UW*F4)/(4*E*I*pow(beta,3));
		return slope;
		break;
	case 3:
		moment = Ma*F1 + (Ra*F2)/(2*beta) - 2*yA*E*I*beta*beta*F3 - thetaA*E*I*beta*F4 - (UW*F3)/(2*beta*beta);
		return moment;
		break;
	case 4:
		shear = Ra*F1 - 2*yA*E*I*pow(beta,3)*F2 - 2*thetaA*E*I*beta*beta*F3 - Ma*beta*F4 - (UW*F2)/(2*beta);
		return shear;
		break;
	default: return 9999999999.9;
	}
}

// Deflection, slope, moment, and shear due to temperature differential
// Return value is dictated by return flag
double temp_defl(double dT, double beta, double t, double gamma, double E, double I, double C1, double C2, double C3, double C4, double C11, double x, int return_type){
	
	// Ra and Ma boundary conditions go to zero for the loading conditions
	// of both ends free
	double Ra, Ma, F1, F2, F3, F4, thetaA, yAT;

	dT = -dT;

	// Boundary for both free ends
	thetaA = (dT*gamma*(C1*C2 + C3*C4 - C2)) / (beta*t*C11);
	yAT = (-(dT)*gamma*(pow(C4, 2) + 2 * C1*C3 - 2 * C3)) / (2 * pow(beta, 2)*t*C11);
	Ra = 0;
	Ma = 0;

	// Boundary with middle of full beam free and separation point fixed
	//thetaA = 0.0;
	//yAT = 0.0;
	//Ra = (dT*gamma * 2 * beta*E*I*-C4) / (t*(2 + C11));
	//Ma = (dT*gamma*E*I*(2 * C1*C1 + C2*C4 - 2 * C1)) / (t*(2 + C11));

	F1 = cosh(beta*x)*cos(beta*x);
	F2 = cosh(beta*x)*sin(beta*x) + sinh(beta*x)*cos(beta*x);
	F3 = sinh(beta*x)*sin(beta*x);
	F4 = cosh(beta*x)*sin(beta*x) - sinh(beta*x)*cos(beta*x);

	double y, slope, moment, shear;
	
	
	switch (return_type)
	{
	case 1:
		y = yAT*F1 + (thetaA*F2) / (2 * beta) + (Ma*F3) / (2 * E*I*pow(beta, 2)) + (Ra*F4) / (4 * E*I*pow(beta, 3)) - ((dT)*gamma*F3) / (2 * t*pow(beta, 2));
		return y; 
		break;
	case 2:
		slope = thetaA*F1 + (Ma*F2) / (2 * E*I*beta) + (Ra*F3) / (2 * E*I*beta*beta) - yAT*beta*F4 - (dT * gamma*F2) / (2 * t*beta);
		return slope; 
		break;
	case 3: 
		moment = Ma*F1 + (Ra*F2) / (2 * beta) - 2 * yAT*E*I*beta*beta*F3 - thetaA*E*I*beta*F4 - (dT*gamma*E*I*(F1 - 1)) / t;
		return moment;
		break;
	case 4: 
		shear = Ra*F1 - 2 * yAT*E*I*pow(beta, 3)*F2 - 2 * thetaA*E*I*beta*beta*F3 - Ma*beta*F4 + (dT*gamma*E*I*beta*F4) / t;
		return shear;
		break;
	default: return 9999999999.9;
	}
}



// Calculate boundary condition factors for system
void C_values(double L, double beta){
	C1 = cosh(beta*L)*cos(beta*L);
	C2 = cosh(beta*L)*sin(beta*L) + sinh(beta*L)*cos(beta*L);
	C3 = sinh(beta*L)*sin(beta*L);
	C4 = cosh(beta*L)*sin(beta*L) - sinh(beta*L)*cos(beta*L);
	C11 = pow(sinh(beta*L), 2) - pow(sin(beta*L), 2);
	C13 = cosh(beta*L)*sinh(beta*L) - cos(beta*L)*sin(beta*L);
}

// Function that calculates a single deflection curve based on user inputs
void type1_Analysis(){
	std::string temp_diff, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue;
	std::string num_points, file_name;
	double dT, h, L, bo, E, UW, gamma, k, x;
	int N;

	std::cout << "Enter temperature differential (deg F): ";
	std::cin >> temp_diff;
	std::cout << temp_diff << " deg F\n";
	std::cout << "Enter beam thickness (in): ";
	std::cin >> beam_thick;
	std::cout << beam_thick << " in\n";
	std::cout << "Enter beam length (in): ";
	std::cin >> beam_leng;
	std::cout << beam_leng << " in\n";
	std::cout << "Enter beam width (in): ";
	std::cin >> beam_width;
	std::cout << beam_width << " in\n";
	std::cout << "Enter beam modulus of elasticity (psi): ";
	std::cin >> Em;
	std::cout << Em << " psi\n";
	std::cout << "Enter beam unit weight (lbs/ft3): ";
	std::cin >> rho;
	std::cout << rho << " lbs/ft3\n";
	std::cout << "Enter the beam coeff. of thermal expansion (1/F x10-5): ";
	std::cin >> cote;
	std::cout << cote << "E-05 /F\n";
	std::cout << "Enter subgrade k-value (pci): ";
	std::cin >> kvalue;
	std::cout << kvalue << " pci\n";
	std::cout << "How many points to generate for plot? ";
	std::cin >> num_points;
	std::cout << num_points << " points will be generated.\n";
	std::cout << "Enter file name for output of points: ";
	std::cin >> file_name;
	std::cout << "Data points will be outputted to " << file_name;

	try{
		dT = stod(temp_diff);
		h = stod(beam_thick);
		L = stod(beam_leng);
		bo = stod(beam_width);
		E = stod(Em);
		UW = stod(rho);
		//Convert UW to lbs/unit volume to have consistent units
		//The Roark derivation uses a density of cross sectional volume
		UW = UW *(1.0 / 1728.0) * (1.0*bo*h);
		gamma = stod(cote);
		//Convert gamma to 1E-5
		gamma = gamma*1E-5;
		k = stod(kvalue);
		N = stoi(num_points);
	}
	catch (const std::invalid_argument ia){
		std::cout << "\nOne of the values entered is incorrect.\nGeneral Error: INVALID ARGUMENT\nSpecific Error: " << ia.what() << "\nTHIS ERROR IS IRRECOVERABLE. PROGRAM TERMINATING...\n";
		exit(EXIT_FAILURE);
	}
	catch (...){
		std::cout << "\nGENERAL EXCEPTION THROWN. UNRECOVERABLE. PROGRAM TERMINATING...\n";
		exit(EXIT_FAILURE);
	}

	std::ofstream output_file1,output_file2,output_file3,output_file4;
	output_file1.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	output_file2.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	output_file3.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	output_file4.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try{
		output_file1.open(file_name+"1");
		output_file2.open(file_name + "2");
		output_file3.open(file_name + "3");
		output_file4.open(file_name + "4");
	}
	catch (std::ofstream::failure e){
		std::cout << "Failure reading creating file. Please make sure the file name does not conflict with an existing name. Also make sure you have write access to the location you are outputting data to.";
		exit(EXIT_FAILURE);
	}
	// Calculate moment of inertia
	double I;
	I = (bo*pow(h, 3)) / 12.0;

	// Calculate beta factor, constant for any location
	double beta;
	beta = pow((bo*k) / (4 * E*I), 1.0 / 4.0);

	// Write header for output file
	/*output_file << "Output file from FuBeam\n";
	output_file << "Temperature differential: " << temp_diff << " deg F\n";
	output_file << "Beam thickness: " << beam_thick << " in\n";
	output_file << "Beam length: " << beam_leng << " in\n";
	output_file << "Beam width: " << beam_width << " in\n";
	output_file << "Elastic Modulus: " << Em << " psi\n";
	output_file << "Unit weight: " << rho << " lbs/ft3\n";
	output_file << "Coeff. of thermal expansion: " << cote << "E-05 /F\n";
	output_file << "Subgrade modulus: " << kvalue << " pci\n";
	output_file << num_points << " points will be generated.\n";*/

	// Try to come up with a close initial guess based simply on the foundational curling
	double L_defl;
	int sep_point_guess = N;
	//C_values(L, beta);
	bool keepGoing = true;
	for (int i = 1; i < N; i++){
		if (keepGoing){
			x = ((double)i / (double)N)*L / 2;
			C_values(2 * x, beta);
			L_defl = temp_defl(dT, beta, h, gamma, E, I, C1, C2, C3, C4, C11, 0, 1) + uni_defl(beta, UW, C1, C2, C3, C4, C11, E, I, 0, 1);
			if (L_defl >= 0){
				sep_point_guess = i;
				keepGoing = false;
			}
		}
	}

	double *ElasticTempAry = new (std::nothrow)double[N];
	double *ElasticUniAry = new (std::nothrow)double[N];
	double *cantAry = new (std::nothrow)double[N];
	double *freeTempAry = new (std::nothrow)double[N];
	double *xAry = new (std::nothrow)double[N];

	double calc_x, cant_x, LCant;
	if (ElasticTempAry == nullptr || ElasticUniAry == nullptr || xAry == nullptr || cantAry == nullptr || freeTempAry == nullptr){
		std::cout << "Memory error in a deflection array initialization. Try again and/or ensure there is sufficient free space for the analysis.";
	}
	else{
		for (int sel = 1; sel <= 4; sel++){
			double sl_in = 0.0;
			for (int i = 0; i < N; i++){
				x = ((double)i / (double)N)*L / 2;
				calc_x = ((double)sep_point_guess / (double)N)*L / 2 - ((double)i / (double)N)*L / 2;
				xAry[i] = x;
				ElasticTempAry[i] = temp_defl(dT, beta, h, gamma, E, I, C1, C2, C3, C4, C11, calc_x, sel);
				ElasticUniAry[i] = uni_defl(beta, UW, C1, C2, C3, C4, C11, E, I, calc_x, sel);
				cantAry[i] = 0.0;
				freeTempAry[i] = 0.0;
				if (i >= sep_point_guess){
					LCant = L / 2 - ((double)sep_point_guess / (double)N)*L / 2;
					cant_x = LCant - (((double)sep_point_guess / (double)N)*L / 2 - ((double)i / (double)N)*L / 2);
					sl_in = temp_defl(dT, beta, h, gamma, E, I, C1, C2, C3, C4, C11, 0, 2);
					cantAry[i] = cant_defl(dT, beta, gamma, UW, LCant, E, I, cant_x, h, sl_in, sel);
					ElasticUniAry[i] = 0;
					ElasticTempAry[i] = 0;
				}
			}
			switch (sel){
			case 1:
				output_file1 << "Point of Separation: " << ((double)sep_point_guess / (double)N)*L / 2 << " inches from middle of beam.\n";
				output_file1 << "Dist. From Middle [in]\t" << "Temp. Deflection [in]\t" << "Load Deflection [in]\t" << "Cantilever Deflection (Load) [in]\t" << "Cantilever Deflection (Temp.) [in]\t" << "Total Deflection\n";

				for (int j = 0; j < N; j++){
					output_file1 << xAry[j] << "\t";
					output_file1 << ElasticTempAry[j] << "\t";
					output_file1 << ElasticUniAry[j] << "\t";
					output_file1 << cantAry[j] << "\t";
					output_file1 << ElasticTempAry[j] + ElasticUniAry[j] + cantAry[j] << "\n";
				}
			case 2:
				output_file2 << "Point of Separation: " << ((double)sep_point_guess / (double)N)*L / 2 << " inches from middle of beam.\n";
				output_file2 << "Dist. From Middle [in]\t" << "Temp. Deflection [in]\t" << "Load Deflection [in]\t" << "Cantilever Deflection (Load) [in]\t" << "Cantilever Deflection (Temp.) [in]\t" << "Total Deflection\n";

				for (int j = 0; j < N; j++){
					output_file2 << xAry[j] << "\t";
					output_file2 << ElasticTempAry[j] << "\t";
					output_file2 << ElasticUniAry[j] << "\t";
					output_file2 << cantAry[j] << "\t";
					output_file2 << ElasticTempAry[j] + ElasticUniAry[j] + cantAry[j] << "\n";
				}
			case 3:
				output_file3 << "Point of Separation: " << ((double)sep_point_guess / (double)N)*L / 2 << " inches from middle of beam.\n";
				output_file3 << "Dist. From Middle [in]\t" << "Temp. Deflection [in]\t" << "Load Deflection [in]\t" << "Cantilever Deflection (Load) [in]\t" << "Cantilever Deflection (Temp.) [in]\t" << "Total Deflection\n";

				for (int j = 0; j < N; j++){
					output_file3 << xAry[j] << "\t";
					output_file3 << ElasticTempAry[j] << "\t";
					output_file3 << ElasticUniAry[j] << "\t";
					output_file3 << cantAry[j] << "\t";
					output_file3 << ElasticTempAry[j] + ElasticUniAry[j] + cantAry[j] << "\n";
				}
			case 4:
				output_file4 << "Point of Separation: " << ((double)sep_point_guess / (double)N)*L / 2 << " inches from middle of beam.\n";
				output_file4 << "Dist. From Middle [in]\t" << "Temp. Deflection [in]\t" << "Load Deflection [in]\t" << "Cantilever Deflection (Load) [in]\t" << "Cantilever Deflection (Temp.) [in]\t" << "Total Deflection\n";

				for (int j = 0; j < N; j++){
					output_file4 << xAry[j] << "\t";
					output_file4 << ElasticTempAry[j] << "\t";
					output_file4 << ElasticUniAry[j] << "\t";
					output_file4 << cantAry[j] << "\t";
					output_file4 << ElasticTempAry[j] + ElasticUniAry[j] + cantAry[j] << "\n";
				}
			default:
				std::cout << "Something went wrong";
			}
		}
	}

	delete[] xAry;
	delete[] ElasticTempAry;
	delete[] ElasticUniAry;
	delete[] cantAry;
	delete[] freeTempAry;
}

int main(){
	std::string typeA_in;
	int typeAnalysis;


	std::cout << "What type of analysis?\n1) Single Curve\n2) Iterative Curve Match\n";
	std::cin >> typeA_in;
	typeAnalysis = stoi(typeA_in);

	if (typeAnalysis == 1){
		type1_Analysis();
	}
	else{
		std::cout << "Invalid analysis method option selected. Program terminating...";
	}
	
}