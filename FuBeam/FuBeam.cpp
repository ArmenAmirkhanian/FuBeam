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

// Reset all array elements to zero
void reset_zero(double *in_Arry, int Num_elements){
	for (int i = 0; i < Num_elements; i++){
		in_Arry[i] = 0;
	}
}

// Creep calculation based on a modified B3 model
double creep_calc(){

}


// Deflection due to catilever condition at the separation point
double cant_defl(double UW, double LCant, double E, double I, double localX){
	double Ra, Rb, Ma, Mb, thetaA, thetaB, yA, yB;

	Rb = UW*LCant;
	Mb = (-UW*pow(LCant, 2)) / 2;
	thetaA = (UW*pow(LCant, 3)) / (6 * E*I);
	yA = (-UW*pow(LCant, 3) * 3 * LCant) / (24 * E*I);

	thetaB = 0;
	yB = 0;
	Ra = 0;
	Ma = 0;

	double y;
	y = yA + thetaA + thetaA*localX + (Ma*pow(localX, 2)) / (2 * E*I) + (Ra*pow(localX, 3)) / (6 * E*I) - (UW*pow(localX, 4)) / (24 * E*I);
	return y;
}

// Deflection due to uniform load and Winkler foundation
double uni_defl(double beta, double UW, double L, double C2, double C3, double C4, double C11, double E, double I, double x){
	// Ra, Ma, and thetaA boundary conditions go to zero for loading conditions
	// The entire equation for uniform load deflection is included in case
	// there is a desire to model different conditions in the future
	double Ra, Ma, thetaA, F1, F2, F3, F4, F5, yA;

	Ra = 0;
	Ma = 0;
	thetaA = 0;

	F1 = cosh(beta*x)*cos(beta*x);
	F2 = cosh(beta*x)*sin(beta*x) + sinh(beta*x)*cos(beta*x);
	F3 = sinh(beta*x)*sin(beta*x);
	F4 = cosh(beta*x)*sin(beta*x) - sinh(beta*x)*cos(beta*x);
	F5 = 1 - cosh(beta*x)*cos(beta*x);

	yA = (UW*(C4*C2 - 2 * pow(C3, 2))) / (4 * E*I*pow(beta, 4)*C11);

	double y;
	y = yA*F1 + (thetaA*F2) / (2 * beta) + (Ma*F3) / (2 * E*I*pow(beta, 2)) + (Ra*F4) / (4 * E*I*pow(beta, 3)) - (UW*F5) / (4 * E*I*pow(beta, 4));
	return y;
}

// Deflection due to temperature differential
double temp_defl(double dT, double beta, double t, double gamma, double E, double I, double C1, double C2, double C3, double C4, double C11, double x){
	
	// Ra and Ma boundary conditions go to zero for the loading conditions
	// The entire equation of temperature deflection is included in case
	// there is a desire to model different conditions in the future
	double Ra, Ma, F1, F2, F3, F4, thetaA, yAT;

	Ra = 0;
	Ma = 0;

	thetaA = (-(dT)*gamma*(C1*C2 + C3*C4 - C2)) / (beta*t*C11);
	yAT = (-(dT)*gamma*(pow(C4, 2) + 2 * C1*C3 - 2 * C3)) / (2 * pow(beta, 2)*t*C11);

	F1 = cosh(beta*x)*cos(beta*x);
	F2 = cosh(beta*x)*sin(beta*x) + sinh(beta*x)*cos(beta*x);
	F3 = sinh(beta*x)*sin(beta*x);
	F4 = cosh(beta*x)*sin(beta*x) - sinh(beta*x)*cos(beta*x);

	double y;
	y = yAT*F1 + (thetaA*F2) / (2 * beta) + (Ma*F3) / (2 * E*I*pow(beta, 2)) + (Ra*F4) / (4 * E*I*pow(beta, 3)) - ((dT)*gamma*F3) / (2 * t*pow(beta, 2));
	return y;
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

	std::ofstream output_file;
	output_file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try{
		output_file.open(file_name);
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

	// Calculate boundary condition factors
	double C1, C2, C3, C4, C11, C13;
	C1 = cosh(beta*L)*cos(beta*L);
	C2 = cosh(beta*L)*sin(beta*L) + sinh(beta*L)*cos(beta*L);
	C3 = sinh(beta*L)*sin(beta*L);
	C4 = cosh(beta*L)*sin(beta*L) - sinh(beta*L)*cos(beta*L);
	C11 = pow(sinh(beta*L), 2) - pow(sin(beta*L), 2);
	C13 = cosh(beta*L)*sinh(beta*L) - cos(beta*L)*sin(beta*L);

	// Write header for output file
	output_file << "Output file from FuBeam\n";
	output_file << "Temperature differential: " << temp_diff << " deg F\n";
	output_file << "Beam thickness: " << beam_thick << " in\n";
	output_file << "Beam length: " << beam_leng << " in\n";
	output_file << "Beam width: " << beam_width << " in\n";
	output_file << "Elastic Modulus: " << Em << " psi\n";
	output_file << "Unit weight: " << rho << " lbs/ft3\n";
	output_file << "Coeff. of thermal expansion: " << cote << "E-05 /F\n";
	output_file << "Subgrade modulus: " << kvalue << " pci\n";
	output_file << num_points << " points will be generated.\n";
	

	double *deflTempAry = new (std::nothrow)double[N];
	double *deflUniAry = new (std::nothrow)double[N];
	double *cantAry = new (std::nothrow)double[N];
	double *xAry = new (std::nothrow)double[N];
	double sep_point, tempX, LCant, offset_defl;
	int sep_flag = 0;

	if (deflTempAry == nullptr || deflUniAry == nullptr || xAry == nullptr || cantAry == nullptr){
		output_file << "Memory error in a deflection array initialization. Try again and/or ensure there is sufficient free space for the analysis.";
	}
	else{
		for (int i = 0; i < N; i++){
			x = ((double)i / (double)N)*L / 2;
			xAry[i] = x;
			deflTempAry[i] = temp_defl(dT, beta, h, gamma, E, I, C1, C2, C3, C4, C11, x);
			deflUniAry[i] = uni_defl(beta, UW, L, C2, C3, C4, C11, E, I, x);
			cantAry[i] = 0;
			if (deflTempAry[i] + deflUniAry[i] > 0.0 && sep_flag == 0){
				sep_point = x;
				sep_flag = 1;
				LCant = L - sep_point;
				offset_defl = deflTempAry[i];
			}
			if (sep_flag == 1){
				tempX = x - sep_point;
				cantAry[i] = cant_defl(UW, LCant, E, I, tempX) - offset_defl;
				deflUniAry[i] = 0;
			}
		}
	}

	output_file << "Point of Separation: " << sep_point << " inches from middle of beam.\n";
	output_file << "Dist. From Middle [in]\t" << "Temp. Deflection [in]\t" << "Load Deflection [in]\t" << "Cantilever Deflection [in]\t" << "Total Deflection\n";

	for (int i = 0; i < N; i++){
		output_file << xAry[i] << "\t";
		output_file << deflTempAry[i] << "\t";
		output_file << deflUniAry[i] << "\t";
		output_file << cantAry[i] << "\t";
		output_file << deflTempAry[i] + deflUniAry[i] + cantAry[i] << "\n";
	}

	delete[] xAry;
	delete[] deflTempAry;
	delete[] deflUniAry;
	delete[] cantAry;
}

// Iterative deflection analysis to match actual beam deflection with equivalent temperature deflection
// Performs calculations in 0.1 deg F increments and uses RMSE to find best fit
// Analysis will calculate over entire 50 deg range then find smallest RMSE
void type2_Analysis(){
	std::string radius, cen_defl, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue;
	std::string num_points, file_name;
	double dT, R, CD, h, L, bo, E, UW, gamma, k, x;
	int N;

	std::cout << "Enter radius of curvature of actual beam (in): ";
	std::cin >> radius;
	std::cout << radius << " in\n";
	std::cout << "Enter deflection at center of actual beam (in): ";
	std::cin >> cen_defl;
	std::cout << cen_defl << " in\n";
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
		R = stod(radius);
		CD = stod(cen_defl);
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

	std::ofstream output_file;
	output_file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try{
		output_file.open(file_name);
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

	// Calculate deflection profile of actual beam
	double *ActDefl = new (std::nothrow)double[N];
	double xR;
	if (ActDefl == nullptr){
		output_file << "Memory error in a deflection array initialization. Try again and/or ensure there is sufficient free space for the analysis.";
	}
	else{
		for (int i = 0; i < N; i++){
			xR = ((double)i / (double)N)*L / 2;
			ActDefl[i] = -(R*sin(acos(xR / R)) - CD - R);
		}
	}
	

	// Calculate boundary condition factors
	double C1, C2, C3, C4, C11, C13;
	C1 = cosh(beta*L)*cos(beta*L);
	C2 = cosh(beta*L)*sin(beta*L) + sinh(beta*L)*cos(beta*L);
	C3 = sinh(beta*L)*sin(beta*L);
	C4 = cosh(beta*L)*sin(beta*L) - sinh(beta*L)*cos(beta*L);
	C11 = pow(sinh(beta*L), 2) - pow(sin(beta*L), 2);
	C13 = cosh(beta*L)*sinh(beta*L) - cos(beta*L)*sin(beta*L);



	// Write header for output file
	output_file << "Output file from FuBeam\n";
	output_file << "Radius of curvature from actual beam: " << radius << " in\n";
	output_file << "Beam thickness: " << beam_thick << " in\n";
	output_file << "Beam length: " << beam_leng << " in\n";
	output_file << "Beam width: " << beam_width << " in\n";
	output_file << "Elastic Modulus: " << Em << " psi\n";
	output_file << "Unit weight: " << rho << " lbs/ft3\n";
	output_file << "Coeff. of thermal expansion: " << cote << "E-05 /F\n";
	output_file << "Subgrade modulus: " << kvalue << " pci\n";
	output_file << num_points << " points will be generated.\n";
	

	double *deflTempAry = new (std::nothrow)double[N];
	double *deflUniAry = new (std::nothrow)double[N];
	double *cantAry = new (std::nothrow)double[N];
	double *xAry = new (std::nothrow)double[N];
	double *RMSE = new (std::nothrow)double[501];
	// Make sure the first value of RMSE is not smallest
	RMSE[0] = 999999;
	double sep_point, tempX, LCant, offset_defl, diff = 0.0;
	int sep_flag = 0;
	int Tf = 1;

	if (deflTempAry == nullptr || deflUniAry == nullptr || xAry == nullptr || cantAry == nullptr || RMSE == nullptr){
		output_file << "Memory error in a deflection array initialization. Try again and/or ensure there is sufficient free space for the analysis.";
	}
	else{
		for (int T = 1; T <= 500; T++)
		{
			reset_zero(deflTempAry, N);
			reset_zero(deflUniAry, N);
			reset_zero(cantAry, N);

			dT = -(double)T*0.1;
			for (int i = 0; i < N; i++){
				x = ((double)i / (double)N)*L / 2;
				xAry[i] = x;
				deflTempAry[i] = temp_defl(dT, beta, h, gamma, E, I, C1, C2, C3, C4, C11, x);
				deflUniAry[i] = uni_defl(beta, UW, L, C2, C3, C4, C11, E, I, x);
				cantAry[i] = 0;
				if (deflTempAry[i] + deflUniAry[i] > 0.0 && sep_flag == 0){
					sep_point = x;
					sep_flag = 1;
					LCant = L - sep_point;
					offset_defl = deflTempAry[i];
				}
				if (sep_flag == 1){
					tempX = x - sep_point;
					cantAry[i] = cant_defl(UW, LCant, E, I, tempX) - offset_defl;
					deflUniAry[i] = 0;
				}
			}
			diff = 0.0;
			for (int j = 0; j < N; j++){
				diff += pow((deflTempAry[j] + deflUniAry[j] + cantAry[j]) - ActDefl[j],2);
			}
			RMSE[T] = sqrt(diff / N);
		}

		double smallestRMSE = RMSE[0];
		for (int j = 1; j <= 500; j++){
			if (smallestRMSE < RMSE[j]){}
			else{ smallestRMSE = RMSE[j]; Tf = j; }
		}
		// Perform final deflection calculation based on point with lowest RMSE
		dT = -(double)Tf*0.1;
		for (int i = 0; i < N; i++){
			x = ((double)i / (double)N)*L / 2;
			xAry[i] = x;
			deflTempAry[i] = temp_defl(dT, beta, h, gamma, E, I, C1, C2, C3, C4, C11, x);
			deflUniAry[i] = uni_defl(beta, UW, L, C2, C3, C4, C11, E, I, x);
			cantAry[i] = 0;
			if (deflTempAry[i] + deflUniAry[i] > 0.0 && sep_flag == 0){
				sep_point = x;
				sep_flag = 1;
				LCant = L - sep_point;
				offset_defl = deflTempAry[i];
			}
			if (sep_flag == 1){
				tempX = x - sep_point;
				cantAry[i] = cant_defl(UW, LCant, E, I, tempX) - offset_defl;
				deflUniAry[i] = 0;
			}
		}
	}

	output_file << "Final temperature differential: -" << (double)Tf*0.1 << " deg F\n";
	output_file << "RMSE of fit: " << RMSE[Tf] << " inches\n";
	output_file << "Point of Separation: " << sep_point << " inches from middle of beam.\n";
	output_file << "Dist. From Middle [in]\t" << "Temp. Deflection [in]\t" << "Load Deflection [in]\t" << "Cantilever Deflection [in]\t" << "Total Deflection [in]\t" << "Actual Deflection [in]\n";

	for (int i = 0; i < N; i++){
		output_file << xAry[i] << "\t";
		output_file << deflTempAry[i] << "\t";
		output_file << deflUniAry[i] << "\t";
		output_file << cantAry[i] << "\t";
		output_file << deflTempAry[i] + deflUniAry[i] + cantAry[i] << "\t";
		output_file << ActDefl[i] << "\n";
	}

	delete[] xAry;
	delete[] deflTempAry;
	delete[] deflUniAry;
	delete[] cantAry;
	delete[] RMSE;
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
	else if (typeAnalysis == 2){
		type2_Analysis();
	}
	else{
		std::cout << "Invalid analysis method option selected. Program terminating...";
	}
	
}