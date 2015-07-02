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

// Function that calculates a single deflection curve based on user inputs
void type1_Analysis(){
	std::string temp_diff, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue;
	std::string num_points, file_name;
	double dT, h, L, bo, E, UW, gamma, k;
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
	//std::cout << "How many points to generate for plot? ";
	//std::cin >> num_points;
	//std::cout << num_points << " points will be generated.\n";
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
		UW = UW *(1.0 / 1728.0) * (1.0*bo*h);
		gamma = stod(cote);
		//Convert gamma to 1E-5
		gamma = gamma*1E-5;
		k = stod(kvalue);

		// Points on each side of separation point
		N = 300;// stoi(num_points);
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
		output_file1.open(file_name + "-Deflection.txt");
		output_file2.open(file_name + "-Slope.txt");
		output_file3.open(file_name + "-Moment.txt");
		output_file4.open(file_name + "-Shear.txt");
	}
	catch (std::ofstream::failure e){
		std::cout << "Failure reading creating file. Please make sure the file name does not conflict with an existing name. Also make sure you have write access to the location you are outputting data to.";
		exit(EXIT_FAILURE);
	}
	// Calculate moment of inertia
	double I;
	I = (bo*pow(h, 3)) / 12.0;

	// Calculate lambda factor, constant for any location
	double lambda;
	lambda = pow((bo*k) / (4 * E*I), 1.0 / 4.0);

	// Calculate temperature moment
	double Mt;
	Mt = (-gamma*E*I*dT) / h;

	// Use Newton-Raphson to find the root. Basically, you're trying to find a length of beam that
	// when placed on the elastic foundation has exactly zero deflection at the ends. Newton-Raphson 
	// is used since the derivative of the deflection function is known analytically.
	
	double tolerance = 1.0E-10; // difference between two predictions
	double epsilon = 1.0E-15; // smallest number to use in a calculation
	int max_iterations = 1000;
	bool solution_found = false;

	double old_x, new_x, elastic_beam_length;
	double fx, fprimeX;
	old_x = 4.0; //This should be reasonably far enough away from 0 so root analysis doesn't fail
	for (int i = 1; i <= max_iterations; i++){
		fx = ((-2 * Mt*lambda*lambda) / (h*k))*((sinh(lambda*old_x) - sin(lambda*old_x)) / (sinh(lambda*old_x) + sin(lambda*old_x))) - (UW*(L/old_x)) / (h*k);
		fprimeX = ((-4 * Mt*lambda*lambda*lambda) / (h*k))*((cosh(lambda*old_x) - cos(lambda*old_x)) / (sinh(lambda*old_x) + sin(lambda*old_x)));
		
		if (abs(fprimeX) < epsilon){
			std::cout << "Epsilon threshold exceeded during Newton-Raphson.";
			break;
		}
		new_x = old_x - fx / fprimeX;

		if (abs(old_x - new_x) / abs(new_x) < tolerance){
			solution_found = true;
			elastic_beam_length = new_x;
			break;
		}

		old_x = new_x;
	}

	// Write headers for output files
	// Should probably be made into functions...
	output_file1 << "Output file from FuBeam\n";
	output_file1 << "Temperature differential: " << temp_diff << " deg F\n";
	output_file1 << "Beam thickness: " << beam_thick << " in\n";
	output_file1 << "Beam length: " << beam_leng << " in\n";
	output_file1 << "Beam width: " << beam_width << " in\n";
	output_file1 << "Elastic Modulus: " << Em << " psi\n";
	output_file1 << "Unit weight: " << rho << " lbs/ft3\n";
	output_file1 << "Coeff. of thermal expansion: " << cote << "E-05 /F\n";
	output_file1 << "Subgrade modulus: " << kvalue << " pci\n";
	output_file1 << "Point of Separation: " << elastic_beam_length/2 << " inches from middle of beam.\n";
	output_file1 << "Dist. From Middle [in]\t" << "El. Temp. Deflection [in]\t" << "El. Self-Weight Deflection [in]\t" << "Cantilever Deflection (Temp. + Self-Weight) [in]\t" << "Total Deflection [in]\n";

	output_file2 << "Output file from FuBeam\n";
	output_file2 << "Temperature differential: " << temp_diff << " deg F\n";
	output_file2 << "Beam thickness: " << beam_thick << " in\n";
	output_file2 << "Beam length: " << beam_leng << " in\n";
	output_file2 << "Beam width: " << beam_width << " in\n";
	output_file2 << "Elastic Modulus: " << Em << " psi\n";
	output_file2 << "Unit weight: " << rho << " lbs/ft3\n";
	output_file2 << "Coeff. of thermal expansion: " << cote << "E-05 /F\n";
	output_file2 << "Subgrade modulus: " << kvalue << " pci\n";
	output_file2 << "Point of Separation: " << elastic_beam_length/2 << " inches from middle of beam.\n";
	output_file2 << "Dist. From Middle [in]\t" << "El. Temp. Slope [rad]\t" << "El. Self-Weight Slope [rad]\t" << "Cantilever Slope (Temp. + Self-Weight) [rad]\t" << "Total Slope [rad]\n";

	output_file3 << "Output file from FuBeam\n";
	output_file3 << "Temperature differential: " << temp_diff << " deg F\n";
	output_file3 << "Beam thickness: " << beam_thick << " in\n";
	output_file3 << "Beam length: " << beam_leng << " in\n";
	output_file3 << "Beam width: " << beam_width << " in\n";
	output_file3 << "Elastic Modulus: " << Em << " psi\n";
	output_file3 << "Unit weight: " << rho << " lbs/ft3\n";
	output_file3 << "Coeff. of thermal expansion: " << cote << "E-05 /F\n";
	output_file3 << "Subgrade modulus: " << kvalue << " pci\n";
	output_file3 << "Point of Separation: " << elastic_beam_length/2 << " inches from middle of beam.\n";
	output_file3 << "Dist. From Middle [in]\t" << "El. Temp. Moment [lb*in]\t" << "El. Self-Weight Moment [lb*in]\t" << "Cantilever Moment (Temp. + Self-Weight) [lb*in]\t" << "Total Moment [lb*in]\n";

	output_file4 << "Output file from FuBeam\n";
	output_file4 << "Temperature differential: " << temp_diff << " deg F\n";
	output_file4 << "Beam thickness: " << beam_thick << " in\n";
	output_file4 << "Beam length: " << beam_leng << " in\n";
	output_file4 << "Beam width: " << beam_width << " in\n";
	output_file4 << "Elastic Modulus: " << Em << " psi\n";
	output_file4 << "Unit weight: " << rho << " lbs/ft3\n";
	output_file4 << "Coeff. of thermal expansion: " << cote << "E-05 /F\n";
	output_file4 << "Subgrade modulus: " << kvalue << " pci\n";
	output_file4 << "Point of Separation: " << elastic_beam_length/2 << " inches from middle of beam.\n";
	output_file4 << "Dist. From Middle [in]\t" << "El. Temp. Shear [lb]\t" << "El. Self-Weight Shear [lb]\t" << "Cantilever Shear (Temp. + Self-Weight) [lb]\t" << "Total Shear [lb]\n";

	// Calculate portion of beam that is on elastic foundation
	double xs; // x-coordinate of the system
	double x; // since we're interested in only half the beam, we do this, follows Hetenyi's notation
	double xp; // this is the remaining portion of the beam, see Hetenyi notation
	double t_defl, w_defl, t_slp, t_mom, t_shr;
	for (int j = 0; j < N; j++){
		xs = ((double)j / (double)N) * elastic_beam_length / 2;
		x = elastic_beam_length/2 + xs;
		xp = elastic_beam_length - x;
		
		t_defl = -((2 * Mt*lambda*lambda) / (h*k*(sinh(lambda*elastic_beam_length) + sin(lambda*elastic_beam_length))))*(sinh(lambda*x)*cos(lambda*xp) - cosh(lambda*x)*sin(lambda*xp) + sinh(lambda*xp)*cos(lambda*x) - cosh(lambda*xp)*sin(lambda*x));
		w_defl = -(UW*(L/elastic_beam_length)) / (h*k);

		t_slp = (-(4 * Mt*lambda*lambda*lambda)*(cosh(lambda*x)*cos(lambda*xp) - cosh(lambda*xp)*cos(lambda*x))) / (h*k*(sinh(lambda*elastic_beam_length) + sin(lambda*elastic_beam_length)));

		t_mom = (Mt / (sinh(lambda*elastic_beam_length) + sin(lambda*elastic_beam_length)))*(sinh(lambda*x)*cos(lambda*xp) + cosh(lambda*x)*sin(lambda*xp) + sinh(lambda*xp)*cos(lambda*x) + cosh(lambda*xp)*sin(lambda*x));

		t_shr = (2 * Mt*lambda)*((sinh(lambda*x)*sin(lambda*xp) - sinh(lambda*xp)*sin(lambda*x)) / (sinh(lambda*elastic_beam_length) + sin(lambda*elastic_beam_length)));

		output_file1 << xs << "\t" << t_defl << "\t" << w_defl << "\t" << "0.000\t" << t_defl + w_defl << "\n";
		output_file2 << xs << "\t" << t_slp << "\t" << "0.000" << "\t" << "0.000\t" << t_slp << "\n";
		output_file3 << xs << "\t" << t_mom << "\t" << "0.000" << "\t" << "0.000\t" << t_mom << "\n";
		output_file4 << xs << "\t" << t_shr << "\t" << "0.000" << "\t" << "0.000\t" << t_shr << "\n";
	}

	// Calculate cantilevered portion of beam
	double LCant = L / 2 - elastic_beam_length / 2;
	double prev_xs = xs;
	double eqslp; // this is the equivalent slope boundary condition
	UW = -1 * UW;
	eqslp = ((4 * Mt*lambda*lambda*lambda) / (h*k))*((cosh(lambda*elastic_beam_length) - cos(lambda*elastic_beam_length)) / (sinh(lambda*elastic_beam_length) + sin(lambda*elastic_beam_length)));
	double c_defl, c_slp, c_mom, c_shr;
	for (int k = 0; k <= N; k++){
		xs = prev_xs + ((double)k / (double)N)*LCant;
		x = ((double)k / (double)N)*LCant;

		c_defl = (-(UW*x*x*x*x) / (24 * E*I) + (UW*LCant*x*x*x) / (6 * E*I) - (gamma*dT*x*x) / (2 * h) + eqslp*x)*-1;
		c_slp = (-(UW*x*x*x) / (6 * E*I) + (UW*LCant*x*x) / (2 * E*I) - (gamma*dT*x) / (h) + eqslp)*-1;
		c_mom = -(UW*x*x) / 2 + UW*LCant*x + Mt;
		c_shr = -UW*x + UW*LCant;

		output_file1 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << c_defl << "\t" << c_defl << "\n";
		output_file2 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << c_slp << "\t" << c_slp << "\n";
		output_file3 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << c_mom << "\t" << c_mom << "\n";
		output_file4 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << c_shr << "\t" << c_shr << "\n";
	}
}

int main(){
	type1_Analysis();	
}