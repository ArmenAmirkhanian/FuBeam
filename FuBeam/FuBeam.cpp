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

// Initialize flags
bool solution_found = false;
bool max_iter_hit = false;

void write_header(std::ofstream& output_file_stream, std::string temp_diff, std::string beam_thick, std::string beam_leng, std::string beam_width, std::string Em, std::string rho, std::string cote, std::string kvalue, double elastic_beam_length){
	output_file_stream << "Output file from FuBeam\n";
	output_file_stream << "Temperature differential or moment: " << temp_diff << " deg F or in*lb\n";
	output_file_stream << "Beam thickness: " << beam_thick << " in\n";
	output_file_stream << "Beam length: " << beam_leng << " in\n";
	output_file_stream << "Beam width: " << beam_width << " in\n";
	output_file_stream << "Elastic Modulus: " << Em << " psi\n";
	output_file_stream << "Unit weight: " << rho << " lbs/ft3\n";
	output_file_stream << "Coeff. of thermal expansion: " << cote << "E-05 /F\n";
	output_file_stream << "Subgrade modulus: " << kvalue << " pci\n";
	output_file_stream << "Point of Separation: " << elastic_beam_length / 2 << " inches from middle of beam.\n";
	output_file_stream << "Dist. From Middle [in]\t" << "El. Temp. Deflection [in]\t" << "El. Self-Weight Deflection [in]\t" << "Cantilever Deflection (Temp. + Self-Weight) [in]\t" << "Total Deflection [in]\n";
}

double newt_rap(double UW, double lambda, double Mt, double L, double k){
	// Use Newton-Raphson to find the root. Basically, you're trying to find a length of beam that
	// when placed on the elastic foundation has exactly zero deflection at the ends. Newton-Raphson 
	// is used since the derivative of the deflection function is known analytically and is smooth.
	double tolerance = 1.0E-10; // difference between two predictions
	double epsilon = 1.0E-15; // smallest number to use in a calculation, roughly one order of magnitude larger than machine epsilon for type double
	int max_iterations = 1000; // this is probably overkill since both f and f' are smooth

	double old_x, new_x, elastic_beam_length;
	double y_dis_sup, y_dis_sup_prime, y_mom_sup, y_mom_sup_prime;
	double fx, fprimeX, MT, ls;
	
	ls = L / 2;
	old_x = ls*0.2; //This should be reasonably far enough away from 0 so root analysis doesn't fail

	for (int i = 1; i <= max_iterations; i++){
		MT = -(UW*((ls - old_x)*(ls - old_x))) / 2 + Mt;
		y_dis_sup = -UW/k;
		y_dis_sup_prime = 0;
		y_mom_sup = -((2*MT*lambda*lambda)/k)*((sinh(lambda*old_x*2)-sin(lambda*old_x*2))/(sinh(lambda*old_x*2)+sin(lambda*old_x*2)));
		y_mom_sup_prime = -((4*MT*lambda*lambda*lambda)/(k))*((cosh(lambda*old_x*2)-cos(lambda*old_x*2))/(sinh(lambda*old_x*2)+sin(lambda*old_x*2)));
		fx = y_dis_sup + y_mom_sup;
		fprimeX = y_dis_sup_prime + y_mom_sup_prime;

		if (abs(fprimeX) < epsilon){
			std::cout << "Epsilon threshold exceeded during Newton-Raphson.";
			elastic_beam_length = -1;
			break;
		}
		new_x = old_x - fx / fprimeX;

		if (abs(old_x - new_x) / abs(new_x) < tolerance){
			solution_found = true;
			if (new_x > ls){
				elastic_beam_length = ls;
			}
			else{
				elastic_beam_length = new_x;
			}
			
			break;
		}

		if (i == max_iterations){
			max_iter_hit = true;
			elastic_beam_length = -1;
		}

		old_x = new_x;
	}

	return elastic_beam_length;
}

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
		// 500 points should be more than adequate
		N = 500;
	}
	catch (const std::invalid_argument ia){
		std::cout << "\nOne of the values entered is incorrect.\nGeneral Error: INVALID ARGUMENT\nSpecific Error: " << ia.what() << "\nTHIS ERROR IS IRRECOVERABLE. PROGRAM TERMINATING...\n";
		exit(EXIT_FAILURE);
	}
	catch (...){
		std::cout << "\nGENERAL EXCEPTION THROWN. UNRECOVERABLE. PROGRAM TERMINATING...\n";
		exit(EXIT_FAILURE);
	}

	std::ofstream output_file1, output_file2, output_file3, output_file4;
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

	// Find separation point
	double elastic_beam_length;
	elastic_beam_length = newt_rap(UW, lambda, Mt, L, k);

	// Calculate prescribed moment on elastic foundation
	// This value is the continuity between the elastic foundation
	// and the cantilevered section
	double MT;
	MT = (-UW*((L/2 - elastic_beam_length)*(L/2 - elastic_beam_length))) / 2 + Mt;

	// Write headers for output files
	write_header(output_file1, temp_diff, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue, elastic_beam_length);
	write_header(output_file2, temp_diff, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue, elastic_beam_length);
	write_header(output_file3, temp_diff, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue, elastic_beam_length);
	write_header(output_file4, temp_diff, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue, elastic_beam_length);

	if (!solution_found){
		output_file1 << "Epsilon threshold exceeded without finding solution!";
		output_file2 << "Epsilon threshold exceeded without finding solution!";
		output_file3 << "Epsilon threshold exceeded without finding solution!";
		output_file4 << "Epsilon threshold exceeded without finding solution!";
	}
	else if (max_iter_hit){
		output_file1 << "Maximum iterations in Newton-Raphson did not yield a root solution.";
		output_file2 << "Maximum iterations in Newton-Raphson did not yield a root solution.";
		output_file3 << "Maximum iterations in Newton-Raphson did not yield a root solution.";
		output_file4 << "Maximum iterations in Newton-Raphson did not yield a root solution.";
	}
	else{
		// Calculate portion of beam that is on elastic foundation
		double xs; // x-coordinate of the system (our system of beam separating from foundation)
		double x; // since we're interested in only half the beam, we do this, follows Hetenyi's notation
		double xp; // this is the remaining portion of the beam, see Hetenyi notation
		double y_dis, y_mom, theta_dis, theta_mom, M_dis, M_mom, Q_dis, Q_mom;

		double le = elastic_beam_length*2;
		
		for (int j = 0; j < N; j++){
			xs = ((double)j / (double)N) * elastic_beam_length;
			x = elastic_beam_length + xs;
			xp = elastic_beam_length*2 - x;

			y_dis = -UW/k;
			y_mom = -((2*MT*lambda*lambda)/k)*(1/(sinh(lambda*le)+sin(lambda*le)))*(sinh(lambda*x)*cos(lambda*xp)-cosh(lambda*x)*sin(lambda*xp)+sinh(lambda*xp)*cos(lambda*x)-cosh(lambda*xp)*sin(lambda*x));

			theta_dis = 0;
			theta_mom = -((4*MT*lambda*lambda*lambda)/(k))*((cosh(lambda*x)*cos(lambda*xp)-cosh(lambda*xp)*cos(lambda*x))/(sinh(lambda*le)+sin(lambda*le)));

			M_dis = 0;
			M_mom = (MT/(sinh(lambda*le)+sin(lambda*le)))*(sinh(lambda*x)*cos(lambda*xp)+cosh(lambda*x)*sin(lambda*xp)+sinh(lambda*xp)*cos(lambda*x)+cosh(lambda*xp)*sin(lambda*x));

			Q_dis = 0;
			Q_mom = (2*MT*lambda)*((sinh(lambda*x)*sin(lambda*xp)-sinh(lambda*xp)*sin(lambda*x))/(sinh(lambda*le)+sin(lambda*le)));

			output_file1 << xs << "\t" << y_mom << "\t" << y_dis << "\t" << "0.000\t" << y_mom + y_dis << "\n";
			output_file2 << xs << "\t" << theta_mom << "\t" << theta_dis << "\t" << "0.000\t" << theta_dis + theta_mom << "\n";
			output_file3 << xs << "\t" << M_mom << "\t" << M_dis << "\t" << "0.000\t" << M_mom + M_dis << "\n";
			output_file4 << xs << "\t" << Q_mom << "\t" << Q_dis << "\t" << "0.000\t" << Q_mom + Q_dis << "\n";
		}

		// Calculate cantilevered portion of beam
		double LCant = L / 2 - elastic_beam_length;
		double prev_xs = xs;
		double eqslp; // this is the equivalent slope boundary condition
		eqslp = ((4*MT*lambda*lambda*lambda)/k)*((cosh(lambda*le)-cos(lambda*le))/(sinh(lambda*le)+sin(lambda*le)));
		double c_def, c_slp, c_mom, prev_c_mom, c_shr;
		prev_c_mom = -10000000;
		for (int k = 0; k <= N; k++){
			xs = prev_xs + ((double)k / (double)N)*LCant;
			x = ((double)k / (double)N)*LCant;

			// M(Lcant) = Mt
			c_shr = -UW*x + UW*LCant;
			c_mom = -(UW*x*x) / 2 + UW*LCant*x - (UW*LCant*LCant) / 2 + Mt;
			c_slp = -(UW*x*x*x) / (6 * E*I) + (UW*LCant*x*x) / (2 * E*I) - (UW*LCant*LCant*x) / (2 * E*I) + (Mt*x) / (E*I) + eqslp;
			c_def = -(UW*x*x*x*x) / (24 * E*I) + (UW*LCant*x*x*x) / (6 * E*I) - (UW*LCant*LCant*x*x) / (4 * E*I) + (Mt*x*x) / (2 * E*I) + eqslp*x;

			output_file1 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << -c_def << "\t" << -c_def << "\n";
			output_file2 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << -c_slp << "\t" << -c_slp << "\n";
			output_file3 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << c_mom << "\t" << c_mom << "\n";
			output_file4 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << -c_shr << "\t" << -c_shr << "\n";

			if (prev_c_mom < c_mom){
				prev_c_mom = c_mom;
			}
		}
		//std::cout << "\nPeak Moment: " << prev_c_mom << "\n";
		std::cout << "\nl_e: " << elastic_beam_length << "\n";
		std::cin >> file_name;
	}

}

void type2_Analysis(){
	std::string app_mom, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue;
	std::string num_points, file_name;
	double aM, dT, h, L, bo, E, UW, gamma, k;
	int N;

	std::cout << "Enter applied moment (in*lb): ";
	std::cin >> app_mom;
	std::cout << app_mom << " in*lb\n";
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
	std::cout << "Enter file name for output of points: ";
	std::cin >> file_name;
	std::cout << "Data points will be outputted to " << file_name;

	try{
		aM = stod(app_mom);
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
		// 500 points should be more than adequate
		N = 500;
	}
	catch (const std::invalid_argument ia){
		std::cout << "\nOne of the values entered is incorrect.\nGeneral Error: INVALID ARGUMENT\nSpecific Error: " << ia.what() << "\nTHIS ERROR IS IRRECOVERABLE. PROGRAM TERMINATING...\n";
		exit(EXIT_FAILURE);
	}
	catch (...){
		std::cout << "\nGENERAL EXCEPTION THROWN. UNRECOVERABLE. PROGRAM TERMINATING...\n";
		exit(EXIT_FAILURE);
	}

	std::ofstream output_file1, output_file2, output_file3, output_file4;
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

	// Calculate temperature differential
	dT = (aM*h) / (-gamma*E*I);

	// Find separation point
	double elastic_beam_length;
	elastic_beam_length = newt_rap(UW, lambda, aM, L, k);

	// Calculate prescribed moment on elastic foundation
	// This value is the continuity between the elastic foundation
	// and the cantilevered section
	double MT;
	MT = (-UW*((L/2 - elastic_beam_length)*(L/2 - elastic_beam_length)) / 2) + aM;

	// Write headers for output files
	write_header(output_file1, app_mom, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue, elastic_beam_length);
	write_header(output_file2, app_mom, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue, elastic_beam_length);
	write_header(output_file3, app_mom, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue, elastic_beam_length);
	write_header(output_file4, app_mom, beam_thick, beam_leng, beam_width, Em, rho, cote, kvalue, elastic_beam_length);

	if (!solution_found){
		output_file1 << "Epsilon threshold exceeded without finding solution!";
		output_file2 << "Epsilon threshold exceeded without finding solution!";
		output_file3 << "Epsilon threshold exceeded without finding solution!";
		output_file4 << "Epsilon threshold exceeded without finding solution!";
	}
	else if (max_iter_hit){
		output_file1 << "Maximum iterations in Newton-Raphson did not yield a root solution.";
		output_file2 << "Maximum iterations in Newton-Raphson did not yield a root solution.";
		output_file3 << "Maximum iterations in Newton-Raphson did not yield a root solution.";
		output_file4 << "Maximum iterations in Newton-Raphson did not yield a root solution.";
	}
	else{
		// Calculate portion of beam that is on elastic foundation
		double xs; // x-coordinate of the system (our system of beam separating from foundation)
		double x; // since we're interested in only half the beam, we do this, follows Hetenyi's notation
		double xp; // this is the remaining portion of the beam, see Hetenyi notation
		double y_dis, y_mom, theta_dis, theta_mom, M_dis, M_mom, Q_dis, Q_mom;

		double le = elastic_beam_length * 2;

		for (int j = 0; j < N; j++){
			xs = ((double)j / (double)N) * elastic_beam_length;
			x = elastic_beam_length + xs;
			xp = elastic_beam_length * 2 - x;

			y_dis = -UW / k;
			y_mom = -((2 * MT*lambda*lambda) / k)*(1 / (sinh(lambda*le) + sin(lambda*le)))*(sinh(lambda*x)*cos(lambda*xp) - cosh(lambda*x)*sin(lambda*xp) + sinh(lambda*xp)*cos(lambda*x) - cosh(lambda*xp)*sin(lambda*x));

			theta_dis = 0;
			theta_mom = -((4 * MT*lambda*lambda*lambda) / (k))*((cosh(lambda*x)*cos(lambda*xp) - cosh(lambda*xp)*cos(lambda*x)) / (sinh(lambda*le) + sin(lambda*le)));

			M_dis = 0;
			M_mom = (MT / (sinh(lambda*le) + sin(lambda*le)))*(sinh(lambda*x)*cos(lambda*xp) + cosh(lambda*x)*sin(lambda*xp) + sinh(lambda*xp)*cos(lambda*x) + cosh(lambda*xp)*sin(lambda*x));

			Q_dis = 0;
			Q_mom = (2 * MT*lambda)*((sinh(lambda*x)*sin(lambda*xp) - sinh(lambda*xp)*sin(lambda*x)) / (sinh(lambda*le) + sin(lambda*le)));

			output_file1 << xs << "\t" << y_mom << "\t" << y_dis << "\t" << "0.000\t" << y_mom + y_dis << "\n";
			output_file2 << xs << "\t" << theta_mom << "\t" << theta_dis << "\t" << "0.000\t" << theta_dis + theta_mom << "\n";
			output_file3 << xs << "\t" << M_mom << "\t" << M_dis << "\t" << "0.000\t" << M_mom + M_dis << "\n";
			output_file4 << xs << "\t" << Q_mom << "\t" << Q_dis << "\t" << "0.000\t" << Q_mom + Q_dis << "\n";
		}

		// Calculate cantilevered portion of beam
		double LCant = L / 2 - elastic_beam_length;
		double prev_xs = xs;
		double eqslp; // this is the equivalent slope boundary condition
		eqslp = ((4 * MT*lambda*lambda*lambda) / k)*((cosh(lambda*le) - cos(lambda*le)) / (sinh(lambda*le) + sin(lambda*le)));
		double c_def, c_slp, c_mom, prev_c_mom, c_shr;
		prev_c_mom = -10000000;
		for (int k = 0; k <= N; k++){
			xs = prev_xs + ((double)k / (double)N)*LCant;
			x = ((double)k / (double)N)*LCant;

			// M(Lcant) = Mt
			c_shr = -UW*x + UW*LCant;
			c_mom = -(UW*x*x) / 2 + UW*LCant*x - (UW*LCant*LCant) / 2 + aM;
			c_slp = -(UW*x*x*x) / (6 * E*I) + (UW*LCant*x*x) / (2 * E*I) - (UW*LCant*LCant*x) / (2 * E*I) + (aM*x) / (E*I) + eqslp;
			c_def = -(UW*x*x*x*x) / (24 * E*I) + (UW*LCant*x*x*x) / (6 * E*I) - (UW*LCant*LCant*x*x) / (4 * E*I) + (aM*x*x) / (2 * E*I) + eqslp*x;

			output_file1 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << -c_def << "\t" << -c_def << "\n";
			output_file2 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << -c_slp << "\t" << -c_slp << "\n";
			output_file3 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << c_mom << "\t" << c_mom << "\n";
			output_file4 << xs << "\t" << "0.000" << "\t" << "0.000" << "\t" << -c_shr << "\t" << -c_shr << "\n";

			if (prev_c_mom < c_mom){
				prev_c_mom = c_mom;
			}
		}
		//std::cout << "\nPeak Moment: " << prev_c_mom << "\n";
		std::cout << "\nl_e: " << elastic_beam_length << "\n";
		std::cin >> file_name;
	}

}

int main(){
	int choice;
	std::cout << "Prescribed temperature differential (1) or moment (2): ";
	std::cin >> choice;
	switch (choice){
	case 1:
		type1_Analysis();
		break;
	case 2:
		type2_Analysis();
		break;
	default:
		break;
	}


}
