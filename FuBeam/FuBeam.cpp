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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <exception>
#include <cmath>

int main(){
	std::string temp_diff, beam_thick, beam_leng, Em, rho, cote, kvalue;
	std::string num_points, file_name;
	double dT, h, L, E, UW, alpha, k;

	std::cout << "Enter temperature differential (deg F): ";
	std::cin >> temp_diff;
	std::cout << temp_diff << " deg F\n";
	std::cout << "Enter beam thickness (in): ";
	std::cin >> beam_thick;
	std::cout << beam_thick << " in\n";
	std::cout << "Enter beam length (in): ";
	std::cin >> beam_leng;
	std::cout << beam_leng << " in\n";
	std::cout << "Enter beam modulus of elasticity (psi): ";
	std::cin >> Em;
	std::cout << Em << " psi\n";
	std::cout << "Enter beam unit weight (lbs/ft3): ";
	std::cin >> rho;
	std::cout << rho << " lbs/ft3\n";
	std::cout << "Enter the beam coeff. of thermal expansion (1/F x10-5): ";
	std::cin >> cote;
	std::cout << cote << "E-05 /F\n";
	std::cout << "Enter subgrad k-value (pci): ";
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
		E = stod(Em);
		UW = stod(rho);
		alpha = stod(cote);
		k = stod(kvalue);
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
	// The underscore for the hardfail option is a bit puzzling. There is no data from MSDN about it.
	output_file.exceptions(std::ofstream::failbit | std::ofstream::badbit | std::ofstream::_Hardfail);
	try{
		output_file.open(file_name);
	}
	catch (std::ofstream::failure e){
		std::cout << "Failure reading creating file. Please make sure the file name does not conflict with an existing name. Also make sure you have write access to the location you are outputting data to.";
		exit(EXIT_FAILURE);
	}
}