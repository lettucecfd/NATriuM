/*
 * FFT.cpp
 *
 *  Created on: Jan 12, 2016
 *      Author: kraemer
 */

#include <algorithm>
#include <fstream>
#include <string>
#include <iostream>
#include "math.h"

#include "fftw3.h"

using namespace std;

int main(int argc, char** argv) {
	cout << "Fourier transform with fftw3" << endl;
	cout << "Usage: fft <dim> <number of points in each direction> <filename_in> <filename_out>" << endl;
	size_t dim = atoi(argv[1]);
	size_t N = atoi(argv[2]);
	string filename_in(argv[3]);
	string filename_out(argv[4]);

	// allocate
	size_t array_size = pow(N,dim);
	fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * array_size);
	fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * array_size);
	int Ns[dim];
	std::fill_n(Ns, dim, N);
	fftw_plan p = fftw_plan_dft(dim, Ns, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	// initialize input (has to be done after creating the plan)
	cout << "Read input..." << endl;
	std::ifstream infile;
	infile.open(filename_in.c_str(), ios::in);
	string line;
	double number;
	if (infile.is_open())
	 {
		size_t i = 0;
	   while ( getline (infile, line) )
	    {
		  if ( i>=array_size ) {
			  break;
		  }
	      cout << stod(line) << '\n';
	      in[i][0] =  stod(line); //Re
	      in[i][1] =  0.00000000; //Im
	      i++;
	    }
	   infile.close();
	  }
	cout << "... done" << endl;
	cout << "Execute Fourier transform..." << endl;

	// execute fourier transform
	fftw_execute(p);
	cout << "done" << endl;


	// write output
	cout << "Write output to " << filename_out << endl;
	std::ofstream outfile;
	outfile.open(filename_out.c_str(), ios::out);
	for (size_t i = 0; i < array_size; i++){
		outfile << i % N << " " << i / N << " " << out[i][0]  << " " <<  out[i][1] << endl;
	}

	// clean up
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
} /** end main */
