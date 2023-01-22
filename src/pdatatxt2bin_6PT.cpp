#include <stdio.h>
#include <memory>
#include <iostream>
#include <cmath>
#include <math.h>
#include <cstdint> // for uint32_t 
#include <chrono> //time
#include <algorithm>
#include <thread>
#include <mutex>
#include <string>


//Variable declarations
uint32_t numrows;
uint32_t numcols;
uint32_t read_step;
float read_PTxy;
float read_PTxz;
float read_PTyz;
float read_PTxx;
float read_PTyy;
float read_PTzz;
char title_line[128];

/* Note that the binary output file is of the following format:

First: numrows (or number of steps) of size uint32_t
Then the following repeating:
PTxy [float] 
PTxz [float]
PTyz [float]
PTxx [float]
PTyy [float]
PTzz [float]

*/


int main(int argc, char*argv[])
{
	//Command Line Argument Interpretation
	// Check the number of parameters
	if (argc != 2) {
		//Notify the user that they need to specify datafile name
		std::cerr << "argc = " << argc << std::endl;
		std::cerr << "Please specify the filename on which to run." << std::endl;
		std::cerr << "(Example of usage: './pdatatxt2bin pressure.data > pressure.bin' )" << std::endl;
		return 1;
	}
	else {
		std::cerr << "Opening file: " << argv[1] << std::endl;
	}
	//End Command Line Argument Interpretation

	//Count lines///////////////////////////////////////////////
	FILE *fp_count;
	int count_lines = 0;
	char chr;
	std::cerr << "Counting variable definitions made. Now opening file." << std::endl;

	fp_count = fopen(argv[1], "r");
	//extract character from file and store in chr
	chr = getc(fp_count);
	std::cerr << "Starting counting loop." << std::endl;
	while (chr != EOF)
	{
		//Count whenever new line is encountered
		if (chr == '\n')
		{
			count_lines = count_lines + 1;
		}
		//take next character from file.
		chr = getc(fp_count);
	}
	fclose(fp_count); //close file.
	numrows = count_lines - 1; //take away one for the title line lammps puts in
	std::cerr << "Found numrows = " << numrows << ", writing to binary file..." << std::endl;
	std::cout.write(reinterpret_cast<char const *>(&numrows), sizeof(uint32_t));
	//Done counting lines///////////////////////////////////////


	//Read lammps txt p tensor data ('r' for read only) and write to bin//////////////////
	FILE *fp;
	fp = fopen(argv[1], "r");
	std::cerr << "Starting conversion..." << std::endl;
	//We need to purge the title line LAMMPS outputs looks like:
	// "# Fix print output for fix fxprint"
	fscanf(fp, "%[^0]", title_line);
	for (uint32_t i = 0; i < numrows; i++)
	{
		fscanf(fp, "%u", &read_step);
		fscanf(fp, "%f", &read_PTxy);
		fscanf(fp, "%f", &read_PTxz);
		fscanf(fp, "%f", &read_PTyz);
		fscanf(fp, "%f", &read_PTxx);
		fscanf(fp, "%f", &read_PTyy);
		fscanf(fp, "%f", &read_PTzz);
		//std::cout.write(reinterpret_cast<char const *>(&read_step), sizeof(uint32_t));
		std::cout.write(reinterpret_cast<char const *>(&read_PTxy), sizeof(float));
		std::cout.write(reinterpret_cast<char const *>(&read_PTxz), sizeof(float));
		std::cout.write(reinterpret_cast<char const *>(&read_PTyz), sizeof(float));
		std::cout.write(reinterpret_cast<char const *>(&read_PTxx), sizeof(float));
		std::cout.write(reinterpret_cast<char const *>(&read_PTyy), sizeof(float));
		std::cout.write(reinterpret_cast<char const *>(&read_PTzz), sizeof(float));
		if (i % 200000 == 0) {
			std::cerr << "Finished scanning line " << i << "." << std::endl;
			std::cerr << "values = " << read_PTxy  << "," << read_PTxz  << ","  << read_PTyz << std::endl;
		}
	}
	//End of binary write////////////////////////////////////////////////////////////////

	return 0;
}
