// Autocorrelation.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
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
#include <cstdlib>
#include "Barrier.hh"
#include <fstream>
#include <vector>
//#include <cublas_v2.h>

std::mutex mutex;


//LM code
class Grid {
	uint32_t _rows;
	uint32_t _columns;
	std::unique_ptr<float[]> data; // data is a pointer
	
public:
	Grid(uint32_t rows, uint32_t columns)
		: _rows{ rows },
		  _columns{ columns },
		  data{ std::make_unique<float[]>(rows * columns) }
	{}
	uint32_t rows() const {
		return _rows;
	}
	uint32_t columns() const {
		return _columns;
	}
	float * operator[](uint32_t row) {
		return row * _columns + data.get();
	}
}; //end LM Code

// Variable Declaration
std::string num_threads_str;
size_t num_threads;
std::string sample_str;
size_t sample;
std::string lagmax_str;
uint32_t lagmax;
uint32_t numrows;
uint32_t numcols = 6;
uint32_t lagrows;
// End of variable declaration

//Fn to read input file to grid
void fn_initialize(Grid & grid)
{
	for (uint32_t row = 0; row <= numrows - 1; row++) //rows
	{
		std::cin.read(reinterpret_cast<char *>(&grid[row][0]), sizeof(float)); //PTxy
		std::cin.read(reinterpret_cast<char *>(&grid[row][1]), sizeof(float)); //PTxz
		std::cin.read(reinterpret_cast<char *>(&grid[row][2]), sizeof(float)); //PTyz
		std::cin.read(reinterpret_cast<char *>(&grid[row][3]), sizeof(float)); //PTxx
		std::cin.read(reinterpret_cast<char *>(&grid[row][4]), sizeof(float)); //PTyy
		std::cin.read(reinterpret_cast<char *>(&grid[row][5]), sizeof(float)); //PTzz
	}
}
//////////////////////////////

//Fn to get average for a cell
float fn_product(Grid & grid, uint32_t row, uint32_t column, uint32_t lag)
{
	float product = (grid[row][column] * grid[row+lag][column]);
	return product / (numrows - lag);
}

//DEBUG FUNCTIONS/////////////////////////////////////////////////////////////////////////////
//Fn to print grid for debug
void fn_debug_printgrid(Grid & grid, uint32_t maxrow)
{
	std::cerr << "Debug: printing first 3 rows of this grid..." << std::endl;
	for (uint32_t row = 0; row < 3; row++) //rows
	{
		std::cerr << grid[row][0] << "     " << grid[row][1] << "     " << grid[row][2] 
				<< "     " << grid[row][3] << "     " << grid[row][4] << "     " << grid[row][5] << std::endl;
	}

	std::cerr << "Debug: printing last 3 rows of this grid.." << std::endl;
	for (uint32_t row = maxrow - 3; row < maxrow; row++) //rows
	{
		std::cerr << grid[row][0] << "     " << grid[row][1] << "     " << grid[row][2] 
				<< "     " << grid[row][3] << "     " << grid[row][4] << "     " << grid[row][5] << std::endl;
	}
}
////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////



//Fn for calculation
void fn_calculate(Grid & datagrid, uint32_t lag, float & sum_PTxy, float & sum_PTxz, float & sum_PTyz
					, float & sum_PTxx, float & sum_PTyy, float & sum_PTzz, int thread_num, size_t & num_threads)
{
	//declare variables
	thread_local float loc_sum_PTxy = 0.0f;
	thread_local float loc_sum_PTxz = 0.0f;
	thread_local float loc_sum_PTyz = 0.0f;
	thread_local float loc_sum_PTxx = 0.0f;
	thread_local float loc_sum_PTyy = 0.0f;
	thread_local float loc_sum_PTzz = 0.0f;
	int startrow;				//starting row for the specific logic core
	int endrow;					//ending row ^
	int lagstodo;				//number of lags we need to do
	int lastlag;				//the last lagstep to run on otherwise it will break
	lagstodo = (numrows-lag);
	lastlag = (numrows-lag);
	//Splitting array into chunks for the threads//////////
	if (thread_num == 0) { 										// For first thread only
		startrow = 0;											// Start at row 0
		endrow = (lagstodo / num_threads);						// Run UPTO just before row: lagstodo/thread
																// With numrows=4m,numthread=24,lag=0: 		0->166,666
																// With numrows=4m,numthread=24,lag=170k:	0->159,583
	}
	else if (thread_num == num_threads - 1) { 					// For last thread only
		startrow = (lagstodo / num_threads) * thread_num;		// Start at rows/thread # of rows from the end
		endrow = lastlag;										// Run UPTO numrows-lag
																// With numrows=4m,numthread=24,lag=0: 		3,833,318->4,000,000
																// With numrows=4m,numthread=24,lag=170k: 	3,670,409->3,830,000
	}
	else { 														// For all other threads
		startrow = (lagstodo / num_threads) * thread_num;		// Start at rows/thread mult. by the thread number
		endrow = (lagstodo / num_threads) * (thread_num + 1);	// End rows/thread rows from start
		if (endrow > (lastlag)) {								// **This line truncates endrow to ensure that fn_product
			endrow = lastlag;									//		is never run with (lag+row)>endrow**
		}														// With numrows=4m,numthread=24,lag=0: 		166,666->333,333
	}															// With numrows=4m,numthread=24,lag=170k: 	166,666->333,333
	//////////////////////////////////////////////////////////////

	loc_sum_PTxy = 0.0f;
	loc_sum_PTxz = 0.0f;
	loc_sum_PTyz = 0.0f;
	loc_sum_PTxx = 0.0f;
	loc_sum_PTyy = 0.0f;
	loc_sum_PTzz = 0.0f;
	
	// 
	for (size_t row = startrow; row < endrow; ++row)						// For 24threads,numrows=4m,numthread=24,lag=500k:
	{																		// Thread 20 has 3.33m->3.5m, exits loop after running row=3,499,999
		loc_sum_PTxy = loc_sum_PTxy + fn_product(datagrid, row, 0, lag);	// Thread 21 has 3.50m->3.5m, does not enter the loop
		loc_sum_PTxz = loc_sum_PTxz + fn_product(datagrid, row, 1, lag);	// Thread 22 has 3.66m->3.5m, does not enter the loop
		loc_sum_PTyz = loc_sum_PTyz + fn_product(datagrid, row, 2, lag);	// Thread 23 has 3.83m->3.5m, does not enter the loop
		loc_sum_PTxx = loc_sum_PTxx + fn_product(datagrid, row, 3, lag);	// Running after 3.50m (w/ lag=500k) breaks fn_product
		loc_sum_PTyy = loc_sum_PTyy + fn_product(datagrid, row, 4, lag);	//
		loc_sum_PTzz = loc_sum_PTzz + fn_product(datagrid, row, 5, lag);	//
	}
	
	//mutex.lock();
	//DEBUG//std::cerr << "Thread " << thread_num << " reporting these sums for lag = " << lag << ":" << std::endl;
	//DEBUG//std::cerr << loc_sum_ppxy << "     " << loc_sum_ppxz << "     " << loc_sum_ppyz
	//DEBUG//		<< "     " << loc_sum_ppxxyy << "     " << loc_sum_ppyyzz << std::endl;
	//mutex.unlock();
	
	// !!!!! Thread-exclusive access to the global sum variables
	mutex.lock();
	sum_PTxy = sum_PTxy + loc_sum_PTxy;
	sum_PTxz = sum_PTxz + loc_sum_PTxz;
	sum_PTyz = sum_PTyz + loc_sum_PTyz;
	sum_PTxx = sum_PTxx + loc_sum_PTxx;
	sum_PTyy = sum_PTyy + loc_sum_PTyy;
	sum_PTzz = sum_PTzz + loc_sum_PTzz;
	mutex.unlock();
/////////////////////////////////////////////////////////////
}


void fn_runlags(Grid & datagrid, Grid & acgrid, float & sum_PTxy, float & sum_PTxz, float & sum_PTyz
					, float & sum_PTxx, float & sum_PTyy, float & sum_PTzz, int thread_num, size_t & num_threads, fsl::Barrier & barrier) 
{
	for (size_t lag = 0; lag <= lagmax; ++lag)
	{
		barrier.wait();
		mutex.lock();		
		sum_PTxy = 0.0f;
		sum_PTxz = 0.0f;
		sum_PTyz = 0.0f;
		sum_PTxx = 0.0f;
		sum_PTyy = 0.0f;
		sum_PTzz = 0.0f;
		mutex.unlock();
		barrier.wait();
		
		//fn_calculate(Grid & datagrid, uint32_t lag, float & sum_ppxy, float & sum_ppxz, float & sum_ppyz
		//				, float & sum_ppxxyy, float & sum_ppyyzz, int thread_num, size_t & num_threads)
		fn_calculate(datagrid, lag, sum_PTxy, sum_PTxz, sum_PTyz, sum_PTxx, sum_PTyy, sum_PTzz, thread_num, num_threads);

		//Run this on one thread only
		barrier.wait();
		mutex.lock();
		if (thread_num==0) {
			//DEBUG//std::cerr << "_____________________________________________________________" << std::endl;
			//DEBUG//std::cerr << "Thread " << thread_num << " reporting these TOTAL sums:" << std::endl;
			//DEBUG//std::cerr << sum_ppxy << "     " << sum_ppxz << "     " << sum_ppyz
			//DEBUG//		<< "     " << sum_ppxxyy << "     " << sum_ppyyzz << std::endl;
			acgrid[lag][0] = sum_PTxy;
			acgrid[lag][1] = sum_PTxz;
			acgrid[lag][2] = sum_PTyz;
			acgrid[lag][3] = sum_PTxx;
			acgrid[lag][4] = sum_PTyy;
			acgrid[lag][5] = sum_PTzz;
			///DEBUG//std::cerr << "ac grid at lag =  " << lag << ",  reads:" << std::endl;
			//DEBUG//std::cerr << acgrid[lag][0] << "     " << acgrid[lag][1] << "     " << acgrid[lag][2]
			//DEBUG//		<< "     " << acgrid[lag][3] << "     " << acgrid[lag][4] << std::endl;
			//DEBUG//std::cerr << "_____________________________________________________________" << std::endl;
			if (lag % 10000 == 0) {
				std::cerr << "Finished lag " << lag << "." << std::endl;
				std::cerr << "values = " << sum_PTxy  << "," << sum_PTxz  << ","  << sum_PTyz  << ","  << sum_PTxx  << "," << sum_PTyy  << "," << sum_PTzz << std::endl;
			}
		}
		mutex.unlock();
		barrier.wait();
		//Resume threads
		
	}
}

////////////////////////////////////////////////////////////////////


/*
//Fn to calculate autocorrelation
void fn_autocorr(Grid & grid, Grid & grid2)
{

}
/////////////////////////////////
*/

/*
// Example kernel definition (matrix add)
__global__ void MatAdd(float A[N][N], float B[N][N], float C[N][N])
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if (i < N && j < N)
        C[i][j] = A[i][j] + B[i][j];
}
*/

int main(int argc, char*argv[])
{
	//Command Line Argument Interpretation
	// Check the number of parameters
	if (argc != 4) {
		//Notify the user that they need to specify number of threads
		std::cerr << "Please specify the number of threads on which to run, as well as the sample step size and number of lag steps" << std::endl;
		std::cerr << "(Example of 16-thread usage: '~/software/DJC_autocorr_threads/autocorrelate_threads_notindep.exe 16 1 200000 < ./pressure.bin > /dev/null' )" << std::endl;
		return 1;
	} else {
		num_threads_str = argv[1];
		num_threads = std::stoi(num_threads_str);
		sample_str = argv[2];
		sample = std::stoi(sample_str);
		lagmax_str = argv[3];
		lagmax = std::stoul(lagmax_str);
		lagrows = lagmax+1;
	}	
	//End Command Line Argument Interpretation
	
	//Read header info
	std::cerr << "Reading header info..." << std::endl;
	std::cin.read(reinterpret_cast<char *>(&numrows), sizeof(uint32_t));
	std::cerr << numrows << " rows found." << std::endl;
	//////////////////
	
	//Setup grids
	std::cerr << "Allocating memory for grids..." << std::endl;
	Grid grid_pdata{ numrows, numcols };
	Grid grid_autocorr{ lagrows, numcols };
	
	std::cerr << "Reading std::cin to datagrid" << std::endl;
	fn_initialize(grid_pdata); //Read binary file via stdin to grid_pdata
	float sum_PTxy = 0.0f;
	float sum_PTxz = 0.0f;
	float sum_PTyz = 0.0f;
	float sum_PTxx = 0.0f;
	float sum_PTyy = 0.0f;
	float sum_PTzz = 0.0f;
	
	//Debug
	std::cerr << "Printing from grid_pdata..." << std::endl;
	fn_debug_printgrid(grid_pdata, numrows);
	
	//MULTITHREADING//////////////////////////////////////////////////////////////////////
	//Vector of threads
	std::cerr << "Creating vector of threads..." << std::endl;
	std::vector<std::thread> threads {};
	
	int thread_num;
	fsl::Barrier barrier(num_threads);
	
	//Spawn threads and run fn_runlags
	//fn_runlags(Grid & datagrid, Grid & acgrid, float & sum_ppxy, float & sum_ppxz, float & sum_ppyz
	//	, float & sum_ppxxyy, float & sum_ppyyzz, int thread_num, size_t & num_threads
	//	, fsl::Barrier & barrier) 
	std::cerr << "Spawning threads and running fn_runlags..." << std::endl;
	for(int i = 0; i < num_threads; ++i){
		thread_num = i;
		threads.emplace_back(fn_runlags, std::ref(grid_pdata), std::ref(grid_autocorr), std::ref(sum_PTxy)
							, std::ref(sum_PTxz), std::ref(sum_PTyz), std::ref(sum_PTxx), std::ref(sum_PTyy)
							, std::ref(sum_PTzz), thread_num, std::ref(num_threads), std::ref(barrier));
	}
	
	//Join threads after stabilized then output
	for (auto & thread: threads) {
		if (thread.joinable()) 
		{
			thread.join();
		}	
	}
	//////////////////////////////////////////////////////////////////////////////////////
	std::cerr << "Processing complete..." << std::endl;
	
	
	/*DEBUG//std::cerr << "Printing lag = 0 pp array..." << std::endl;
	//fn_product(Grid & grid, uint32_t row, uint32_t column, uint32_t lag)
	// Print first 3 of lag = 0 pp array
	for (int row = 0; row < 3; ++row) {
		float prpxy = fn_product(grid_pdata, row, 0, 0) * (numrows-0);
		float prpxz = fn_product(grid_pdata, row, 1, 0) * (numrows-0);
		float prpyz = fn_product(grid_pdata, row, 2, 0) * (numrows-0);
		float prpxxyy = fn_product(grid_pdata, row, 3, 0) * (numrows-0);
		float prpyyzz = fn_product(grid_pdata, row, 4, 0) * (numrows-0);
		std::cerr << row << "     "  << prpxy << "     "  << prpxz << "     "  << prpyz << "     "  << prpxxyy << "     "  << prpyyzz << "     " << std::endl;
	}
	// Print last 3 of lag = 0 pp array
	for (int row = numrows - 3; row < numrows; ++row) {
		float prpxy = fn_product(grid_pdata, row, 0, 0) * (numrows-0);
		float prpxz = fn_product(grid_pdata, row, 1, 0) * (numrows-0);
		float prpyz = fn_product(grid_pdata, row, 2, 0) * (numrows-0);
		float prpxxyy = fn_product(grid_pdata, row, 3, 0) * (numrows-0);
		float prpyyzz = fn_product(grid_pdata, row, 4, 0) * (numrows-0);
		std::cerr << row << "     "  << prpxy << "     "  << prpxz << "     "  << prpyz << "     "  << prpxxyy << "     "  << prpyyzz << "     " << std::endl;
	}
	//Write to file
	//std::cerr << "Printing result autocorr array for debug..." << std::endl;
	//fn_debug_printgrid(grid_autocorr, lagrows);
	/////////////////
	*/

	std::ofstream myfile;
	myfile.open ("autocorr_sample" + sample_str + ",6PT.csv");
	for(int lag = 0; lag < lagrows; ++lag){
		int lagsample = lag * sample;
		myfile << lagsample << "," << grid_autocorr[lag][0] << "," << grid_autocorr[lag][1]
				<< "," << grid_autocorr[lag][2] << "," << grid_autocorr[lag][3]
				<< "," << grid_autocorr[lag][4] << "," << grid_autocorr[lag][5] << std::endl;
	}
	myfile.close();
	
	
	return 0;
}

