# autocorrelation_multithreaded
A multi-threaded autocorrelation calculator written in C++



##########################################################################################
#######  readme for autocorr_threads.exe  ################################################
##########################################################################################

This executable runs autocorrelation on the binary file provided by pdatatxt2bin.

Results are outputted to autocorr.csv for processing.


##########################################################################################
##################  CMAKE  ###############################################################
##########################################################################################

module load cmake
module load gcc/7
cd src									# Go into src dir
mkdir build								# Make build dir
cd build								# Move into build dir
cmake -DCMAKE_BUILD_TYPE=Release ..		# Produce makefiles in current folder using
										# 		CMakeLists.txt from 'src' dir.
										# Compile with 'release' flag
make									# Make the executable


##########################################################################################
##################  USAGE  ###############################################################
##########################################################################################

"./autocorrelate_threads.exe 16 1 200000 < ./pressure.bin > /dev/null"
			Where: 	
				-pressure.bin is the binary file produced by pdatatxt2bin
				-16 is the number of threads to use (more always better--this is always
														highly parallelizable process)
				-1 is the sample stepper. 1 means we look at every data point. 10 would 
														mean we take every 10 data points
				-200000 is the number of lagsteps to run
														
														
##########################################################################################
##############  SPECIAL NOTES  ###########################################################
##########################################################################################




##########################################################################################
##################  INPUT  ###############################################################
##########################################################################################

The bin file used here is a the binary file created by pdatatxt2bin.
It is of the following format:
	First:
		numrows (or number of steps) [uint32_t]
	Then the following repeating:
		PTxy [float] 
		PTxz [float]
		PTyz [float]
		PTxx [float]
		PTyy [float]
		PTzz [float]
		
Note that:
	sizeof(float) = 4 bytes
	sizeof(uint32_t) = 4 bytes

	Thus it should be approx of size: 
		24 bytes x {numsteps}
		
		Ex: 4.0m steps = 96MB

							
##########################################################################################
##################  OUTPUT  ##############################################################
##########################################################################################

The resulting autocorr.csv is formatted in the following manner:
	lag, PTxy*PTxy, PTxz*PTxz, PTyz*PTyz, PTxx*PTxx, PTyy*PTyy, PTzz*PTzz

