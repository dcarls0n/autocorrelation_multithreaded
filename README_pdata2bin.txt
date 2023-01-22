
##########################################################################################
#######  readme for pdatatxt2bin.exe  ####################################################
##########################################################################################

This executable changes the plaintext pressure.data file outputted from lammps into a binary '.bin' file.


##########################################################################################
##################  CMAKE  ###############################################################
##########################################################################################

cd src									# Go into src dir
mkdir build								# Make build dir
cd build								# Move into build dir
cmake -DCMAKE_BUILD_TYPE=Release ..		# Produce makefiles in current folder using CMakeLists.txt from 'src' dir
make									# Make the executable


##########################################################################################
##################  USAGE  ###############################################################
##########################################################################################

"./pdatatxt2bin pressure.data > pressure.bin"
			Where: pressure.data is the lammps prepared data file



##########################################################################################
##############  DATAFILE PREP  ###########################################################
##########################################################################################

Please ensure the 'pressure.data' file is prepared in the following manner:

	-Formatted in the following manner:
			step,PTxy,PTxz,PTyz,PTxx,PTyy,PTzz
	-Duplicate rows are removed. index should go from 1 to X (X being the total number of steps), 
			increasing by 1 with each row (with no repeats)
	-Any rows including comments or text are removed
			Commonly: "# Fix print output for fix fxprint" is found for each restart. This can be removed with sed via:
							"sed -i '# Fix print output for fix fxprint' ./pressure.data"
							

##########################################################################################
##################  OUTPUT  ##############################################################
##########################################################################################

The bin file created here is a binary file of the following format:
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
