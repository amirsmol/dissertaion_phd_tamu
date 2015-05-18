# ============================================================================
# Name        : Makefile
# Author      : amir sohrabi
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean run

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran
exec_file=    bin\run.exe 
module_files= src\fem_functions_and_parameters.f90 src/linearsolvers.f90 src/femgeometry.f90 src/materialbehavior.f90 src/femlibs.f90   
driver_file = driver\elec_mech_fem.f90
compiled_module_files=fem_functions_and_parameters.mod fem_geometry.mod
output_files=*.mod *.exe *.csv *.out *.txt fort.* *.vtu *.pdf *.eps *.gnu .\gnuplot\*.gnu .\excel\*.csv .\vtu\*.vtu
library_files=c:\lapack\liblapack.a c:\lapack\librefblas.a
gnu_plot_command_file=excel\gnu_plotter.txt
#flags=-O1 -Wall -g -fcheck=all
flags=-malign-double 

all: 	
	del $(output_files) $(exec_file) ;
	$(FORTRAN_COMPILER)  $(flags) -o\
	$(exec_file) \
	$(module_files) $(driver_file)	$(library_files) 	
	$(exec_file)
#	gnuplot $(gnu_plot_command_file)  


#all: 
#	del $(output_files) $(exec_file)	
#	$(FORTRAN_COMPILER) $(flags) -o\
#	$(exec_file) \
#	$(module_files) $(driver_file) $(library_files)
#	$(exec_file)
##	gnuplot $(gnu_plot_command_file)  

fullall: 
	del $(output_files) $(exec_file)	
	$(FORTRAN_COMPILER)  $(flags) -o\
	$(exec_file) \
	$(module_files) $(driver_file)	$(library_files) 
	$(exec_file)
#	gnuplot $(gnu_pl

clean:
	del $(output_files) $(exec_file) ;
	
	
run: 
	del $(output_files) $(exec_file)
	$(FORTRAN_COMPILER) $(flags)  -o\
	$(exec_file) \
	$(module_files) $(driver_file) $(library_files)
	$(exec_file)
	gnuplot $(gnu_plot_command_file)
	
plot: 
	gnuplot $(gnu_plot_command_file)
	C:\Program Files (x86)\Google\Chrome\Application\chrome.exe electric_field_vs_polarization.pdf 	

build: 
	del $(output_files) $(exec_file)	
	$(FORTRAN_COMPILER)  $(flags) -o\
	$(exec_file) \
	$(module_files) $(driver_file)	$(library_files) 	 