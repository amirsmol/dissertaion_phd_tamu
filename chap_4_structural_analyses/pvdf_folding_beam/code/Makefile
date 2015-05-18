# ============================================================================
# Name        : Makefile
# Author      : amir sohrabi
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean run

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = ifort
exec_file=    bin/run.exe
module_files= src/linearsolvers.f90 src/femgeometry.f90 src/materialbehavior.f90 src/femlibs.f90 
driver_file = src/fem_piezeo_3d__visco_polarization_switching.f90
compiled_module_files=fem_geometry.mod fem_libs.mod linearsolvers.mod material_behavior.mod
output_files=*.mod *.exe *.csv *.out *.txt fort.* *.vtu *.pdf *.eps *.gnu
#library_files=$(home)/lapack/liblapack.a $(home)/lapack/libblas.a $(home)/lapack/libtmglib.a
library_files=-mkl -lmkl_lapack95_lp64
gnu_plot_command_file=src/gnu_plotter.txt
flags=-O1 -Wall -g -fcheck=all

all: 
	$(FORTRAN_COMPILER) -o\
	$(exec_file) \
	$(module_files) $(driver_file) $(library_files)
	./$(exec_file)
#	gnuplot $(gnu_plot_command_file)  
	


fullall: 
	del $(output_files) $(exec_file)	
	$(FORTRAN_COMPILER)  $(flags) -o\
	$(exec_file) \
	$(module_files) $(driver_file)	$(library_files) 
	$(exec_file)
#	gnuplot $(gnu_pl



clean:
	rm $(output_files) $(exec_file) ;
	
	
run: 
	rm $(output_files) $(exec_file)
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