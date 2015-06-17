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
FORTRAN_COMPILER_intel = ifort

exec_file= ./bin/run.exe
module_files= src/linearsolvers.f90 src/femgeometry.f90 src/materialbehavior.f90 src/femlibs.f90 
driver_file = driver/elec_mech_fem.f90
compiled_module_files=fem_geometry.mod fem_libs.mod linearsolvers.mod material_behavior.mod
output_files=*.mod *.txt *.vtu fort.* ./gnuplot/*.gnu ./excel/*.csv
output_files_heavy=*.mod *.txt


library_files_intel=-mkl -lmkl_lapack95_lp64


library_files=c:/lapack/liblapack.a c:/lapack/librefblas.a
gnu_plot_command_file = ./gnuplot/gnu_plotter_afc_frequency.txt

#flags=-O1 -Wall -g -fcheck=all
flags=-O3   

all: 
#   del $(output_files) $(exec_file)
	$(FORTRAN_COMPILER) $(flags) -o\
	$(exec_file) \
	$(module_files) $(driver_file) $(library_files)
	$(exec_file)	
# 	gnuplot $(gnu_plot_command_file) 


all_intel: 
	$(FORTRAN_COMPILER_intel) $(flags) -o\
	$(exec_file) \
	$(module_files) $(driver_file) $(library_files_intel)
	$(exec_file)
	 gnuplot $(gnu_plot_command_file)  

light: 
	rm $(output_files_heavy) $(exec_file) ;

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