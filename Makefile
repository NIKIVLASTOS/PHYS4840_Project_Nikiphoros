#==========================                                          
# Fixed Makefile                                                     
#==========================

FC = gfortran                                                        # Fortran compiler
FFLAGS = -O2 -Wall                                                   # Compiler flags optimize and show all warnings

MODS = potential_functions.f90 matrix_tools.f90 normalization.f90   # Module source files
MAIN = main_solver.f90                                              # Main program file
EXEC = schrodinger_solver                                           # Executable name

all: $(EXEC)                                                        # build the executable

$(EXEC): $(MODS) $(MAIN)                                            # compile and link all source files
	$(FC) $(FFLAGS) -c $(MODS)                                       # Compile modules to object files
	$(FC) $(FFLAGS) -o $(EXEC) *.o $(MAIN)                           # Link object files with main to create executable

run: $(EXEC)                                                        # run  program
	./$(EXEC)                                                       # Execute fortran program
	python3 plot_results.py                                         # run the Python script to plot results

clean:                                                              # to clean up build artifacts
	rm -f $(EXEC) *.o *.mod                                         # Remove executable, object files, and module files
	rm -rf __pycache__ plots                                       # Remove python cache and plots folder
