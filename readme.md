# Project 7 - Noah Egger

## Program Goal

This program solves the time-dependent Schrodinger equation and evolves the initial gaussian wave function in time using the Crank-Nicolson method. Two separate potentials are consider: one where V(x) = 0, another where V(x) = 1/2kx^2.

Namelist files are provided that specify the following information: the length of the potential region in which the wave function resides (-length to length), the number of sampling point across -length to length, the number of time steps, the size of the time steps, the width of the initial gaussian wave function, the central maximum of the initial wave function, and the "spring" constant value k that determines the strength of the harmonic oscillator potential. Moreover, the namelist files provide file names to store desired output data to. Note: if a namelist file is not used, default values are provided.

The program works by discretizing the wave function into a series of spatial vectors specificed in time. 
After constructing the Hamiltonian as a matrix, we split up the hamiltonian using the Crank-Nicolson method which divides the hamiltonian into half implicit and half explicit part. Rearranging this equation, we obtain a matrix multiplication problem for evolving the wave function to the next step. Ignoring complex variables, we work with purely reals. This is accomplished by re-writing our main equation as two coupled equations and combining them into a "super-matrix" where the diagonal "elements" are the identity matrix and the off-diagonal "elements" contain multiples of the hamiltonian matrix. Therefore if the wave function is discretized as a (complex) vector of length N, we replace it with a purely real vector of length 2N and purely real matrices of size 2N x 2N. Evolving the wave function to the next step in time becomes a matrix multiplication problem where we multiply the inverse of the left "super-matrix" by the right "super-matrix" times the initial wave function. From the wave function as a function of space and time, we may calculate the probability density as a function of space and time as well by taking the sum of squared magnitudes of the real and imaginary parts of the wave fucntion. Having found the density, we may calculate and store the expectation values for position, width, and normalization all as a function of time. Results are ultimately written to two files: a density file and a time file. For the former, the probability density is written as a 2D matrix indicating the value of the density in position and time. For the latter, said expectation values are written along with the corresponding moment in time along with an analytic values for the width as a function of time for comparison later. Moreover, analytic values for the position as a function of time are included in an additional output file for the case of a harmonic potential used. 


## Directions for Usage
Ensure all files are contained within your direction. Navigate to your directory through terminal and type `make`. Press "enter".  This will compile all files and create an executable called `schrodinger`. Type `./schrodinger file_name.namelist` into your terminal, with "file_name" being the name of the namelist file you wish to use. If you do not wish to use a namelist file, ignore everything after the executable. Press "enter". The results for a particular set of parameters specified in chosen namelist will be written to a file that the user specifies within the namelist file. If the default values are used, the program will print results to a file named `density_results.dat` and `time_results.dat`.If the `zero_potential.namelist` file is used, the file will print to the same file names. If the `harmonic_potential.namelist` file is used, results will be printed to a files named `density_results_harmonic.dat` and `time_results_harmonic.dat`. If the `harmonic_potential_2.namelist` file is used, results will be printed to a files named `density_results_harmonic_2.dat` and `time_results_harmonic_2.dat`. These files are read by `project7_analysis.ipynb`. To visualize results, open  `project7_analysis.ipynb` and click "run all". Note, some plots will only appear if the user has run the program with all provided namelist files. These plots include visualizations for snapshots of the probability density, width v.s. time, position v.s. time, normalization v.s time, as well as numerical v.s. analytic comparisons in width and position as a function of time, all for the two distinct potentials and a modification to the harmonic potential.



## Program Contents

`read_write.f90` reads the input file (or default values), and writes results to new files

`quantum.f90` contains the subroutines which creates the sampling locations, constructs the initial wave function at t=0, constructs the matrix that evolves the wave function, performs the evolution of the wave function, and calculates the expectation values. 

`hamiltonian.f90` constructs the kinetic energy, the potential energy, and thus the hamiltonian in its matrix form, as well as performing the operations which fill arrays for the normalization, position, and gaussian width at each time step. 

`linear_algebra.f90` contains the subroutines for the main mathematics employed in this program, namely matrix multiplication, inversion, and a system solver.

`types.f90` contains the argument types in the program   

`assignment_07.ipynb` plots results of program

`makefile` compiles all .f90 files and creates the executable.    

`main.f90` contains the main calls to run the program.  

`zero_potential.namelist` zero potential namelist file 

`harmonic_potential.namelist` harmonic potential namelist file

`harmonic_potential_2.namelist` harmonic potential namelist file but with different initial width vallue