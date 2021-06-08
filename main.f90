! Program: schrodinger
! By: Noah Egger
!-----------------------------------------------------------------------------
! This program solves the time dependent schrodinger equation using the Crank -
! Nicolson method with varying potentials. Moreover, the program receives 
! specifications from a namelist file, or default values, that dictate the values
! for the length of the potential region, the number of sampling points in space,
! the number of points in time, the step size in time, the initial width of the
! wave functions distribution, the central maximum location, and the vallue of
! the harmonic oscillator spring constant. The output of this program will be
! files that contain the expectation values for position, width, and normalization
! as a function of time, the probability density as a function of space and time,
! and the analytic position and widths as a function of time.
!-----------------------------------------------------------------------------
program schrodinger 

use types
use read_write, only : read_input, write_time_evolution, write_expectation_values, &
& write_analytic_position
use quantum, only : sample_box, construct_initial_wavefunction, construct_time_evolution_matrix, &
    evolve_wave_function, expectation_values
use hamiltonian, only : kinetic_energy, construct_hamiltonian, sigma_array_analytic, &
& calculate_harmonic_v

implicit none

real(dp) :: length, delta_t, width, center, k_oscillator
integer :: n_points, n_steps
character(len=1024) :: time_file, density_file, position_file
real(dp), allocatable :: x_vector(:) !will be of size n_points.
real(dp), allocatable :: wave_function(:)! will be of size 2*n_points.
real(dp), allocatable :: evolution_matrix(:,:) !will be of size 2*n_points by 2*n_points
real(dp), allocatable :: time_wave_function(:,:) !will be of size n_points by n_steps + 1 (the +1 is so that you can store the t=0 value)
real(dp), allocatable :: norm(:), postion(:), sigma(:), sigma_analytic(:) !all of size n_steps + 1 
real(dp), allocatable :: ke_diag(:), ke_offdiag(:), hamiltonian(:,:)
real(dp), allocatable :: position(:), sigma_array(:)
real(dp), allocatable :: potential(:)

call read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator, &
    & time_file, density_file, position_file)

! Time evolution of the initial wave function. 
allocate(x_vector(1:n_points))
allocate(evolution_matrix(1:2*n_points,1:2*n_points))
allocate(norm(1:n_steps+1))
allocate(position(1:n_steps+1))
allocate(sigma_array(1:n_steps+1))
allocate(sigma_analytic(1:n_steps+1))
allocate(potential(1:n_points))

call sample_box(length, n_points, x_vector) 

call construct_initial_wavefunction(x_vector, width, center, n_points, wave_function)

call kinetic_energy(n_points, length, ke_diag, ke_offdiag )

call calculate_harmonic_v(x_vector, n_points, k_oscillator, potential)

call construct_hamiltonian(n_points, potential, ke_diag, ke_offdiag, hamiltonian)

call construct_time_evolution_matrix(delta_t, n_points, hamiltonian, evolution_matrix)

call evolve_wave_function(wave_function, n_points, n_steps, evolution_matrix, time_wave_function)

call expectation_values(time_wave_function, n_steps, length, n_points, x_vector, norm, position, sigma_array)

call write_time_evolution(density_file, x_vector, n_steps, n_points, time_wave_function)

call sigma_array_analytic(delta_t, n_steps, width, sigma_analytic)

if(k_oscillator /=0) then
    call write_analytic_position(position_file, delta_t, center, n_steps, k_oscillator)
endif

call write_expectation_values(time_file, sigma_analytic, position, norm, sigma_array, n_steps, delta_t)

end program schrodinger