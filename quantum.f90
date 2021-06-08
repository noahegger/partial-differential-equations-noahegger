!-----------------------------------------------------------------------
!Module: quantum
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This module is responsible for creating the array that holds the sampling
!! points in space, the initial wave function array, the time evolution matrix,
!! the array that holds the wave function at all points in time and space,
!! as well as the arrays that hold the expectation values in position and
!! width as a function of time.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! sample_box
!! construct_initial_wavefunction
!! evolve_wave_function
!! construct_time_evolution_matrix
!! expectation_values
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module quantum

use types
use linear_algebra, only : invert_matrix, right_supermatrix, left_supermatrix, &
& eq_solver, matrix_mult
use hamiltonian, only : normal_array, position_array, sigma_array

implicit none
private
public sample_box, construct_initial_wavefunction, evolve_wave_function, construct_time_evolution_matrix, expectation_values

contains

!-----------------------------------------------------------------------
!! Subroutine: sample_box
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine constructs the 1D array of position sampling points,
!! calculated from user input 'length'.
!!----------------------------------------------------------------------
!! Input:
!!
!! length       real        half length of potential
!! n_points     integer     number of sampling points
!!----------------------------------------------------------------------
!! Output:
!!
!! x_vector(:)  real        1D array of position sampling points from -L to L
!!----------------------------------------------------------------------
subroutine sample_box(length, n_points, x_vector)
    implicit none
    integer, intent(in) :: n_points
    real(dp), intent(in) :: length
    real(dp), allocatable, intent(out) :: x_vector(:)
    real(dp) :: delta_x
    integer :: i

    if(allocated(x_vector)) deallocate(x_vector)
    allocate(x_vector(1:n_points))

    ! Step size in position
    delta_x = 2._dp*length/(n_points - 1)

    do i = 1, n_points
        if (i == 1) then
            x_vector(i) = -length
        else 
            x_vector(i) = x_vector(i - 1) + delta_x
        end if
    enddo
end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: construct_initial_wavefunction
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine evaluate the initial wave function at each sampling
!! point in space and stores the values in an array. The second half is
!! imaginary, which we fill with zeros. 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector         real        1D array containing the sampling points between -L and +L
!! n_points         integer     number of sampling points
!! center           real        central location of the gaussian distribution
!! width            real        width of the starting wave function gaussian distribution
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! wave_function    real        array containing the initial wave function (x,t=0)
!!----------------------------------------------------------------------
subroutine construct_initial_wavefunction(x_vector, width, center, n_points, wave_function)
    implicit none
    real(dp), intent(in) :: width, center, x_vector(:)
    real(dp), allocatable, intent(out) :: wave_function(:)
    integer, intent(in) :: n_points
    real(dp) :: coeff, denominator
    integer :: i, j

    if(allocated(wave_function)) deallocate(wave_function)
    allocate(wave_function(1:2*n_points))

    ! Break terms up to clean computation
    coeff = (2*pi*width**2)**(-0.25)
    denominator = -1/(4*width**2)

    ! Construct wave function array evaluated at each x position for first half
    do i = 1, n_points
        wave_function(i) = coeff*exp(denominator*(x_vector(i) - center)**2)
    enddo

    ! Set imaginary second half of wave function to 0 
    do j = n_points + 1, 2*n_points 
        wave_function(j) = 0
    enddo

end subroutine construct_initial_wavefunction

!-----------------------------------------------------------------------
!! Subroutine: construct_time_evolution_matrix
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine calls the routines needed to construct the left and 
!! right "super-matrices" as well as invert the left matrix and multiply
!! the two together. As a result, we obtain the time_evolution_matrix that,
!! when multiplied by the initial wave function, produces the wave function
!! at the next time step.
!!----------------------------------------------------------------------
!! Input:
!!
!! n_points             integer     number of sampling points
!! delta_t              real        size of time step
!! hamiltonian_mat      real        2D array representing the hamiltonian
!!----------------------------------------------------------------------
!! Output:
!!
!! evolution_matrix     real        2D array used for evolving wave function in time
!!----------------------------------------------------------------------

subroutine construct_time_evolution_matrix(delta_t, n_points, hamiltonian_mat, evolution_matrix)
    implicit none
    real(dp), intent(in) :: hamiltonian_mat(:,:)
    real(dp), allocatable, intent(out) :: evolution_matrix(:,:)
    real(dp), intent(in) :: delta_t
    integer, intent(in) :: n_points
    real(dp), allocatable :: left_inv(:,:), left(:,:), right(:,:)
    integer, allocatable :: identity(:,:)
    integer :: i, limit

    if(allocated(evolution_matrix)) deallocate(evolution_matrix)

    allocate(identity(1:n_points, 1:n_points))
    allocate(evolution_matrix(1:2*n_points, 1:2*n_points))
    allocate(right(1:2*n_points, 1:2*n_points))
    allocate(left(1:2*n_points, 1:2*n_points))
    allocate(left_inv(1:2*n_points, 1:2*n_points))

    ! Build identity matrix
    do i = 1, n_points
        identity(i,i) = 1
    enddo

    ! Build left "super-matrix"
    call left_supermatrix(delta_t, n_points, identity, left, hamiltonian_mat)

    ! Invert the left "super-matrix"
    call invert_matrix(left, left_inv)

    ! Build the right "super-matrix"
    call right_supermatrix(delta_t, n_points, identity, right, hamiltonian_mat)

    ! Multiply these matrices in order to get the evolution matrix
    call matrix_mult(left_inv, right, evolution_matrix)

end subroutine construct_time_evolution_matrix

!-----------------------------------------------------------------------
!! Subroutine: evolve_wave_function
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine utilizes the initial wave function at t=0 and multiplies
!! the evolution matrix in order to evolve the wave function to the next
!! time step. Each wave function snapshot in time is stored as a column as 
!! the subroutine gradually builds the 2D time wave function.
!!----------------------------------------------------------------------
!! Input:
!!
!! wave_function        real        1D array containing wave function defined at t=0
!! n_steps              integer     number of time steps
!! n_points             integer     number of smapling points
!! evolution_matrix     real        2D array containing matrix needed for evolving the wave function      
!!----------------------------------------------------------------------
!! Output:
!!
!! time_wave_function   real        2D array containing the wave function (x,t)
!!----------------------------------------------------------------------
subroutine evolve_wave_function(wave_function, n_points, n_steps, evolution_matrix, time_wave_function)
    implicit none
    real(dp), intent(in) :: wave_function(:), evolution_matrix(:,:)
    real(dp), allocatable, intent(out) :: time_wave_function(:,:)
    real(dp), allocatable :: initial_vector(:), final_vector(:)
    integer, intent(in) :: n_steps, n_points
    integer :: i, j

    if(allocated(time_wave_function)) deallocate(time_wave_function)
    allocate(time_wave_function(1:2*n_points, 1:n_steps+1))
    allocate(initial_vector(1:2*n_points))
    allocate(final_vector(1:2*n_points))

    ! First column is wave function at t = 0
    time_wave_function(:,1) = wave_function

    do i = 1, n_steps
        initial_vector = time_wave_function(:,i)

        ! Evolve wave function across time steps 
        ! and store value at subsequent positions
        call eq_solver(initial_vector, evolution_matrix, final_vector)
        time_wave_function(:,i+1) = final_vector
    enddo

end subroutine evolve_wave_function

!-----------------------------------------------------------------------
!! Subroutine: expectation_values
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroute performs the relevant calls to other subroutines in order
!! to construct the expectation value arrays for position and sigma.
!!----------------------------------------------------------------------
!! Input:
!!
!! time_wave_function       real        2D array containing wave function (x,t)
!! length                   real        half length of potential region
!! n_steps                  integer     number of time n_steps
!! n_points                 integer     number of sampling points in space
!! x_vector                 real        1D array containing sampling points
!!----------------------------------------------------------------------
!! Output:
!!
!! norm                     real        1D array containing normalization constant at each time step
!! position                 real        1D array containing position in space at each time step 
!! sigma                    real        1D array containing the width of the gaussian distribution at each time step
!!----------------------------------------------------------------------
subroutine expectation_values(time_wave_function, n_steps, length, n_points, x_vector, norm, position, sigma)
    implicit none
    real(dp), intent(in) :: time_wave_function(:,:), x_vector(:), length
    real(dp), allocatable, intent(out) :: norm(:), position(:), sigma(:)
    real(dp), allocatable :: prob_dens(:,:)
    integer, intent(in) :: n_steps, n_points
    integer :: i, j

    allocate(prob_dens(1:n_points, 1:n_steps+1))
    ! Construct probability density 
    do j = 1, n_steps + 1
        do i = 1, n_points
            prob_dens(i,j) = time_wave_function(i,j)**2 + time_wave_function(i+n_points,j)**2
        enddo
    enddo

    ! Fill normalization constant array
    call normal_array(prob_dens, time_wave_function, n_steps, length, n_points, norm)

    ! Fill expectation value in position array
    call position_array(prob_dens, time_wave_function, n_steps, length, n_points, x_vector, norm, position)

    ! Fill expectation value of sigma array
    call sigma_array(prob_dens, time_wave_function, n_steps, length, n_points, x_vector, norm, position, sigma)

end subroutine expectation_values
  
end module quantum