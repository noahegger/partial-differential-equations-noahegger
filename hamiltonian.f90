!!-----------------------------------------------------------------------
!Module: hamiltonian
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module houses the subroutines necessary for constructing the kinetic
!! energy matrix as well as the potential energy matrix. One must construct
!! both the diagonal and off-diagonal terms for both matrices, as well
!! as for the distinct potentials. More specifically, the potential energy 
!! will be a diagonal matrix and the kinetic energy will be a tridiagonal matrix, 
!! thus making the hamiltonian a tridiagonal matrix as well. Additionally, this 
!! module holds subroutines to compute the normalization constant of the wave
!! function for each time step, as well as subroutines for computing the standard 
!! deviation in position both numerically and analytically.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_harmonic_v
!! kinetic_energy
!! construct_hamiltonian
!! normal_array
!! position_array
!! sigma_array
!! sigma_array_analytic
!! 
!!----------------------------------------------------------------------
module hamiltonian
use types
implicit none

! Since more than one subroutine in this module will use the value for
! hbar and mass, it would be a good idea to define them here as 
! parameters

real(dp), parameter :: h_bar = 1
real(dp), parameter :: mass = 1


private
public :: calculate_harmonic_v, kinetic_energy, construct_hamiltonian, normal_array, &
& position_array, sigma_array, sigma_array_analytic

contains

!-----------------------------------------------------------------------
!! Subroutine: calculate_harmonic_v
!-----------------------------------------------------------------------
!! By: Noah 
!!
!! This subroutine constructs the diagonal elements for the harmonic 
!! oscillator potential defined as V = 0.5*k*x^2 where h is hbar and 
!! x is the position we are sampling for the potential. This 1-D array is 
!! then used to construct the full matrix. 
!!----------------------------------------------------------------------
!! Input:
!! 
!! x_vector(:)      real        array containing all sampling points between -L and +L
!! n_points         integer     number of sampling points
!! k_oscillator     real        "spring" constant for harmonic potential
!-----------------------------------------------------------------------
!! Output:
!!
!! potential(:)     real        array containing the potential evaluated at corresponding x
!!
!-----------------------------------------------------------------------

subroutine calculate_harmonic_v(x_vector, n_points, k_oscillator, potential)
    implicit none
    real(dp), intent(out), allocatable :: potential(:)
    real(dp), intent(in) :: x_vector(:), k_oscillator
    integer, intent(in) :: n_points
    integer :: i
	if(allocated(potential)) deallocate(potential)
    allocate(potential(1:n_points))

    ! Fill the potential array which contains the diagonal values
    do i = 1, n_points
        potential(i) = 0.5*k_oscillator*x_vector(i)**2
    enddo
end subroutine calculate_harmonic_v

!-----------------------------------------------------------------------
!! Subroutine: kinetic_energy
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine constructs and fills the kinetic energy diagonal and
!! off diagonal arrays for use in the qm_solver.f90 module, in which the 
!! hamiltonians are constructed and eigenvectors and eigenvalues are found. 
!!----------------------------------------------------------------------
!! Input:
!!
!! length           real        half length of the potential well
!! n_points         integer     number of sampling points
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! ke_diag          real        array containing kinetic energy diagonal terms
!! ke_offdiag       real        array containing kinetic energy off diagonal terms
!!
!-----------------------------------------------------------------------

subroutine kinetic_energy(n_points, length, ke_diag, ke_offdiag)
    implicit none
    real(dp), intent(in) :: length
    integer, intent(in) :: n_points
    real(dp), allocatable, intent(out) :: ke_diag(:), ke_offdiag(:)
    real(dp) :: dx
    integer :: i, j
    if(allocated(ke_diag)) deallocate(ke_diag)
    if(allocated(ke_offdiag)) deallocate(ke_offdiag)
    allocate(ke_diag(1:n_points))
    allocate(ke_offdiag(1:n_points-1))

    ! Step size
    dx = 2._dp*length/(n_points-1._dp)

    ! Construct kinetic energy diagonal
    do i = 1, n_points
        ke_diag(i) = (h_bar**2)/mass/(dx**2)
    enddo

    ! Construct kinetic energy off diagonal
    do j = 1, n_points-1
        ke_offdiag(j) = -0.5_dp*(h_bar**2)/mass/(dx**2)
    enddo

    ! Arrays will called by another subroutine to solve TDSE and find
    ! eigenvalues and eigenvectors 

end subroutine kinetic_energy

!-----------------------------------------------------------------------
!! Subroutine: construct hamiltonian
!-----------------------------------------------------------------------
!! By: Noah 
!!
!! This subroutine combines the potential and kinetic energy arrays to
!! construct the hamiltonian matrix. 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector(:)         real        1D array containing sampling locations
!! ke_diag (:)         real        1D of kinetic energy diagonal terms
!! ke_offdiag(:)       real        1D array of kinetic energy off-diagonal terms
!! n_points(:)         integer     number of sampling points
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! hamiltonian(:)     real        2D hamiltonian matrix
!!
!-----------------------------------------------------------------------

subroutine construct_hamiltonian(n_points, potential, ke_diag, ke_offdiag, hamiltonian)
    implicit none
    real(dp), intent(in) :: ke_diag(:), ke_offdiag(:), potential(:)
    integer, intent(in) :: n_points
    real(dp), intent(out), allocatable :: hamiltonian(:,:)
    integer :: i, j, k
    if(allocated(hamiltonian)) deallocate(hamiltonian)
    allocate(hamiltonian(1:n_points, 1:n_points))

    ! Construct diagonal terms
    do i = 1, n_points
        hamiltonian(i,i) = ke_diag(i) + potential(i)
    enddo

    ! Construct upper off-diagonal terms
    do j = 1, n_points-1
        hamiltonian(j,j+1) = ke_offdiag(1)
    enddo

    ! Construct lower off-diagonal terms
    do k = 2, n_points
        hamiltonian(k,k-1) = ke_offdiag(1)
    enddo

end subroutine construct_hamiltonian

!-----------------------------------------------------------------------
!! Subroutine: normal_array
!-----------------------------------------------------------------------
!! By: Noah 
!!
!! This subroutine constructs a 1D array that contains the normalization
!! constant at each time step for the wave function. This is computed 
!! using the probability density at each point in space and time.
!!----------------------------------------------------------------------
!! Input:
!!
!! time_wave_function(:,:)       real        2-D array containing wave function(x,t)
!! length                        real        half length of the potential well
!! n_points                      integer     number of sampling points
!! n_steps                       integer     number of time steps
!! prob_dens                     real        2D array containing probability density(x,t) 
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! norm         real        1D array containing normalization constant for each time step
!!
!-----------------------------------------------------------------------

subroutine normal_array(prob_dens, time_wave_function, n_steps, length, n_points, norm)
    implicit none
    real(dp), intent(in) :: time_wave_function(:,:), length, prob_dens(:,:)
    integer, intent(in) :: n_points, n_steps
    real(dp), intent(out), allocatable :: norm(:)
    real(dp) :: sum, dx
    integer :: i, j

    if(allocated(norm)) deallocate(norm)
    allocate(norm(1:n_steps+1))

    ! Step size
    dx = 2._dp*length/(n_points-1._dp)

    ! j indexes columns

    do j = 1, n_steps + 1
        sum = 0._dp

        ! i indexes rows
        ! For a given column, sum all the rows to get the normalization
        ! for a given time step
        do i = 1, n_points
            sum = sum + prob_dens(i,j)*dx
        enddo

        ! End results is an array of normalization constants for each time step in j
        norm(j) = sum
    enddo
end subroutine normal_array

!-----------------------------------------------------------------------
!! Subroutine: position_array
!-----------------------------------------------------------------------
!! By: Noah 
!!
!! This subroutine constucts a 1D array to store the position values for
!! the wafe function at each time step. The position is calculated as 
!! the normalized expectation value at each time step. 
!!----------------------------------------------------------------------
!! Input:
!!
!! length               real        half length of potential region
!! n_points             integer     number of sampling 
!! n_steps              integer     number of time steps 
!! time_wave_function   real        2D array containing the wave function at each position and time step
!! normalization        real        1D array containing the normalization constant for each time step  
!! prob_dens            real        array containing the probability density(x,t) 
!! x_vector             real        array containing the sampling locations
!! 
!!----------------------------------------------------------------------
!! Output:
!!
!! position     real        1D array containing expectation value of position at each time step 
!!
!!----------------------------------------------------------------------
subroutine position_array(prob_dens, time_wave_function, n_steps, length, &
    & n_points, x_vector, norm, position)
    implicit none
    real(dp), intent(in) :: time_wave_function(:,:), x_vector(:), length, prob_dens(:,:)
    real(dp), intent(in) :: norm(:)
    integer, intent(in) :: n_steps, n_points
    real(dp), intent(out), allocatable :: position(:)
    real(dp), allocatable :: position_numerator(:)
    real(dp) :: sum, dx
    integer :: i, j, k 

    if(allocated(position)) deallocate(position)
    allocate(position(1:n_steps+1))

    allocate(position_numerator(1:n_steps+1))

    dx = 2._dp*length/(n_points - 1._dp)

    do i = 1, n_steps + 1
        sum = 0._dp

        do j = 1, n_points
            sum = sum + x_vector(j)*prob_dens(j,i)
        enddo

        position_numerator(i) = sum*dx
    enddo

    do k = 1, n_steps + 1
        position(k) = position_numerator(k)/norm(k)
    enddo

end subroutine position_array

!-----------------------------------------------------------------------
!! Subroutine: sigma_array
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine constructs a 1D array containing the standard deviation,
!! otherwise known as the width of the distribution, at each time step. 
!! This is constructed from sqrt(<x^2> - <x>^2). Since we already have
!! <x>, all that is needed is to construct <x^2> and compute the width.
!!----------------------------------------------------------------------
!! Input:
!!
!! length               real            half length of potential region
!! n_points             integer         number of sampling points across the region
!! n_steps              integer         number of time steps 
!! time_wave_function   real            2D array containing the wave function(x,t)
!! norm                 real            1D array containing the normalization constant for each time step  
!! prob_dist            real            2D array containing the probability density (x,t) 
!! x_vector             real            1D array of sampling locations
!! position             real            1D array containing position at each time step 
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! sigma      real        1D array containing the width of the gaussian distribution at each time step
!!
!!----------------------------------------------------------------------
subroutine sigma_array(prob_dens, time_wave_function, n_steps, length, n_points, &
& x_vector, norm, position, sigma)
    implicit none
    real(dp), intent(in) :: time_wave_function(:,:), x_vector(:), length, position(:)
    real(dp), intent(in) :: norm(:), prob_dens(:,:)
    integer, intent(in) :: n_steps, n_points
    real(dp), allocatable, intent(out) :: sigma(:)
    real(dp), allocatable :: position_squared(:), position_squared_numerator(:)
    real(dp) :: sum, dx
    integer :: i, j, k, l

    if(allocated(sigma)) deallocate(sigma)


    allocate(position_squared(1:n_steps+1))
    allocate(position_squared_numerator(1:n_steps+1))
    allocate(sigma(1:n_steps+1))

    ! Step size
     dx = 2._dp*length/(n_points-1._dp)

    ! Construct <x^2> = sum_{j} x_j^2*rho_j*dx/(sum_{i} rho_i*dx)
    do i = 1,n_steps + 1
        sum = 0._dp

        do j = 1, n_points

        sum = sum + (x_vector(j)**2)*prob_dens(j,i)
        enddo

        ! Put numerator together separately
        position_squared_numerator(i) = sum*dx
    enddo
 
        do k = 1, n_steps + 1
            ! Divide by normalization to get <x^2>
            position_squared(k) = position_squared_numerator(k)/norm(k)
        enddo

    do i = 1, n_steps + 1
        ! Obtain the standard deviation
        sigma(i) = sqrt(position_squared(i)-position(i)**2)
    enddo


end subroutine sigma_array

!-----------------------------------------------------------------------
!! Subroutine: sigma_array_analytic
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine computes the standard deviation of the position 
!! analytically. Each element of the 1D corresponds to the value of the 
!! wdith at that time step. Ultimately, it will be used as comparison to 
!! our numerical calculation. 
!!----------------------------------------------------------------------
!! Input:
!!
!! n_steps      integer         number of time steps
!! delta_t      real            time step size
!! width        real            width of the gaussian distribution at t=0
!! 
!!----------------------------------------------------------------------
!! Output:
!!
!! sigma_analytic   real        1D array containing analytic width of gaussian at each time step 
!!
!!----------------------------------------------------------------------
subroutine sigma_array_analytic(delta_t, n_steps, width, sigma_analytic)
    implicit none
    real(dp), intent(in) :: delta_t, width
    integer, intent(in) :: n_steps
    real(dp), allocatable, intent(out) :: sigma_analytic(:)
    integer :: i
    real(dp) :: time

    allocate(sigma_analytic(1:n_steps+1))

    time = 0

    do i = 1,n_steps + 1

        sigma_analytic(i) = sqrt((width**2) + (h_bar**4)*(time**2)/(4*(mass**2)*width**2))

        time = time + delta_t
    enddo

end subroutine sigma_array_analytic


end module hamiltonian