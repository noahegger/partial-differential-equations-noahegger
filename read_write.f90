!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This module is responsible for reading a user inputted namelist file
!! that then assigns relevant parameters for the calculation. Moreover, 
!! this module contains routines to write the probability density as a
!! function of space and time to a file, the expectation values for each
!! point in time along with the time increment, and also the analytic 
!! poisition as a function of time.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! write_time_evolution
!! write_expectation_values
!! write_analytic_position
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module read_write
use types

implicit none

private
public :: read_input, write_time_evolution, write_expectation_values, &
& write_analytic_position

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Noah
!!
!! Receives user provided namelist file to specificy parameters used in 
!! the program calculations. If namelist file is not provided, default 
!! values are given.
!!----------------------------------------------------------------------
!! Input:
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! length       real            half length of potential region
!! n_points     integer         number of sampling points 
!! n_steps      integer         number of time steps 
!! delta_t      real            time step size
!! width        real            width of the starting wave function gaussian distribution
!! center       real            central maximum of the initial wave function gaussian distribution
!! k_oscillator real            "spring" constant for harmonic oscillator potential
!!
!!----------------------------------------------------------------------
subroutine read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator &
    , time_file, density_file, position_file)
    implicit none
    real(dp), intent(out) :: length, delta_t, width, center, k_oscillator
    integer, intent(out) :: n_points, n_steps
    character(len=*) :: time_file, density_file, position_file
    character(len=1024) :: namelist_file 
    logical :: file_exists
    integer :: unit, ierror, n_arguments

    ! Subsection names of namelist file
    ! names of parameters
    namelist /integration/ length, n_points, n_steps, delta_t
    namelist /wave_function/ width, center
    namelist /oscillator/ k_oscillator
    namelist /output/ time_file, density_file, position_file

    ! Default values
    length = 5._dp
    n_points = 100
    n_steps = 100
    delta_t = 0.05_dp
    width = 0.5_dp
    center = 0._dp
    k_oscillator = 0.0_dp
    time_file = 'time_results.dat'
    density_file = 'density_results.dat'

    n_arguments = command_argument_count()
    if (n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        inquire(file=trim(namelist_file), exist = file_exists)
        if (file_exists) then
            open(newunit = unit, file=namelist_file)
            read(unit, nml = integration, iostat = ierror)
            if(ierror /= 0) then
                print*, 'error reading integration namelist'
                stop
            endif
            read(unit, nml = wave_function, iostat = ierror)
            if(ierror /= 0) then
                print*, 'error reading wave_function namelist'
                stop
            endif
            read(unit, nml = oscillator, iostat = ierror)
            if(ierror /= 0) then
                print*, 'error reading oscillator namelist'
                stop
            endif
            read(unit, nml = output, iostat = ierror)
            if(ierror /= 0) then
                print*, 'error reading output namelist'
                stop
            endif
            close(unit)
        else
    ! If namelist doesn't exist, end program
            print*, 'Argument, ', trim(namelist_file)
            print*, 'does not exist. Ending program'
            stop
        endif
    ! If there isnt exactly 0 or 1 arguments after executable is typed into 
    ! terminal by user, print error and stop program. 
    else if(n_arguments /= 0 ) then
        print*, 'Incorrect number of arguments'
        print*, 'Program takes either 0 or 1 arguments'
        print*, 'See documentation in README.md'
        stop
    endif

end subroutine read_input


!-----------------------------------------------------------------------
!! Subroutine: write_time_evolution
!-----------------------------------------------------------------------
!! By: Noah 
!!
!! This subroutine writes the probability density as a function of space
!! and time to a file.
!!----------------------------------------------------------------------
!! Input:
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! n_steps              integer         number of time steps 
!! n_points             integer         number of sampling points
!! x_vector             real            1D array containing sampling points 
!! time_wave_function   real            2D array containing wave function (x,t)\
!! density_file         string          file name for probability density
!!----------------------------------------------------------------------
subroutine write_time_evolution(density_file, x_vector, n_steps, n_points, time_wave_function)
    implicit none
    character(len=*), intent(in) :: density_file
    integer, intent(in) :: n_steps, n_points
    real(dp), intent(in) :: x_vector(:), time_wave_function(:,:)
    real(dp), allocatable :: prob_dens(:,:)
    integer :: unit, i, j

    allocate(prob_dens(1:n_points,1:n_steps + 1))
    !This subroutine should write to the density_file file the probability 
    !density at different times. The first LINE should contain the sample 
    !points along the x axis.

    !The successive lines should contain the probability density at  
    !different time steps.

    do j = 1, n_steps + 1
        do i = 1, n_points
            prob_dens(i,j) = time_wave_function(i,j)**2 + time_wave_function(i+n_points, j)**2
        enddo
    enddo

    open(newunit=unit,file=trim(density_file))
    ! First line is sample points
    write(unit, *) x_vector

    ! subsequent lines are the probability density at different time steps
    ! time is incremented column to column
    ! space is incremented row to row
    do i=1,n_steps+1
        write(unit,*) prob_dens(:,i)
    enddo
    close(unit)
    
end subroutine write_time_evolution

!-----------------------------------------------------------------------
!! Subroutine: write_expectation_values
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine writes to the time_file file the expectation values 
!! as a function time. The first COLUMN contains the times at which the 
!! wave function was calculated .The successive columns contain the 
!! expectation values (normalization, position, width) at the respective 
!! times.
!!----------------------------------------------------------------------
!! Input:
!!
!! n_steps          integer         number of time steps 
!! delta_t          real            size of the time step
!! time_file        string          name of the file results are written to
!! sigma_analytic   real            1D array containing analytic width at each time step
!! position         real            1D array containing expectationv value of position at each time step
!! norm             real            1D array containing normalization constant at each time step
!! width_array      real            1D array containing numerical width at each time step
!!----------------------------------------------------------------------
!! Output:
!!
!!----------------------------------------------------------------------
subroutine write_expectation_values(time_file, sigma_analytic, position, norm, sigma_array, n_steps, delta_t)
    implicit none
    character(len=*), intent(in) :: time_file
    real(dp), intent(in) :: sigma_analytic(:), position(:), norm(:)
    real(dp), intent(in) :: sigma_array(:), delta_t
    real(dp), allocatable :: time(:)
    integer, intent(in) :: n_steps
    integer :: i, unit1

    if(allocated(time)) deallocate(time)
    allocate(time(1:n_steps+1))
    ! Create array to store time values at each step
    do i = 1, n_steps+1
        if (i == 1) then
            time(i) = 0._dp
        else 
            time(i) = time(i - 1) + delta_t
        end if
    enddo

    ! Create file to write normalization, width, position for each time step to.
    open(newunit=unit1, file=trim(time_file))
    write(unit1, '(4a28)') 'time', 'position', 'normalization', 'sigma_array'! 'sigma_analytic'
    do i = 1,n_steps+1
        write(unit1, *) time(i), position(i), norm(i), sigma_array(i), sigma_analytic(i)
    enddo
    close(unit1)
    
end subroutine write_expectation_values

!-----------------------------------------------------------------------
!! Subroutine: write_expectation_values
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine writes to a file the analytic position at each time 
!! step. 
!!----------------------------------------------------------------------
!! Input:
!!
!! n_steps          integer         number of time steps
!! delta_t          real            size of the time steps f
!! center           real            central maximum of the initial wave function gaussian distribution
!! position_file    string          file name
!! k_oscillator     real            spring constant
!!----------------------------------------------------------------------
!! Output:
!!
!!----------------------------------------------------------------------
subroutine write_analytic_position(position_file, delta_t, center, n_steps, k_oscillator)
    implicit none
    character(len=*), intent(in) :: position_file
    real(dp), intent(in) :: delta_t, center, k_oscillator
    integer, intent(in) :: n_steps
    real(dp), allocatable :: time(:)
    integer :: i, unit

    if(allocated(time)) deallocate(time)
    allocate(time(1:n_steps+1))
    open(newunit=unit, file=trim(position_file))
    write(unit, '(4a28)') 'time', 'analytic position'

    do i=1,n_steps+1
        if (i == 1) then
            time(i) = 0._dp
        else 
            time(i) = time(i - 1) + delta_t
        end if
        write(unit, *) time(i), center*cos(time(i)*sqrt(k_oscillator))
    enddo
    close(unit)
    
end subroutine write_analytic_position

end module read_write