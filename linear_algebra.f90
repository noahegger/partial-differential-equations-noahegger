!-----------------------------------------------------------------------
!Module: linear_algebra
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module computes the inverse of the alpha matrix constructed within
!! the nuclear model module, then solves the system of equations to obtain
!! the parameters needed for the binding energy calculations. Included
!! are checks to guarantee the matrix is square as well as the row and column
!! sizes match the b vector and paramater vector sizes, respectively.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_linear_system
!! test_array_sizes
!! invert_matrix
!! ludcmp
!! lubksb
!! eq_solver
!!
!!----------------------------------------------------------------------
module linear_algebra
use types
implicit none

real(dp), parameter :: h_bar = 1
real(dp), parameter :: mass = 1

private
public :: left_supermatrix, right_supermatrix, invert_matrix, &
& matrix_mult, eq_solver
contains

!-----------------------------------------------------------------------
!! Subroutine: left_supermatrix
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine constructs the 2Nx2N "super-matrix" used for evolving
!! the wave function in time. The diagonal "terms" contain the identity
!! matrix while the off-diagonal "terms" contain multiples of the hamiltonian.
!!----------------------------------------------------------------------
!! Input:
!!
!! delta_t             real       time step size
!! n_points            integer    number of sampling points 
!! hamiltonian_mat     real(dp)   2D array containing the hamiltonian matrix 
!! identity            real(dp)   2D array containing the identity matrix
!! 
!!----------------------------------------------------------------------
!! Output:
!!
!! a_left   real    left "super-matrix"
!!
!!----------------------------------------------------------------------
subroutine left_supermatrix(delta_t, n_points, identity, a_left, hamiltonian_mat)
    implicit none
    integer, intent(in) :: identity(:,:)
    real(dp), intent(in) :: delta_t, hamiltonian_mat(:,:)
    integer, intent(in) :: n_points
    real(dp), allocatable, intent(out) :: a_left(:,:)
    integer :: i, j, k

    allocate(a_left(1:2*n_points,1:2*n_points))

    ! Construct top half
    do i = 1, n_points

        ! Top left is identity
        do j = 1, n_points
            a_left(i,j) = identity(i,j)
        enddo

        ! Top right is -delta_t/2 * hamiltonian
        do k = n_points + 1,2*n_points
            a_left(i,k) = (-delta_t/(2*h_bar) )*hamiltonian_mat(i,k-n_points)
        enddo
    enddo

    ! Construct bottom half
    do i = n_points+1, 2*n_points

        ! Botto left is delta_t/2 *hamiltonian
        do j = 1,n_points
            a_left(i,j) = (delta_t/(2*h_bar) )*hamiltonian_mat(i-n_points,j)
        enddo

        ! Bottom right is identity matrix
        do k = n_points+1, 2*n_points
            a_left(i,k) = identity(i-n_points,k-n_points)
        enddo
    enddo
end subroutine left_supermatrix

!-----------------------------------------------------------------------
!! Subroutine: right_supermatrix
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine constructs the 2Nx2N "super-matrix" used for evolving
!! the wave function in time. The diagonal "terms" contain the identity
!! matrix while the off-diagonal "terms" contain multiples of the hamiltonian.
!!----------------------------------------------------------------------
!! Input:
!!
!! delta_t             real       time step size
!! n_points            integer    number of sampling points 
!! hamiltonian_mat     real(dp)   2D array containing the hamiltonian matrix 
!! identity            real(dp)   2D array containing the identity matrix
!! 
!!----------------------------------------------------------------------
!! Output:
!!
!! a_right     real    right "super-matrix"
!!
!!----------------------------------------------------------------------
subroutine right_supermatrix(delta_t, n_points, identity, a_right, hamiltonian_mat)
    implicit none
    integer, intent(in) :: identity(:,:)
    real(dp), intent(in) :: delta_t, hamiltonian_mat(:,:)
    integer, intent(in) :: n_points
    real(dp), allocatable, intent(out) :: a_right(:,:)
    integer :: i, j, k

    allocate(a_right(1:2*n_points,1:2*n_points))

    ! Construct top half
    do i = 1, n_points

        ! Top left is identity
        do j = 1, n_points
            a_right(i,j) = identity(i,j)
        enddo

        ! Top right is delta_t/2 * hamiltonian
        do k = n_points + 1,2*n_points
            a_right(i,k) = (delta_t/(2*h_bar))*hamiltonian_mat(i,k-n_points)
        enddo
    enddo

    ! Construct bottom half
    do i = n_points+1, 2*n_points

        ! Botto left is -delta_t/2 *hamiltonian
        do j = 1,n_points
            a_right(i,j) = (-delta_t/(2*h_bar) )*hamiltonian_mat(i-n_points,j)
        enddo

        ! Bottom right is identity matrix
        do k = n_points+1, 2*n_points
            a_right(i,k) = identity(i-n_points,k-n_points)
        enddo
    enddo
    
end subroutine right_supermatrix

!-----------------------------------------------------------------------
!! Subroutine: invert_matrix
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Given a non singular matrix $a$, returns its inverse $a^{-1}$
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real    2D array containing the $a$ matrix
!!----------------------------------------------------------------------
!! Output:
!!
!! a_inv    real    2D array with the $a^{-1}$ matrix
!-----------------------------------------------------------------------
subroutine invert_matrix(a, a_inv)
    implicit none
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: a_inv(:,:)
    real(dp), allocatable :: a_work(:,:)
    integer :: shape_a(1:2), n, i
    real(dp) :: d
    integer, allocatable :: indx(:)

    ! This I'll give you for free. It's the LU decomposition we discussed in
    ! class with the following back-substitution
    
    allocate(a_work,mold=a)
    shape_a = shape(a)
    n = shape_a(1)
    allocate(indx(1:n))
    
    ! ludcmp destroys the input matrix a. In order to preserve a we will copy
    ! it into a work array that will be used in ludcmp
    a_work = a
    call ludcmp(a_work,indx,d)
    
    ! We construct a matrix that has orthogonal unit vectors as columns
    a_inv = 0._dp
    do i=1,n
        a_inv(i,i) = 1._dp
    enddo

    ! And then feed each column to the back-substitution routine
    do i = 1,n
        call lubksb(a_work,indx,a_inv(:,i))
    enddo
    ! This results in a_inv being the inverse of a
end subroutine invert_matrix

!-----------------------------------------------------------------------
!! Subroutine: matrix_mult
!-----------------------------------------------------------------------
!! By: Noah
!!
!! This subroutine performs matrix multiplication.
!!----------------------------------------------------------------------
!! Input:
!!
!! m1       real(dp)    left matrix to be multiplied
!! m2       real(dp)    right matrix to be multipled
!!----------------------------------------------------------------------
!! Output:
!! 
!! result      real(dp)    result of the matrix multiplication
!!
!-----------------------------------------------------------------------
subroutine matrix_mult(m1, m2, result)
implicit none
    real(dp), intent(in) :: m1(:,:), m2(:,:)
    real(dp), allocatable, intent(out) ::  result(:,:)
    integer :: col2, row1, col1, row2, i, j, k
    integer, allocatable :: m1_size(:), m2_size(:)

    ! allocate arrays to contain row/column number for input matrices
    allocate(m1_size(1:2))
    allocate(m2_size(1:2))

    ! Fill arrays with row and column number
    m1_size = shape(m1)
    m2_size = shape(m2)

    ! Make variables based on row and column numbers for both matrices 
    col1 = m1_size(2)
    col2 = m2_size(2)
    row1 = m1_size(1)
    row2 = m2_size(1)

    ! Ensure column count of left matrix matches row count of 
    ! right matrix. 
    if(col1 /= row2) then
        print*, 'Matrices wrong size for multiplication'
        print*, 'Matrix_mult subroutine in linear_algebra.f90'
        stop
    endif

    ! allocate matrix to be solution to matrix multiplication. 
    allocate(result(1:row1, 1:col2))

! Do loop to perform matrix multiplication.
result = 0

do i = 1, col2
    do j = 1, row1
        do k = 1, col1
            result(j, i) = result(j, i) + m1(j, k)*m2(k, i)
        enddo
    enddo
enddo

end subroutine matrix_mult

! The subroutines below were taken from numerical recipes and were slightly
! modified to work with double precision reals.

! Notice how much harder it is to understand what a code does when 
! explicit informative names are not used for the different variables
! and processes 

!-----------------------------------------------------------------------
!! Subroutine: ludcmp
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Adapted from numerical recipes subroutine.
!! Performs LU decomposition on a non singular matrix $a$.
!! The original $a$ matrix is destroyed as the LU decomposition is returned
!! in the same array
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real        2D array containing the $a$ matrix
!!----------------------------------------------------------------------
!! Output:
!!
!! a        real        2D array with LU decomposition of the $a$ matrix
!! indx     integer     1D array that records the row permutation effected by the partial pivoting
!! d        real        +1 or -1 depending on whether the number of row interchanges was even or odd, respectively
!-----------------------------------------------------------------------
subroutine ludcmp(a, indx, d)
    implicit none
    real(dp), intent(inout) :: a(:,:)
    integer, intent(out) :: indx(:)
    real(dp), intent(out) :: d
    integer :: n,i,imax,j,k
    real(dp) aamax,dum,sum
    real(dp), allocatable :: vv(:)
    n = size(indx)
    allocate(vv(1:n))
    d=1._dp
    do i=1,n
        aamax=0._dp
        do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0._dp) then
            print *, 'singular matrix in ludcmp'
            stop
        endif
        vv(i)=1._dp/aamax
    enddo

    do j=1,n
        do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
        enddo
        aamax=0._dp
        do i=j,n
            sum=a(i,j)
            do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
                imax=i
                aamax=dum
            endif
        enddo
        if (j.ne.imax)then
            do k=1,n
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
            enddo
            d=-d
            vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0._dp) a(j,j) = tiny(1._sp)
        if(j.ne.n)then
            dum=1._dp/a(j,j)
            do i=j+1,n
                a(i,j)=a(i,j)*dum
            enddo
        endif
    enddo
end subroutine ludcmp

!-----------------------------------------------------------------------
!! Subroutine: lubksb
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Adapted from numerical recipes subroutine.
!!
!! Performs back-substitution after a LU decomposition in order to solve the
!! linear system of equations $a \cdot x = b$. The $b$ vector is given in the b
!! array (which is destroyed) and the solution $x$ is returned in its place
!!----------------------------------------------------------------------
!! Input:
!!
!! a        real        2D array containing the LU decomposition $a$ matrix (as returned by ludecomp)
!! indx     integer     1D array with the record of the row permutation effected by the partial pivoting (as returned by ludecomp)
!! b        real        1D array containing the $b$ vector
!!----------------------------------------------------------------------
!! Output:
!! b        real        1D array containing the $x$ vector
!-----------------------------------------------------------------------
subroutine lubksb(a, indx, b)
    implicit none
    real(dp), intent(in) :: a(:,:)
    integer, intent(in) :: indx(:)
    real(dp), intent(inout) :: b(:)

    integer :: n
    integer :: i,ii,j,ll
    real(dp) :: sum

    n = size(b)
    ii=0
    do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            enddo
        else if (sum.ne.0.) then
            ii=i
        endif
        b(i)=sum
    enddo
    do i=n,1,-1
        sum=b(i)
        do j=i+1,n
            sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
    enddo
end subroutine lubksb

!-----------------------------------------------------------------------
!! Subroutine: eq_solver
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine uses a_inverse and b_vector to solve the system
!! of equations for x_vector. i.e $A x_vec = b_vec$ --> $x_vec = A^-1 b_vec$
!!----------------------------------------------------------------------
!! Input:
!!
!! a_inverse  real      2D array containing the inverse of the matrix "a".
!! b_vector   real      1D array containing the $b$ vector
!!----------------------------------------------------------------------
!! Output:
!! x_vector        real        1D array containing the $x$ vector, BE parameters
!-----------------------------------------------------------------------
subroutine eq_solver(b_vector, a_inverse, x_vector)
    implicit none
    real(dp), intent(in) :: a_inverse(:,:), b_vector(:)
    real(dp), intent(out) ::  x_vector(:)
    integer :: b_size, i, j
    integer, allocatable :: a_size(:)

    allocate(a_size(1:2))

    b_size=size(b_vector)
    a_size=shape(a_inverse)
    ! Do loop to perform matrix multiplication.
    ! i indexes through rows of a_inverse.
    do i = 1, a_size(2) 
        x_vector(i) = 0
        ! j indexes through columns of a_inverse and elements of b_vector
        do j = 1, b_size 
            ! To ensure efficiency, we choose a given i
            ! then loop through the column elements j of a_inverse
            ! for that given i multiped by the jth element of b_vector.
            ! Once the loop for the given i is finished, we move onto
            ! the next row and loop through j once again.
            x_vector(i) = x_vector(i) + a_inverse(i,j)*b_vector(j) 
        enddo
     
    enddo

end subroutine eq_solver

end module linear_algebra