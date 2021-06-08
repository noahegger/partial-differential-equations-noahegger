!----------------------------------------------------------------------
!Module: types
!---------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Defines the integer parameters sp, dp and qp to
!! be used as kinds to define real variables as single precision 
!! (32-bit), double precision (64-bit) and quadruple precision 
!! (128-bit) (depending on the machine and compiler, this last one may
!! not always be available).
!!
!! The intrinsic module iso_fortran_env (Fortran 2008 and later) is
!! used.
!!
!! More types may be defined later (i.e. larger integers)
!!
!! Fundamental constants (i.e. &pi;, e, ...) may also be defined here if
!! desired. This is a good place to define this type of constants as
!! all modules and the main program will (in principle) use this module
!---------------------------------------------------------------------
module types

use iso_fortran_env

implicit none

integer, parameter :: sp = REAL32 !< single precision kind
integer, parameter :: dp = REAL64 !< double precision kind
integer, parameter :: qp = REAL128!< quadruple precision kind
real(dp), parameter :: pi=acos(-1._dp)!< &pi; = 3.141592...

end module types
