!*****************************************************************************************
!> author: Jacob Williams
!
!  Real and integer kinds used by the DAQP package.
!
!  The real kind is selected by a preprocessor flag:
!  `-DREAL32`, `-DREAL64` (the default), or `-DREAL128`.

    module daqp_kinds

    use iso_fortran_env, only: real32, real64, real128, int32

    implicit none

    private

#ifdef REAL32
    integer, parameter, public :: daqp_wp = real32   !! real kind used by the package [4 bytes]
#elif REAL128
    integer, parameter, public :: daqp_wp = real128  !! real kind used by the package [16 bytes]
#else
    integer, parameter, public :: daqp_wp = real64   !! real kind used by the package [8 bytes]
#endif

    integer, parameter, public :: daqp_ip = int32    !! integer kind used by the package

    end module daqp_kinds
!*****************************************************************************************
