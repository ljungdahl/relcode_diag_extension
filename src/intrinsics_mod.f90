! This module imports intrinsic values, ie precision intrinsics, for use
! throughout the program
module intrinsics_mod
    use, intrinsic :: iso_fortran_env, only: REAL64, REAL32
    implicit none

    private

    public dp, sp

    integer, parameter :: dp = REAL64
    integer, parameter :: sp = REAL32

end module intrinsics_mod