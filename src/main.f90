#include "logger.h"

program main
    use intrinsics_mod, only: dp
    implicit none

    complex(dp) :: zomplex
    complex(dp), pointer :: pointer_3x3(:, :)

    integer :: N

    N = 3

    zomplex = (1.0_dp, 2.0_dp)

    LOG_WRITE zomplex

    call test_in_place_reshape(pointer_3x3)

    LOG_WRITE "after in place reshape, supposedly destroyed"
    do j = 1, N
        do i = 1, N
            LOG_WRITE i, j, pointer_3x3(i,j)
        end do
    end do


end program main

! -------------------------------------------------

subroutine test_in_place_reshape(pointer_3x3)
    ! This routine tests the use of C_F_POINTER to appropriately
    ! reshape an array in place using a pointer.
    use intrinsics_mod, only: dp
    use, intrinsic :: iso_c_binding, only: C_LOC, C_F_POINTER
    implicit none

    complex(dp), pointer, intent(out) :: pointer_3x3(:, :)


    complex(dp), allocatable, target :: raw_array(:)

    integer :: N, i, j, idx

    N = 3
    allocate(raw_array(N*N))

    do j = 1, N
        do i = 1, N
            idx = j+N*(i-1)
            raw_array(idx) = cmplx(idx, j, kind=dp)
        end do
    end do

    do i=1,size(raw_array)
        LOG_WRITE i, raw_array(i)
    end do

    call C_F_POINTER(C_LOC(raw_array), pointer_3x3, (/3, 3/))

    do j = 1, N
        do i = 1, N
            LOG_WRITE i, j, pointer_3x3(i,j)
        end do
    end do

    deallocate(raw_array)


end subroutine test_in_place_reshape