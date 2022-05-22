#include "logger.h"
#define MAX_CHAR_LENGTH 1024
module types_mod
    use intrinsics_mod, only : dp
    implicit none

    private

    public HartreeFock_t, data_from_relcode_dir, Diag_t, &
            Matrix_Parameters_t, Channel_Indices_t, Bsplines_t, Channel_Index_t

    character(len = *), parameter :: data_from_relcode_dir = "./data_from_relcode/"

    ! -----

    type :: Diag_t
        integer :: system_size
        complex(dp), allocatable, dimension(:) :: eigenvalues
        complex(dp), allocatable, dimension(:) :: matrix_elements
        complex(dp), allocatable, dimension(:, :) :: eigenvectors
    end type Diag_t

    ! -----
    ! Below are types "copied" from the relcode program.
    ! -----

    type :: HartreeFock_t
        integer :: kappa_list_size, size_eig, size_rhs
        complex(dp), allocatable, dimension(:, :, :) :: eigenstates
        complex(dp), allocatable, dimension(:, :) :: eigenvalues
    end type HartreeFock_t

    ! -----

    type Bsplines_t
        integer :: size1_large, size2_large
        integer :: size1_small, size2_small
        integer :: radial_pos_vector_size
        complex*16, allocatable, dimension(:, :) :: B_large
        !complex*16, allocatable, dimension(:, :) :: dB_large
        complex*16, allocatable, dimension(:, :) :: B_small
        !complex*16, allocatable, dimension(:, :) :: dB_small
        complex*16, allocatable, dimension(:) :: radial_pos_vector
    end type Bsplines_t

    ! -----

    type Matrix_Parameters_t
        integer :: startorb, orbmin, system_size
        integer, dimension(:, :), allocatable :: qn_indices
        complex*16, dimension(:), allocatable :: energies
    end type Matrix_Parameters_t

    ! -----

    type Channel_Index_t
        integer :: hole_kappa, hole_l, hole_j
        integer :: final_kappa
        integer :: start_idx, end_idx
		complex(dp) :: hole_binding_energy
        character(len = 1024) :: print_string
    contains
        procedure, private :: channel_index_write_unformatted
        procedure, private :: channel_index_read_unformatted
        generic :: write(unformatted) => channel_index_write_unformatted
        generic :: read(unformatted) => channel_index_read_unformatted
    end type Channel_Index_t

    ! -----

    type Channel_Indices_t
        integer :: num_channels
        type(Channel_Index_t), dimension(:), allocatable :: channel
    contains
        procedure, private :: channel_indices_write_unformatted
        procedure, private :: channel_indices_read_unformatted
        generic :: write(unformatted) => channel_indices_write_unformatted
        generic :: read(unformatted) => channel_indices_read_unformatted
    end type Channel_Indices_t

    ! -----

contains

    ! ============================================================================

    subroutine channel_index_write_unformatted(dtv, unit, iostat, iomsg)
        class(Channel_Index_t), intent(in) :: dtv
        integer, intent(in) :: unit
        integer, intent(out) :: iostat
        character(len = *), intent(inout) :: iomsg

        write(unit, iostat = iostat, iomsg = iomsg) dtv%hole_kappa, dtv%hole_l, dtv%hole_j
        write(unit, iostat = iostat, iomsg = iomsg) dtv%final_kappa
        write(unit, iostat = iostat, iomsg = iomsg) dtv%start_idx, dtv%end_idx
		write(unit, iostat = iostat, iomsg = iomsg) dtv%hole_binding_energy
        write(unit, iostat = iostat, iomsg = iomsg) dtv%print_string

    end subroutine channel_index_write_unformatted

    ! ============================================================================

    subroutine channel_index_read_unformatted(dtv, unit, iostat, iomsg)
        class(Channel_Index_t), intent(inout) :: dtv
        integer, intent(in) :: unit
        integer, intent(out) :: iostat
        character(len = *), intent(inout) :: iomsg

        read(unit, iostat = iostat, iomsg = iomsg) dtv%hole_kappa, dtv%hole_l, dtv%hole_j
        read(unit, iostat = iostat, iomsg = iomsg) dtv%final_kappa
        read(unit, iostat = iostat, iomsg = iomsg) dtv%start_idx, dtv%end_idx
		read(unit, iostat = iostat, iomsg = iomsg) dtv%hole_binding_energy
        read(unit, iostat = iostat, iomsg = iomsg) dtv%print_string

    end subroutine channel_index_read_unformatted

    ! ============================================================================

    subroutine channel_indices_write_unformatted(dtv, unit, iostat, iomsg)
        class(Channel_Indices_t), intent(in) :: dtv
        integer, intent(in) :: unit
        integer, intent(out) :: iostat
        character(len = *), intent(inout) :: iomsg

        integer :: i

        write(unit, iostat = iostat, iomsg = iomsg) dtv%num_channels

        do i = 1, dtv%num_channels
            write(unit, iostat = iostat, iomsg = iomsg) dtv%channel(i)
        end do

    end subroutine channel_indices_write_unformatted

    ! ============================================================================

    subroutine channel_indices_read_unformatted(dtv, unit, iostat, iomsg)
        class(Channel_Indices_t), intent(inout) :: dtv
        integer, intent(in) :: unit
        integer, intent(out) :: iostat
        character(len = *), intent(inout) :: iomsg

        integer :: i

        read(unit, iostat = iostat, iomsg = iomsg) dtv%num_channels
        if(allocated(dtv%channel)) then
            LOG_FATAL("Channel_Indices_t%channel already allocated when trying to read binary!")
        end if
        allocate(dtv%channel(dtv%num_channels))
        do i = 1, dtv%num_channels
            read(unit, iostat = iostat, iomsg = iomsg) dtv%channel(i)
        end do

    end subroutine channel_indices_read_unformatted

end module types_mod