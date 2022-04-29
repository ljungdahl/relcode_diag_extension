#include "logger.h"
! This program reads in binary (mostly) data from a relcode diagonalisation run.
! The idea is the following:
! With this one can form the perturbed wave function per ionisation channel to get, f. ex, the partial cross sections.
! In particular we can remove certain states (resonances) from the diagonalisation spectrum to investigate the
! contribution of these states to f.ex the partial cross section.
program main
    use intrinsics_mod, only : dp

    use in_data_mod, only : In_Metadata_t, read_metadata, check_if_subdir_exists, &
            read_complex_array_binary_rank1, read_complex_array_binary_rank2, read_complex_array_binary_rank3, &
            read_integer_array_binary_rank2

    use types_mod, only : HartreeFock_t, Diag_t, Matrix_Parameters_t, Channel_Indices_t, Bsplines_t

    use perturbed_wavefunction_mod, only: Perturbed_Wavefunction_t, init_pwf, &
            compute_bspline_coefficients

    use photons_mod, only : Photons_t, get_xuv_photons_linspace

    implicit none

    type(In_Metadata_t) :: diag_eigenvalues_metadata, diag_matrix_elements_metadata
    type(Diag_t) :: diag
    type(In_Metadata_t) :: hf_eigenstates_metadata, hf_eigenvalues_metadata ! NOTE(anton): Not sure if eigenvalues are used.
    type(HartreeFock_t) :: hf
    type(In_Metadata_t) :: matrix_params_qn_indices_metadata
    type(Matrix_Parameters_t) :: matrix_params
    type(In_Metadata_t) :: bspl_large_metadata ! Large component
    type(In_Metadata_t) :: bspl_small_metadata ! Small component
    type(In_Metadata_t) :: bspl_radial_pos_vec_metadata ! Vector of radial coords in which the bsplines are defined.
    type(Bsplines_t)    :: bsplines
    type(Channel_Indices_t) :: channel_indices ! This is read using type derived I/O, ie overload of the read/write operators!

    ! We have a list of perturbed wavefunctions, one for each channel, so that w ecan compute
    ! partial probability currents (rates) etc.
    type(Perturbed_Wavefunction_t), allocatable, dimension(:) :: list_of_pwfs

    type(Photons_t) :: photons
    real(dp) :: omega_XUV_au
    character(len = *), parameter :: filename_channel_indices = "./data_from_relcode/channel_indices.bin"

    integer :: i, j, number_of_channels, number_of_photons

    ! #################################################################################################################
    ! # Read data from relcode program
    ! #################################################################################################################
    ! NOTE(anton): The routines for reading data from relcoded assumes the data has been copied to
    ! a subdirectory to the directory from which this program is run, called "data_from_relcode"
    ! TODO(anton): All of this could be done much more easier and compact if we used "derived type I/O" which is
    ! supported by modern Fortran standards. But we would need to define the write/read operations on each derived type,
    ! both in this program and in relcode proper.
    call check_if_subdir_exists()

    ! The metadata files describe the binary data from relcode and is necessary to properly read the binary information
    ! into arrays in this program.
    ! Diag data
    call read_metadata("diag_eigenvalues", diag_eigenvalues_metadata)
    call read_metadata("diag_matrix_elements", diag_matrix_elements_metadata)

    if(diag_eigenvalues_metadata%sizes(1) /= diag_matrix_elements_metadata%sizes(1)) then
        LOG_WRITE "diag_eigenvalues size = ", diag_eigenvalues_metadata%sizes(1)
        LOG_WRITE "diag_matrix_elements size = ", diag_matrix_elements_metadata%sizes(1)
        LOG_FATAL("Eigenvalues and matrix elements read from relcode should have the same dimensions!")
    end if
    diag%system_size = diag_matrix_elements_metadata%sizes(1)

    allocate(diag%eigenvalues(diag_eigenvalues_metadata%sizes(1)))
    call read_complex_array_binary_rank1("diag_eigenvalues", diag_eigenvalues_metadata, diag%eigenvalues)

    allocate(diag%matrix_elements(diag_matrix_elements_metadata%sizes(1)))
    call read_complex_array_binary_rank1("diag_matrix_elements", diag_matrix_elements_metadata, diag%matrix_elements)

    ! Hartree-Fock data
    call read_metadata("hf_eigenstates", hf_eigenstates_metadata)
    call read_metadata("hf_eigenvalues", hf_eigenvalues_metadata)

    allocate(hf%eigenstates(hf_eigenstates_metadata%sizes(1), &
            hf_eigenstates_metadata%sizes(2), hf_eigenstates_metadata%sizes(3)) &
            )
    call read_complex_array_binary_rank3("hf_eigenstates", hf_eigenstates_metadata, hf%eigenstates)

    allocate(hf%eigenvalues(hf_eigenvalues_metadata%sizes(1), hf_eigenvalues_metadata%sizes(2)))
    call read_complex_array_binary_rank2("hf_eigenvalues", hf_eigenvalues_metadata, hf%eigenvalues)

    ! This is what the sizes are when written out from relcode.
    hf%size_rhs        = hf_eigenstates_metadata%sizes(1)
    hf%size_eig        = hf_eigenstates_metadata%sizes(2)
    hf%kappa_list_size = hf_eigenstates_metadata%sizes(3)

    ! Matrix params - "quantum number indices". This is a map for indices to be used when getting
    ! hartree-fock basis states.
    call read_metadata("matrix_params_qn_indices", matrix_params_qn_indices_metadata)

    allocate(matrix_params%qn_indices(matrix_params_qn_indices_metadata%sizes(1), &
            matrix_params_qn_indices_metadata%sizes(2))&
            )
    call read_integer_array_binary_rank2("matrix_params_qn_indices", &
            matrix_params_qn_indices_metadata, matrix_params%qn_indices)

    ! Bspline for the large and small radial components
    call read_metadata("bspl_large", bspl_large_metadata)
    bsplines%size1_large = bspl_large_metadata%sizes(1)
    bsplines%size2_large = bspl_large_metadata%sizes(2)
    allocate(bsplines%B_large(bsplines%size1_large, bsplines%size2_large))
    call read_complex_array_binary_rank2("bspl_large", bspl_large_metadata, bsplines%B_large)

    call read_metadata("bspl_small", bspl_small_metadata)
    bsplines%size1_small = bspl_small_metadata%sizes(1)
    bsplines%size2_small = bspl_small_metadata%sizes(2)
    allocate(bsplines%B_small(bsplines%size1_small, bsplines%size2_small))
    call read_complex_array_binary_rank2("bspl_small", bspl_small_metadata, bsplines%B_small)

    call read_metadata("bspl_radial_pos_vec", bspl_radial_pos_vec_metadata)
    bsplines%radial_pos_vector_size = bspl_radial_pos_vec_metadata%sizes(1)
    allocate(bsplines%radial_pos_vector(bsplines%radial_pos_vector_size))
    call read_complex_array_binary_rank1("bspl_radial_pos_vec", &
            bspl_radial_pos_vec_metadata, bsplines%radial_pos_vector)

    ! Channel indices:
    ! Here we actual use a working implementation of derived type file I/O!
    ! We can just read into the channel_indices instance. The compiler translates this for us
    ! according to our implementation in the derived type. See definitions of Channel_Indices_t and
    ! Channel_Index_t.
    open(123, file = filename_channel_indices, form = "unformatted")
    read(123) channel_indices
    close(123)

    !        ! Test derived type file I/O
    !        do i=1,channel_indices%num_channels
    !            LOG_WRITE channel_indices%channel(i)%hole_kappa, channel_indices%channel(i)%final_kappa, &
    !                    channel_indices%channel(i)%start_idx, channel_indices%channel(i)%end_idx
    !        end do

    ! #################################################################################################################
    ! #
    ! #################################################################################################################
    ! When we know how many channels we have we can allocate the list of perturbed wave functions,
    ! and initialise the arrays.
    number_of_channels = channel_indices%num_channels
    allocate(list_of_pwfs(number_of_channels))

    do i = 1,number_of_channels
        ! The number of bspline coefficients are hf%size_rhs so we pass this.
        call init_pwf(list_of_pwfs(i), channel_indices%channel(i), hf)
    end do


    ! We should decide on some list of XUV photon energies here.
    ! TODO(anton): Make the photon energies set by external parameters (ie config file).
    call get_xuv_photons_linspace(photons)

    ! For each photon we now have to form the perturbed wave function for each channel.
    ! From the perturbed wave function we compute necessary information, such as the
    ! ionisation rate.
    number_of_photons = photons%size

    LOG_WRITE "Computing for ",number_of_photons," XUV photons,"
    LOG_WRITE "starting from ",photons%start_eV," eV to ",photons%end_eV, " eV"
    do i =1, number_of_photons
        omega_XUV_au = photons%list_au(i)

        ! TODO(anton): Also we would need some way of choosing what resonances to exclude!
        ! Before implementing this we need to make sure that the we can get proper
        ! partial ionisation rates with the complete spectrum.

        ! Here we form the bspline coefficients for all channels, for a particular photon energy.
        call compute_bspline_coefficients(number_of_channels, diag, matrix_params, &
                channel_indices, hf, bsplines, omega_XUV_au, list_of_pwfs)

        ! TODO(anton): Implement calculation of ionisation rate. This will just be
        ! the same thing as in get_amp.f90 in relcode.

        LOG_WRITE "Done with photon ", i, " out of ", number_of_photons
    end do


end program main


! Test for in place reshaping:
!program main
!    use intrinsics_mod, only: dp
!    implicit none
!
!    complex(dp) :: zomplex
!    complex(dp), pointer :: pointer_3x3(:, :)
!
!    interface
!        subroutine test_in_place_reshape(ptr)
!            use intrinsics_mod, only: dp
!            complex(dp), pointer, intent(out) :: ptr(:,:)
!        end subroutine test_in_place_reshape
!    end interface
!
!    integer :: N, i, j
!
!    N = 3
!
!    call test_in_place_reshape(pointer_3x3)
!
!    LOG_WRITE "after in place reshape, supposedly destroyed"
!    if(associated(pointer_3x3)) then
!    do j = 1, N
!        do i = 1, N
!            LOG_WRITE i, j, pointer_3x3(i,j)
!        end do
!    end do
!    else
!        LOG_WRITE "pointer not associated"
!    end if
!
!
!end program main
!
!! -------------------------------------------------
!
!subroutine test_in_place_reshape(pointer_3x3)
!    ! This routine tests the use of C_F_POINTER to appropriately
!    ! reshape an array in place using a pointer.
!    use intrinsics_mod, only: dp
!    use, intrinsic :: iso_c_binding, only: C_LOC, C_F_POINTER
!    implicit none
!
!    complex(dp), pointer, intent(out) :: pointer_3x3(:, :)
!
!
!    complex(dp), allocatable, target :: raw_array(:)
!
!    integer :: N, i, j, idx
!
!    N = 3
!    allocate(raw_array(N*N))
!
!    do j = 1, N
!        do i = 1, N
!            idx = j+N*(i-1)
!            raw_array(idx) = cmplx(idx, j, kind=dp)
!        end do
!    end do
!
!    do i=1,size(raw_array)
!        LOG_WRITE i, raw_array(i)
!    end do
!
!    call C_F_POINTER(C_LOC(raw_array), pointer_3x3, (/3, 3/))
!
!    do j = 1, N
!        do i = 1, N
!            LOG_WRITE i, j, pointer_3x3(i,j)
!        end do
!    end do
!
!    deallocate(raw_array)
!    pointer_3x3 => null()
!
!
!end subroutine test_in_place_reshape