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
            read_integer_array_binary_rank2, Input_Parameters_t, read_input_parameters

    use output_data_mod, only : Output_Data_t
    use types_mod, only : HartreeFock_t, Diag_t, Matrix_Parameters_t, Channel_Indices_t, Bsplines_t

    use perturbed_wavefunction_mod, only : Perturbed_Wavefunction_t, init_pwf, &
            compute_bspline_coefficients

    use photons_mod, only : Photons_t, get_xuv_photons_linspace

    use ionisation_mod, only : calculate_ionisation_rate_and_amp_and_phase, test_total_absorption

    implicit none

    type(Input_Parameters_t) :: input_parameters
    type(Output_Data_t) :: output_data

    type(In_Metadata_t) :: diag_eigenvalues_metadata, diag_matrix_elements_metadata, diag_eigenvectors_metadata
    type(Diag_t) :: diag
    type(In_Metadata_t) :: hf_eigenstates_metadata, hf_eigenvalues_metadata ! NOTE(anton): Not sure if eigenvalues are used.
    type(HartreeFock_t) :: hf
    type(In_Metadata_t) :: matrix_params_qn_indices_metadata
    type(Matrix_Parameters_t) :: matrix_params
    type(In_Metadata_t) :: bspl_large_metadata ! Large component
    type(In_Metadata_t) :: bspl_small_metadata ! Small component
    type(In_Metadata_t) :: bspl_radial_pos_vec_metadata ! Vector of radial coords in which the bsplines are defined.
    type(Bsplines_t) :: bsplines
    type(Channel_Indices_t) :: channel_indices ! This is read using type derived I/O, ie overload of the read/write operators!

    ! We have a list of perturbed wavefunctions, one for each channel, so that w ecan compute
    ! partial probability currents (rates) etc.
    type(Perturbed_Wavefunction_t), allocatable, dimension(:) :: list_of_pwfs

    type(Photons_t) :: photons
    real(dp) :: omega_XUV_au
    complex(dp), dimension(:), allocatable :: rhoF_tmp
    character(len = *), parameter :: filename_channel_indices = "./data_from_relcode/channel_indices.bin"
    character(len = *), parameter :: output_pcur_filename = "omega_and_all_channels_pcur.dat"
    character(len = *), parameter :: output_amp_filename = "omega_and_all_channels_amp.dat"
    character(len = *), parameter :: output_phaseF_filename = "omega_and_all_channels_phaseF.dat"
    character(len = *), parameter :: output_phaseG_filename = "omega_and_all_channels_phaseG.dat"
    character(len = *), parameter :: pert_wave_re_filename = "./pert_waves/pwf_real"
    character(len = *), parameter :: pert_wave_im_filename = "./pert_waves/pwf_imag"
    character(len = *), parameter :: total_cs_filename = "total_absorption_cs.dat"
    character(len = 1024) :: pwf_re_final_fname, pwf_im_final_fname
    character(len = 24) :: pcur_and_omega_fmt_str

    integer :: i, j, number_of_channels, number_of_photons, m, n
    integer :: omega_index, pcur_out_file_id, pwf_re_file_id, pwf_im_file_id, amp_out_file_id, phaseF_out_file_id
    integer :: total_cs_file_id, index_to_remove, phaseG_out_file_id

    ! #################################################################################################################
    ! # Read input parameters
    ! #################################################################################################################
    call read_input_parameters(input_parameters)

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
    call read_metadata("diag_eigenvectors", diag_eigenvectors_metadata)

    if(diag_eigenvalues_metadata%sizes(1) /= diag_matrix_elements_metadata%sizes(1) .or. &
        diag_eigenvectors_metadata%sizes(1) /= diag_eigenvalues_metadata%sizes(1) .or. &
        diag_eigenvectors_metadata%sizes(2) /= diag_eigenvalues_metadata%sizes(1)) then
        LOG_WRITE "diag_eigenvectors_size1 = ", diag_eigenvectors_metadata%sizes(1), &
                "diag_eigenvectors_size2 = ", diag_eigenvectors_metadata%sizes(2)
        LOG_WRITE "diag_eigenvalues size = ", diag_eigenvalues_metadata%sizes(1)
        LOG_WRITE "diag_matrix_elements size = ", diag_matrix_elements_metadata%sizes(1)
        LOG_FATAL("Eigenvectors/values and matrix elements read from relcode should have the same dimensions!")
    end if
    diag%system_size = diag_matrix_elements_metadata%sizes(1)

    allocate(diag%eigenvectors(diag%system_size, diag%system_size))
    call read_complex_array_binary_rank2("diag_eigenvectors", diag_eigenvectors_metadata, diag%eigenvectors)

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
    hf%size_rhs = hf_eigenstates_metadata%sizes(1)
    hf%size_eig = hf_eigenstates_metadata%sizes(2)
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

!    LOG_WRITE "bsplines%size1_large =", bsplines%size1_large
!    LOG_WRITE "bsplines%size2_large =", bsplines%size2_large
!    LOG_WRITE "bsplines%size1_small =", bsplines%size1_small
!    LOG_WRITE "bsplines%size2_small =", bsplines%size2_small

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
    ! Check that we're not trying to remove an index that is out of bounds of the diag data.
    do i = 1, input_parameters%number_of_indices_to_remove
        index_to_remove = input_parameters%indices_to_remove(i)
        if(index_to_remove > diag%system_size) then
            LOG_WRITE "diag%system_size = ", diag%system_size
            LOG_WRITE "index_to_remove = ", index_to_remove
            LOG_FATAL("We can't remove an index larger than size of diag data!")
        end if
    end do

    ! We should decide on some list of XUV photon energies here.
    call get_xuv_photons_linspace(input_parameters, photons)

    ! For each photon we now have to form the perturbed wave function for each channel.
    ! From the perturbed wave function we compute necessary information, such as the
    ! ionisation rate.
    number_of_photons = photons%size

    ! When we know how many channels we have we can allocate the list of perturbed wave functions,
    ! and initialise the arrays.
    number_of_channels = channel_indices%num_channels
    write(*,*) "Channels read from relcode data:"
	do i = 1, number_of_channels
        write(*,*) trim(channel_indices%channel(i)%print_string), ", with hole binding energy: ", channel_indices%channel(i)%hole_binding_energy, " au"
	end do
    allocate(list_of_pwfs(number_of_channels))

    ! Then we also can allocate the output pcur array.
    ! First column is photon energies, then we have number_of_channels columns.
    allocate(output_data%pcur(number_of_photons, number_of_channels + 1))
    allocate(output_data%amp(number_of_photons, number_of_channels + 1))
    allocate(output_data%phaseF(number_of_photons, number_of_channels + 1))
    allocate(output_data%phaseG(number_of_photons, number_of_channels + 1))
    output_data%pcur = (0.0_dp, 0.0_dp)
    output_data%amp = (0.0_dp, 0.0_dp)
    output_data%phaseF = (0.0_dp, 0.0_dp)
    output_data%phaseG = (0.0_dp, 0.0_dp)

    ! Create the format string for writing output. 
    write(pcur_and_omega_fmt_str, '(a,i0,a7)')"(", size(output_data%pcur, 2), "e20.12)"
    LOG_WRITE "trim(pcur_and_omega_fmt_str) = ", trim(pcur_and_omega_fmt_str)
    
    do i = 1, number_of_channels
        ! The number of bspline coefficients are hf%size_rhs so we pass this.
        ! We also allocate the ionisation rate arrays (member var of perturbed wave function type)
        ! to have number_of_photons elements.
        call init_pwf(list_of_pwfs(i), channel_indices%channel(i), hf, number_of_photons)
    end do

    ! allocate output data for perturbed wave functions, real imag part.
    !allocate(output_data%pwf_real(bsplines%radial_pos_vector_size, number_of_channels))
    !allocate(output_data%pwf_imag(bsplines%radial_pos_vector_size, number_of_channels))
    !allocate(rhoF_tmp(number_of_channels))
    ! File is open for writing here. Will be written to in every energy step.
    pcur_out_file_id = 35
    open(pcur_out_file_id, file = trim(output_pcur_filename), action = "write")
    LOG_WRITE "Opened file for writing ", trim(output_pcur_filename)

    amp_out_file_id = pcur_out_file_id + 1
    open(amp_out_file_id, file = trim(output_amp_filename), action = "write")
    LOG_WRITE "Opened file for writing ", trim(output_amp_filename)
    
    phaseF_out_file_id = amp_out_file_id + 1
    open(phaseF_out_file_id, file = trim(output_phaseF_filename), action = "write")
    LOG_WRITE "Opened file for writing ", trim(output_phaseF_filename)
    
    phaseG_out_file_id = phaseF_out_file_id + 1
    open(phaseG_out_file_id, file = trim(output_phaseG_filename), action = "write")
    LOG_WRITE "Opened file for writing ", trim(output_phaseG_filename)
    
    total_cs_file_id = 77
    open(total_cs_file_id, file=trim(total_cs_filename), action="write")
    LOG_WRITE "Opened file for writing ", trim(total_cs_filename)

    LOG_WRITE "Computing for ", number_of_photons, " XUV photons,"
    LOG_WRITE "starting from ", photons%start_eV, " eV to ", photons%end_eV, " eV"
    LOG_WRITE "with step size ", photons%step_eV, " eV."
    do i = 1, number_of_photons
        omega_XUV_au = photons%list_au(i)

        ! NB:
        ! The "input_parameters" structure contains the indices used to remove states from the spectrum.

        ! To verify we got diagonalisation data correctly we compute total absorption cross section here also.
        ! This subroutine writes the total cross section to file for each omegaXUV.
        call test_total_absorption(input_parameters, diag, omega_XUV_au, size(diag%matrix_elements), total_cs_file_id)

        ! Here we form the bspline coefficients for all channels, for a particular photon energy.

        call compute_bspline_coefficients(input_parameters, number_of_channels, diag, matrix_params, &
                channel_indices, hf, bsplines, omega_XUV_au, list_of_pwfs)

        omega_index = i
        call calculate_ionisation_rate_and_amp_and_phase(input_parameters%Z, &
                bsplines, omega_index, omega_XUV_au, &
                number_of_channels, channel_indices, list_of_pwfs, output_data)

        ! Note that we can use the same output format string for pcur, amp, phaseF/G
        write(pcur_out_file_id,   trim(pcur_and_omega_fmt_str)) real(output_data%pcur(i, :))
        write(amp_out_file_id,    trim(pcur_and_omega_fmt_str)) output_data%amp(i, :)
        write(phaseF_out_file_id, trim(pcur_and_omega_fmt_str)) output_data%phaseF(i, :)
        write(phaseG_out_file_id, trim(pcur_and_omega_fmt_str)) output_data%phaseG(i, :)
        ! Temp writing out perturbed wave function for debugging
!        if(i == 500) then
!            pwf_re_file_id = 44
!            pwf_im_file_id = pwf_re_file_id + 1
!            write(pwf_re_final_fname, '(a,a,i0,a)') pert_wave_re_filename, "_", i, ".dat"
!            write(pwf_im_final_fname, '(a,a,i0,a)') pert_wave_im_filename, "_", i, ".dat"
!            open(pwf_re_file_id, file = trim(pwf_re_final_fname), action = "write")
!            LOG_WRITE "Opened file for writing ", trim(pwf_re_final_fname)
!            open(pwf_im_file_id, file = trim(pwf_im_final_fname), action = "write")
!            LOG_WRITE "Opened file for writing ", trim(pwf_im_final_fname)
!
!            do m = 1, bsplines%radial_pos_vector_size
!                do j = 1, number_of_channels
!                    rhoF_tmp(j) = sum(bsplines%B_large(m, 2:(bsplines%size2_large - 2) + 1) * &
!                            list_of_pwfs(j)%bspline_coefficients(1:bsplines%size2_large - 2))
!                end do
!                ! The same fmt string works since we have radial coord instead of omega now.
!                write(pwf_re_file_id, pcur_and_omega_fmt_str) real(bsplines%radial_pos_vector(m)), real(rhoF_tmp)
!                write(pwf_im_file_id, pcur_and_omega_fmt_str) real(bsplines%radial_pos_vector(m)), imag(rhoF_tmp)
!            end do
!
!            close(pwf_re_file_id)
!            close(pwf_im_file_id)
!            LOG_WRITE "Wrote to file ", trim(pwf_re_final_fname)
!            LOG_WRITE "Wrote to file ", trim(pwf_im_final_fname)
!        end if

        LOG_WRITE "Done with photon ", i, " out of ", number_of_photons
    end do

    LOG_WRITE "Done with all photons."

    close(pcur_out_file_id)
    LOG_WRITE "Wrote to file ", trim(output_pcur_filename)
    close(amp_out_file_id)
    LOG_WRITE "Wrote to file ", trim(output_amp_filename)
    close(phaseF_out_file_id)
    LOG_WRITE "Wrote to file ", trim(output_phaseF_filename)
    close(phaseG_out_file_id)
    LOG_WRITE "Wrote to file ", trim(output_phaseG_filename)

    close(total_cs_file_id)
    LOG_WRITE "Wrote to file ",  trim(total_cs_filename)

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