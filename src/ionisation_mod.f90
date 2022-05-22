#include "logger.h"
module ionisation_mod
    use intrinsics_mod, only : dp
    use utility_mod, only : l_from_kappa
    use types_mod, only : Channel_Index_t, HartreeFock_t, Matrix_Parameters_t, Diag_t, Bsplines_t, &
            Channel_Indices_t
    use perturbed_wavefunction_mod, only : Perturbed_Wavefunction_t
    use output_data_mod, only : Output_Data_t
    use in_data_mod, only : Input_Parameters_t
    implicit none

    private

    public calculate_ionisation_rate_and_amp_and_phase, test_total_absorption

    ! -----

    logical :: are_breakpoints_set ! so we only set breakpoints once.
    integer, parameter :: number_of_breakpoints = 5
    integer, parameter :: pick_breakpoint_index = 3

    integer, dimension(number_of_breakpoints) :: breakpoint_indices
    ! -----

contains

    ! ============================================================================

    ! This routine sets up indices used for breakpoints. Compare with get_amp.f90 in relcode.
    subroutine set_breakpoint_indices(bsplines)
        type(Bsplines_t), intent(in) :: bsplines

        real(dp) :: factor
        integer :: k, index_step, N, i

        index_step = 75

        are_breakpoints_set = .false.

        factor = 1.1_dp/3.0_dp

        N = bsplines%radial_pos_vector_size
        associate (pos_vec => bsplines%radial_pos_vector)
            k = 1
            do while(real(pos_vec(k)) < factor*real(pos_vec(N)))
                k = k + 1
                if (k > N) then
                    LOG_FATAL("while loop in set_breakpoint_indices() failed, k > N.")
                end if
            end do
        end associate

        if((k + (number_of_breakpoints-1)*index_step) > N) then
            LOG_FATAL("breakpoint indices computed to be larger than length of radial pos vec array!")
        end if

        do i = 1,number_of_breakpoints
            breakpoint_indices(i) = k + (i-1)*index_step
        end do

        are_breakpoints_set = .true.

    end subroutine set_breakpoint_indices


    ! ============================================================================

    ! The get_amp.f90 function in relcode also computes amplitude and phase (sans Coulomb part)
    ! of the perturbed wave function. Here we only compute ionisation rates per channel.
    subroutine calculate_ionisation_rate_and_amp_and_phase(Z, &
        bsplines, omega_index, omega_XUV_au, &
        num_channels, channel_indices, list_of_pwfs, output_data)
        use utility_mod, only : l_from_kappa
        
        real(dp), intent(in) :: Z
        type(Bsplines_t), intent(in) :: bsplines
        integer, intent(in) :: omega_index
        real(dp), intent(in) :: omega_XUV_au
        integer, intent(in) :: num_channels
        type(Channel_Indices_t), intent(in) :: channel_indices
        type(Perturbed_Wavefunction_t), dimension(num_channels), intent(in) :: list_of_pwfs
        type(Output_Data_t), intent(inout) :: output_data

        real(dp) :: radial_coordinate, amp, phaseF, phaseG
        complex(dp) :: i_imag, rhoF, rhoG
		! For Coulomb functions
		complex(dp) :: clmb1, clmb2, dclmb1, dclmb2, clmb_tot, sig, kinetic_energy
        integer :: i, out_channel_index, pos_index
        integer :: size1, size2, num_coeffs
        integer :: l_hole, hole_kappa, l_coulomb, l_final

        ! NOTE(anton):
        ! "size1" is kind of a misnomer here, it's just to conform to get_amp.f90 in relcode.
        ! I usually use size1 and size2 to refer to what dimension (what rank index) in an array.
        ! Here it means number of bsplines for small and large component respectively.
        ! (assuming small component has more elements than large, else this breaks!)
        size2 = bsplines%size2_small-2 ! In get_amp.f90 the bsplines have sizes size2+2 and size1+2 respectively.
        size1 = bsplines%size2_large-2

        i_imag = (0.0_dp, 1.0_dp)

        ! Functions and vars local to this module for figuring out breakpoint indices.
        if (.not.are_breakpoints_set) then
            call set_breakpoint_indices(bsplines)
        end if

        pos_index = breakpoint_indices(pick_breakpoint_index)

        output_data%pcur(omega_index, 1)   = omega_XUV_au
        output_data%amp(omega_index, 1)    = omega_XUV_au
        output_data%phaseF(omega_index, 1) = omega_XUV_au
        output_data%phaseG(omega_index, 1) = omega_XUV_au
        
        ! For each channel we compute the ioniosation rate at breakpoint with index pos_index.
        do i=1,num_channels
            num_coeffs = list_of_pwfs(i)%num_bspl_coeffs
            out_channel_index = i+1 ! For use in output data, first col is photon energies so we add 1.

            rhoF = sum(bsplines%B_large(pos_index, 2:size1+1)*&
                    list_of_pwfs(i)%bspline_coefficients(1:size1))

            rhoG = sum(bsplines%B_small(pos_index, 2:size2+1)*&
                    list_of_pwfs(i)%bspline_coefficients(size1+1:num_coeffs))

            output_data%pcur(omega_index, out_channel_index) = &
                    2.0_dp*i_imag*(conjg(rhoF)*rhoG - conjg(rhoG)*rhoF)
					
			! Compute coulomb function part of this perturbed wave function, at this coordinate
			kinetic_energy = omega_XUV_au + channel_indices%channel(i)%hole_binding_energy
            l_hole = channel_indices%channel(i)%hole_l
            hole_kappa = channel_indices%channel(i)%hole_kappa
            l_final = l_from_kappa(channel_indices%channel(i)%final_kappa)
            ! If we go down in l we compute a different coulomb function
            if(l_hole > l_final) then
                l_coulomb = l_hole-1
            else
                l_coulomb = l_hole+1
            end if
            
            radial_coordinate = bsplines%radial_pos_vector(pos_index)
            call coul_gen_cc_rel(Z, kinetic_energy, radial_coordinate, l_coulomb, .false., &
                clmb1, clmb2, dclmb1, dclmb2, sig)
            
            clmb_tot = (clmb2 + i_imag*clmb1)             

            amp    = (abs(rhoF) + abs(rhoG)) / abs(clmb_tot)
            phaseF = atan2(aimag(rhoF), real(rhoF))-atan2(aimag(clmb_tot),real(clmb_tot))
            phaseG = atan2(aimag(-i_imag*rhoG), real(-i_imag*rhoG)) - atan2(aimag(clmb_tot), real(clmb_tot))
            
            output_data%amp(omega_index, out_channel_index)    = amp
            output_data%phaseF(omega_index, out_channel_index) = phaseF
            output_data%phaseG(omega_index, out_channel_index) = phaseG
                
! Debug prints.
!            LOG_WRITE "channel i=",i,"bsplines%radial_pos_vector(pos_index)=",bsplines%radial_pos_vector(pos_index)
!            LOG_WRITE "sum(bsplines%B_large(pos_index, 2:size1+1))=", sum(bsplines%B_large(pos_index, 2:size1+1))
!            LOG_WRITE "sum(list_of_pwfs(i)%bspline_coefficients(1:size1))=",sum(list_of_pwfs(i)%bspline_coefficients(1:size1))
!            LOG_WRITE "sum(bsplines%B_small(pos_index, 2:size2+1))=",sum(bsplines%B_small(pos_index, 2:size2+1))
!            LOG_WRITE "sum(list_of_pwfs(i)%bspline_coefficients(size1+1:num_coeffs))=",sum(list_of_pwfs(i)%bspline_coefficients(size1+1:num_coeffs))
!            LOG_WRITE "rhoF=",rhoF,"rhoG=",rhoG,"pcur=",output_data%pcur(omega_index, channel_index)
!            LOG_WRITE " "
        end do


    end subroutine calculate_ionisation_rate_and_amp_and_phase

    ! ============================================================================

    ! This subroutine is just so we can confirm that the diagonalisation data we load in is correct.
    subroutine test_total_absorption(input_parameters, diag, omega, array_size, total_cs_file_id)
        type(Input_Parameters_t), intent(in) :: input_parameters
        type(Diag_t), intent(in) :: diag
        real(dp), intent(in) :: omega ! a.u.
        integer, intent(in) :: array_size, total_cs_file_id

        integer :: remove_index, i
        complex(dp), dimension(array_size) :: M
        complex(dp) :: imag_i
        real(dp) :: au_to_Mbarn, imag_term, c_au, C, out_Mbarn
        real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

        imag_i = (0.0_dp, 1.0_dp)
        ! Since the dipole transitions elements are defined with an overall extra i-factor in Fortran,
        ! we just multiply by -i here again to remove that.
        ! What is called "matrix element" M here is M_k = sum_i x_k^i * dipole_element_i,
        ! where x_k^i are the coefficients in the array representing eigenvector x_k after diagonalisation.
        M = -imag_i*diag%matrix_elements

        ! Set matrix element to zero at indices for those states that are to be removed.
        do i = 1, input_parameters%number_of_indices_to_remove
            remove_index = input_parameters%indices_to_remove(i)
            M(remove_index) = cmplx(0.0_dp, 0.0_dp)
        end do

        associate(E => diag%eigenvalues)
            ! NOTE(anton): Since the matrix is constructed with an energy diagonal of
            ! e_s - e_a (excited - ground), we need to subtract the photon energy here.
            ! Compare to ie Jimmy's lic with a denominator (e_a - e_s + Omega)
            imag_term = imag(sum(M * M / (E - omega)))
            c_au = 137.035999074_dp
            au_to_Mbarn = (0.529177210903_dp) * (0.529177210903_dp) * 100_dp

            ! NOTE(anton): There used to be an unexplained overall factor -1 here,
            ! but it comes from a definition of the dipole element with an extra i-factor in Fortran.
            ! This resulted in a -1 factor when taking the square of M above.
            ! Now the M are loaded with that i-factor removed, however!
            C = (4.0_dp*pi/3.0_dp)*(1.0_dp/c_au)*au_to_Mbarn

            out_Mbarn = C*imag_term*omega
        end associate

        write(total_cs_file_id, '(1e20.12)') out_Mbarn

    end subroutine test_total_absorption

    ! ============================================================================
end module ionisation_mod