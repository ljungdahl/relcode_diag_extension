#include "logger.h"
module perturbed_wavefunction_mod
    use intrinsics_mod, only : dp
    use types_mod, only : Channel_Index_t, HartreeFock_t, Matrix_Parameters_t, Diag_t, Bsplines_t, &
            Channel_Indices_t
    use in_data_mod, only : Input_Parameters_t
    implicit none

    private

    public Perturbed_Wavefunction_t, init_pwf, compute_bspline_coefficients

    ! -----

    type :: Perturbed_Wavefunction_t
        integer :: hole_kappa, final_kappa
        ! Note that the Bspline coefficients array contain coefficients both for large
        ! and for the small components. First we have bspl size1_large number of elements
        ! for the large component, then size_rhs-size1_large elements for the small component.
        ! (check get_amp.f90 in relcode program to see this).
        ! All in all this means that we have size_rhs coefficients (and we get size_rhs from the
        ! Hartree_Fock_t structure).
        complex(dp), allocatable, dimension(:) :: bspline_coefficients
        integer :: num_bspl_coeffs
        ! We also store the ionisation rate per photon energy omegaXUV here.
        complex(dp), allocatable, dimension(:) :: rate
    end type Perturbed_Wavefunction_t

    ! -----

contains

    ! ============================================================================

    subroutine init_pwf(pwf, channel, hf, num_photons)
        type(Perturbed_Wavefunction_t), intent(out) :: pwf
        type(Channel_Index_t), intent(in) :: channel
        type(HartreeFock_t), intent(in) :: hf
        integer, intent(in) :: num_photons

        integer :: num_bspl_coeffs

        num_bspl_coeffs = hf%size_rhs
        pwf%num_bspl_coeffs = num_bspl_coeffs

        pwf%hole_kappa = channel%hole_kappa
        pwf%final_kappa = channel%final_kappa

        if(allocated(pwf%bspline_coefficients)) then
            LOG_FATAL("pwf%bspline_coefficients already allocated!")
        else
            allocate(pwf%bspline_coefficients(num_bspl_coeffs))
        end if

        if(allocated(pwf%rate)) then
            LOG_FATAL("pwf%rate already allocated")
        else
            allocate(pwf%rate(num_photons))
        end if

    end subroutine init_pwf

    ! ============================================================================

    subroutine compute_bspline_coefficients(input_params, num_channels, diag, matrix_params, &
            chan_idxs, hf, bsplines, omega, list_of_pwfs)
        type(Input_Parameters_t), intent(in) :: input_params
        integer, intent(in) :: num_channels
        type(Diag_t), intent(in) :: diag
        type(Matrix_Parameters_t), intent(in) :: matrix_params
        type(Channel_Indices_t), intent(in) :: chan_idxs
        type(HartreeFock_t), intent(in) :: hf
        type(Bsplines_t), intent(in) :: bsplines
        real(dp), intent(in) :: omega ! XUV photon energy in atomic units
        type(Perturbed_Wavefunction_t), dimension(num_channels), intent(inout) :: list_of_pwfs

        !complex(dp), dimension(diag%system_size) :: pertmat
        complex(dp), dimension(size(diag%matrix_elements)) :: M
        complex(dp) :: pertmat

        integer :: i, j, is, ks, current_channel_index, remove_index
        ! This routine forms all the Bspline coefficients for the perturbed wave function, by
        ! multiplying and summing up all the "matrix element" coefficients from diagonalisation.
        ! The same procedure is done in relcode after solution of the linear system of equations.
        ! (see conform_solutions_to_pert_array.f90 in relcode repository)
        ! A big difference there is that the photon energy is "baked in" into the solution coefficients.
        ! Here we need to form the energy denominator explicity using diag eigenvalues and the provided omega.

        ! First we need to zero perturbed wave function bspline coefficients.
        do i = 1, num_channels
            list_of_pwfs(i)%bspline_coefficients = (0.0_dp, 0.0_dp)
        end do

        M = diag%matrix_elements
        do i = 1, input_params%number_of_indices_to_remove
            remove_index = input_params%indices_to_remove(i)
            M(remove_index) = cmplx(0.0_dp, 0.0_dp)
        end do
        ! Starting channel
        current_channel_index = 1

!        do i = 1, diag%system_size
!            ! Here we form the diag coefficents c_k with correct energy denominator.
!            ! pertmat = matrix_elements(i)/(eigenvalues(i) - omega), with
!            ! matrix_elements(i) = sum(psi_i*dipole_elements), where psi_i is eigenvector i
!            ! and eigenvalues(i) = (\epsilon_a-\epsilon_s)_i after diagonalisation.
!
!        end do

        do i = 1, diag%system_size
            ! These indices are used to get the proper hf coefficients.
            ! (Ie how we get coefficients for different quantum numbers in Jimmy's system)
            is = matrix_params%qn_indices(i, 1)
            ks = matrix_params%qn_indices(i, 2)

            ! TODO(anton): Understand this again... basis change?
            ! "pertmat" notation is from my old krylov version (perturbed_wf_construction.f90)
            pertmat = sum( diag%eigenvectors(i, :)*M / (diag%eigenvalues - dcmplx(omega, 0.0_dp)))

            ! Check if we have a new channel, if so increment channel index.
            if(i > chan_idxs%channel(current_channel_index)%end_idx) then
                current_channel_index = current_channel_index + 1
                if(current_channel_index > num_channels) then
                    LOG_FATAL("current_channel_index > num_channels")
                end if
            end if

            ! For each channel:
            ! Sum up the c_i multiplied on all the hf coefficients.
            ! Reminder: HF coefficients are the coefficients for expressing HF orbitals in
            ! Bspline basis.
            list_of_pwfs(current_channel_index)%bspline_coefficients(:) = &
                    list_of_pwfs(current_channel_index)%bspline_coefficients(:) + &
                            hf%eigenstates(:, is, ks) * pertmat

        end do

    end subroutine

    ! ============================================================================

end module perturbed_wavefunction_mod
