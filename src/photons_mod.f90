#include "logger.h"

module photons_mod
    use intrinsics_mod, only : dp
    use in_data_mod, only : Input_Parameters_t
    implicit none

    private

    public Photons_t, get_xuv_photons_linspace

!
!    integer, parameter :: g_num_photons = 1000
!    real(dp), parameter :: g_start_eV = 30.0_dp
!    real(dp), parameter :: g_end_eV = 60.0_dp
    real(dp), parameter :: g_eV_per_Hartree = 27.211396641308_dp

    ! -----

    type :: Photons_t
        integer :: size ! Number of photons
        real(dp) :: start_eV, end_eV, step_eV
        real(dp), allocatable, dimension(:) :: list_au
    end type Photons_t

    ! ----

contains

    subroutine get_xuv_photons_linspace(input_params, photons)
        type(Input_Parameters_t), intent(in) :: input_params
        type(Photons_t), intent(out) :: photons

        real(dp) :: delta_eV, energy_au, start_au, delta_au
        integer :: N, i

        N = input_params%number_of_photons

        if(allocated(photons%list_au)) then
            LOG_FATAL("Photons list already allocated!")
        else
            allocate(photons%list_au(N))
        end if


        photons%size = N
        photons%start_eV = input_params%start_eV
        photons%end_eV = input_params%end_eV
        delta_eV = (photons%end_eV-photons%start_eV)/real(N, kind=dp)
        photons%step_eV = delta_eV
        photons%list_au = 0.0_dp

        start_au = photons%start_eV/g_eV_per_Hartree
        delta_au = delta_eV/g_eV_per_Hartree

        LOG_WRITE "Creating photon linspace with delta ", delta_eV, " eV"

        do i=1,N
            energy_au = start_au + real(i-1, kind=dp)*delta_au
            photons%list_au(i) = energy_au
            ! Debug print
            !LOG_WRITE "photons%list_au(i) = ",photons%list_au(i)
        end do

    end subroutine get_xuv_photons_linspace

end module photons_mod