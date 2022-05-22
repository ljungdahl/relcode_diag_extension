module utility_mod
    use intrinsics_mod, only : dp
    implicit none

    private

    public l_from_kappa, j_from_kappa

contains

    ! ============================================================================

    integer function l_from_kappa(kappa)
        integer, intent(in) :: kappa
        !integer :: l_from_kappa

        if(kappa < 0) then
            l_from_kappa = (-kappa) - 1
        else
            l_from_kappa = kappa
        end if

    end function l_from_kappa

    ! ============================================================================

    real(dp) function j_from_kappa(kappa)
        integer, intent(in) :: kappa
        real(dp) :: kappa_dbl
        !real(dp) :: j_from_kappa

        kappa_dbl = dble(kappa)

        if(kappa < 0) then
            j_from_kappa = (-2.0_dp*kappa_dbl)-1
            else
            j_from_kappa = 2.0_dp*kappa_dbl-1
        end if

        j_from_kappa = j_from_kappa*0.5_dp

    end function j_from_kappa

    ! ============================================================================


end module utility_mod