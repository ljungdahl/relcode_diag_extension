module types_mod
    use intrinsics_mod, only : dp, sp
    implicit none

    private

    public Perturbed_Wavefunction_t, In_Metadata_t


    ! -----

    type :: Perturbed_Wavefunction_t

    end type Perturbed_Wavefunction_t

    ! -----

    type :: In_Metadata_t

    end type In_Metadata_t

    ! -----

    type :: In_Data_t
        type(In_Metadata_t) :: meta

        ! We load the raw binary data into a one-dimensional array.
        ! Then we can use the connected metadata to
        complex(dp), allocatable, dimension(:) :: raw
    end type In_Data_t

contains


end module types_mod