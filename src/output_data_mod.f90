#include "logger.h"

#define MAX_CHAR_LENGTH 1024

module output_data_mod
    use intrinsics_mod, only : dp
    implicit none

    private

    public Output_Data_t

    type Output_Data_t
        ! Since Python scripts that handle ionisation rate data from relcode assumes a certain layout,
        ! we use the same layout here for writing out rates (ie as in pcur_all.dat).
        ! That is we have number_of_photons number of rows, and then the first column is photon energy.
        ! The following columns are for each ioniosation channel.
        ! Note that it's not exactly as in pcur_all, since we now have all channels in one file (ie not several pert-folders).
        complex(dp), allocatable, dimension(:,:) :: pcur
        complex(dp), allocatable, dimension(:,:) :: pwf_real, pwf_imag
        real(dp), allocatable, dimension(:, :) :: amp, phaseF, phaseG
    end type Output_Data_t

    contains


end module output_data_mod