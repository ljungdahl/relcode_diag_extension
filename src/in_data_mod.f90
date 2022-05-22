#include "logger.h"

#define MAX_CHAR_LENGTH 1024
#define MAX_RANK 7
#define MAX_NUM_PARAMETERS 10

module in_data_mod
    use intrinsics_mod, only : dp, sp
    use types_mod, only : data_from_relcode_dir
    implicit none

    private

    public In_Metadata_t, read_metadata, check_if_subdir_exists, &
            read_complex_array_binary_rank1, read_complex_array_binary_rank2, read_complex_array_binary_rank3, &
            read_integer_array_binary_rank2, read_input_parameters, Input_Parameters_t

    character(len = *), parameter :: input_params_filename = "diag_ext.input"
    character(len = *), parameter :: bin_metadata_arrsize_format_str = &
            "(i0, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0)"

    ! -----

    type :: In_Metadata_t
        integer :: number_of_non_comment_lines
        character(len = MAX_CHAR_LENGTH), allocatable, dimension(:) :: non_comment_lines
        character(len = MAX_CHAR_LENGTH) :: filename
        integer, allocatable, dimension(:) :: sizes
        integer :: rank
        logical :: is_read
    end type In_Metadata_t

    ! -----

    type :: Input_Parameters_t
        integer :: number_of_photons, number_of_indices_to_remove
        real(dp) :: start_eV, end_eV
        real(dp) :: Z
        integer, allocatable, dimension(:) :: indices_to_remove
    end type Input_Parameters_t

contains
    ! ================================================================================================

    subroutine read_input_parameters(in_params)
        type(Input_Parameters_t), intent(out) :: in_params

        integer :: file_id, stat, i
        character(len = MAX_CHAR_LENGTH) :: comment_line

        file_id = 102
        open(file_id, file = trim(input_params_filename), status = 'old', action = 'read', iostat = stat)

        if(stat /= 0) then
            write(*, *) "FATAL ERROR! Cannot open file ", trim(input_params_filename)
            call exit(1)
        end if

        read(file_id, '(a)') comment_line
        read(file_id, *) in_params%Z
        read(file_id, '(a)') comment_line
        read(file_id, *) in_params%number_of_photons
        read(file_id, '(a)') comment_line
        read(file_id, *) in_params%start_eV
        read(file_id, '(a)') comment_line
        read(file_id, *) in_params%end_eV
        read(file_id, '(a)') comment_line
        read(file_id, *) in_params%number_of_indices_to_remove

        if(in_params%number_of_indices_to_remove > 0) then
            read(file_id, '(a)') comment_line
            allocate(in_params%indices_to_remove(in_params%number_of_indices_to_remove))
            do i = 1, in_params%number_of_indices_to_remove
                read(file_id, *) in_params%indices_to_remove(i)
            end do
        end if

        close(file_id)

        LOG_WRITE "Read input parameters from file ", trim(input_params_filename)
        LOG_WRITE "Atom Z = ", in_params%Z
        LOG_WRITE "number_of_photons = ", in_params%number_of_photons
        LOG_WRITE "start_eV = ", in_params%start_eV
        LOG_WRITE "end_eV = ", in_params%end_eV
        LOG_WRITE "number of indices to remove =", in_params%number_of_indices_to_remove
        if(in_params%number_of_indices_to_remove > 0) then
            LOG_WRITE "indices to remove:"
            do i = 1, in_params%number_of_indices_to_remove
                LOG_WRITE in_params%indices_to_remove(i)
            end do
        end if

        LOG_WRITE " "

    end subroutine read_input_parameters

    ! ================================================================================================

    subroutine read_metadata(filename_prefix, meta_data)
        character(len = *), intent(in) :: filename_prefix
        type(In_Metadata_t), intent(out) :: meta_data

        character(len = *), parameter :: suffix = ".metadata"
        character(len = MAX_CHAR_LENGTH) :: filename, line
        character(len = MAX_CHAR_LENGTH), dimension(MAX_NUM_PARAMETERS) :: param_string_list
        integer :: file_id, line_counter, param_counter, stat, num_params, i
        integer :: in_data_rank
        integer, dimension(MAX_RANK) :: sizes

        meta_data%is_read = .false.
        sizes = 0

        filename = trim(data_from_relcode_dir) // trim(filename_prefix) // trim(suffix)

        meta_data%filename = trim(filename)

        ! Parse the metadata file
        file_id = 102
        open(file_id, file = trim(filename), status = 'old', action = 'read', iostat = stat)

        if(stat /= 0) then
            write(*, *) "FATAL ERROR! Cannot open file ", trim(filename)
            LOG_FATAL("")
        end if

        line_counter = 1
        param_counter = 0

        stat = 0
        do while(stat == 0)

            read(file_id, '(a)', iostat = stat) line

            if (stat /= 0) cycle

            if (line(1:1) == "#" .or. line(1:2) == " #") then
                !LOG_WRITE "Line ", line_counter, " is a comment"
            else

                param_counter = param_counter + 1
                param_string_list(param_counter) = trim(line)

            end if
            line_counter = line_counter + 1
        end do

        close(file_id)

        num_params = param_counter

        if(num_params .lt. 1) then
            LOG_FATAL("FATAL ERROR: Didn't ready any parameters!")
        end if

        ! Translate the parsed information to the metadata structure
        meta_data%number_of_non_comment_lines = num_params
        allocate(meta_data%non_comment_lines(num_params))

        do i = 1, num_params
            meta_data%non_comment_lines(i) = trim(param_string_list(i))
            !LOG_WRITE trim(meta_data%non_comment_lines(i))
        end do

        ! The first non_comment line should always be the line with sizes for all dimensions
        ! Read the size of each dimension and save here.
        read(meta_data%non_comment_lines(1), *) &
                sizes(1), sizes(2), sizes(3), sizes(4), sizes(5), sizes(6), sizes(7)

        !        do i = 1, size(sizes)
        !            LOG_WRITE i, sizes(i)
        !        end do

        ! Whenever we encounter a zero in the sizes we know we've got the max rank.
        in_data_rank = 0
        do i = 1, size(sizes)
            if(sizes(i) == 0) exit
            in_data_rank = i
        end do

        if(in_data_rank == 0) then
            LOG_FATAL("FATAL ERROR! Data to be read must have at least rank 1, got zero!")
        end if
        meta_data%rank = in_data_rank

        ! Store the dimension sizes in the meta_data structure.
        allocate(meta_data%sizes(in_data_rank))
        do i = 1, in_data_rank
            meta_data%sizes(i) = sizes(i)
        end do

        ! Now we're done with reading metadata. We can add more things above here.
        meta_data%is_read = .true.

    end subroutine read_metadata

    ! ================================================================================================

    subroutine read_complex_array_binary_rank1(filename_prefix, meta_data, array_to_be_filled)
        use intrinsics_mod, only : dp
        use types_mod, only : data_from_relcode_dir

        character(len = *), intent(in) :: filename_prefix
        type(In_Metadata_t), intent(in) :: meta_data
        complex(dp), dimension(meta_data%sizes(1)), intent(out) :: array_to_be_filled

        character(len = MAX_CHAR_LENGTH) :: filename
        integer :: file_id

        file_id = 45

        filename = trim(data_from_relcode_dir) // trim(filename_prefix) // ".bin"

        open(file_id, file = trim(filename), form = "unformatted")
        read(file_id) array_to_be_filled
        close(file_id)

    end subroutine read_complex_array_binary_rank1

    ! ================================================================================================

    subroutine read_complex_array_binary_rank2(filename_prefix, meta_data, array_to_be_filled)
        use intrinsics_mod, only : dp
        use types_mod, only : data_from_relcode_dir

        character(len = *), intent(in) :: filename_prefix
        type(In_Metadata_t), intent(in) :: meta_data
        complex(dp), dimension(meta_data%sizes(1), meta_data%sizes(2)), intent(out) :: array_to_be_filled

        character(len = MAX_CHAR_LENGTH) :: filename
        integer :: file_id

        file_id = 45

        filename = trim(data_from_relcode_dir) // trim(filename_prefix) // ".bin"

        open(file_id, file = trim(filename), form = "unformatted")
        read(file_id) array_to_be_filled
        close(file_id)

    end subroutine read_complex_array_binary_rank2

    ! ================================================================================================

    subroutine read_complex_array_binary_rank3(filename_prefix, meta_data, array_to_be_filled)
        use intrinsics_mod, only : dp
        use types_mod, only : data_from_relcode_dir

        character(len = *) :: filename_prefix
        type(In_Metadata_t), intent(in) :: meta_data
        complex(dp), dimension(meta_data%sizes(1), meta_data%sizes(2), meta_data%sizes(3)), intent(out) :: array_to_be_filled

        character(len = MAX_CHAR_LENGTH) :: filename
        integer :: file_id

        file_id = 45

        filename = trim(data_from_relcode_dir) // trim(filename_prefix) // ".bin"

        open(file_id, file = trim(filename), form = "unformatted")
        read(file_id) array_to_be_filled
        close(file_id)

    end subroutine read_complex_array_binary_rank3

    ! ================================================================================================

    subroutine read_integer_array_binary_rank2(filename_prefix, meta_data, array_to_be_filled)
        use intrinsics_mod, only : dp
        use types_mod, only : data_from_relcode_dir

        character(len = *), intent(in) :: filename_prefix
        type(In_Metadata_t), intent(in) :: meta_data
        integer, dimension(meta_data%sizes(1), meta_data%sizes(2)), intent(out) :: array_to_be_filled

        character(len = MAX_CHAR_LENGTH) :: filename
        integer :: file_id

        file_id = 45

        filename = trim(data_from_relcode_dir) // trim(filename_prefix) // ".bin"

        open(file_id, file = trim(filename), form = "unformatted")
        read(file_id) array_to_be_filled
        close(file_id)

    end subroutine read_integer_array_binary_rank2

    ! ================================================================================================
    subroutine check_if_subdir_exists()
        ! If data_from_relcode_dir does not exist we throw an error, since there is no data to be read.
        use types_mod, only : data_from_relcode_dir
        implicit none

        character(len = MAX_CHAR_LENGTH) :: dir_path_no_ending_slash
        integer :: dir_str_len
        logical :: dir_exists

        dir_str_len = len(trim(data_from_relcode_dir))
        dir_path_no_ending_slash = trim(data_from_relcode_dir)
        dir_path_no_ending_slash = dir_path_no_ending_slash(1:dir_str_len - 1)
        !LOG_WRITE trim(dir_path_no_ending_slash)
        inquire(file = trim(dir_path_no_ending_slash), exist = dir_exists)

        if(.not. dir_exists) then
            LOG_WRITE "FATAL ERROR: Could not find subdirector with data:", trim(dir_path_no_ending_slash)
            LOG_WRITE "Make sure it exists and the data from relcode is present!"
            LOG_FATAL("")
        end if

    end subroutine check_if_subdir_exists

    ! ================================================================================================

end module in_data_mod


! NOTE(anton): The following was an attempt to have some general system of loading data from file and then
! reshaping a rank 1 array into whatever rank it should be. For now this is on hold and we just write the
! code needed to load in the data straight up. We can revisit this if the number of structures that has to be
! transfered between the programs are too many to handle explicitly.
!module in_data_mod
!    use intrinsics_mod, only : dp, sp
!    implicit none
!
!    private
!
!    public In_Metadata_t, In_Data_t
!
!
!    ! -----
!
!    type :: In_Metadata_t
!        integer :: number_of_non_comment_lines
!        character(len = MAX_CHAR_LENGTH), allocatable, dimension(:) :: non_comment_lines
!        logical :: is_read
!    end type In_Metadata_t
!
!    ! -----
!
!    type :: In_Data_t
!        type(In_Metadata_t) :: meta
!        ! We load the raw binary data into a one-dimensional array.
!        ! Then we can use the connected metadata to
!        complex(dp), allocatable, dimension(:) :: raw
!        integer :: rank ! The actual rank of the data, ie this should be the same as the
!        ! rank of the pointer which will be in-place reshaped pointing to raw.
!    contains ! Type methods
!        procedure, public :: deallocate => in_data_deallocate
!    end type In_Data_t
!
!    ! -----
!
!
!contains
!
!    subroutine in_data_deallocate(this)
!        class(In_Data_t), intent(out) :: this
!
!        if(allocated(this%raw)) then
!            deallocate(this%raw)
!        end if
!
!        if(allocated(this%meta%non_comment_lines)) then
!            deallocate(this%meta%non_comment_lines)
!        end if
!
!    end subroutine in_data_deallocate
!
!
!end module in_data_mod