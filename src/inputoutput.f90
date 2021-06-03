module inputoutput
  use parallel
  use communication
  implicit none

  private
  integer,parameter :: len_max = 512

  integer :: id_inputfile = 16
  character(64) :: name_inputfile = './inp'
  integer :: id_input_log = 17
  character(64) :: name_input_log   = './inp_log.out'

  integer :: id_outputfile_sequential


  public :: init_inputoutput, &
            read_basic_input, &
            read_vector_input, &
            read_matrix_input, &
            fin_inputoutput


  interface read_basic_input
     module procedure read_basic_input_character
     module procedure read_basic_input_integer
     module procedure read_basic_input_real8
     module procedure read_basic_input_logical
  end interface read_basic_input

  interface read_vector_input
     module procedure read_vector_input_integer
     module procedure read_vector_input_real8
  end interface read_vector_input

  interface read_matrix_input
     module procedure read_matrix_input_integer
     module procedure read_matrix_input_real8
  end interface

contains
!-------------------------------------------------------------------------------
  subroutine init_inputoutput
    implicit none

    if(if_root_global)then
      open(id_inputfile,file=name_inputfile)
      open(id_input_log,file=name_input_log)
    end if

    sequential = 20

  end subroutine init_inputoutput
!-------------------------------------------------------------------------------
  subroutine fin_inputoutput
    implicit none

    if(if_root_global)then
      close(id_inputfile)
      close(id_input_log)
    end if

  end subroutine fin_inputoutput
!-------------------------------------------------------------------------------
  subroutine read_basic_input_character(name, val, val_default, if_default)
    implicit none
    character(*),intent(in) :: name
    character(*),intent(out) :: val
    character(*),intent(in),optional :: val_default
    logical,intent(out),optional :: if_default
    character(len_max) :: val_t
    logical :: if_found
    integer :: index_equal, length_trimed
    
    if(if_root_global)then
      if(present(val_default))val = val_default
      call lookup_input(name, if_found, val_t)
      if(if_found)then
        index_equal= index(val_t,'=')
        length_trimed = len_trim(val_t)
        val = trim(adjustl(val_t(index_equal+1:length_trimed)))
      end if
    end if
    call comm_bcast(val)

    if(present(if_default))then
      if_default = .not.if_found
      call comm_bcast(if_default)
    end if

  end subroutine read_basic_input_character
!-------------------------------------------------------------------------------
  subroutine read_basic_input_integer(name, val, val_default, if_default)
    implicit none
    character(*),intent(in) :: name
    integer,intent(out) :: val
    integer,intent(in),optional :: val_default
    logical,intent(out),optional :: if_default
    character(len_max) :: val_t
    logical :: if_found
    integer :: index_equal, length_trimed
    
    if(if_root_global)then
      if(present(val_default))val = val_default
      call lookup_input(name, if_found, val_t)
      if(if_found)then
        index_equal= index(val_t,'=')
        length_trimed = len_trim(val_t)
        val_t = trim(val_t(index_equal+1:length_trimed))
        read(val_t,*) val
      end if
    end if
    call comm_bcast(val)

    if(present(if_default))then
      if_default = .not.if_found
      call comm_bcast(if_default)
    end if

  end subroutine read_basic_input_integer
!-------------------------------------------------------------------------------
  subroutine read_basic_input_real8(name, val, val_default, if_default)
    implicit none
    character(*),intent(in) :: name
    real(8),intent(out) :: val
    real(8),intent(in),optional :: val_default
    logical,intent(out),optional :: if_default
    character(len_max) :: val_t
    logical :: if_found
    integer :: index_equal, length_trimed
    
    if(if_root_global)then
      if(present(val_default))val = val_default
      call lookup_input(name, if_found, val_t)
      if(if_found)then
        index_equal= index(val_t,'=')
        length_trimed = len_trim(val_t)
        val_t = trim(val_t(index_equal+1:length_trimed))
        read(val_t,*) val
      end if
    end if
    call comm_bcast(val)

    if(present(if_default))then
      if_default = .not.if_found
      call comm_bcast(if_default)
    end if

  end subroutine read_basic_input_real8
!-------------------------------------------------------------------------------
  subroutine read_basic_input_logical(name, val, val_default, if_default)
    implicit none
    character(*),intent(in) :: name
    logical,intent(out) :: val
    logical,intent(in),optional :: val_default
    logical,intent(out),optional :: if_default
    character(len_max) :: val_t
    character(len=3) :: val3
    character(len=3),parameter :: cyes = 'yes', cy = 'y', cno = 'no', cn = 'n'
    logical :: if_found, if_error
    integer :: index_equal, length_trimed
    
    if_error =.false.

    if(if_root_global)then
      if(present(val_default))val = val_default
      call lookup_input(name, if_found, val_t)
      if(if_found)then
        index_equal= index(val_t,'=')
        length_trimed = len_trim(val_t)
        val_t = trim(adjustl(val_t(index_equal+1:length_trimed)))
        if(trim(val_t) == 'y' .or. trim(val_t) == 'yes')then
          val = .true.
        else if(trim(val_t) == 'n' .or. trim(val_t) == 'no')then
          val = .false.
        else
          if_error =.true.
        end if
      end if
    end if
    call comm_bcast(if_error)

!Error check
    if(if_error)then
      message(1) = 'Error: Invalid input for '//trim(name)//'.'
      call error_finalize(message(1))
    end if


    call comm_bcast(val)

    if(present(if_default))then
      if_default = .not.if_found
      call comm_bcast(if_default)
    end if

  end subroutine read_basic_input_logical
!-------------------------------------------------------------------------------
  subroutine read_vector_input_integer(name, val, val_default, if_default)
    implicit none
    character(*),intent(in) :: name
    integer,intent(out) :: val(:)
    integer,intent(in),optional :: val_default(:)
    logical,intent(out),optional :: if_default
    character(len_max) :: val_t
    logical :: if_found
    integer :: index_equal, length_trimed
    
    if(if_root_global)then
      if(present(val_default))val = val_default
      call lookup_input(name, if_found)
      if(if_found)then
        read(id_inputfile, *) val(:)
      end if
    end if
    call comm_bcast(val)

    if(present(if_default))then
      if_default = .not.if_found
      call comm_bcast(if_default)
    end if

  end subroutine read_vector_input_integer
!-------------------------------------------------------------------------------
  subroutine read_vector_input_real8(name, val, val_default, if_default)
    implicit none
    character(*),intent(in) :: name
    real(8),intent(out) :: val(:)
    real(8),intent(in),optional :: val_default(:)
    logical,intent(out),optional :: if_default
    character(len_max) :: val_t
    logical :: if_found
    integer :: index_equal, length_trimed
    
    if(if_root_global)then
      if(present(val_default))val = val_default
      call lookup_input(name, if_found)
      if(if_found)then
        read(id_inputfile, *) val(:)
      end if
    end if
    call comm_bcast(val)

    if(present(if_default))then
      if_default = .not.if_found
      call comm_bcast(if_default)
    end if

  end subroutine read_vector_input_real8
!-------------------------------------------------------------------------------
  subroutine read_matrix_input_integer(name, val, val_default, if_default)
    implicit none
    character(*),intent(in) :: name
    integer,intent(out) :: val(:,:)
    integer,intent(in),optional :: val_default(:,:)
    logical,intent(out),optional :: if_default
    character(len_max) :: val_t
    logical :: if_found
    integer :: index_equal, length_trimed
    integer :: i_lbound, i_ubound, i
    
    if(if_root_global)then
      i_lbound = lbound(val,1); i_ubound = ubound(val,1)

      if(present(val_default))val = val_default
      call lookup_input(name, if_found)
      if(if_found)then
        do i = i_lbound, i_ubound
          read(id_inputfile, *) val(i,:)
        end do
      end if
    end if
    call comm_bcast(val)

    if(present(if_default))then
      if_default = .not.if_found
      call comm_bcast(if_default)
    end if

  end subroutine read_matrix_input_integer
!-------------------------------------------------------------------------------
  subroutine read_matrix_input_real8(name, val, val_default, if_default)
    implicit none
    character(*),intent(in) :: name
    real(8),intent(out) :: val(:,:)
    real(8),intent(in),optional :: val_default(:,:)
    logical,intent(out),optional :: if_default
    character(len_max) :: val_t
    logical :: if_found
    integer :: index_equal, length_trimed
    integer :: i_lbound, i_ubound, i
    
    if(if_root_global)then
      i_lbound = lbound(val,1); i_ubound = ubound(val,1)

      if(present(val_default))val = val_default
      call lookup_input(name, if_found)
      if(if_found)then
        do i = i_lbound, i_ubound
          read(id_inputfile, *) val(i,:)
        end do
      end if
    end if
    call comm_bcast(val)

    if(present(if_default))then
      if_default = .not.if_found
      call comm_bcast(if_default)
    end if

  end subroutine read_matrix_input_real8

!-------------------------------------------------------------------------------
  subroutine lookup_input(char_in, if_found, char_out)
    implicit none
    character(*),intent(in) :: char_in
    logical,intent(out) :: if_found
    character(len_max),intent(out),optional :: char_out
    character(len_max) :: char_t
    integer :: istat, len_t, index_equal

    if(if_root_global)then
      if_found = .false.
      len_t = len(char_in)
      rewind(id_inputfile)

      do 
        read(id_inputfile, '(a)', iostat=istat) char_t
        if(istat < 0)exit
        index_equal = index(char_t,'=')
        if(index_equal == 0)then
          if(char_in == trim(char_t))then
            if_found = .true.
            if(present(char_out))char_out = char_t
            return
          end if
        else
          if(char_in == trim(char_t(1:index_equal-1)))then
            if_found = .true.
            if(present(char_out))char_out = char_t
            return
          end if
        end if

      end do

    end if

  end subroutine lookup_input


!-------------------------------------------------------------------------------
end module inputoutput
