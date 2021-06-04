module electron_dynamics
  use parallel
  use communication
  use inputoutput
  use electronic_system
  use laser
  use constants
  implicit none
  private

  public :: calc_electron_dynamics

! time propagation parameter
  real(8) :: total_propagation_time, time_step
  integer :: num_time_step

  contains
!-------------------------------------------------------------------------------
    subroutine calc_electron_dynamics
      implicit none
      integer :: it

      call init_elec_system
      call init_electron_dynamics


      do it = 0, num_time_step

      end do

    end subroutine calc_electron_dynamics
!-------------------------------------------------------------------------------
    subroutine init_electron_dynamics
      implicit none
      real(8) :: total_propagation_time_fs

      call read_basic_input('total_propagation_time_fs' &
        ,total_propagation_time_fs,val_default = -1.0d0)
      total_propagation_time = total_propagation_time_fs*fs
      call read_basic_input('time_step', time_step,val_default = -1.0d0)

      if(if_root_global)then
        write(*,"(A,2x,e26.16e3)")'time_step (initial)=',time_step
        num_time_step = aint(total_propagation_time/time_step)+1
        num_time_step = max(num_time_step, 0)
        time_step = total_propagation_time/num_time_step

        write(*,"(A,2x,e26.16e3)")'time_step (refined)=',time_step
        write(*,"(A,2x,I9)")'num_time_step            =',num_time_step
      end if
      call comm_bcast(time_step)

      call init_laser

    end subroutine init_electron_dynamics
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


end module electron_dynamics
