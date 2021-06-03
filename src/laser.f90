module laser
  use parallel
  use communication
  use inputoutput
  implicit none
  private

  logical :: if_pulse_1
! pulse 1
  real(8) :: E1_x, E1_y, pulse_width_1, pulse_center_1
  real(8) :: omega_1, phi_CEP_1_x, phi_CEP_1_y

  contains
!-------------------------------------------------------------------------------
    subroutine init_laser
      implicit none
      real(8) :: E1_x_Vm, E1_y_vm, pulse_width_1_fs, pulse_center_1_fs
      real(8) :: omega_1_ev, phi_CEP_1_x_2pi, phi_CEP_1_y_2pi


    end subroutine init_laser
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  

end module laser
