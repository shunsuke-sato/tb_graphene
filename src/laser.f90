module laser
  use inputoutput
  use math
  use constants
  implicit none
  private

  public :: init_laser, calc_vector_potential_time

! laser 1
  logical :: if_pulse_1

  real(8) :: E0_pulse_1_x, E0_pulse_1_y, pulse_width_1
  real(8) :: pulse_center_1, time_delay_1
  real(8) :: omega_1, phi_CEP_1_x, phi_CEP_1_y

  contains
!-------------------------------------------------------------------------------
    subroutine init_laser
      implicit none
      real(8) :: E0_pulse_1_x_Vm, E0_pulse_1_y_Vm, pulse_width_1_fs
      real(8) :: pulse_center_1_fs, time_delay_1_fs
      real(8) :: omega_1_ev, phi_CEP_1_x_2pi, phi_CEP_1_y_2pi


      call read_basic_input('E0_pulse_1_x_Vm',E0_pulse_1_x_Vm,val_default = 0d0)
      call read_basic_input('E0_pulse_1_y_Vm',E0_pulse_1_y_Vm,val_default = 0d0)

      call read_basic_input('pulse_width_1_fs',pulse_width_1_fs,val_default = 0d0)
      call read_basic_input('pulse_center_1_fs',pulse_center_1_fs &
                             ,val_default = pulse_width_1_fs*0.5d0)
      call read_basic_input('time_delay_1_fs',time_delay_1_fs,val_default = 0d0)

      call read_basic_input('omega_1_ev',omega_1_ev,val_default = 0d0)

      call read_basic_input('phi_CEP_1_x_2pi',phi_CEP_1_x_2pi,val_default = 0d0)
      call read_basic_input('phi_CEP_1_y_2pi',phi_CEP_1_y_2pi,val_default = 0d0)

      E0_pulse_1_x = E0_pulse_1_x_Vm*ev/(angstrom*1d10)
      E0_pulse_1_y = E0_pulse_1_y_Vm*ev/(angstrom*1d10)
      pulse_width_1 = pulse_width_1_fs*fs
      pulse_center_1 = pulse_center_1_fs*fs
      time_delay_1 = time_delay_1_fs*fs
      omega_1        = omega_1_ev*ev
      phi_CEP_1_x    = phi_CEP_1_x_2pi*2d0*pi 
      phi_CEP_1_y    = phi_CEP_1_y_2pi*2d0*pi 

      if_pulse_1 = .true.
      if(abs(E0_pulse_1_x) + abs(E0_pulse_1_y) == 0) if_pulse_1 = .false.


    end subroutine init_laser
!-------------------------------------------------------------------------------
    subroutine calc_vector_potential_time(tt, Act_x, Act_y)
      implicit none
      real(8),intent(in) :: tt ! time
      real(8),intent(out) :: Act_x, Act_y
      real(8) :: xx

      Act_x = 0d0
      Act_y = 0d0

      if(if_pulse_1)then

        xx = tt-pulse_center_1
        if(abs(xx)<0.5d0*pulse_width_1)then
          Act_x = -(E0_pulse_1_x/omega_1)*sin(omega_1*xx+phi_CEP_1_x) &
            *cos(pi*xx/pulse_width_1)**4

          Act_y = -(E0_pulse_1_y/omega_1)*sin(omega_1*xx+phi_CEP_1_y) &
            *cos(pi*xx/pulse_width_1)**4
        end if
      end if


    end subroutine calc_vector_potential_time
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  

end module laser
