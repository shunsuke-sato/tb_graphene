module laser
  use inputoutput
  use math
  use constants
  implicit none
  private

  public :: init_laser, &
            calc_vector_potential_time, &
            calc_electric_field_time

! laser 1
  logical :: if_pulse_1

  real(8) :: E0_pulse_1_x, E0_pulse_1_y, pulse_width_1
  real(8) :: E0_pulse_1, pol_theta_1
  real(8) :: pulse_center_1, time_delay_1
  real(8) :: omega_1, phi_CEP_1_x, phi_CEP_1_y

! laser 2
  logical :: if_pulse_2

  real(8) :: E0_pulse_2_x, E0_pulse_2_y, pulse_width_2
  real(8) :: E0_pulse_2, pol_theta_2
  real(8) :: pulse_center_2, time_delay_2
  real(8) :: omega_2, phi_CEP_2_x, phi_CEP_2_y

! impulse
  logical :: if_impulse_1
  real(8) :: A0_impulse_1_x, A0_impulse_1_y
  real(8) :: impulse_center_1

! dc-field
  logical :: if_dc_field
  real(8) :: E0_dc_x, E0_dc_y


  contains
!-------------------------------------------------------------------------------
    subroutine init_laser
      implicit none
! laser 1
      real(8) :: E0_pulse_1_x_Vm, E0_pulse_1_y_Vm, pulse_width_1_fs
      real(8) :: E0_pulse_1_Vm, pol_theta_1_2pi
      real(8) :: pulse_center_1_fs, time_delay_1_fs
      real(8) :: omega_1_ev, phi_CEP_1_x_2pi, phi_CEP_1_y_2pi
! laser 2
      real(8) :: E0_pulse_2_x_Vm, E0_pulse_2_y_Vm, pulse_width_2_fs
      real(8) :: E0_pulse_2_Vm, pol_theta_2_2pi
      real(8) :: pulse_center_2_fs, time_delay_2_fs
      real(8) :: omega_2_ev, phi_CEP_2_x_2pi, phi_CEP_2_y_2pi
! impulse 1
      real(8) :: impulse_center_1_fs
! dc-field
      real(8) :: E0_dc_x_Vm, E0_dc_y_Vm

! laser 1
      call read_basic_input('E0_pulse_1_x_Vm',E0_pulse_1_x_Vm,val_default = 0d0)
      call read_basic_input('E0_pulse_1_y_Vm',E0_pulse_1_y_Vm,val_default = 0d0)

      call read_basic_input('E0_pulse_1_Vm',E0_pulse_1_Vm,val_default = 0d0)
      call read_basic_input('pol_theta_1_2pi',pol_theta_1_2pi,val_default = 0d0)

      call read_basic_input('pulse_width_1_fs',pulse_width_1_fs,val_default = 0d0)
      call read_basic_input('pulse_center_1_fs',pulse_center_1_fs &
                             ,val_default = pulse_width_1_fs*0.5d0)
      call read_basic_input('time_delay_1_fs',time_delay_1_fs,val_default = 0d0)

      call read_basic_input('omega_1_ev',omega_1_ev,val_default = 1d0)

      call read_basic_input('phi_CEP_1_x_2pi',phi_CEP_1_x_2pi,val_default = 0d0)
      call read_basic_input('phi_CEP_1_y_2pi',phi_CEP_1_y_2pi,val_default = 0d0)


      E0_pulse_1_x = E0_pulse_1_x_Vm*ev/(angstrom*1d10)
      E0_pulse_1_y = E0_pulse_1_y_Vm*ev/(angstrom*1d10)

      if(abs(E0_pulse_1_x) + abs(E0_pulse_1_y) == 0)then
        E0_pulse_1 = E0_pulse_1_Vm*ev/(angstrom*1d10)
        pol_theta_1 = pol_theta_1_2pi*2d0*pi

        E0_pulse_1_x = cos(pol_theta_1)*E0_pulse_1
        E0_pulse_1_y = sin(pol_theta_1)*E0_pulse_1
      end if

      pulse_width_1 = pulse_width_1_fs*fs
      pulse_center_1 = pulse_center_1_fs*fs
      time_delay_1 = time_delay_1_fs*fs
      omega_1        = omega_1_ev*ev
      phi_CEP_1_x    = phi_CEP_1_x_2pi*2d0*pi 
      phi_CEP_1_y    = phi_CEP_1_y_2pi*2d0*pi 

      if_pulse_1 = .true.
      if(abs(E0_pulse_1_x) + abs(E0_pulse_1_y) == 0) if_pulse_1 = .false.

! laser 2
      call read_basic_input('E0_pulse_2_x_Vm',E0_pulse_2_x_Vm,val_default = 0d0)
      call read_basic_input('E0_pulse_2_y_Vm',E0_pulse_2_y_Vm,val_default = 0d0)

      call read_basic_input('E0_pulse_2_Vm',E0_pulse_2_Vm,val_default = 0d0)
      call read_basic_input('pol_theta_2_2pi',pol_theta_2_2pi,val_default = 0d0)

      call read_basic_input('pulse_width_2_fs',pulse_width_2_fs,val_default = 0d0)
      call read_basic_input('pulse_center_2_fs',pulse_center_2_fs &
                             ,val_default = pulse_width_2_fs*0.5d0)
      call read_basic_input('time_delay_2_fs',time_delay_2_fs,val_default = 0d0)

      call read_basic_input('omega_2_ev',omega_2_ev,val_default = 1d0)

      call read_basic_input('phi_CEP_2_x_2pi',phi_CEP_2_x_2pi,val_default = 0d0)
      call read_basic_input('phi_CEP_2_y_2pi',phi_CEP_2_y_2pi,val_default = 0d0)

      E0_pulse_2_x = E0_pulse_2_x_Vm*ev/(angstrom*1d10)
      E0_pulse_2_y = E0_pulse_2_y_Vm*ev/(angstrom*1d10)

      if(abs(E0_pulse_2_x) + abs(E0_pulse_2_y) == 0)then
        E0_pulse_2 = E0_pulse_2_Vm*ev/(angstrom*1d10)
        pol_theta_2 = pol_theta_2_2pi*2d0*pi

        E0_pulse_2_x = cos(pol_theta_2)*E0_pulse_2
        E0_pulse_2_y = sin(pol_theta_2)*E0_pulse_2
      end if


      pulse_width_2 = pulse_width_2_fs*fs
      pulse_center_2 = pulse_center_2_fs*fs
      time_delay_2 = time_delay_2_fs*fs
      omega_2        = omega_2_ev*ev
      phi_CEP_2_x    = phi_CEP_2_x_2pi*2d0*pi 
      phi_CEP_2_y    = phi_CEP_2_y_2pi*2d0*pi 

      if_pulse_2 = .true.
      if(abs(E0_pulse_2_x) + abs(E0_pulse_2_y) == 0) if_pulse_2 = .false.

! impulse 1
      call read_basic_input('A0_impulse_1_x',A0_impulse_1_x,val_default = 0d0)
      call read_basic_input('A0_impulse_1_y',A0_impulse_1_y,val_default = 0d0)
      call read_basic_input('impulse_center_1_fs',impulse_center_1_fs,val_default = 0d0)
      impulse_center_1 = impulse_center_1_fs*fs

      if_impulse_1 = .true.
      if(abs(A0_impulse_1_x) + abs(A0_impulse_1_y) == 0) if_impulse_1 = .false.

! dc-field
      call read_basic_input('E0_dc_x_Vm',E0_dc_x_Vm,val_default = 0d0)
      call read_basic_input('E0_dc_y_Vm',E0_dc_y_Vm,val_default = 0d0)

      E0_dc_x = E0_dc_x_Vm*ev/(angstrom*1d10)
      E0_dc_y = E0_dc_y_Vm*ev/(angstrom*1d10)

      if_dc_field = .true.
      if(abs(E0_dc_x) + abs(E0_dc_y) == 0) if_dc_field = .false.

    end subroutine init_laser
!-------------------------------------------------------------------------------
    subroutine calc_vector_potential_time(tt, Act_x, Act_y)
      implicit none
      real(8),intent(in) :: tt ! time
      real(8),intent(out) :: Act_x, Act_y
      real(8) :: xx

      Act_x = 0d0
      Act_y = 0d0

! laser 1
      if(if_pulse_1)then

        xx = tt-pulse_center_1 - time_delay_1
        if(abs(xx)<0.5d0*pulse_width_1)then
          Act_x = Act_x -(E0_pulse_1_x/omega_1)*sin(omega_1*xx+phi_CEP_1_x) &
            *cos(pi*xx/pulse_width_1)**4

          Act_y = Act_y -(E0_pulse_1_y/omega_1)*sin(omega_1*xx+phi_CEP_1_y) &
            *cos(pi*xx/pulse_width_1)**4
        end if
      end if

! laser 2
      if(if_pulse_2)then

        xx = tt-pulse_center_1 - time_delay_1 - time_delay_2
        if(abs(xx)<0.5d0*pulse_width_2)then
          Act_x = Act_x -(E0_pulse_2_x/omega_2)*sin(omega_2*xx+phi_CEP_2_x) &
            *cos(pi*xx/pulse_width_2)**4

          Act_y = Act_y -(E0_pulse_2_y/omega_2)*sin(omega_2*xx+phi_CEP_2_y) &
            *cos(pi*xx/pulse_width_2)**4
        end if
      end if

! impulse 1
      if(if_impulse_1)then
        if(tt>= impulse_center_1)then
          Act_x = Act_x + A0_impulse_1_x
          Act_y = Act_y + A0_impulse_1_y
        end if
      end if

! dc-field
      if(if_dc_field)then
        Act_x = Act_x - E0_dc_x*tt
        Act_y = Act_y - E0_dc_y*tt
      end if

    end subroutine calc_vector_potential_time
!-------------------------------------------------------------------------------
    subroutine calc_electric_field_time(tt, dt, Et_x, Et_y)
      implicit none
      real(8),intent(in) :: tt, dt ! time
      real(8),intent(out) :: Et_x, Et_y
      real(8) :: Act_x_p, Act_y_p
      real(8) :: Act_x_n, Act_y_n
      real(8) :: xx

      xx = tt + 0.5d0*dt
      call calc_vector_potential_time(xx, Act_x_p, Act_y_p)

      xx = tt - 0.5d0*dt
      call calc_vector_potential_time(xx, Act_x_n, Act_y_n)

      Et_x = -(Act_x_p - Act_x_n)/dt
      Et_y = -(Act_y_p - Act_y_n)/dt

    end subroutine calc_electric_field_time
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  

end module laser
