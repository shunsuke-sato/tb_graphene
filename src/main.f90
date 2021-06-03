program main
  use parallel
  use inputoutput
  use electronic_system
  implicit none

  call init_parallel
  call init_inputoutput
  write(*,*)'myrank=',comm_id_global

  call calc_electron_dynamics

  call fin_inputoutput
  call fin_parallel

end program main
