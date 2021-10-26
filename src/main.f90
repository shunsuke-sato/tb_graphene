program main
  use parallel
  use inputoutput
  use electron_dynamics
  implicit none

  call init_parallel
  call init_inputoutput
  if(if_root_global)then
     write(*,"(A)")"!! Warning!! This simaulation is done with the single-band approximation, which is physically incorrect!!"
  end if
!  write(*,*)'myrank=',comm_id_global

  call calc_electron_dynamics

  if(if_root_global)then
     write(*,"(A)")"!! Warning!! This simaulation is done with the single-band approximation, which is physically incorrect!!"
  end if  
  call fin_inputoutput
  call fin_parallel

end program main
