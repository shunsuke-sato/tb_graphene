program main
  use parallel
  implicit none

  call init_parallel
  call init_inputoutput
  write(*,*)'myrank=',comm_id_global

  call fin_inputoutput
  call fin_parallel

end program main
