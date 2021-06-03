program main
  use parallel
  implicit none

  call init_parallel

  write(*,*)'myrank=',comm_id_global


  call fin_parallel

end program main
