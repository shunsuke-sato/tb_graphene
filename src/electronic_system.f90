module electronic_system
  use parallel
  use communication
  use math
  implicit none
  private

! parameters are given by Rev. Mod. Phys. 81 109, (2009).

! lattice vectors
  real(8) :: a0 ! C-C bond length
  real(8) :: a_l(2,2), b_l(2,2)

! nearest-enighbor vectors and hopping
  real(8) :: delta_l(2,3), t_hop

! k-grids
  integer :: nk1, nk2, nk, nk_s, nk_e
  real(8),allocatable :: kxy(:,:)



end module electronic_system

