module constants
  implicit none

  real(8),parameter :: fs = 1d0/0.024189d0
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: angstrom = 1d0/0.52917721067d0
  real(8),parameter :: clight = 137.035999139d0

contains

  function Fermi_Dirac_distribution(eps,e_Fermi, kbT) result(y)
    real(8),intent(in) :: eps, e_Fermi, kbT
    real(8) :: y

    y = 0d0
    if(kbT <= 0d0)then
      if(eps <= e_Fermi)y=1d0
    else
      y = 1d0/(exp((eps - e_Fermi)/kbT) + 1d0)
    end if
    
  end function Fermi_Dirac_distribution

end module constants
