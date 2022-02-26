module math
  implicit none
  private

  public :: vdCorput_sequence

  real(8),parameter,public :: pi = 3.141592653589793d0
  complex(8),parameter,public :: zI=(0.d0,1.d0)

  contains
!-----------------------------------------
    subroutine vdCorput_sequence(n_in,nbase,vdc_res)
      implicit none
      integer,intent(in) :: n_in, nbase
      real(8),intent(out) :: vdc_res
      integer :: n
      real(8) :: bk
      
      n = n_in
      bk = 1d0/nbase
      vdc_res = 0d0
      
      do while(n>0)
        vdc_res = vdc_res + mod(n,nbase)*bk
        n = n/nbase
        bk = bk/nbase
        
      end do
      
    end subroutine vdCorput_sequence

end module math
