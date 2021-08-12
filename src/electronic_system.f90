module electronic_system
  use parallel
  use communication
  use inputoutput
  use math
  use constants
  implicit none
  private

  public :: init_elec_system &
           ,calc_commutator_2x2_gra_density_matrix &
           ,calc_eigv_2x2_gra &
           ,zfk_tb &
           ,calc_current_elec_system

! parameters are given by Rev. Mod. Phys. 81 109, (2009).

! lattice vectors
  real(8) :: a0_l ! C-C bond length
  real(8) :: a_l(2,2), b_l(2,2)

! nearest-enighbor vectors and hopping
  real(8),public :: delta_l(2,3), t_hop

! k-grids
  integer,public :: nk1, nk2, nk, nk_start, nk_end, nk_average, nk_remainder
  real(8),allocatable,public :: kx0(:),ky0(:),kx(:),ky(:)
  integer,allocatable,public :: ik_table(:,:)
  real(8) :: dkxy

! density matrix
  complex(8),allocatable,public :: zrho_dm(:,:,:)
  real(8),public :: t1_relax, t2_relax, t1_relax_fs, t2_relax_fs
  real(8),public :: mu_F, mu_F_ev
  real(8),public :: kbT, kbT_K, doping_per_site
  logical,public :: if_chemical_potential_is_given


  real(8),allocatable :: eps_dist(:,:),occ_dist(:,:)
  contains
!-------------------------------------------------------------------------------
    subroutine init_elec_system
      implicit none
      integer :: ik1, ik2, ik
      integer :: id_file
      logical :: if_default_mu_F, if_default_doping_per_site

! set prameters

      a0_l = 1.42d0*angstrom

      a_l(1,1) = 3d0*a0_l/2d0
      a_l(2,1) = sqrt(3d0)*a0_l/2d0

      a_l(1,2) = 3d0*a0_l/2d0
      a_l(2,2) = -sqrt(3d0)*a0_l/2d0

      b_l(1,1) = 2d0*pi/(3d0*a0_l)
      b_l(2,1) = 2d0*pi*sqrt(3d0)/(3d0*a0_l)

      b_l(1,2) = 2d0*pi/(3d0*a0_l)
      b_l(2,2) = -2d0*pi*sqrt(3d0)/(3d0*a0_l)


      delta_l(1,1) = a0_l/2d0
      delta_l(2,1) = a0_l*sqrt(3d0)/2d0

      delta_l(1,2) = a0_l/2d0
      delta_l(2,2) = -a0_l*sqrt(3d0)/2d0

      delta_l(1,3) = -a0_l
      delta_l(2,3) = 0d0

      t_hop = 2.8d0*ev

! read prameters
      call read_basic_input('nk1',nk1,val_default = -1)
      call read_basic_input('nk2',nk2,val_default = -1)

      call read_basic_input('T1_relax_fs',T1_relax_fs,val_default = -1d0)
      call read_basic_input('T2_relax_fs',T2_relax_fs,val_default = -1d0)

      T1_relax = T1_relax_fs*fs
      T2_relax = T2_relax_fs*fs

      call read_basic_input('mu_F_ev',mu_F_ev,val_default = 0d0 &
        ,if_default = if_default_mu_F)
      mu_F = mu_F_ev*ev
      if_chemical_potential_is_given = .not. if_default_mu_F

      call read_basic_input('doping_per_site',doping_per_site,val_default = 0d0 &
        ,if_default = if_default_doping_per_site)

      if((.not.if_default_mu_F) .and. (.not.if_default_doping_per_site))then
        message(1)= &
          'Error message: mu_F_ev and doping_per_site cannot be used simultaneously.'
        call error_finalize(message(1))
      end if

      call read_basic_input('kbT_K',kbT_K,val_default = 0d0)
      kbT  = kbT_K/11604.505d0*ev

! init grids and parallelization
      nk = nk1*nk2

      nk_average = nk/comm_nproc_global
      nk_remainder = mod(nk,comm_nproc_global)
      if(comm_id_global+1 <= nk_remainder)then
        nk_start = 1 + comm_id_global*(nk_average+1)
        nk_end   = nk_start + (nk_average + 1) -1
      else
        nk_start = 1 + nk_remainder*(nk_average+1) + nk_average*(comm_id_global - nk_remainder)
        nk_end    = nk_start + nk_average  -1
      end if


      allocate(kx(nk),ky(nk))
      allocate(kx0(nk),ky0(nk))    
      allocate(zrho_dm(2,2,nk_start:nk_end))
      allocate(ik_table(0:nk1-1,0:nk2-1))
      allocate(eps_dist(2,nk),occ_dist(2,nk))

      dkxy = abs(b_l(1,1)*b_l(2,2)-b_l(2,1)*b_l(1,2))/nk
      
      ! initialize k-grids
      ik = 0
      do ik1 = 0, nk1-1
        do ik2 = 0, nk2-1
          ik = ik + 1
          
          kx0(ik) = b_l(1,1)*ik1/dble(nk1) + b_l(1,2)*ik2/dble(nk2)
          ky0(ik) = b_l(2,1)*ik1/dble(nk1) + b_l(2,2)*ik2/dble(nk2)
 
          ik_table(ik1,ik2) = ik
          
        end do
      end do

      kx = kx0
      ky = ky0

      call initialize_density_matrix
      if(if_root_global)then
        call get_newfile_id(id_file)
        open(id_file, file='electronic_structure.out')
        do ik1 = 0, nk1-1
          do ik2 = 1, nk2-1
            write(id_file,"(999e26.16e3)")kx(ik_table(ik1,ik2)) &
                                         ,ky(ik_table(ik1,ik2)) &
                                         ,eps_dist(:,ik_table(ik1,ik2)) &
                                         ,occ_dist(:,ik_table(ik1,ik2))
          end do
          write(id_file,*)
        end do
        close(id_file)
      end if

      
    end subroutine init_elec_system
!-------------------------------------------------------------------------------
    subroutine initialize_density_matrix
      implicit none
      integer :: ik, iter
      real(8) :: kx_t, ky_t
      complex(8) :: zfk_t, zalpha
      complex(8) :: zeigv(2,2)
      real(8) :: eig(2), occ(2)
      real(8) :: mu_F_min, mu_F_max, pop_per_site, pop_per_site_ref

      eps_dist = 0d0
      occ_dist = 0d0

      do ik = nk_start, nk_end
        kx_t = kx(ik)
        ky_t = ky(ik)

        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t
        call calc_eigv_2x2_gra(zalpha, zeigv, eig)
        occ(1) = Fermi_Dirac_distribution(eig(1), mu_F, kbT)
        occ(2) = Fermi_Dirac_distribution(eig(2), mu_F, kbT)

        zrho_dm(1,1, ik) = occ(1)*abs(zeigv(1,1))**2 &
                         + occ(2)*abs(zeigv(1,2))**2

        zrho_dm(2,1, ik) = occ(1)*zeigv(2,1)*conjg(zeigv(1,1)) &
                         + occ(2)*zeigv(2,2)*conjg(zeigv(1,2))

        zrho_dm(2,2, ik) = occ(1)*abs(zeigv(2,1))**2 &
                         + occ(2)*abs(zeigv(2,2))**2

        zrho_dm(1,2, ik) = conjg(zrho_dm(2,1, ik))


        eps_dist(:,ik) = eig(:)
        occ_dist(:,ik) = occ(:)

      end do

      call comm_allreduce(eps_dist)
      call comm_allreduce(occ_dist)

      if(if_chemical_potential_is_given)return
      pop_per_site_ref = 1d0 + doping_per_site

! defining chemical potential by bi-section method
      if(doping_per_site>1d0)then
        message(1)= &
          'Error message: doping_per_site is too large (>1).'
        call error_finalize(message(1))
      else if(doping_per_site<-1d0)then
        message(1)= &
          'Error message: doping_per_site is too small (<-1).'
        call error_finalize(message(1))
      end if


      pop_per_site = 0.5d0*2d0*sum(occ_dist)/nk
      call comm_bcast(pop_per_site)
      if(pop_per_site > pop_per_site_ref)then
        mu_F_max = mu_F
        mu_F = minval(eps_dist)
        call comm_bcast(mu_F)

        occ_dist = 0d0
        do ik = nk_start, nk_end
          occ_dist(1,ik) = Fermi_Dirac_distribution(eps_dist(1,ik), mu_F, kbT)
          occ_dist(2,ik) = Fermi_Dirac_distribution(eps_dist(2,ik), mu_F, kbT)
        end do
        call comm_allreduce(occ_dist)
        pop_per_site = 0.5d0*2d0*sum(occ_dist)/nk
        call comm_bcast(pop_per_site)
        if(pop_per_site > pop_per_site_ref)then
          message(1)= &
            'Error message: something wrong happans in the doping parameters (see initialize_density_matrix routine).'
          message(2)= 'pop_per_site > doping_per_site'
          call error_finalize(message(1:2))
        else
          mu_F_min = mu_F
        end if

      else
        mu_F_min = mu_F
        mu_F = maxval(eps_dist)
        call comm_bcast(mu_F)

        occ_dist = 0d0
        do ik = nk_start, nk_end
          occ_dist(1,ik) = Fermi_Dirac_distribution(eps_dist(1,ik), mu_F, kbT)
          occ_dist(2,ik) = Fermi_Dirac_distribution(eps_dist(2,ik), mu_F, kbT)
        end do
        call comm_allreduce(occ_dist)
        pop_per_site = 0.5d0*2d0*sum(occ_dist)/nk
        call comm_bcast(pop_per_site)
        if(pop_per_site <= pop_per_site_ref)then
          message(1)= &
            'Error message: something wrong happans in the doping parameters (see initialize_density_matrix routine).'
          message(2)= 'pop_per_site <= doping_per_site'
          call error_finalize(message(1:2))
        else
          mu_F_max = mu_F
        end if
      end if


! serch minimum
      iter = 0
      do
        iter = iter + 1
        mu_F = 0.5d0*(mu_F_max + mu_F_min)
        occ_dist = 0d0
        do ik = nk_start, nk_end
          occ_dist(1,ik) = Fermi_Dirac_distribution(eps_dist(1,ik), mu_F, kbT)
          occ_dist(2,ik) = Fermi_Dirac_distribution(eps_dist(2,ik), mu_F, kbT)
        end do
        call comm_allreduce(occ_dist)
        pop_per_site = 0.5d0*2d0*sum(occ_dist)/nk
        call comm_bcast(pop_per_site)
        if(pop_per_site > pop_per_site_ref)then
          mu_F_max = mu_F
        else
          mu_F_min = mu_F
        end if

        if(abs(mu_F_max-mu_F_min) < 1d-10)exit
        if(iter > 1000)then
          message(1)= &
            'Error message: chemical potential search does not converge.'
          call error_finalize(message(1))
        end if

      end do

      mu_F = 0.5d0*(mu_F_max + mu_F_min)

      if(if_root_global)then
        write(*,"(A,2x,I7)")'Chemical potential is computed with the iteration of',iter
        write(*,"(A,2x,e26.16e3,2x,A)")'The computed chemical potential is ',mu_F,"a.u."
        write(*,"(A,2x,e26.16e3,2x,A)")'The computed chemical potential is ',mu_F/ev*1d3,"meV."
      end if


      eps_dist = 0d0
      occ_dist = 0d0

      do ik = nk_start, nk_end
        kx_t = kx(ik)
        ky_t = ky(ik)

        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t
        call calc_eigv_2x2_gra(zalpha, zeigv, eig)
        occ(1) = Fermi_Dirac_distribution(eig(1), mu_F, kbT)
        occ(2) = Fermi_Dirac_distribution(eig(2), mu_F, kbT)

        zrho_dm(1,1, ik) = occ(1)*abs(zeigv(1,1))**2 &
                         + occ(2)*abs(zeigv(1,2))**2

        zrho_dm(2,1, ik) = occ(1)*zeigv(2,1)*conjg(zeigv(1,1)) &
                         + occ(2)*zeigv(2,2)*conjg(zeigv(1,2))

        zrho_dm(2,2, ik) = occ(1)*abs(zeigv(2,1))**2 &
                         + occ(2)*abs(zeigv(2,2))**2

        zrho_dm(1,2, ik) = conjg(zrho_dm(2,1, ik))


        eps_dist(:,ik) = eig(:)
        occ_dist(:,ik) = occ(:)

      end do

      call comm_allreduce(eps_dist)
      call comm_allreduce(occ_dist)

    end subroutine initialize_density_matrix
!-------------------------------------------------------------------------------
    subroutine calc_current_elec_system(jxy_t)
      implicit none
      real(8),intent(out) :: jxy_t(2)
      complex(8) :: zjx_op(2,2), zjy_op(2,2)
      complex(8) :: zdfk_dk(2), zmat_tmp(2,2)
      real(8) :: kx_t, ky_t
      integer :: ik

      zjx_op = 0d0
      zjy_op = 0d0

      jxy_t = 0d0
      do ik = nk_start, nk_end
        kx_t = kx(ik)
        ky_t = ky(ik)

        zdfk_dk(:) = zdfk_dk_tb(kx_t, ky_t)
        zjx_op(1,2) = t_hop*zdfk_dk(1) ! *(-1)*(-1)
        zjx_op(2,1) = conjg(zjx_op(2,1))

        zjy_op(1,2) = t_hop*zdfk_dk(2) ! *(-1)*(-1)
        zjy_op(2,1) = conjg(zjy_op(2,1))

        zmat_tmp = matmul(zjx_op,zrho_dm(:,:,ik))
        jxy_t(1) = jxy_t(1) + zmat_tmp(1,1)+ zmat_tmp(2,2)

        zmat_tmp = matmul(zjy_op,zrho_dm(:,:,ik))
        jxy_t(2) = jxy_t(2) + zmat_tmp(1,1)+ zmat_tmp(2,2)

      end do

      call comm_allreduce(jxy_t)
      jxy_t = 2d0*jxy_t*dkxy/(2d0*pi)**2

    end subroutine calc_current_elec_system
!-------------------------------------------------------------------------------
! diagonalize a special matrix
! ( 0,       za)
! (conjg(za), 0)
    subroutine calc_eigv_2x2_gra(zalpha, zeigv, eig)
      implicit none
      complex(8),intent(in) :: zalpha
      complex(8),intent(out) :: zeigv(2,2)
      real(8),intent(out) :: eig(2)
      real(8) :: aa, phi

      aa =abs(zalpha)
      if(aa /= 0d0)then
        phi = -aimag(log(zalpha))
        eig(1) = aa
        eig(2) = -aa
        zeigv(1,1)=1d0/sqrt(2d0)
        zeigv(2,1)=exp(zi*phi)/sqrt(2d0)
        zeigv(1,2)=1d0/sqrt(2d0)
        zeigv(2,2)=-exp(zi*phi)/sqrt(2d0)
      else
        eig(1) = 0d0
        eig(2) = 0d0
        zeigv(1,1)=1d0/sqrt(2d0)
        zeigv(2,1)=1d0/sqrt(2d0)
        zeigv(1,2)=1d0/sqrt(2d0)
        zeigv(2,2)=-1d0/sqrt(2d0)
      end if
      

    end subroutine calc_eigv_2x2_gra
!-------------------------------------------------------------------------------
    function zfk_tb(kx_t, ky_t) result(zfk_result)
      implicit none
      real(8),intent(in) :: kx_t, ky_t
      complex(8) :: zfk_result

      zfk_result = exp(zi*(delta_l(1,1)*kx_t + delta_l(2,1)*ky_t)) &
                 + exp(zi*(delta_l(1,2)*kx_t + delta_l(2,2)*ky_t)) &
                 + exp(zi*(delta_l(1,3)*kx_t + delta_l(2,3)*ky_t)) 

    end function zfk_tb
!-------------------------------------------------------------------------------
    function zdfk_dk_tb(kx_t, ky_t) result(zdfk_dk_result)
      implicit none
      real(8),intent(in) :: kx_t, ky_t
      complex(8) :: zdfk_dk_result(2)

      zdfk_dk_result(:) = zi*( &
        delta_l(:,1)*exp(zi*(delta_l(1,1)*kx_t + delta_l(2,1)*ky_t)) &
      + delta_l(:,2)*exp(zi*(delta_l(1,2)*kx_t + delta_l(2,2)*ky_t)) &
      + delta_l(:,3)*exp(zi*(delta_l(1,3)*kx_t + delta_l(2,3)*ky_t)))

    end function zdfk_dk_tb
!-------------------------------------------------------------------------------
! diagonalize a special matrix
! H = ( 0,       za),  compute [H, rho]
!     (conjg(za), 0)
    subroutine calc_commutator_2x2_gra_density_matrix(zalpha, zdrho_dm_t, zcomm_out)
      implicit none
      complex(8),intent(in) :: zalpha, zdrho_dm_t(2,2)
      complex(8),intent(out) :: zcomm_out(2,2)
      real(8) :: rho22_rho11
      complex(8) :: zs

      zs = zalpha*zdrho_dm_t(2,1)-conjg(zalpha)*zdrho_dm_t(1,2)
      rho22_rho11 = zdrho_dm_t(2,2)-zdrho_dm_t(1,1)

      zcomm_out(1,1) = zs
      zcomm_out(2,1) = -conjg(zalpha)*rho22_rho11
      zcomm_out(1,2) = zalpha*rho22_rho11
      zcomm_out(2,2) = -zs


    end subroutine calc_commutator_2x2_gra_density_matrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    


end module electronic_system

