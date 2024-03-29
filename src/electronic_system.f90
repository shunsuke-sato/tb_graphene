module electronic_system
  use mpi
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
           ,calc_current_elec_system &
           ,calc_energy_elec_system &
           ,calc_carrier_distribution &
           ,integrate_transition_current_density &
           ,output_transition_current_density


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

  logical,public :: if_output_kspace_distribution

! transition current density
  logical,public :: if_calc_transition_current_density
  complex(8),allocatable :: zjxy_tcd(:,:,:)
  real(8) :: omega_tcd

  contains
!-------------------------------------------------------------------------------
    subroutine init_elec_system
      implicit none
      integer :: ik1, ik2, ik
      integer :: id_file
      logical :: if_default_mu_F, if_default_doping_per_site
      real(8) :: Eelec_t
      real(8),allocatable :: kx_g(:),ky_g(:)
      real(8),allocatable :: eps_dist_g(:,:),occ_dist_g(:,:)
      real(8) :: omega_tcd_ev

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


      allocate(kx(nk_start:nk_end),ky(nk_start:nk_end))
      allocate(kx0(nk_start:nk_end),ky0(nk_start:nk_end))    
      allocate(zrho_dm(2,2,nk_start:nk_end))
      allocate(ik_table(0:nk1-1,0:nk2-1))
      allocate(eps_dist(2,nk_start:nk_end),occ_dist(2,nk_start:nk_end))

      dkxy = abs(b_l(1,1)*b_l(2,2)-b_l(2,1)*b_l(1,2))/nk
      
      ! initialize k-grids
      ik = 0
      do ik1 = 0, nk1-1
        do ik2 = 0, nk2-1
          ik = ik + 1

          if(ik >= nk_start .and. ik <= nk_end)then
             kx0(ik) = b_l(1,1)*ik1/dble(nk1) + b_l(1,2)*ik2/dble(nk2)
             ky0(ik) = b_l(2,1)*ik1/dble(nk1) + b_l(2,2)*ik2/dble(nk2)
          end if
 
          ik_table(ik1,ik2) = ik
          
        end do
      end do

      kx = kx0
      ky = ky0

      call initialize_density_matrix

      call read_basic_input('if_output_kspace_distribution' &
           ,if_output_kspace_distribution,val_default = .true.)

      if(if_output_kspace_distribution)then

         allocate(kx_g(nk),ky_g(nk))
         allocate(eps_dist_g(2,nk),occ_dist_g(2,nk))
         kx_g = 0d0; ky_g = 0d0
         kx_g(nk_start:nk_end) = kx(nk_start:nk_end)
         ky_g(nk_start:nk_end) = ky(nk_start:nk_end)
         eps_dist_g = 0d0; occ_dist_g = 0d0
         eps_dist_g(:,nk_start:nk_end) = eps_dist(:,nk_start:nk_end)
         occ_dist_g(:,nk_start:nk_end) = occ_dist(:,nk_start:nk_end)

         call comm_allreduce(kx_g)
         call comm_allreduce(ky_g)
         call comm_allreduce(eps_dist_g)
         call comm_allreduce(occ_dist_g)

         if(if_root_global)then
            call get_newfile_id(id_file)
            open(id_file, file='electronic_structure.out')
            do ik1 = 0, nk1-1
               do ik2 = 0, nk2-1
                  write(id_file,"(999e26.16e3)")kx_g(ik_table(ik1,ik2)) &
                                               ,ky_g(ik_table(ik1,ik2)) &
                                               ,eps_dist_g(:,ik_table(ik1,ik2)) &
                                               ,occ_dist_g(:,ik_table(ik1,ik2))
               end do
               write(id_file,*)
            end do
            close(id_file)
         end if
      end if

      call calc_energy_elec_system(Eelec_t)
      if(if_root_global)then
        write(*,"(A,2x,e26.16e3,2x,A)")'The initial electronic energy is ',Eelec_t,"a.u."
        write(*,"(A,2x,e26.16e3,2x,A)")'The initial electronic energy is ',Eelec_t/ev*1d3,"meV."
      end if

      call read_basic_input('if_calc_transition_current_density' &
           ,if_calc_transition_current_density,val_default = .false.)

      call read_basic_input('omega_tcd_ev' &
           ,omega_tcd_ev,val_default = 0d0)
      omega_tcd = omega_tcd_ev*ev

      if(if_calc_transition_current_density)then
        allocate(zjxy_tcd(0:nk1-1,0:nk2-1,2))
        zjxy_tcd = 0d0
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
      call comm_allreduce(pop_per_site)
      if(pop_per_site > pop_per_site_ref)then
        mu_F_max = mu_F
        mu_F = minval(eps_dist)
        call comm_allreduce(mu_F, method = MPI_MIN)

        occ_dist = 0d0
        do ik = nk_start, nk_end
          occ_dist(1,ik) = Fermi_Dirac_distribution(eps_dist(1,ik), mu_F, kbT)
          occ_dist(2,ik) = Fermi_Dirac_distribution(eps_dist(2,ik), mu_F, kbT)
        end do

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
        call comm_allreduce(mu_F, method = MPI_MAX)

        occ_dist = 0d0
        do ik = nk_start, nk_end
          occ_dist(1,ik) = Fermi_Dirac_distribution(eps_dist(1,ik), mu_F, kbT)
          occ_dist(2,ik) = Fermi_Dirac_distribution(eps_dist(2,ik), mu_F, kbT)
        end do

        pop_per_site = 0.5d0*2d0*sum(occ_dist)/nk
        call comm_allreduce(pop_per_site)
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

        pop_per_site = 0.5d0*2d0*sum(occ_dist)/nk
        call comm_allreduce(pop_per_site)
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

    end subroutine initialize_density_matrix
!-------------------------------------------------------------------------------
    subroutine calc_current_elec_system(jxy_t, jxy_intra_c_t, jxy_intra_v_t)
      implicit none
      real(8),intent(out) :: jxy_t(2),jxy_intra_c_t(2),jxy_intra_v_t(2)
      complex(8) :: zjx_op(2,2), zjy_op(2,2)
      complex(8) :: zdfk_dk(2), zmat_tmp(2,2)
      complex(8) :: zfk_t, zalpha
      complex(8) :: zeigv(2,2), zrho_dm_occ(2,2)
      real(8) :: eig(2), occ(2)
      real(8) :: kx_t, ky_t
      integer :: ik

      zjx_op = 0d0
      zjy_op = 0d0

      jxy_t = 0d0
      jxy_intra_c_t=0d0
      jxy_intra_v_t=0d0
      do ik = nk_start, nk_end
        kx_t = kx(ik)
        ky_t = ky(ik)

        zdfk_dk(:) = zdfk_dk_tb(kx_t, ky_t)
        zjx_op(1,2) = t_hop*zdfk_dk(1) ! *(-1)*(-1)
        zjx_op(2,1) = conjg(zjx_op(1,2))

        zjy_op(1,2) = t_hop*zdfk_dk(2) ! *(-1)*(-1)
        zjy_op(2,1) = conjg(zjy_op(1,2))

        zmat_tmp = matmul(zjx_op,zrho_dm(:,:,ik))
        jxy_t(1) = jxy_t(1) + zmat_tmp(1,1)+ zmat_tmp(2,2)

        zmat_tmp = matmul(zjy_op,zrho_dm(:,:,ik))
        jxy_t(2) = jxy_t(2) + zmat_tmp(1,1)+ zmat_tmp(2,2)

        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t
! H = ( 0,       za)
!     (conjg(za), 0)        
        call calc_eigv_2x2_gra(zalpha, zeigv, eig)

        zrho_dm_occ = matmul( transpose(conjg(zeigv)),  matmul(zrho_dm(:,:,ik), zeigv))
        occ(1) = zrho_dm_occ(1,1)
        occ(2) = zrho_dm_occ(2,2)

        jxy_intra_c_t(:) = jxy_intra_c_t(:) &
          - occ(1)*abs(t_hop)*zdfk_dk(:)*conjg(zfk_t)/abs(zfk_t)
        jxy_intra_v_t(:) = jxy_intra_v_t(:) &
          + occ(2)*abs(t_hop)*zdfk_dk(:)*conjg(zfk_t)/abs(zfk_t)

      end do

      call comm_allreduce(jxy_t)
      jxy_t = 2d0*jxy_t*dkxy/(2d0*pi)**2

      call comm_allreduce(jxy_intra_c_t)
      jxy_intra_c_t = 2d0*jxy_intra_c_t*dkxy/(2d0*pi)**2

      call comm_allreduce(jxy_intra_v_t)
      jxy_intra_v_t = 2d0*jxy_intra_v_t*dkxy/(2d0*pi)**2

    end subroutine calc_current_elec_system
!-------------------------------------------------------------------------------
    subroutine calc_energy_elec_system(Eelec_t)
      implicit none
      real(8),intent(out) :: Eelec_t
      complex(8) :: zfk_t, zalpha
      real(8) :: kx_t, ky_t
      integer :: ik

      Eelec_t = 0d0

      do ik = nk_start, nk_end
        kx_t = kx(ik)
        ky_t = ky(ik)

        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t

! H = ( 0,       za)
!     (conjg(za), 0)        

! Tr[H*zrho_dm]
        Eelec_t = Eelec_t + zalpha*zrho_dm(2,1,ik) &
                          + conjg(zalpha)*zrho_dm(1,2,ik)

!! debug
!        zham = 0d0
!        zham(1,1) = 0d0; zham(1,2) = zalpha
!        zham(2,1) = conjg(zalpha); zham(2,2) = 0d0
!
!        zmat = matmul(zham,zrho_dm(:,:,ik))
!        Eelec_t = Eelec_t + zmat(1,1)+zmat(2,2)

      end do

      call comm_allreduce(Eelec_t)
      Eelec_t = 2d0*Eelec_t*dkxy/(2d0*pi)**2

    end subroutine calc_energy_elec_system
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
    subroutine calc_carrier_distribution(cfilename_t)
      implicit none
      character(*),intent(in),optional :: cfilename_t
      integer :: ik, ik1, ik2, id_file_t
      real(8) :: kx_t, ky_t
      complex(8) :: zfk_t, zalpha
      complex(8) :: zeigv(2,2)
      real(8) :: eig(2), occ(2), occ_t(2)
      real(8) :: kx_g(nk),ky_g(nk)
      real(8) :: eps_dist_g(2,nk)
      real(8) :: occ_dist_g(2,nk)
      real(8) :: docc_dist_g(2,nk)

      kx_g = 0d0; kx_g = 0d0
      eps_dist_g = 0d0
      occ_dist_g = 0d0
      docc_dist_g = 0d0

      do ik = nk_start, nk_end
        kx_t = kx(ik)
        ky_t = ky(ik)

        kx_g(ik) = kx(ik)
        ky_g(ik) = ky(ik)

        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t
        call calc_eigv_2x2_gra(zalpha, zeigv, eig)
        occ_t(1) = Fermi_Dirac_distribution(eig(1), mu_F, kbT)
        occ_t(2) = Fermi_Dirac_distribution(eig(2), mu_F, kbT)


        occ(1) = sum(conjg(zeigv(:,1))* matmul(zrho_dm(:,:,ik),zeigv(:,1)))
        occ(2) = sum(conjg(zeigv(:,2))* matmul(zrho_dm(:,:,ik),zeigv(:,2)))

        eps_dist_g(:,ik) = eig(:)
        occ_dist_g(:,ik) = occ(:)

        docc_dist_g(:,ik) = occ(:)-occ_t(:)

      end do

      call comm_allreduce(kx_g)
      call comm_allreduce(ky_g)
      call comm_allreduce(eps_dist_g)
      call comm_allreduce(occ_dist_g)
      call comm_allreduce(docc_dist_g)

      if(if_root_global)then
         call get_newfile_id(id_file_t)
         open(id_file_t,file=trim(cfilename_t))
         write(id_file_t,"(A)")"#kx, ky, eps(1:2), occ(1:2), docc(1:2)"
         do ik1 = 0, nk1-1
            do ik2 = 0, nk2-1
               write(id_file_t,"(999e26.16e3)")kx_g(ik_table(ik1,ik2)) &
                                            ,ky_g(ik_table(ik1,ik2)) &
                                            ,eps_dist_g(:,ik_table(ik1,ik2)) &
                                            ,occ_dist_g(:,ik_table(ik1,ik2)) &
                                            ,docc_dist_g(:,ik_table(ik1,ik2))
            end do
            write(id_file_t,*)
         end do
         close(id_file_t)
      end if
      
    end subroutine calc_carrier_distribution
!-------------------------------------------------------------------------------
    subroutine integrate_transition_current_density(time_t)
      implicit none
      real(8),intent(in) :: time_t
      real(8) :: kx_t, ky_t, k1_t, k2_t
      real(8) :: x1,x2
      integer :: ik
      complex(8) :: zjx_op(2,2), zjy_op(2,2)
      complex(8) :: zdfk_dk(2), zmat_tmp(2,2)
      real(8) :: jxy_t(2)
      real(8) :: Binv(2,2), detB
      integer :: n1, n2, n1_p, n2_p
      real(8) :: dx1, dx2
      real(8) :: ww(0:1,0:1), ss

      detB = b_l(1,1)*b_l(2,2)-b_l(1,2)*b_l(2,1)
      Binv(1,1) = b_l(2,2)
      Binv(2,2) = b_l(1,1)
      Binv(1,2) = -b_l(1,2)
      Binv(2,1) = -b_l(2,1)
      Binv = Binv/detB

      do ik = nk_start, nk_end
        kx_t = kx(ik)
        ky_t = ky(ik)

        zdfk_dk(:) = zdfk_dk_tb(kx_t, ky_t)
        zjx_op(1,2) = t_hop*zdfk_dk(1) ! *(-1)*(-1)
        zjx_op(2,1) = conjg(zjx_op(1,2))

        zjy_op(1,2) = t_hop*zdfk_dk(2) ! *(-1)*(-1)
        zjy_op(2,1) = conjg(zjy_op(1,2))

        zmat_tmp = matmul(zjx_op,zrho_dm(:,:,ik))
        jxy_t(1) = zmat_tmp(1,1)+ zmat_tmp(2,2)

        zmat_tmp = matmul(zjy_op,zrho_dm(:,:,ik))
        jxy_t(2) = zmat_tmp(1,1)+ zmat_tmp(2,2)

        k1_t = Binv(1,1)*kx_t +Binv(1,2)*ky_t
        k2_t = Binv(2,1)*kx_t +Binv(2,2)*ky_t

        x1 = k1_t - int(k1_t)
        x2 = k2_t - int(k2_t)

        if(x1 < 0d0)x1 = x1+1d0
        if(x2 < 0d0)x2 = x2+1d0

        x1 = x1*nk1
        x2 = x2*nk2

        n1 = x1
        n2 = x2

        n1_p = n1 +1
        n2_p = n2 +1

        ww(0,0) = 1d0/(dble(x1-n1)**2 + dble(x2-n2)**2 + 1d-6)
        ww(1,0) = 1d0/(dble(x1-n1_p)**2 + dble(x2-n2)**2 + 1d-6)
        ww(0,1) = 1d0/(dble(x1-n1)**2 + dble(x2-n2_p)**2 + 1d-6)
        ww(1,1) = 1d0/(dble(x1-n1_p)**2 + dble(x2-n2_p)**2 + 1d-6)
        ss = sum(ww)
        ww = ww/ss


        if(n1 >= nk1)n1 = n1-nk1
        if(n2 >= nk2)n2 = n2-nk2
        if(n1_p >= nk1)n1_p = n1_p-nk1
        if(n2_p >= nk2)n2_p = n2_p-nk2

        zjxy_tcd(n1,n2,:) = zjxy_tcd(n1,n2,:) + ww(0,0)*jxy_t(:)*exp(zi*omega_tcd*time_t)
        zjxy_tcd(n1_p,n2,:) = zjxy_tcd(n1_p,n2,:) + ww(1,0)*jxy_t(:)*exp(zi*omega_tcd*time_t)
        zjxy_tcd(n1,n2_p,:) = zjxy_tcd(n1,n2_p,:) + ww(0,1)*jxy_t(:)*exp(zi*omega_tcd*time_t)
        zjxy_tcd(n1_p,n2_p,:) = zjxy_tcd(n1_p,n2_p,:) + ww(1,1)*jxy_t(:)*exp(zi*omega_tcd*time_t)

        
      end do


    end subroutine integrate_transition_current_density
!-------------------------------------------------------------------------------
    subroutine output_transition_current_density
      implicit none
      integer ::ik1, ik2
      real(8) :: kx_t, ky_t
      integer :: id_file_t

      call comm_allreduce(zjxy_tcd)
      zjxy_tcd = zjxy_tcd*2d0/(2d0*pi)**2

      if(if_root_global)then
        call get_newfile_id(id_file_t)
        open(id_file_t,file="transition_current_density.out")
        write(id_file_t,"(A,2x,999e26.16e3)")"# omega_tcd=",omega_tcd

        do ik1 = 0, nk1-1
          do ik2 = 0, nk2-1

            kx_t = b_l(1,1)*ik1/dble(nk1) + b_l(1,2)*ik2/dble(nk2)
            ky_t = b_l(2,1)*ik1/dble(nk1) + b_l(2,2)*ik2/dble(nk2)

            write(id_file_t,"(999e26.16e3)")kx_t, ky_t, zjxy_tcd(ik1,ik2,:)
          end do
          write(id_file_t,*)
        end do
        close(id_file_t)
      end if

    end subroutine output_transition_current_density
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    


end module electronic_system

