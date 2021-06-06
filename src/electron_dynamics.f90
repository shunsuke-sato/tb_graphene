module electron_dynamics
  use parallel
  use communication
  use inputoutput
  use electronic_system
  use laser
  use math
  use constants
  implicit none
  private

  public :: calc_electron_dynamics

! time propagation parameter
  real(8) :: total_propagation_time, time_step
  integer :: num_time_step

  contains
!-------------------------------------------------------------------------------
    subroutine calc_electron_dynamics
      implicit none
      integer :: it, id_file_current
      real(8) :: jxy_t(2)
      real(8) :: tt, Act_x, Act_y

      call init_elec_system
      call init_electron_dynamics
      tt = 0d0
      call calc_vector_potential_time(tt, Act_x, Act_y)
      kx = kx0 + Act_x; ky = ky0 + Act_y
      call calc_current_elec_system(jxy_t)

      if(if_root_global)then
        call get_newfile_id(id_file_current)
        open(id_file_current, file='jt_act.out')
        write(id_file_current, "(A)")"# time, jx, jy, Acx, Acy"
        write(id_file_current,"(999e26.16e3)")tt,jxy_t(:),act_x,act_y
      end if

      do it = 0, num_time_step
        call dt_evolve_elec(it)
        tt = (it+1)*time_step
        call calc_vector_potential_time(tt, Act_x, Act_y)
        kx = kx0 + Act_x; ky = ky0 + Act_y
        call calc_current_elec_system(jxy_t)
        if(if_root_global)write(id_file_current,"(999e26.16e3)")tt,jxy_t(:),act_x,act_y
      end do

      if(if_root_global)then
        close(id_file_current)
      end if

    end subroutine calc_electron_dynamics
!-------------------------------------------------------------------------------
    subroutine init_electron_dynamics
      implicit none
      real(8) :: total_propagation_time_fs

      call read_basic_input('total_propagation_time_fs' &
        ,total_propagation_time_fs,val_default = -1.0d0)
      total_propagation_time = total_propagation_time_fs*fs
      call read_basic_input('time_step', time_step,val_default = -1.0d0)

      if(if_root_global)then
        write(*,"(A,2x,e26.16e3)")'time_step (initial)=',time_step
        num_time_step = aint(total_propagation_time/time_step)+1
        num_time_step = max(num_time_step, 0)
        time_step = total_propagation_time/num_time_step

        write(*,"(A,2x,e26.16e3)")'time_step (refined)=',time_step
        write(*,"(A,2x,I9)")'num_time_step            =',num_time_step
      end if
      call comm_bcast(time_step)
      call comm_bcast(num_time_step)

      call init_laser

    end subroutine init_electron_dynamics
!-------------------------------------------------------------------------------
! propagate dnesity matrix with RK4
    subroutine dt_evolve_elec(it)
      implicit none
      integer,intent(in) :: it
      integer :: ik
      real(8) :: tt, Act_x, Act_y
      complex(8) :: zrho_dm_t(2,2), zrho_dm_RK4(2,2,4)
      complex(8) :: zLrho_dm_t(2,2), zrho_col(2,2)
      real(8) :: kx_t, ky_t
      complex(8) :: zfk_t, zalpha
      complex(8) :: zeigv(2,2)
      real(8) :: eig(2), occ(2)

      do ik = nk_start, nk_end
! RK1
        tt = it*time_step
        zrho_dm_t(:,:) = zrho_dm(:,:,ik)

        call calc_vector_potential_time(tt, Act_x, Act_y)
        kx_t = kx0(ik) + Act_x
        ky_t = ky0(ik) + Act_y
        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t
        call calc_commutator_2x2_gra_density_matrix(zalpha, zrho_dm_t, zLrho_dm_t)
        zLrho_dm_t = -zi*zLrho_dm_t

        call calc_eigv_2x2_gra(zalpha, zeigv, eig)
        occ(1) = Fermi_Dirac_distribution(eig(1), mu_F, kbT)
        occ(2) = Fermi_Dirac_distribution(eig(2), mu_F, kbT)

! transform to instantaneous eigen basis
        zrho_col = matmul( transpose(conjg(zeigv)),  matmul(zrho_dm_t, zeigv))
        zrho_col(1,1) = -(zrho_col(1,1)-occ(1))/T1_relax
        zrho_col(1,2) = -zrho_col(1,2)/T2_relax
        zrho_col(2,1) = -zrho_col(2,1)/T2_relax
        zrho_col(2,2) = -(zrho_col(2,2)-occ(2))/T1_relax
! transform back to site basis
        zrho_col = matmul(zeigv,  matmul(zrho_col, transpose(conjg(zeigv))))

        zrho_dm_RK4(:,:,1) = zLrho_dm_t + zrho_col
! RK2
        tt = it*time_step +0.5d0*time_step
        zrho_dm_t(:,:) = zrho_dm(:,:,ik) + zrho_dm_RK4(:,:,1)*0.5d0*time_step

        call calc_vector_potential_time(tt, Act_x, Act_y)
        kx_t = kx0(ik) + Act_x
        ky_t = ky0(ik) + Act_y
        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t
        call calc_commutator_2x2_gra_density_matrix(zalpha, zrho_dm_t, zLrho_dm_t)
        zLrho_dm_t = -zi*zLrho_dm_t

        call calc_eigv_2x2_gra(zalpha, zeigv, eig)
        occ(1) = Fermi_Dirac_distribution(eig(1), mu_F, kbT)
        occ(2) = Fermi_Dirac_distribution(eig(2), mu_F, kbT)
! transform to instantaneous eigen basis
        zrho_col = matmul( transpose(conjg(zeigv)),  matmul(zrho_dm_t, zeigv))
        zrho_col(1,1) = -(zrho_col(1,1)-occ(1))/T1_relax
        zrho_col(1,2) = -zrho_col(1,2)/T2_relax
        zrho_col(2,1) = -zrho_col(2,1)/T2_relax
        zrho_col(2,2) = -(zrho_col(2,2)-occ(2))/T1_relax
! transform back to site basis
        zrho_col = matmul(zeigv,  matmul(zrho_col, transpose(conjg(zeigv))))

        zrho_dm_RK4(:,:,2) = zLrho_dm_t + zrho_col

! RK3
        tt = it*time_step +0.5d0*time_step
        zrho_dm_t(:,:) = zrho_dm(:,:,ik) + zrho_dm_RK4(:,:,2)*0.5d0*time_step

        call calc_vector_potential_time(tt, Act_x, Act_y)
        kx_t = kx0(ik) + Act_x
        ky_t = ky0(ik) + Act_y
        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t
        call calc_commutator_2x2_gra_density_matrix(zalpha, zrho_dm_t, zLrho_dm_t)
        zLrho_dm_t = -zi*zLrho_dm_t

        call calc_eigv_2x2_gra(zalpha, zeigv, eig)
        occ(1) = Fermi_Dirac_distribution(eig(1), mu_F, kbT)
        occ(2) = Fermi_Dirac_distribution(eig(2), mu_F, kbT)
! transform to instantaneous eigen basis
        zrho_col = matmul( transpose(conjg(zeigv)),  matmul(zrho_dm_t, zeigv))
        zrho_col(1,1) = -(zrho_col(1,1)-occ(1))/T1_relax
        zrho_col(1,2) = -zrho_col(1,2)/T2_relax
        zrho_col(2,1) = -zrho_col(2,1)/T2_relax
        zrho_col(2,2) = -(zrho_col(2,2)-occ(2))/T1_relax
! transform back to site basis
        zrho_col = matmul(zeigv,  matmul(zrho_col, transpose(conjg(zeigv))))

        zrho_dm_RK4(:,:,3) = zLrho_dm_t + zrho_col

! RK4
        tt = it*time_step + time_step
        zrho_dm_t(:,:) = zrho_dm(:,:,ik) + zrho_dm_RK4(:,:,3)*time_step

        call calc_vector_potential_time(tt, Act_x, Act_y)
        kx_t = kx0(ik) + Act_x
        ky_t = ky0(ik) + Act_y
        zfk_t = zfk_tb(kx_t, ky_t)
        zalpha = -t_hop*zfk_t
        call calc_commutator_2x2_gra_density_matrix(zalpha, zrho_dm_t, zLrho_dm_t)
        zLrho_dm_t = -zi*zLrho_dm_t

        call calc_eigv_2x2_gra(zalpha, zeigv, eig)
        occ(1) = Fermi_Dirac_distribution(eig(1), mu_F, kbT)
        occ(2) = Fermi_Dirac_distribution(eig(2), mu_F, kbT)
! transform to instantaneous eigen basis
        zrho_col = matmul( transpose(conjg(zeigv)),  matmul(zrho_dm_t, zeigv))
        zrho_col(1,1) = -(zrho_col(1,1)-occ(1))/T1_relax
        zrho_col(1,2) = -zrho_col(1,2)/T2_relax
        zrho_col(2,1) = -zrho_col(2,1)/T2_relax
        zrho_col(2,2) = -(zrho_col(2,2)-occ(2))/T1_relax
! transform back to site basis
        zrho_col = matmul(zeigv,  matmul(zrho_col, transpose(conjg(zeigv))))

        zrho_dm_RK4(:,:,4) = zLrho_dm_t + zrho_col

! RK final
        zrho_dm(:,:,ik) = zrho_dm(:,:,ik) + time_step/6d0*(&
          zrho_dm_RK4(:,:,1) + 2d0*zrho_dm_RK4(:,:,2) &
         +2d0*zrho_dm_RK4(:,:,3) + zrho_dm_RK4(:,:,4))

      end do


    end subroutine dt_evolve_elec
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


end module electron_dynamics
