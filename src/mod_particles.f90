!=====================================================
! The subroutine initialize_particles calculates particles position xp, vpar and mu
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_particles
  use mod_precision
  use mod_config
  use mod_grid
  use mod_cic
  use mod_random

  implicit none

  type :: particles_t
     integer                :: np=0
     real(dp), allocatable  :: x(:), y(:), z(:)
     real(dp), allocatable  :: vpar(:), mu(:)
     real(dp), allocatable  :: ekin_initial(:), ekin_final(:)
     logical, allocatable   :: active(:)
  end type particles_t

contains

  subroutine allocate_particles(part,np)
    type(particles_t), intent(inout)  :: part
    integer, intent(in)               :: np
    part%np=np
    allocate(part%x(np), part%y(np), part%z(np), part%vpar(np), part%mu(np))
    allocate(part%ekin_initial(np), part%ekin_final(np), part%active(np))
    part%ekin_initial=0.0_dp; part%ekin_final=0.0_dp; part%active=.true.
  end subroutine allocate_particles

  subroutine initialize_particles(cfg,grid,part)
    type(config_t), intent(in)        :: cfg
    type(grid_t), intent(in)          :: grid
    type(particles_t), intent(inout)  :: part
    integer                           :: ip
    real(dp)                          :: sigma_v, rx, ry, rz, vperp1, vperp2, vperp2_tot, eps
    type(field_sample_t)              :: fs
    logical                           :: inside

    call allocate_particles(part,cfg%npart)

    ! sigma_v = sqrt(kb * T / m)
    sigma_v = sqrt(cfg%kb*cfg%temperature/cfg%m_particle)
    eps = 1.0e-12_dp
    do ip=1,part%np
       rx=uniform_random(); ry=uniform_random(); rz=uniform_random()
       ! initial position xp = xmin + random*(xmax - xmin)
       part%x(ip)=cfg%xmin + rx*(cfg%xmax-cfg%xmin)
       part%y(ip)=cfg%ymin + ry*(cfg%ymax-cfg%ymin)
       part%z(ip)=cfg%zmin + rz*(cfg%zmax-cfg%zmin)
       ! if xp > xmax, xp = xmax - value*(xmax-xmin)
       if (part%x(ip) >= grid%xmax) part%x(ip)=grid%xmax-eps*(grid%xmax-grid%xmin)
       if (part%y(ip) >= grid%ymax) part%y(ip)=grid%ymax-eps*(grid%ymax-grid%ymin)
       if (part%z(ip) >= grid%zmax) part%z(ip)=grid%zmax-eps*(grid%zmax-grid%zmin)

       ! vpar = sigma_v * gaussian random
       part%vpar(ip)=sigma_v*gaussian_random()

       ! vperp1,2 = sigma_v*gaussian
       vperp1=sigma_v*gaussian_random(); vperp2=sigma_v*gaussian_random()

       ! vperp2_tot = vperp1**2 + vperp2**2
       vperp2_tot=vperp1*vperp1 + vperp2*vperp2

       ! to obtain Bmag in the particle position and calculate mu = 0.5 m vperp2_tot / Bmag
       call cic_sample_fields(grid,part%x(ip),part%y(ip),part%z(ip),fs,inside)
       if (.not. inside .or. fs%bmag <= 0.0_dp) then
          part%active(ip)=.false.; part%mu(ip)=0.0_dp
       else
          part%mu(ip)=0.5_dp*cfg%m_particle*vperp2_tot/fs%bmag
       end if
    end do
  end subroutine initialize_particles
end module mod_particles
