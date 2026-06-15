!=====================================================
! The subroutine compute_energies calculates the kinetic energy
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_diagnostics
  use mod_precision
  use mod_config
  use mod_grid
  use mod_particles
  use mod_cic
  implicit none
contains
  subroutine compute_energies(cfg,grid,part,initial)
    type(config_t), intent(in) :: cfg
    type(grid_t), intent(in) :: grid
    type(particles_t), intent(inout) :: part
    logical, intent(in) :: initial
    integer :: ip
    real(dp) :: ekin
    type(field_sample_t) :: fs
    logical :: inside
    do ip=1,part%np
       if (.not. part%active(ip)) then
          ekin=-1.0_dp
       else
          call cic_sample_fields(grid,part%x(ip),part%y(ip),part%z(ip),fs,inside)
          if (.not. inside) then
             ekin=-1.0_dp
          else
             ekin=0.5_dp*cfg%m_particle*part%vpar(ip)**2 + part%mu(ip)*fs%bmag
          end if
       end if
       if (initial) then
          part%ekin_initial(ip)=ekin
       else
          part%ekin_final(ip)=ekin
       end if
    end do
  end subroutine compute_energies

  function count_active(part) result(nactive)
    type(particles_t), intent(in) :: part
    integer :: nactive, ip
    nactive=0
    do ip=1,part%np
       if (part%active(ip)) nactive=nactive+1
    end do
  end function count_active
end module mod_diagnostics
