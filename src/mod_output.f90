!=====================================================
! The module contains different outputs
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_output
  use mod_precision
  use mod_config
  use mod_particles
  implicit none

contains

  subroutine ensure_output_dir(cfg)
    type(config_t), intent(in) :: cfg
    call execute_command_line('mkdir -p '//trim(cfg%output_dir))
  end subroutine ensure_output_dir

  subroutine write_outputs(cfg,part)
    type(config_t), intent(in) :: cfg
    type(particles_t), intent(in) :: part
    call ensure_output_dir(cfg)
    call write_particle_table(cfg,part)
    call write_energy_table(cfg,part)
    call write_energy_histogram(cfg,part)
  end subroutine write_outputs

  subroutine write_particle_table(cfg,part)
    type(config_t), intent(in) :: cfg
    type(particles_t), intent(in) :: part
    character(len=1024) :: filename
    integer :: unit, ip
    filename=trim(cfg%output_dir)//'/'//trim(cfg%output_prefix)//'_particles_final.dat'
    open(newunit=unit,file=trim(filename),status='replace',action='write')
    write(unit,'(A)') '# id x y z vpar mu active'
    do ip=1,part%np
       write(unit,'(I10,1X,5(ES24.15,1X),L1)') ip,part%x(ip),part%y(ip),part%z(ip),part%vpar(ip),part%mu(ip),part%active(ip)
    end do
    close(unit)
  end subroutine write_particle_table

  subroutine write_energy_table(cfg,part)
    type(config_t), intent(in) :: cfg
    type(particles_t), intent(in) :: part
    character(len=1024) :: filename
    integer :: unit, ip
    filename=trim(cfg%output_dir)//'/'//trim(cfg%output_prefix)//'_energy_table.dat'
    open(newunit=unit,file=trim(filename),status='replace',action='write')
    write(unit,'(A)') '# id E_initial E_final active'
    do ip=1,part%np
       write(unit,'(I10,1X,2(ES24.15,1X),L1)') ip,part%ekin_initial(ip),part%ekin_final(ip),part%active(ip)
    end do
    close(unit)
  end subroutine write_energy_table

  subroutine write_energy_histogram(cfg,part)
    type(config_t), intent(in) :: cfg
    type(particles_t), intent(in) :: part
    character(len=1024) :: filename
    integer :: unit, ip, ibin, nbins
    integer, allocatable :: hi(:), hf(:)
    real(dp) :: emin, emax, de, ecenter, ei, ef
    logical :: first
    nbins=cfg%nbins_energy
    allocate(hi(nbins),hf(nbins)); hi=0; hf=0
    first=.true.; emin=0.0_dp; emax=1.0_dp
    do ip=1,part%np
       ei=part%ekin_initial(ip); ef=part%ekin_final(ip)
       if (ei >= 0.0_dp .and. ef >= 0.0_dp) then
          if (first) then
             emin=min(ei,ef); emax=max(ei,ef); first=.false.
          else
             emin=min(emin,min(ei,ef)); emax=max(emax,max(ei,ef))
          end if
       end if
    end do
    if (emax <= emin) then
       if (emin == 0.0_dp) then; emin=-0.5_dp; emax=0.5_dp
       else; emin=emin-0.5_dp*abs(emin); emax=emax+0.5_dp*abs(emax)
       end if
    end if
    de=(emax-emin)/real(nbins,dp)
    do ip=1,part%np
       ei=part%ekin_initial(ip); ef=part%ekin_final(ip)
       if (ei >= 0.0_dp) then
          ibin=int((ei-emin)/de)+1; ibin=min(max(ibin,1),nbins); hi(ibin)=hi(ibin)+1
       end if
       if (ef >= 0.0_dp) then
          ibin=int((ef-emin)/de)+1; ibin=min(max(ibin,1),nbins); hf(ibin)=hf(ibin)+1
       end if
    end do
    filename=trim(cfg%output_dir)//'/'//trim(cfg%output_prefix)//'_energy_histogram.dat'
    open(newunit=unit,file=trim(filename),status='replace',action='write')
    write(unit,'(A)') '# E_center N_initial N_final'
    do ibin=1,nbins
       ecenter=emin+(real(ibin,dp)-0.5_dp)*de
       write(unit,'(ES24.15,1X,I12,1X,I12)') ecenter,hi(ibin),hf(ibin)
    end do
    close(unit)
    deallocate(hi,hf)
  end subroutine write_energy_histogram
end module mod_output
