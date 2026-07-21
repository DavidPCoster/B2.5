!> core-edge actor: ETS <--> SOLPS/ITER
!!
!! \brief This implements a coupling between ETS and SOLPS-ITER
!!
!! \author David Coster (David.Coster@ipp.mpg.de)
!!
!! \copyright Copyright David Coster (Max Planck Institute for Plasma Physics).  
!!
!! \license{License}
!!  This work is licensed under the EUROPEAN UNION PUBLIC LICENCE v. 1.2.
!!  (https://joinup.ec.europa.eu/sites/default/files/custom-page/attachment/2020-03/EUPL-1.2%20EN.txt)
!!
!! @param[in]  core_profiles_in               [ids_core_profiles] providing the input core_profiles IDS
!! @param[in]  core_transport_in              [ids_core_transport] providing the input transport IDS
!! @param[in]  equilibrium_in                 [ids_equilibrium] providing the input equilibrium IDS
!! @param[in]  transport_solver_numerics_in   [ids_transport_solver_numerics] providing the input transport_solver_numerics IDS
!! @param[out] transport_solver_numerics_out  [ids_transport_solver_numerics] providing the output transport_solver_numerics IDS
!! @param[in]  codeparam                      [ids_parameters_input] providing the code parameters
!! @param[out] outputFlag                     [int] providing the return code (0 for OK)
!! @param[out] diagnosticInfo                 [pointer to char(*)] diagnostic information
!!
!! \todo implement the ability to set the type of boundary condition
!!
!! \todo implement a minimal physics based edge model
!!

module core_edge

  use b2mod_types , B2R8 => R8
  implicit none
  
  CHARACTER(len=255) :: old_path
  CHARACTER(len=255), save :: new_path, solps_directory
  integer, allocatable, save :: direction(:)
  integer, allocatable, save :: b2etsmap(:,:)
  real (kind=B2R8), allocatable, save :: bc_ce_na(:), bc_ec_na(:)
  real (kind=B2R8), allocatable, save :: bc_ce_ti(:)
  real (kind=B2R8), save :: bc_ce_te, bc_ec_te, bc_ec_ti
  real (kind=B2R8), allocatable :: te_ave_iy(:), ti_ave_iy(:), ne_ave_iy(:), na_ave_iy(:,:), vol_ave_iy(:)
  real (kind=B2R8), allocatable :: map_vol(:, :)

contains
  subroutine b2mn_ets(core_profiles_in, core_transport_in, equilibrium_in,  &
     transport_solver_numerics_in, transport_solver_numerics_out,  &
     codeparam, outputFlag, diagnosticInfo)

    use ids_schemas             ! IGNORE
    use ids_types               ! IGNORE
    use ids_routines            ! IGNORE
    use b2mod_constants
    use b2mod_types , B2R8 => R8
    use b2mod_b2cmpa
    use b2mod_boundary_namelist
    use b2mod_ual_io
    use b2mod_geo
    use b2mod_indirect
    use b2mod_plasma
    use b2mod_main

    implicit none

    type (ids_core_profiles)             :: core_profiles_in
    type (ids_core_transport)            :: core_transport_in
    type (ids_equilibrium)               :: equilibrium_in
    type (ids_transport_solver_numerics) :: transport_solver_numerics_in, transport_solver_numerics_out
    type (ids_parameters_input)          :: codeparam
    integer (ids_int),        intent(out)  :: outputFlag
    character(len=:), pointer, intent(out) :: diagnosticInfo

    integer, save :: niter
    real (kind=B2R8), save :: rxf_ce, rxf_ec

    real (kind=B2R8), allocatable :: tsn_sqrt_norm_vol(:)

    integer, save :: imodel, ieq, irho, irho_bc
    real (kind=B2R8) :: surface_area, ion_energy_flux, bc_rho, d_rho
    real (kind=B2R8) :: volave, neave, naave, teave, tiave, rxf

    real (kind=B2R8) :: navol_ave, centroid_r, centroid_z, tmp_vol_1, tmp_vol_2
    integer :: iy_core, centroid_n
    character (383) :: command
    integer :: EXITSTAT, CMDSTAT
    character (128) :: CMDMSG

    integer iloop, nrho, ix, iy, is0, is1, is2, inuc
    INTEGER :: return_status
    logical, save :: first, firstpass
    logical :: fileexist
    data first /.true./, firstpass /.true./

    call prgini ('ids_b2mn')
    outputFlag = 0
    nullify(diagnosticInfo)
    
    call getcwd(old_path)
    write(*,*) 'Found CWD to be ', trim(old_path)

! initialize on the first call only

    first_loop: if(first) then
! Note: we don't have access to nx, ny, ns until the call to b2mn_init below       
       allocate(direction(0:1024))
       direction=1
       niter=10
       rxf_ce=1.0_B2R8
       rxf_ec=1.0_B2R8
       solps_directory=''
       call assign_code_parameters(codeparam%parameters_value, return_status)
       solps_directory =  adjustl(solps_directory)
       write(*,*) 'solps_directory = ', trim(solps_directory)
       if(solps_directory.eq.'') then
          write(*,*) 'solps_directory not set -- using environment variable SOLPSPATH'
          call getenv("SOLPSPATH", new_path)
          write(*,*) 'SOLPSPATH = ', trim(new_path)
       else
          if(solps_directory(1:1).eq.'/') then
             new_path=solps_directory
          else
             call getenv("SOLPSTOP", new_path)
             new_path = trim(new_path) // '/' // trim(solps_directory)
          endif
       endif
       command = '../../../../bin/setup_SOLPS-ITER ' // new_path
       write(*,*) 'Executing the command ', trim(command)
       CMDMSG = ''
       CALL EXECUTE_COMMAND_LINE(command, .true., EXITSTAT, CMDSTAT, CMDMSG)
       write(*,*) 'EXITSTAT = ', EXITSTAT
       write(*,*) 'CMDSTAT = ', CMDSTAT
       write(*,*) 'CMDMSG = ', trim(CMDMSG)
!       call chdir(trim(new_path))
!       write(*,*) 'Set CWD to be ', trim(new_path)      
       call status_logger('Started the SOLPS initialization')
       open(99, file='Core-Edge.coupling')
       write(99,'(a)') 'Core-Edge data transfers'
       close(99)
       call b2mn_init
! Note: we now have access to nx, ny, ns
       first=.false.
       write(*,*) 'nx, ny, ns = ', nx, ny, ns
       allocate(te_ave_iy(-1:ny), ti_ave_iy(-1:ny), ne_ave_iy(-1:ny), na_ave_iy(-1:ny, 0:ns-1), vol_ave_iy(-1:ny))
       allocate(map_vol(-1:ny, 2))
       allocate(bc_ce_na(0:ns-1), bc_ec_na(0:ns-1), bc_ce_ti(0:ns-1))
       bc_ce_na = 0.0_B2r8
       bc_ec_na = 0.0_B2r8
       bc_ce_ti = 0.0_B2r8
       bc_ec_ti = 0.0_B2r8
       bc_ce_te = 0.0_B2r8
       bc_ec_te = 0.0_B2r8
       inquire(file='Core-Edge.dat',exist=fileexist)
       if(fileexist) then
          open(99, file='Core-Edge.dat')
          read(99,*) bc_ce_na, bc_ec_na, bc_ce_ti, bc_ec_ti, bc_ce_te, bc_ec_te
          close(99)
          firstpass=.false.
       endif
       write(*,*) 'From assign_code_parameters'
       write(*,*) 'NITER = ', niter
       write(*,*) 'DIRECTION = ', direction(1:ns)
       ! note that only one direction is coded (so far): fluxes from core; values to core
       write(*,*) '*******************************************************************************'
       write(*,*) '*** The only mapping implemented so far is fluxes from core, values to core ***'
       write(*,*) '*******************************************************************************'
       ! set up the mapping
       write(*,*) 'B2 ns = ', ns
       allocate(b2etsmap(-1:ns-1,0:4))
       ! b2etsmap(,0) == 1 if main ion, 2 if impurity
       ! b2etsmap(,1) == position in core_profiles
       ! b2etsmap(,2) == 0 if main ion, == charge state if impurity
       ! b2etsmap(,3) equation in transport_solver_numerics for density
       ! b2etsmap(,4) equation in transport_solver_numerics for temperature
       b2etsmap=-1    ! nothing found yet
       write(*,*) 'B2'
       do is0 = 0, ns-1
          write(*,'(i3,4(1x,f8.3))') is0, zn(is0), zamin(is0), zamax(is0), am(is0)
       end do
       call xertst(associated(core_profiles_in%profiles_1d), &
            'core_profiles_in%profiles_1d) must be associated')
       call xertst(size(core_profiles_in%profiles_1d) .eq. 1, &
            'core_profiles_in%profiles_1d length must be 1')
       if(associated(core_profiles_in%profiles_1d(1)%ion)) then
          write(*,*) 'CORE_PROFILES'
          do is1=1, SIZE(core_profiles_in%profiles_1d(1)%ion)
             associate (ion => core_profiles_in%profiles_1d(1)%ion(is1))
               if(ion%multiple_states_flag .le. 0) then
                  if(size(ion%element) .ne. 1) then
                     write(*,*) 'error: have only coded the case with 1 element in ion'
                     stop
                  endif
                  do is0=0,ns-1
                     if(abs(zn(is0)-ion%element(1)%z_n).lt.1e-1_B2R8.and.  &
                          abs(am(is0)-ion%element(1)%a).lt.1e-1_B2R8.and.  &
                          abs(zamin(is0)-ion%z_ion).lt.1e-1_B2R8.and.  &
                          abs(zamax(is0)-ion%z_ion).lt.1e-1_B2R8) then
                        if(b2etsmap(is0,1).lt.0) then
                           b2etsmap(is0,0)=1
                           b2etsmap(is0,1)=is1
                           b2etsmap(is0,2)=0
                        else
                           write(*,*) 'Overlap found ', is1, b2etsmap(is0,1)
                        endif
                     endif
                  enddo
                  write(*,'(i3,4(1x,f8.3))') is1, ion%element(1)%z_n, &
                       ion%z_ion, &
                       ion%z_ion, &
                       ion%element(1)%a
               else
                  if(size(ion%element) .ne. 1) then
                     write(*,*) 'error: have only coded the case with 1 element in ion'
                     stop
                  endif
                  do is2 = 1, size(ion%state)
                     associate (state => core_profiles_in%profiles_1d(1)%ion(is1)%state(is2))
                       do is0=0,ns-1
                          if(abs(zn(is0)-ion%element(1)%z_n).lt.1e-1_B2R8.and.  &
                               abs(am(is0)-ion%element(1)%a).lt.1e-1_B2R8.and.  &
                               abs(zamin(is0)-state%z_min).lt.1e-1_B2R8.and.  &
                               abs(zamax(is0)-state%z_max).lt.1e-1_B2R8) then
                             if(b2etsmap(is0,1).lt.0) then
                                b2etsmap(is0,0)=2
                                b2etsmap(is0,1)=is1
                                b2etsmap(is0,2)=is2
                             else
                                write(*,*) 'Overlap found ', is1, b2etsmap(is0,1)
                             endif
                          endif
                       enddo
                       write(*,'(i3,4(1x,f8.3))') is1, ion%element(1)%z_n, &
                            state%z_min, &
                            state%z_max, &
                            ion%element(1)%a
                     end associate
                  end do
               endif
             end associate
          enddo
       endif
! find the mappings to equations in transport_solver_numerics
       call xertst(associated(transport_solver_numerics_in%solver_1d), &
            'transport_solver_numerics_in%profiles_1d) must be associated')
       call xertst(size(transport_solver_numerics_in%solver_1d) .eq. 1, &
            'transport_solver_numerics_in%profiles_1d length must be 1')
       if (associated(transport_solver_numerics_in%solver_1d(1)%equation)) then
          do ieq = 1, size(transport_solver_numerics_in%solver_1d(1)%equation)
             associate (eq => transport_solver_numerics_in%solver_1d(1)%equation(ieq))
               write(*,*) ieq, trim(eq%primary_quantity%identifier%name(1)), &
                    eq%primary_quantity%ion_index, eq%primary_quantity%state_index, &
                    eq%primary_quantity%neutral_index
               select case (eq%primary_quantity%identifier%index)
               case (1)   ! psi erquation
                  write(*,*) 'Psi equation skipped'
               case (2)   ! density equation
                  if(eq%primary_quantity%ion_index .eq. 0) then   ! electrons
                     b2etsmap(-1,3) = ieq
                  else if(eq%primary_quantity%ion_index .gt. 0) then  ! ions
                     if (eq%primary_quantity%state_index .le. 0) then ! charge state within the ion
                        is1 = find_ets_ion_in_b2etsmap(eq%primary_quantity%ion_index, &
                             b2etsmap, ns)
                        if (b2etsmap(is1,1) .eq. eq%primary_quantity%ion_index) then
                           b2etsmap(is1,3) = ieq
                        else
                           write(*,*) 'Not equal', b2etsmap(is1,1), eq%primary_quantity%ion_index
                        endif
                     else
                        is1 = find_ets_ion_and_state_in_b2etsmap(eq%primary_quantity%ion_index, &
                             eq%primary_quantity%state_index, b2etsmap, ns)
                        if ((b2etsmap(is1,1) .eq. eq%primary_quantity%ion_index) .and. &
                             (b2etsmap(is1,2) .eq. eq%primary_quantity%state_index)) then
                           b2etsmap(is1,3) = ieq
                        else
                           write(*,*) 'Not equal', b2etsmap(is1,1), eq%primary_quantity%ion_index
                        endif
                     endif
                  endif
               case (3)   ! temperature equation
                  if(eq%primary_quantity%ion_index .eq. 0) then   ! electrons
                     b2etsmap(-1,4) = ieq
                  else if(eq%primary_quantity%ion_index .gt. 0) then  ! ions
                     if (eq%primary_quantity%state_index .le. 0) then ! charge state within the ion
                        is1 = find_ets_ion_in_b2etsmap(eq%primary_quantity%ion_index, &
                             b2etsmap, ns)
                        if (b2etsmap(is1,1) .eq. eq%primary_quantity%ion_index) then
                           b2etsmap(is1,4) = ieq
                        else
                           write(*,*) 'Not equal', b2etsmap(is1,1), eq%primary_quantity%ion_index
                        endif
                     else
                        is1 = find_ets_ion_and_state_in_b2etsmap(eq%primary_quantity%ion_index, &
                             eq%primary_quantity%state_index, b2etsmap, ns)
                        if ((b2etsmap(is1,1) .eq. eq%primary_quantity%ion_index) .and. &
                             (b2etsmap(is1,2) .eq. eq%primary_quantity%state_index)) then
                           b2etsmap(is1,4) = ieq
                        else
                           write(*,*) 'Not equal', b2etsmap(is1,1), eq%primary_quantity%ion_index
                        endif
                     endif
                  endif
               case (4)   ! toroidal rotation equation
                  write(*,*) 'Toroidal rotation equation skipped'
               case (5)   ! poloidal rotation equation
                  write(*,*) 'Poloidal rotation equation skipped'
               end select
             end associate
          end do
       end if
       write(*,*) 'B2-ETS Species Map'
       do is0 = -1, ns-1
          write(*,'(6(i3,1x))') is0, b2etsmap(is0,:)
       enddo

       ! find the model associated with the transport_solver
       call xertst(associated(core_transport_in%model), &
            'core_transport_in%model) must be associated')
       call xertst(size(core_transport_in%model) .ge. 1, &
            'core_transport_in%profiles_1d length must be 1')
       do imodel = 1, size(core_transport_in%model)
          if(core_transport_in%model(imodel)%identifier%index .eq. 2) then     ! ToDo: replace hard-coded 2 with the right logic
             exit
          endif
       enddo
       call xertst(imodel .le. size(core_transport_in%model), &
            'transport_solver not found in core-transport%models')
       
       call xertst(associated(core_transport_in%model(imodel)%profiles_1d), &
            'core_transport_in%model(imodel)%profiles_1d) must be associated')
       call xertst(size(core_transport_in%model(imodel)%profiles_1d) .eq. 1, &
            'core_transport_in%model(imodel)%profiles_1d length must be 1')

       call status_logger('Ended the SOLPS initialization')
    else
!       call chdir(trim(new_path))
!       write(*,*) 'Set CWD to be ', trim(new_path)
    endif first_loop
    
    ! map the information from core_transport_in onto the B2 data structures (CORE -> EDGE)
    
    call status_logger('Mapping the core data onto the edge')
    if(firstpass) then
       rxf = 1.0_B2r8
    else
       rxf = rxf_ce
    endif
    open(99, file='Core-Edge.coupling', position='append')
    associate (model => core_transport_in%model(imodel)%profiles_1d(1))
      ! electron power flux
      ieq = b2etsmap(-1,4)
      write(*,*) trim(transport_solver_numerics_in%solver_1d(1)%equation(ieq)%primary_quantity%identifier%name(1))
      bc_rho = transport_solver_numerics_in%solver_1d(1)%equation(ieq)%boundary_condition(2)%position
      write(*,*) bc_rho
      irho_bc = 1
      d_rho = abs(model%grid_flux%rho_tor_norm(irho_bc) - bc_rho)
      nrho=size(model%grid_flux%rho_tor_norm)
      do irho = 2, nrho
         if (abs(model%grid_flux%rho_tor_norm(irho) - bc_rho) .lt. d_rho) then
            irho_bc = irho
            d_rho = abs(model%grid_flux%rho_tor_norm(irho_bc) - bc_rho)
         endif
      end do
      write(*,*) 'Closest rho_tor_norm to ', bc_rho, ' is for position ', irho_bc, &
           ' with value ', model%grid_flux%rho_tor_norm(irho_bc)
      surface_area = model%grid_flux%surface(irho_bc)
      bcene(1) = 8   ! ToDo single null geometry assumed
      bc_ce_te = (1.0_b2r8 - rxf) * bc_ce_te + rxf * model%electrons%energy%flux(irho_bc) * surface_area
      enepar(1,1)=bc_ce_te
      write(*,'(a,1p,1(1x,g15.5))') 'Core->Edge: FLUX_EE', bc_ce_te
      write(99,'(a,1p,1(1x,g15.5))') 'Core->Edge: FLUX_EE', bc_ce_te
      ! ion power flux
      bceni(1)=8
      ion_energy_flux = 0.0_B2R8
      do is1 = 1, size(model%ion)
         is0 = find_ets_ion_in_b2etsmap(is1, b2etsmap, ns)
         if (is0 .ge. 0) then
            ieq = b2etsmap(is0,4)
            if (ieq .gt. 0) then
               write(*,*) trim(transport_solver_numerics_in%solver_1d(1)%equation(ieq)%primary_quantity%identifier%name(1))
               call xertst(transport_solver_numerics_in%solver_1d(1)%equation(ieq)%boundary_condition(2)%position .eq. bc_rho,  &
                    'Case of differing b.c. position not handled')
               ion_energy_flux = ion_energy_flux + model%ion(is1)%energy%flux(irho_bc)   ! ToDo: what about impurities?
            else
               write(*,*) 'is0 = ', is0, ' does not match an equation'
            end if
         else
            write(*,*) 'No match for is1 = ', is1
         end if
      end do
      bc_ce_ti(:) = (1.0_b2r8 - rxf) * bc_ce_ti(:) + rxf * ion_energy_flux * surface_area
      enipar(1,1) = bc_ce_ti(1)
      write(*,'(a,1p,1(1x,g15.5))') 'Core->Edge: FLUX_IE', bc_ce_ti(1)
      write(99,'(a,1p,1(1x,g15.5))') 'Core->Edge: FLUX_IE', bc_ce_ti(1)
      do is0 = 0, ns-1
         if(b2etsmap(is0,0).eq.1) then ! main ion
            is1 = b2etsmap(is0,1)
            ieq = b2etsmap(is0,3)
            if (ieq .gt. 0) then
               write(*,*) trim(transport_solver_numerics_in%solver_1d(1)%equation(ieq)%primary_quantity%identifier%name(1))
               call xertst(transport_solver_numerics_in%solver_1d(1)%equation(ieq)%boundary_condition(2)%position .eq. bc_rho,  &
                    'Case of differing b.c. position not handled')
               bccon(is0,1)=13
               bc_ce_na(is0) = (1.0_b2r8 - rxf) * bc_ce_na(is0) + rxf * model%ion(is1)%particles%flux(irho_bc) * surface_area
               conpar(is0,1,1) = bc_ce_na(is0)
               call xertst(conpar(is0,1,2) .gt. 0.0_R8, 'conpar(,,2) must be > 0')
               write(*,'(a,i3,1x,1p,g15.5)') 'Core->Edge: FLUX_NI', is0, bc_ce_na(is0)
               write(99,'(a,i3,1x,1p,g15.5)') 'Core->Edge: FLUX_NI', is0, bc_ce_na(is0)
            end if
         else if(b2etsmap(is0,0) .eq. 2) then ! impurity
            is1 = b2etsmap(is0,1)
            is2 = b2etsmap(is0,2)
            ieq = b2etsmap(is0,3)
            if (ieq .gt. 0) then
               write(*,*) trim(transport_solver_numerics_in%solver_1d(1)%equation(ieq)%primary_quantity%identifier%name(1))
               call xertst(transport_solver_numerics_in%solver_1d(1)%equation(ieq)%boundary_condition(2)%position .eq. bc_rho,  &
                    'Case of differing b.c. position not handled')
               bc_ce_na(is0) = (1.0_b2r8 - rxf) * bc_ce_na(is0) + &
                    rxf * model%ion(is1)%state(is2)%particles%flux(irho_bc) * surface_area
               if (.false.) then
                  bccon(is0,1)=13
                  conpar(is0,1,1) = bc_ce_na(is0)
                  call xertst(conpar(is0,1,2) .gt. 0.0_R8, 'conpar(,,2) must be > 0')
               else
                  bccon(is0,1)=8
                  conpar(is0,1,1) = bc_ce_na(is0)
               endif
               write(*,'(a,i3,1x,1p,g15.5)') 'Core->Edge: FLUX_NI', is0, bc_ce_na(is0)
               write(99,'(a,i3,1x,1p,g15.5)') 'Core->Edge: FLUX_NI', is0, bc_ce_na(is0)
            end if
         elseif(b2etsmap(is0,0).eq.3) then
            stop 'Case for neutrals not coded yet'
         else
            write(*,*) 'Case not coded yet', b2etsmap(is0,0)
         endif
      enddo
    end associate
    write(99,'(a)') 'Core->Edge: End'
    close(99)
    call write_b2mod_boundary_namelist
    
    ! advance the edge plasma
    call status_logger('Calling the edge code')
    call b2mn_step (niter)
    call status_logger('Finished the edge code')
    
    ! save the output to the edge cpo
    call status_logger('Transferring the edge data to the EDGE CPO')
! ToDo    allocate(edge_out(1))
! ToDo    call write_cpo(edge_out(1))
    
    ! map the information from the B2 data structures onto transport_solver_numerics_out (EDGE -> CORE)
    call status_logger('Mapping the edge data onto the core')
    if(firstpass) then
       rxf = 1.0_B2r8
    else
       rxf = rxf_ec
    endif
    call ids_copy(transport_solver_numerics_in, transport_solver_numerics_out)
    
    open(99, file='Core-Edge.coupling', position='append')

    volave=0.0_B2R8
    neave=0.0_B2R8
    teave=0.0_B2R8
    centroid_r = 0.0
    centroid_z = 0.0
    centroid_n = 0
    do ix = -1, nx
       iy=FindBottomRealCell(ny,ix)
       if(iy.eq.ny+1) cycle
       if(mod(region(ix,iy,0),4).eq.1) then
          centroid_r = centroid_r + cr(ix,iy)
          centroid_z = centroid_z + cz(ix,iy)
          centroid_n = centroid_n + 1
          volave=volave+vol(bottomix(ix,iy),bottomiy(ix,iy))
          neave=neave+ne(bottomix(ix,iy),bottomiy(ix,iy))* &
               & vol(bottomix(ix,iy),bottomiy(ix,iy))
          teave=teave+te(bottomix(ix,iy),bottomiy(ix,iy))* &
               &  ne(bottomix(ix,iy),bottomiy(ix,iy))* &
               & vol(bottomix(ix,iy),bottomiy(ix,iy))
       endif
    enddo
    centroid_r = centroid_r / centroid_n
    centroid_z = centroid_z / centroid_n
    write(*,*) 'centroid_r = ', centroid_r, 'centroid_z = ', centroid_z
    

    te_ave_iy = 0.0_R8
    ti_ave_iy = 0.0_R8
    ne_ave_iy = 0.0_R8
    na_ave_iy = 0.0_R8
    vol_ave_iy = 0.0_R8
    map_vol = 0.0_R8
    iy_core = -2

    do ix = -1, nx
       do iy = -1, ny
          if(mod(region(ix,iy,0),4).eq.1) then
             iy_core = max(iy_core, iy)
             tmp_vol_1 = triangle_volume(centroid_r, centroid_z, crx(ix,iy,0), cry(ix,iy,0), crx(ix,iy,1), cry(ix,iy,1))
             map_vol(iy,1) = map_vol(iy,1) + tmp_vol_1
             tmp_vol_2 = triangle_volume(centroid_r, centroid_z, crx(ix,iy,2), cry(ix,iy,2), crx(ix,iy,3), cry(ix,iy,3))
             map_vol(iy,2) = map_vol(iy,2) + tmp_vol_2
!!!             write(*,*) 'VOLUME-ELEMENTS', ix, iy, tmp_vol_1, tmp_vol_2                  
             te_ave_iy(iy) = te_ave_iy(iy) + te(ix,iy) * ne(ix,iy) * vol(ix,iy)
             ne_ave_iy(iy) = ne_ave_iy(iy) + ne(ix,iy) * vol(ix,iy)
             vol_ave_iy(iy) = vol_ave_iy(iy) + vol(ix,iy)
             do is0 = 0, ns-1
                ti_ave_iy(iy) = ti_ave_iy(iy) + ti(ix,iy) * na(ix,iy,is0) * vol(ix,iy)
                na_ave_iy(iy,is0) = na_ave_iy(iy,is0) + na(ix,iy,is0) * vol(ix,iy)
             enddo
          endif
       enddo
    enddo
    write(*,'(a3,1000(1x,a16))') 'iy', 'map_vol(iy,1)', 'map_vol(iy,2)', 'te_ave_iy(iy)', 'ti_ave_iy(iy)', &
         'ne_ave_iy(iy)', ('na_ave_iy(iy,is)', is0=0,ns-1)
    do iy = -1, iy_core
       te_ave_iy(iy) = te_ave_iy(iy) / ne_ave_iy(iy) / ev
       ne_ave_iy(iy) = ne_ave_iy(iy) / vol_ave_iy(iy)
       navol_ave = 0.0_R8
       do is0 = 0, ns-1
          navol_ave = navol_ave + na_ave_iy(iy,is0)
          na_ave_iy(iy,is0) = na_ave_iy(iy,is0) / vol_ave_iy(iy)
       enddo
       ti_ave_iy(iy) = ti_ave_iy(iy) / navol_ave / ev
       write(*,'(i3,1p,1000(1x,g16.6))') iy, map_vol(iy,1), map_vol(iy,2), te_ave_iy(iy), ti_ave_iy(iy), &
            ne_ave_iy(iy), (na_ave_iy(iy,is0), is0=0,ns-1)
! switch map_vol to normalised volume and sqrt(normalised volume)
       map_vol(iy,1) = 0.5_R8 * (map_vol(iy,1) + map_vol(iy,2)) / map_vol(iy_core,2)
       map_vol(iy,2) = sqrt(map_vol(iy,1))
    enddo

    volave=0.0_B2R8
    neave=0.0_B2R8
    teave=0.0_B2R8
    do ix = -1, nx
       iy=FindBottomRealCell(ny,ix)
       if(iy.eq.ny+1) cycle
       if(mod(region(ix,iy,0),4).eq.1) then
          volave=volave+vol(bottomix(ix,iy),bottomiy(ix,iy))
          neave=neave+ne(bottomix(ix,iy),bottomiy(ix,iy))* &
               & vol(bottomix(ix,iy),bottomiy(ix,iy))
          teave=teave+te(bottomix(ix,iy),bottomiy(ix,iy))* &
               &  ne(bottomix(ix,iy),bottomiy(ix,iy))* &
               & vol(bottomix(ix,iy),bottomiy(ix,iy))
       endif
    enddo
    teave=teave/neave/ev
    neave=neave/volave
    write(*,*) 'teave, neave = ', teave, neave
    bc_ec_te = (1.0_b2r8 - rxf) * bc_ec_te + rxf * teave
    write(*,*) b2etsmap(-1,4)
    nrho=size(transport_solver_numerics_out%solver_1d(1)%grid%volume)
    allocate(tsn_sqrt_norm_vol(nrho))
    tsn_sqrt_norm_vol = sqrt(transport_solver_numerics_out%solver_1d(1)%grid%volume / &
         transport_solver_numerics_out%solver_1d(1)%grid%volume(nrho))
    associate (eq => transport_solver_numerics_out%solver_1d(1)%equation(b2etsmap(-1,4)))
      write(*,*) trim(eq%primary_quantity%identifier%name(1))
      if (eq%computation_mode%index .ne. 2) then   ! predictive
         write(*,*) 'computation mode not predictive'
      else
         allocate(eq%boundary_condition(2)%type%name(1), eq%boundary_condition(2)%type%description(1))
         eq%boundary_condition(2)%type%name = 'value'
         eq%boundary_condition(2)%type%index = 1
         eq%boundary_condition(2)%type%description = 'Boundary condition is the value of the equations primary quantity'
         write(*,*) 'Old BC', eq%boundary_condition(2)%value(1)
         eq%boundary_condition(2)%value(1) = bc_ec_te
         write(*,*) 'New BC', eq%boundary_condition(2)%value(1)
         if ((eq%boundary_condition(2)%position .eq. bc_rho) .and. (bc_rho .lt. 1.0_R8)) then
            call interpolate_edge(nrho, irho_bc, iy_core, te_ave_iy, map_vol(:,2), &
                 eq%primary_quantity%profile, tsn_sqrt_norm_vol)
            write(*,*) 'DPC: debug interpolated edge', eq%boundary_condition(2)%value(1), te_ave_iy(-1), &
                 eq%primary_quantity%profile(irho_bc:nrho)
         end if
      endif
    end associate
    write(*,'(a,1p,1(1x,g15.5))') 'Core<-Edge: TE ', bc_ec_te
    write(99,'(a,1p,1(1x,g15.5))') 'Core<-Edge: TE ', bc_ec_te
! calculate the ion temperature b.c. == note B2 has a common ion temperature
    tiave = 0.0_B2R8
    naave = 0.0_B2R8
    do is0 = 0, ns-1
       do ix = -1, nx
          iy=FindBottomRealCell(ny,ix)
          if(iy.eq.ny+1) cycle
          if(mod(region(ix,iy,0),4).eq.1) then
             naave = naave + na(bottomix(ix,iy),bottomiy(ix,iy),is0)*  &
                  &       vol(bottomix(ix,iy),bottomiy(ix,iy))
             tiave = tiave + ti(bottomix(ix,iy),bottomiy(ix,iy))*  &
                  &       na(bottomix(ix,iy),bottomiy(ix,iy),is0)*  &
                  &       vol(bottomix(ix,iy),bottomiy(ix,iy))
          endif
       enddo
    end do
    tiave=tiave/naave/ev
    bc_ec_ti = (1.0_b2r8 - rxf) * bc_ec_ti + rxf * tiave
    write(*,'(a,1p,1(1x,g15.5))') 'Core<-Edge: TI ', bc_ec_ti
    write(99,'(a,1p,1(1x,g15.5))') 'Core<-Edge: TI ', bc_ec_ti
! loop over B2 species specifying core b.c. for matching species       
    do is0 = 0, ns-1
       if(b2etsmap(is0,0) .ge. 1) then
          naave = 0.0_B2R8
          do ix = -1, nx
             iy=FindBottomRealCell(ny,ix)
             if(iy.eq.ny+1) cycle
             if(mod(region(ix,iy,0),4).eq.1) then
                naave = naave+na(bottomix(ix,iy),bottomiy(ix,iy),is0)*  &
                     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
             endif
          enddo
          naave = naave/volave
          write(*,*) 'is0, naave = ', is0, naave
          bc_ec_na(is0) = (1.0_b2r8 - rxf) * bc_ec_na(is0) + rxf * naave
          ! ion density
          if (b2etsmap(is0,3) .gt. 0) then
             ieq = b2etsmap(is0,3)
             associate (eq => transport_solver_numerics_out%solver_1d(1)%equation(ieq))
               write(*,*) trim(eq%primary_quantity%identifier%name(1)), eq%primary_quantity%ion_index
               call xertst(eq%primary_quantity%ion_index .eq. b2etsmap(is0,1),   &
                    'ion_index in transport_solver_numerics does not match b2etsmap value')
               if (eq%computation_mode%index .ne. 2) then   ! predictive
                  write(*,*) 'computation mode not predictive'
               else
                  allocate(eq%boundary_condition(2)%type%name(1), eq%boundary_condition(2)%type%description(1))
                  eq%boundary_condition(2)%type%name = 'value'
                  eq%boundary_condition(2)%type%index = 1
                  eq%boundary_condition(2)%type%description = 'Boundary condition is the value of the equations primary quantity'
                  write(*,'(a,i3,1x,1p,g15.5)') 'Core<-Edge: IS, NI ', is0, bc_ec_na(is0)
                  write(99,'(a,i3,1x,1p,g15.5)') 'Core<-Edge: IS, NI ', is0, bc_ec_na(is0)
                  write(*,*) 'Old BC', eq%boundary_condition(2)%value(1)
                  eq%boundary_condition(2)%value(1) = bc_ec_na(is0)
                  write(*,*) 'New BC', eq%boundary_condition(2)%value(1)
                  if ((eq%boundary_condition(2)%position .eq. bc_rho) .and. (bc_rho .lt. 1.0_R8)) then
                     call interpolate_edge(nrho, irho_bc, iy_core, na_ave_iy(:,is0), map_vol(:,2), &
                          eq%primary_quantity%profile, tsn_sqrt_norm_vol)
                  end if
               endif
             end associate
          end if
          ! ion temperature
          if (b2etsmap(is0,4) .gt. 0) then
             associate (eq => transport_solver_numerics_out%solver_1d(1)%equation(b2etsmap(is0,4)))
               write(*,*) trim(eq%primary_quantity%identifier%name(1)), eq%primary_quantity%ion_index
               call xertst(eq%primary_quantity%ion_index .eq. b2etsmap(is0,1),   &
                    'ion_index in transport_solver_numerics does not match b2etsmap value')
               if (eq%computation_mode%index .ne. 2) then   ! predictive
                  write(*,*) 'computation mode not predictive'
               else
                  allocate(eq%boundary_condition(2)%type%name(1), eq%boundary_condition(2)%type%description(1))
                  eq%boundary_condition(2)%type%name = 'value'
                  eq%boundary_condition(2)%type%index = 1
                  eq%boundary_condition(2)%type%description = 'Boundary condition is the value of the equations primary quantity'
                  write(*,*) 'Old BC', eq%boundary_condition(2)%value(1)
                  eq%boundary_condition(2)%value(1) = bc_ec_ti
                  write(*,*) 'New BC', eq%boundary_condition(2)%value(1)
               endif
               if ((eq%boundary_condition(2)%position .eq. bc_rho) .and. (bc_rho .lt. 1.0_R8)) then
                  call interpolate_edge(nrho, irho_bc, iy_core, ti_ave_iy, map_vol(:,2), &
                       eq%primary_quantity%profile, tsn_sqrt_norm_vol)
               end if
             end associate
          end if
       else
          write(*,*) 'Skipped B2 species ', is0
       endif
    enddo
    deallocate(tsn_sqrt_norm_vol)

    write(99,'(a)') 'Core<-Edge: End'
    close(99)
    open(99, file='Core-Edge.dat')
    write(99,*) bc_ce_na, bc_ec_na, bc_ce_ti, bc_ec_ti, bc_ce_te, bc_ec_te
    close(99)
    firstpass=.false.

    call status_logger('Returning to the core code')
!    call chdir(trim(old_path))
!    write(*,*) 'Reset CWD to be ', trim(old_path)
    call prgend()

    return

  CONTAINS

    integer function find_ets_ion_in_b2etsmap(i0, b2etsmap, ns)
      implicit none
      integer b2etsmap(-1:ns-1, 0:4)
      integer i0, ns, is

      find_ets_ion_in_b2etsmap = -1
      do is = 0, ns-1
         if (b2etsmap(is, 1) .eq. i0) then
            find_ets_ion_in_b2etsmap = is
            return
         end if
      enddo
      return
    end function find_ets_ion_in_b2etsmap

    integer function find_ets_ion_and_state_in_b2etsmap(i0, i1, b2etsmap, ns)
      implicit none
      integer b2etsmap(-1:ns-1, 0:4)
      integer i0, i1, ns, is

      find_ets_ion_and_state_in_b2etsmap = -1
      do is = 0, ns-1
         if ((b2etsmap(is, 1) .eq. i0) .and. (b2etsmap(is, 2) .eq. i1)) then
            find_ets_ion_and_state_in_b2etsmap = is
            return
         end if
      enddo
      return
    end function find_ets_ion_and_state_in_b2etsmap

    subroutine interpolate_edge(nrho, irho_bc, iy_core, q_ave_iy, map_vol, &
         profile, sqrt_norm_vol)
      implicit none
      integer :: nrho, irho_bc, iy_core
      REAL (kind=B2R8), DIMENSION(-1:iy_core), INTENT(IN) :: q_ave_iy, map_vol
      REAL (kind=B2R8), DIMENSION(nrho), INTENT(IN) :: sqrt_norm_vol
      REAL (kind=B2R8), DIMENSION(nrho), INTENT(INOUT) :: profile

      call L3interp(q_ave_iy, map_vol, size(q_ave_iy), &
           profile(irho_bc+1:nrho), sqrt_norm_vol(irho_bc+1:nrho), nrho-irho_bc)

    end subroutine interpolate_edge

    SUBROUTINE L3interp(y_in,x_in,nr_in,y_out,x_out,nr_out)
!!$  Abramowitz & Stegun, Eqs 25.2.1 & 2 on p 878 for general form
!!$      equidistant is Eq 25.2.13 on p 879 for comparison and checking
!!$  4 pt general form coded by B Scott in 200x
!!$  added capability for decreasing abcissa in 2010
!!$  it works best if there is no extrapolation, ie, fixed end pts
!!$  a buffer region also works and is best for periodic cases
!!$      on which I usually use a double cycle even though it is overkill
!!$      because then you don't have to think about deformation or
!!$      resolution

      IMPLICIT NONE

      INTEGER :: nr_in,nr_out
      REAL (kind=B2R8), DIMENSION(nr_in), INTENT(IN) :: y_in,x_in
      REAL (kind=B2R8), DIMENSION(nr_out), INTENT(IN) :: x_out
      REAL (kind=B2R8), DIMENSION(nr_out), INTENT(OUT) :: y_out

      REAL (kind=B2R8) :: x,aintm,aint0,aint1,aint2,xm,x0,x1,x2
      INTEGER :: j,jm,j0,j1,j2
      INTEGER :: jstart,jfirst,jlast,jstep

      IF (x_in(nr_in) > x_in(1)) THEN
         jstart=3
         jfirst=1
         jlast=nr_out
         jstep=1
      ELSE
         jstart=nr_out-2
         jfirst=nr_out
         jlast=1
         jstep=-1
      END IF

      j1=jstart
      DO j=jfirst,jlast,jstep
         x=x_out(j)
         DO WHILE (x >= x_in(j1) .AND. j1 < nr_in-1 .AND. j1 > 2) 
            j1=j1+jstep
         END DO
         j2=j1+jstep
         j0=j1-jstep
         jm=j1-2*jstep

!...  extrapolate inside or outside

         x2=x_in(j2)
         x1=x_in(j1)
         x0=x_in(j0)
         xm=x_in(jm)

         aintm=(x-x0)*(x-x1)*(x-x2)/((xm-x0)*(xm-x1)*(xm-x2))
         aint0=(x-xm)*(x-x1)*(x-x2)/((x0-xm)*(x0-x1)*(x0-x2))
         aint1=(x-xm)*(x-x0)*(x-x2)/((x1-xm)*(x1-x0)*(x1-x2))
         aint2=(x-xm)*(x-x0)*(x-x1)/((x2-xm)*(x2-x0)*(x2-x1))

         y_out(j)=aintm*y_in(jm)+aint0*y_in(j0) &
              +aint1*y_in(j1)+aint2*y_in(j2)

      END DO

    END SUBROUTINE L3interp


    SUBROUTINE assign_code_parameters(codeparam_string, return_status)

      !-----------------------------------------------------------------------
      ! calls the XML parser for the code parameters and assign the
      ! resulting values to the corresponding variables
      !-----------------------------------------------------------------------

      USE xml2eg_mdl      ! IGNORE
      
      IMPLICIT NONE
      
      character(len=132), pointer :: codeparam_string(:)
      integer                     :: return_status
      type(type_xml2eg_document)  :: doc
      logical                     :: errorflag
      integer                     :: nion, nimp, iion, iimp, ndir

      call xml2eg_parse_memory  (codeparam_string, doc, warning_level=0)

      return_status = 0      ! no error

      call xml2eg_get (doc,  'coupling/niter', niter, errorflag)
      call xml2eg_get (doc,  'coupling/rxf_ce', rxf_ce, errorflag)
      call xml2eg_get (doc,  'coupling/rxf_ec', rxf_ec, errorflag)
!      call xml2eg_get (doc,  'coupling/direction', direction, errorflag)  !! multiple
      call xml2eg_get (doc,  'solps/directory', solps_directory, errorflag)

      call xml2eg_free_doc(doc)
      
      return

    end SUBROUTINE assign_code_parameters

  end subroutine b2mn_ets
  
  subroutine b2mn_ets_finalize (outputFlag, diagnosticInfo)
    use ids_types               ! IGNORE
    use b2mod_main
    implicit none
    integer (ids_int),        intent(out)  :: outputFlag
    character(len=:), pointer, intent(out) :: diagnosticInfo

    outputFlag = 0
    nullify(diagnosticInfo)

    write(*,*) 'b2mn_ets_finalize'
!    call chdir(trim(new_path))
!    write(*,*) 'Set CWD to be ', trim(new_path)
    call status_logger('Started the SOLPS finalization')
    call b2mn_fin
    deallocate(te_ave_iy, ti_ave_iy, ne_ave_iy, na_ave_iy, vol_ave_iy, map_vol)
    deallocate(direction)
    deallocate(b2etsmap)
    deallocate(bc_ce_na, bc_ec_na, bc_ce_ti)
    call status_logger('Completed the SOLPS finalization, returning to the core code')
!    call chdir(trim(old_path))
!    write(*,*) 'Reset CWD to be ', trim(old_path)
  end subroutine b2mn_ets_finalize

  subroutine status_logger(message)

    character (len=*) :: message
    real cputime
    character*8 date
    character*10 time
    character*5 zone
    logical, save :: first
    data first /.true./

    if(first) then
       open(99, file='Core-Edge.status')
       first=.false.
    else
       open(99, file='Core-Edge.status', position='append')
    endif
    call cpu_time(cputime)
    call date_and_time(date,time,zone)
    write(99,'(f10.3,4(1x,a))')  cputime, date(1:4)//'-'//date(5:6)//'-'//date(7:8), time(1:2)//':'//time(3:4)//':'//time(5:10), &
         & zone, message
    close(99)

  end subroutine status_logger

  function triangle_volume(x1, y1, x2, y2, x3, y3)
    use b2mod_constants
    implicit none
    real (kind=B2R8) :: x1, y1, x2, y2, x3, y3, triangle_volume
    triangle_volume =  abs(0.5_R8 * (x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))) * &
         (x1 + x2 + x3) / 3.0_R8 * 2.0_R8 * pi
  end function triangle_volume

end module core_edge
