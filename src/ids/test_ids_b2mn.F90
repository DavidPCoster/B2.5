!> Wrapper for testing ids_b2mn
!!
!! \brief The input and output IDS's are specified via the command line
!!
!! To use
!!
!!    compile the code with, e.g., 
!!
!!      <tt>make obj/$ITM_INTEL_OBJECTCODE/test_source_fusion.exe</tt>
!!
!!    run the code with, e.g., 
!!
!!      <tt>obj/$ITM_INTEL_OBJECTCODE/test_source_fusion.exe -u "imas:hdf5?path=/afs/eufus.eu/user/g/g2dpc/public/imasdb/iter_g2tjohns/3/134173/17716" -t 210  </tt>
!!
!! Supported options
!!   \li <b>-u</b> URI
!!   \li <b>-t, -T</b> time (0.0)
!!   \li <b>-i</b> iterations (1)
!!
!! \author David Coster (David.Coster@ipp.mpg.de)
!!
program test_ids_b2mn

  use ids_schemas                           ! IGNORE
  use ids_routines                          ! IGNORE
  use xml_file_reader, only: fill_param     ! IGNORE
  use core_edge

  implicit none

  type (ids_core_profiles)             :: core_profiles_in
  type (ids_core_transport)            :: core_transport_in
  type (ids_equilibrium)               :: equilibrium_in
  type (ids_transport_solver_numerics) :: transport_solver_numerics_in, transport_solver_numerics_out
  type (ids_parameters_input)          :: codeparam
  integer (ids_int)                    :: user_out_outputFlag
  character(len=:), pointer            :: user_out_diagnosticInfo

  real(ids_real)                       :: time=0.0_ids_real
  integer                              :: count, icount, idx, ids_status
  character(len=256)                   :: arg, uri
  character(:), allocatable            :: retmsg
  integer                              :: iterations=1, i, iter

  count = command_argument_count()

  icount = 1
  uri = ''

  do while (icount .le. count)
     call get_command_argument(icount, arg)
     select case(arg)
     case("-u")
        if(icount .eq. count) then
           stop 'error specifying -u'
        else
           call get_command_argument(icount+1, uri)
        end if
     case("-T", "-t")
        if(icount .eq. count) then
           stop 'error specifying -T'
        else
           call get_command_argument(icount+1, arg)
           read(arg,*) time
        end if
     case("-i")
        if(icount .eq. count) then
           stop 'error specifying -i'
        else
           call get_command_argument(icount+1, arg)
           read(arg,*) iterations
        end if
     case default
        write(*,*) 'Option ', trim(arg), ' not recognized'
        stop 'error sepcifying option'
     end select
     icount = icount + 2
  end do

  if (uri .eq. "") then
     write(*,*) 'ERROR: no URI specified'
     stop
  endif
  write(*,"(a,a,' @ ',g13.6)") "Input case: ", trim(uri), time

  call imas_open(uri, OPEN_PULSE, idx, ids_status, retmsg)
  if (ids_status .ne. 0) then
     write(*,*) 'ERROR received from imas_open, STATUS=', ids_status
     write(*,'(a)') trim(retmsg)
     stop
  else
     write(*,*) 'imas_open OK'
  end if

  call ids_get_slice (idx, 'core_profiles', core_profiles_in, time, closest_interp, ids_status)
  if (ids_status .ne. 0) then
     write(*,*)'ERROR received from ids_get_slices reading core_profiles, STATUS=', ids_status
     stop
  end if

  call ids_get_slice (idx, 'core_transport', core_transport_in, time, closest_interp, ids_status)
  if (ids_status .ne. 0) then
     write(*,*)'ERROR received from ids_get_slices reading core_transport, STATUS=', ids_status
     stop
  end if

  call ids_get_slice (idx, 'equilibrium', equilibrium_in, time, closest_interp, ids_status)
  if (ids_status .ne. 0) then
     write(*,*)'ERROR received from ids_get_slices reading equilibrium, STATUS=', ids_status
     stop
  end if

  call ids_get_slice (idx, 'transport_solver_numerics', transport_solver_numerics_in, time, closest_interp, ids_status)
  if (ids_status .ne. 0) then
     write(*,*)'ERROR received from ids_get_slices reading transport_solver_numerics, STATUS=', ids_status
     stop
  end if
  call imas_close (idx)

  ! get XML code parameters
  CALL FILL_PARAM (codeparam%parameters_value, codeparam%schema,          codeparam%parameters_default, &
       'xml/core-edge.xml',  'xsd/core-edge.xsd', 'xml/core-edge.xml')

  do iter = 1, iterations
     call b2mn_ets(core_profiles_in, core_transport_in, equilibrium_in,  &
     transport_solver_numerics_in, transport_solver_numerics_out,  &
     codeparam, user_out_outputFlag, user_out_diagnosticInfo)
     write(*,*) 'Returned from b2mn_ets'
     if(associated(user_out_diagnosticInfo)) then
        if(user_out_outputFlag.ne.0) then
           write(*,*) 'Error reported: ', user_out_outputFlag, trim(user_out_diagnosticInfo)
        endif
        deallocate(user_out_diagnosticInfo)
     else
        if(user_out_outputFlag.ne.0) then
           write(*,*) 'Error reported: ', user_out_outputFlag
        endif
     endif

     write(*,*)
     write(*,*) 'Time = ', transport_solver_numerics_out%time
     write(*,*)

! write as HDF5
     call imas_open ('imas:hdf5?path=./testdb', FORCE_CREATE_PULSE, idx, ids_status)
     write(*,*) 'IDX = ', idx
     write(*,*) 'ids_STATUS = ', ids_status
     if(ids_status.ne.0) then
        write(0,*) 'Failure opening IDS!'
        stop 1
     endif
     write(*,*) 'Calling ids_put'
     call ids_put(idx, 'transport_solver_numerics', transport_solver_numerics_out, ids_status)
     write(*,*) 'ids_STATUS = ', ids_status
     if(ids_status.ne.0) then
        write(0,*) 'Failure writing IDS!'
        stop 1
     endif
     write(*,*) 'ids_put done'
     call imas_close  (idx)
     call ids_deallocate(transport_solver_numerics_out)
     write(*,*)
  enddo

  call b2mn_ets_finalize(user_out_outputFlag, user_out_diagnosticInfo)
  write(*,*) 'Returned from b2mn_ets (finalize)'
  if(associated(user_out_diagnosticInfo)) then
     if(user_out_outputFlag.ne.0) then
        write(*,*) 'Error reported: ', user_out_outputFlag, trim(user_out_diagnosticInfo)
     endif
     deallocate(user_out_diagnosticInfo)
  else
     if(user_out_outputFlag.ne.0) then
        write(*,*) 'Error reported: ', user_out_outputFlag
     endif
  endif
  call ids_deallocate(core_profiles_in)
  call ids_deallocate(core_transport_in)
  call ids_deallocate(equilibrium_in)
  call ids_deallocate(transport_solver_numerics_in)
  deallocate (codeparam%parameters_value, codeparam%schema, codeparam%parameters_default)
!!! NYI  call ets_b2mn_finalize()

end program test_ids_b2mn
