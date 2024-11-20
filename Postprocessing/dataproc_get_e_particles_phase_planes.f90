
!--------------------------------------
!
MODULE ParticleTracing

!  INTEGER, PARAMETER :: maximal_N_cells_tracing = 50
!  INTEGER, PARAMETER :: maximal_N_part_to_trace_in_cell = 20

  LOGICAL save_e_vz
  LOGICAL save_e_irrEx
  LOGICAL save_e_irrEy
  LOGICAL save_e_extEz
  LOGICAL save_e_solEx
  LOGICAL save_e_solEy
  LOGICAL save_e_solEz
  LOGICAL save_e_solBx
  LOGICAL save_e_solBy
  LOGICAL save_e_solBz
  LOGICAL save_e_extBx
  LOGICAL save_e_extBy
  LOGICAL save_e_extBz

!  INTEGER Tcntr_e_tracing_start
!  INTEGER Tcntr_e_tracing_end
!  INTEGER Tcntr_e_tracing_step
!  INTEGER T_cntr_save_traced_electrons

  INTEGER e_record_length

!  INTEGER N_of_cells_e_tracing

!  TYPE ParticleTracingSetup
!     INTEGER icell
!     INTEGER jcell
!     INTEGER N_to_trace
!     LOGICAL choose_negative_VX
!     LOGICAL choose_positive_VX
!     LOGICAL choose_negative_VY
!     LOGICAL choose_positive_VY
!  END TYPE ParticleTracingSetup

!  TYPE(ParticleTracingSetup) cell_e_tracing(1:maximal_N_cells_tracing)

END MODULE ParticleTracing

!------------------------------------
!
program get_traced_epp_at_snapshot_times

  use ParticleTracing

  implicit none

  logical exists

  character(1) buf
  integer ios
  integer itmp

  integer N_of_all_snaps
  integer, allocatable :: Tcntr_snapshot(:)
  integer alloc_err
  real rtmp

  logical first_record

  integer Tcntr

  integer current_snapshot
  integer i

  real t_ns
  integer N_particles_this_step

  integer k

  character(20) epp_filename   ! _NNNN_traced_epp.dat
!                                ----x----I----x----I

  integer tag
  real rvector(2:18)
  integer pos

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE
 
  inquire( file = 'electron_particles_vst.dat', exist = exists )
  if (.not.exists) then
     print '("file electron_particles_vst.dat not found, program terminated")'
     stop
  end if

  inquire( file = '_snapmoments.dat', exist = exists )
  if (.not.exists) then
     print '("file _snapmoments.dat not found, program terminated")'
     stop
  end if

! find number of snaphot moments
  open (41, file = '_snapmoments.dat')
!                 "--****-----*******.*****----********----*----*----*----*----*"
!  WRITE (41, '(" number       time(ns)       T_cntr    vdf  pp  ioniz icbo ecbo ")')
  read (41, '(A1)') buf
!     DO i = 1, N_of_all_snaps
  do
     read (41, '(2x,i5,5x,f13.5,4x,i8,4x,i1,4x,i1,4x,L1,4x,L1,4x,L1)', iostat = ios) itmp
     if (ios.ne.0) exit
     N_of_all_snaps = itmp
!        WRITE (41, '(2x,i4,5x,f13.5,4x,i8,4x,i1,4x,i1,4x,L1,4x,L1,4x,L1)') &
!             & i, &
!             & Tcntr_snapshot(i) * 1.0d9 * delta_t_s, &
!             & Tcntr_snapshot(i), &
!             & save_evdf_snapshot(i), &
!             & save_pp_snapshot(i), &
!             & save_ionization_rates_2d(i), &
!             & save_ions_collided_with_bo(i), &
!             & save_e_collided_with_bo(i)
  end do
  close (41, status = 'KEEP')

  allocate(Tcntr_snapshot(N_of_all_snaps), stat = alloc_err)

! read snapshot moments
  open (41, file = '_snapmoments.dat')
  read (41, '(A1)') buf
  do i = 1, N_of_all_snaps
     read (41, '(2x,i5,5x,f13.5,4x,i8,4x,i1,4x,i1,4x,L1,4x,L1,4x,L1)') itmp, rtmp, Tcntr_snapshot(i)
  end do
  close (41, status = 'KEEP')

  open (9, file = 'electron_particles_vst.dat', access = 'stream')

  read (9) save_e_vz
  read (9) save_e_irrEx
  read (9) save_e_irrEy
  read (9) save_e_extEz
  read (9) save_e_solEx
  read (9) save_e_solEy
  read (9) save_e_solEz
  read (9) save_e_solBx
  read (9) save_e_solBy
  read (9) save_e_solBz
  read (9) save_e_extBx
  read (9) save_e_extBy
  read (9) save_e_extBz

  e_record_length = 5
  if (save_e_vz) e_record_length = e_record_length+1
  if (save_e_irrEx) e_record_length = e_record_length+1
  if (save_e_irrEy) e_record_length = e_record_length+1
  if (save_e_extEz) e_record_length = e_record_length+1
  if (save_e_solEx) e_record_length = e_record_length+1
  if (save_e_solEy) e_record_length = e_record_length+1
  if (save_e_solEz) e_record_length = e_record_length+1
  if (save_e_solBx) e_record_length = e_record_length+1
  if (save_e_solBy) e_record_length = e_record_length+1
  if (save_e_solBz) e_record_length = e_record_length+1
  if (save_e_extBx) e_record_length = e_record_length+1
  if (save_e_extBy) e_record_length = e_record_length+1
  if (save_e_extBz) e_record_length = e_record_length+1

  first_record = .true.

  do

     read (9, iostat = ios) Tcntr
     if (ios.ne.0) exit

     if (first_record) then
! find the first snapshot which is at or after the time of the first record
        current_snapshot = -1
        do i = 1, N_of_all_snaps
           if (Tcntr_snapshot(i).ge.Tcntr) then
! assume the simplest case when particles were save at each step, 
! do not worry that there may be no snapshot at the time of the first record
              current_snapshot = i
              exit
           end if
        end do
     end if

     first_record = .false.

     read (9) t_ns
     read (9) N_particles_this_step

     if (Tcntr_snapshot(current_snapshot).ne.Tcntr) then

        do k = 1, N_particles_this_step
           read (9) itmp
           do i = 2, e_record_length
              read (9) rtmp
           end do
        end do

     else   !###   if (Tcntr_snapshot(current_snapshot).ne.Tcntr) then

! create filename
        epp_filename = '_NNNN_traced_epp.dat'
!                              ----x----I----x----I
        epp_filename(2:5) = convert_int_to_txt_string(current_snapshot, 4)
        open (10, file = epp_filename)

        do k = 1, N_particles_this_step
           read (9) tag
           read (9) rvector(2)  ! x_m
           read (9) rvector(3)  ! y_m
           read (9) rvector(4)  ! vx_ms
           read (9) rvector(5)  ! vy_ms
           pos=6
           if (save_e_vz) then
              read (9) rvector(pos)  ! vz_ms
              pos = pos+1
           end if
           if (save_e_irrEx) then
              read (9) rvector(pos)  ! irrEx_Vm
              pos = pos+1
           end if
           if (save_e_irrEy) then
              read (9) rvector(pos)  ! irrEy_Vm
              pos = pos+1
           end if
           if (save_e_extEz) then
              read (9) rvector(pos)  ! extEz_Vm
              pos = pos+1
           end if
           if (save_e_solEx) then
              read (9) rvector(pos)  ! solEx_Vm
              pos = pos+1
           end if
           if (save_e_solEy) then
              read (9) rvector(pos)  ! solEy_Vm
              pos = pos+1
           end if
           if (save_e_solEz) then
              read (9) rvector(pos)  ! solEz_Vm
              pos = pos+1
           end if
           if (save_e_solBx) then
              read (9) rvector(pos)  ! solBx_T
              pos = pos+1
           end if
           if (save_e_solBy) then
              read (9) rvector(pos)  ! solBy_T
              pos = pos+1
           end if
           if (save_e_solBz) then
              read (9) rvector(pos)  ! solBz_T
              pos = pos+1
           end if
           if (save_e_extBx) then
              read (9) rvector(pos)  ! extBx_T
              pos = pos+1
           end if
           if (save_e_extBy) then
              read (9) rvector(pos)  ! extBy_T
              pos = pos+1
           end if
           if (save_e_extBz) then
              read (9) rvector(pos)  ! extBz_T
              pos = pos+1
           end if

           write (10, '(2x,i9,17(2x,e16.9))') &
             & tag, &
             & rvector(2:e_record_length)

        end do   !###   do k = 1, N_particles_this_step

        close (10, status = 'keep')
        print '("created file ",A20)', epp_filename

        current_snapshot = current_snapshot+1
        if (current_snapshot.gt.N_of_all_snaps) exit

     end if   !###   if (Tcntr_snapshot(current_snapshot).ne.Tcntr) then
         
  end do   !###   do

  close (9, status = 'keep')
  print '("finished")'

end program get_traced_epp_at_snapshot_times

!-----------------------------------
! creates a string of length "length_of_string" out of an integer number "int_number"
!
function convert_int_to_txt_string(int_number, length_of_string)

  implicit none

  integer int_number

  integer length_of_string
  character*(length_of_string) convert_int_to_txt_string

  character(5) format_string
  character(2) length2_txt
  character(1) length1_txt

  character*(length_of_string) number_txt
  
  integer blanks_number
  integer i

! produce format string
  if ((length_of_string.gt.0).and.(length_of_string.lt.10)) then
     write (length1_txt, '(i1)') length_of_string
     format_string = '(iN) '
     format_string(3:3) = length1_txt
  else if ((length_of_string.ge.10).and.(length_of_string.lt.100)) then
     write (length2_txt, '(i2)') length_of_string
     format_string = '(iNN)'
     format_string(3:4) = length2_txt
  else
     print *, "ERROR in CONVERT_INT_TO_TXT_STRING:"
     print *, "incorrect string length requested: ", length_of_string
     stop
  end if

  WRITE (number_txt, format_string) int_number
  number_txt = ADJUSTL(TRIM(number_txt))
  blanks_number = length_of_string - LEN_TRIM(number_txt)
  number_txt = ADJUSTR(number_txt)
  do i = 1, blanks_number
     number_txt(i:i) = '0'
  end do

  convert_int_to_txt_string = number_txt

end function convert_int_to_txt_string
