
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
program get_separate_e_history_vst

  use ParticleTracing

  implicit none

  logical exists
  logical first_record
  integer step_counter

  integer ios

  integer Tcntr
  real t_ns
  integer N_particles_this_step

  integer k

  integer tag
  real rvector(2:18)
  integer pos

  character(21) particle_filename   ! traced_e_NNNN_vst.dat
!                                     ----x----I----x----I-

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
  step_counter = 0

  do

     read (9, iostat = ios) Tcntr
     if (ios.ne.0) exit
     step_counter = step_counter+1
     read (9) t_ns
     read (9) N_particles_this_step

     print '("step ",i9," N_particles ",i4)', Tcntr, N_particles_this_step

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
! create filename
        particle_filename = 'traced_e_NNNN_vst.dat'
!                            ----x----I----x----I-
        particle_filename(10:13) = convert_int_to_txt_string(tag/10000, 4)

        if (first_record) then
           open (10, file = particle_filename, status = 'replace')
        else
           inquire (file = particle_filename, exist = exists)
           if (.not.exists) then
              print '("something wrong, a file for particle with tag ",i9," should have been created already")'
              stop
           end if
           open (10, file = particle_filename, position = 'append')
        end if

        write (10, '(2x,i9,2x,f13.3,17(2x,e16.9))') &
             & Tcntr, &
             & t_ns, &
             & rvector(2:e_record_length)

        close (10, status = 'keep')
        
     end do   !###   do k = 1, N_particles_this_step
     first_record = .false.
  end do   !###   do

  close (9, status = 'keep')
  print '("finished, total number of steps processed ",i9)', step_counter

end program get_separate_e_history_vst

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
