
!--------------------

SUBROUTINE INITIATE_AVERAGED_SNAPSHOTS

  USE ParallelOperationValues
  USE AvgSnapshots
  USE CurrentProblemValues, ONLY : delta_t_s, N_subcycles, Start_T_cntr, string_length, delta_x_m
  USE Checkpoints, ONLY : use_checkpoint
  USE MCCollisions, ONLY : N_neutral_spec, collision_e_neutral, en_collisions_turned_off
  USE mod_print, ONLY: print_parser_error
  USE IonParticles, ONLY: N_spec
  USE mod_print, ONLY: print_message

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  LOGICAL exists
  CHARACTER (1) buf

  INTEGER ios
  INTEGER saveflagi(1:39)           ! integer flags used to set values of logical flags controlling saving of data files

  INTEGER N_of_snap_groups          ! number of sets of snapshots, read from file
  INTEGER i, n, p, count

  REAL(8) Rqst_snap_start_ns        ! requested start of current set of snapshots [ns], read from file
  REAL(8) Rqst_snap_finish_ns       ! requested finish of current set of snapshots [ns], read from file 
  INTEGER Rqst_n_of_snaps           ! requested number of snapshots in current set, read from file
  
! temporary arrays
  INTEGER, ALLOCATABLE :: timestep_begin(:)   ! moments when average data collection for a snapshot begins
  INTEGER, ALLOCATABLE :: timestep_end(:)     ! moments when average data collection for a snapshot ends and the snapshot is saved
  INTEGER ALLOC_ERR

  INTEGER T1, T2, T2prev
  INTEGER large_step
  INTEGER time_begin, time_end
  CHARACTER(LEN=string_length) :: message, plane_name, file_name
  INTEGER :: idx, species

  INTEGER :: i_dummy

! default values ensure that if init_avgsnapshots.dat is not found 
! procedures COLLECT_DATA_FOR_AVERAGED_SNAPSHOT and CREATE_AVERAGED_SNAPSHOT do nothing
  N_of_all_avgsnaps = 0
  current_avgsnap = 1
  save_avg_data = .FALSE.
  time_begin = 0
  time_end = 0
!  avg_data_collection_offset = -1
  
  num_plane_x_locations = 0
  num_plane_y_locations = 0
  save_collision_freq_ee = .FALSE.
  avg_flux_and_history = .FALSE.
! read / write the data file 
  INQUIRE (FILE = 'init_avgsnapshots.dat', EXIST = exists)

  IF (.NOT.exists) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### File init_avgsnapshots.dat not found. Time-averaged snapshots will not be created. ###")'
   !   IF (num_plane_x_locations>0 .OR. num_plane_y_locations>0) THEN
   !    WRITE( message, '(A,I10,A,I10,A)'), "You required either cut planes in X (NX) or Y (NY) directions to measure the flux through it.& 
   !                                         Got NX = ",num_plane_x_locations," and NY = ",num_plane_y_locations,". &
   !                                         You must have an init_avgsnapshots.dat file to define the time average period for the flux."//achar(10)
   !    CALL print_parser_error(message)      
   !   ENDIF 
     RETURN
  END IF

  IF (Rank_of_process.EQ.0) PRINT '("### File init_avgsnapshots.dat is found. Reading the data file... ###")'

  OPEN (9, FILE = 'init_avgsnapshots.dat')

  saveflagi = 0

  READ(9, '(A1)') buf  !--- save 2D maps of the following TIME-AVERAGED parameters? (1=yes, 0=no)
  READ(9, '(A1)') buf  !-----F----EX----EY--JXsum--JYsum--JZsum--e-n collision frequencies   (type flag values below)
  READ(9, *, iostat = ios) saveflagi(1:6) , saveflagi(39)

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat first set of flags is probably incomplete, missing flag(s) set to zero ###")', ios
     BACKSPACE(9)
  END IF

  READ(9, '(A1)') buf  !----Ne----JXe---JYe---JZe---VXe---VYe---VZe---WXe---WYe---WZe---TXe---TYe---TZe---QXe---QYe---QZe  (type flag values below)
  READ(9, *, iostat = ios) saveflagi(7:22)

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat second set of flags is probably incomplete, missing flag(s) set to zero ###")', ios
     BACKSPACE(9)
  END IF

  READ(9, '(A1)') buf  !----Ni----JXi---JYi---JZi---VXi---VYi---VZi---WXi---WYi---WZi---TXi---TYi---TZi---QXi---QYi---QZi  (type flag values below)
  READ (9, *, iostat = ios) saveflagi(23:38)

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat third set of flags is probably incomplete, missing flag(s) set to zero ###")', ios
     BACKSPACE(9)
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: saveflagi(1:6,39)     = ",7(2x,i2))', saveflagi(1:6) , saveflagi(39)
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: saveflagi(7:22)       = ",16(2x,i2))', saveflagi(7:22)
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: saveflagi(23:38)      = ",16(2x,i2))', saveflagi(7:22)
  END IF

  DO i = 1, 39
     IF (saveflagi(i).GT.0) save_avg_data(i) = .TRUE.
  END DO

  IF (en_collisions_turned_off) save_avg_data(39) = .FALSE.

! if electron-neutral collision frequencies are requested, check that at least one active collisional process has save_collfreq_2d == .TRUE.
  IF (save_avg_data(39)) THEN
     count = 0
     DO n = 1, N_neutral_spec
        DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
           IF (collision_e_neutral(n)%colproc_info(p)%save_collfreq_2d) count = count + 1
        END DO
     END DO
     IF (count.EQ.0) THEN
        save_avg_data(39) = .FALSE.
        IF (Rank_of_process.EQ.0) PRINT '("### WARNING :: Saving e-n collision frequencies was requested in init_avgsnapshots.dat but not set in init_neutral_AAAAAA.dat so it is turned off ###")'
     END IF
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: save_avg_data(1:6,39) = ",7(3x,L1))', save_avg_data(1:6) , save_avg_data(39)
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: save_avg_data(7:22)   = ",16(3x,L1))', save_avg_data(7:22)
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: save_avg_data(23:38)  = ",16(3x,L1))', save_avg_data(7:22)
  END IF

  READ (9, '(A1)') buf !--- number of groups of snapshots (>0, no averaged snapshots if <=0), type below
  READ (9, *, iostat=ios) N_of_snap_groups

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat cannot read the number of snapshot groups, set it to zero ###")', ios
     N_of_snap_groups = 0
  END IF

  IF (N_of_snap_groups.LE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### ### Time-averaged snapshots will not be created. ### ###")'
     CLOSE (9, STATUS = 'KEEP')
     RETURN
  END IF

  READ (9, '(A1)') buf !--- For each group, type below start (ns), finish (ns), number of snapshots (>=0)

  T2prev = -1

  ALLOCATE(timestep_begin(1:99999), STAT = ALLOC_ERR)
  ALLOCATE(timestep_end(1:99999), STAT = ALLOC_ERR)
     
  DO i = 1, N_of_snap_groups
! read the parameters of current set of snapshot from the data file
     READ (9, *, iostat = ios) Rqst_snap_start_ns, Rqst_snap_finish_ns, Rqst_n_of_snaps

     IF (ios.NE.0) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat while reading snapshot group ",i3,", skip ###")', ios, i
        CYCLE
     END IF

! try the next group of snapshots if the current group snapshot number is zero
     IF (Rqst_n_of_snaps.LT.1) CYCLE
     IF (Rqst_snap_start_ns.GE.Rqst_snap_finish_ns) CYCLE

! timestep when averaging for the first snapshot in the group begins 
     T1 = Rqst_snap_start_ns / (delta_t_s * 1.0d9)

! timestep when averaging for the last snapshot in the group ends
     T2 = Rqst_snap_finish_ns / (delta_t_s * 1.0d9)

     T1 = (T1 / N_subcycles) * N_subcycles     ! + avg_data_collection_offset
     T2 = (T2 / N_subcycles) * N_subcycles     ! + avg_data_collection_offset

     T1 = MAX(T1, T2prev)

     large_step = (T2 - T1) / Rqst_n_of_snaps
     large_step = (large_step / N_subcycles) * N_subcycles
     IF (large_step.LE.N_subcycles) THEN
      WRITE( message,'(A,I0,A,I0,A)') "Time step for averaged snapshots is less then the ion cycle. Reduce number of snapshots. No average snapshots will be generated: large_step =  ",large_step," vs N_subcycles",N_subcycles,achar(10)
      CALL print_message(message)      
      CYCLE
     END IF

! for all possible snapshots of the current set
     DO n = 1, Rqst_n_of_snaps  !Fact_n_of_snaps
! Calculate and save the snapshot moment in the temporary array
        time_begin = T1 + (n - 1) * large_step
        time_end = time_begin + large_step - 1 !N_subcycles
        IF (time_begin.GE.T2) EXIT
        IF (time_end.GT.T2) EXIT  !time_end = T2

        N_of_all_avgsnaps = N_of_all_avgsnaps + 1
        timestep_begin(N_of_all_avgsnaps) =  time_begin
        timestep_end(N_of_all_avgsnaps) = time_end
     END DO        ! end of cycle over snapshots in one set

     IF (N_of_all_avgsnaps.LE.0) CYCLE
     T2prev = timestep_end(N_of_all_avgsnaps)

  END DO           ! end of cycle over sets of snapshots     

  CLOSE (9, STATUS = 'KEEP')

  IF (N_of_all_avgsnaps.EQ.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### Time-averaged snapshots will NOT be created ... ###")'
     WRITE( message,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') "Debugging Large step = ",large_step," T1 = ",T1," and T2 = ",T2," time_begin = ",time_begin," time_end = ",time_end,achar(10)
     CALL print_message(message)
! cleanup
     DEALLOCATE(timestep_begin, STAT = ALLOC_ERR)
     DEALLOCATE(timestep_end, STAT = ALLOC_ERR)
     RETURN
  END IF 

! if we are here, snapshots will be created ...

! allocate the array of moments of snapshots
  ALLOCATE(avgsnapshot(1:N_of_all_avgsnaps), STAT=ALLOC_ERR)

! move the calculated snapshot moments from the temporary array to the allocated array 
  DO n = 1, N_of_all_avgsnaps
     avgsnapshot(n)%T_cntr_begin = timestep_begin(n)
     avgsnapshot(n)%T_cntr_end = timestep_end(n)
  END DO

  DEALLOCATE(timestep_begin, STAT = ALLOC_ERR)
  DEALLOCATE(timestep_end, STAT = ALLOC_ERR)

  IF (Rank_of_process.EQ.0) THEN 
     PRINT '("### The program will create ",i4," averaged snapshots ###")', N_of_all_avgsnaps

! write moments of snapshot creation into the file
     OPEN (41, FILE = '_avg_snapmoments.dat', STATUS = 'REPLACE')
     WRITE (41, '(" number   start_time(ns)   end_time(ns)   start_T_cntr   end_T_cntr   N_of_avg_points,e/i")')
     DO n = 1, N_of_all_avgsnaps
        WRITE (41, '(2x,i5,2x,2(2x,f13.5),2x,2(2x,i9),2(4x,i6))') &
             & n, &
             & avgsnapshot(n)%T_cntr_begin * 1.0d9 * delta_t_s, &
             & avgsnapshot(n)%T_cntr_end * 1.0d9 * delta_t_s, &
             & avgsnapshot(n)%T_cntr_begin, &
             & avgsnapshot(n)%T_cntr_end, &
             & avgsnapshot(n)%T_cntr_end - avgsnapshot(n)%T_cntr_begin + 1, &
             & (avgsnapshot(n)%T_cntr_end - avgsnapshot(n)%T_cntr_begin + 1) / N_subcycles
     END DO
     CLOSE (41, STATUS = 'KEEP')
  END IF

! consistency check 1
  DO n = 1, N_of_all_avgsnaps
     IF (avgsnapshot(n)%T_cntr_begin.GE.avgsnapshot(n)%T_cntr_end) THEN
        IF (Rank_of_process.EQ.0) PRINT '("ERROR-1 in INITIATE_AVERAGED_SNAPSHOTS, snapshot "i4," begins at ",i9," and ends at ",i9)', &
             & n, avgsnapshot(n)%T_cntr_begin, avgsnapshot(n)%T_cntr_end
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF
  END DO

  IF (Rank_of_process.EQ.0) PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: average snapshot timing passed consistency check 1")'

! consistency check 2
  DO n = 1, N_of_all_avgsnaps-1
     IF (avgsnapshot(n)%T_cntr_end.GE.avgsnapshot(n+1)%T_cntr_begin) THEN
        IF (Rank_of_process.EQ.0) PRINT '("ERROR-2 in INITIATE_AVERAGED_SNAPSHOTS, snapshot ",i4,2x,i9,2x,i9," overlaps with snapshot "i4,2x,i9,2x,i9)', &
             & n,   avgsnapshot(n)%T_cntr_begin,   avgsnapshot(n)%T_cntr_end, &
             & n+1, avgsnapshot(n+1)%T_cntr_begin, avgsnapshot(n+1)%T_cntr_end
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF
  END DO

  IF (Rank_of_process.EQ.0) PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: average snapshot timing passed consistency check 2")'

! note, since we are here, N_of_all_avgsnaps is not zero
! overwrite default value of the first snapshot number if the system is initialized using a checkpoint
  IF (use_checkpoint.EQ.1) THEN 
      ! Check if we should due an avg snapshot given the current starting time
   IF ( avgsnapshot(N_of_all_avgsnaps)%T_cntr_begin < Start_T_cntr ) THEN ! If we enter here it means we already have all our snapshots
      current_avgsnap = N_of_all_avgsnaps +1 ! We artifically set this index to prevent from generating an average snapshot  
      IF (Rank_of_process.EQ.0) PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: the simulation time is greater than the starting time for averaged snapshots. &
                                       All snapshots have been generated already")'
   ELSE ! It means we still we have some snapshots to do. 
      DO n = 1, N_of_all_avgsnaps
         IF (avgsnapshot(n)%T_cntr_begin.GE.Start_T_cntr) THEN
            current_avgsnap = n
            IF (Rank_of_process.EQ.0) PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: adjusted number of the first snapshot ",i4)', current_avgsnap
            EXIT
         END IF
      END DO
   ENDIF
  END IF

  CALL ADDITIONAL_PARAMETERS_AVG_SNAPHOTS

   IF (Rank_of_process==0) THEN
      IF ( use_checkpoint/=1 ) THEN ! Fresh start 
         DO idx=1,num_plane_x_locations+num_plane_y_locations
            WRITE( file_name,'(A,I3.3,A)') "Flux_through_plane_",idx,".dat"
            OPEN (21, FILE = file_name, STATUS = 'REPLACE')      
            IF ( idx<=num_plane_x_locations ) THEN
               WRITE( plane_name,'(A,ES15.7,A)') "X = ",plane_x_cuts_location(idx)*delta_x_m ," [m]"
            ELSE 
               WRITE( plane_name,'(A,ES15.7,A)') "Y = ",plane_y_cuts_location(idx-num_plane_x_locations)*delta_x_m," [m]"
            END IF
            WRITE (21,'(A)' ) plane_name
            WRITE( 21,'(A)', ADVANCE='no' ) "Time_counter,Time[s],electrons[s-1]"
            DO species=2,N_spec+1
               WRITE (21, '(A,I1,A)', ADVANCE='no' ) ",ions_species_",idx,"[s-1]"
            END DO       
            CLOSE (21, STATUS = 'KEEP')
         ENDDO 
      ELSE
         DO idx=1,num_plane_x_locations+num_plane_y_locations
            exists = .FALSE.
            WRITE( file_name,'(A,I3.3,A)') "Flux_through_plane_",idx,".dat"
            INQUIRE ( FILE = file_name, EXIST = exists)
            IF (exists) THEN                                                       
               OPEN (21, FILE = file_name, STATUS = 'OLD')         
               ! Skip header
               DO n=1,3
                  READ(21, '(A)') 
               END DO
               ! Read data
               DO 
                  READ (21, '(2x,i10,2x,ES18.10)', iostat = ios) i_dummy
                  IF (ios.NE.0) EXIT
                  IF (i_dummy.GE.Start_T_cntr) EXIT
               END DO
               BACKSPACE(21) ! Backspace is restored. I can safely eliminate last poitn because it should be recomputed as checkpoint automatically restarts at start of average window
               ENDFILE 21         
            ELSE       
               OPEN (21, FILE = file_name, STATUS = 'REPLACE')      
               IF ( idx<=num_plane_x_locations ) THEN
                  WRITE( plane_name,'(A,ES15.7,A)') "X = ",plane_x_cuts_location(idx)*delta_x_m ," [m]"
               ELSE 
                  WRITE( plane_name,'(A,ES15.7,A)') "Y = ",plane_y_cuts_location(idx-num_plane_x_locations)*delta_x_m," [m]"
               END IF
               WRITE (21,'(A)' ) plane_name
               WRITE( 21,'(A)', ADVANCE='no' ) "Time_counter,Time[s],electrons[s-1]"
               DO species=2,N_spec+1
                  WRITE (21, '(A,I1,A)', ADVANCE='no' ) ",ions_species_",idx,"[s-1]"
               END DO       
            END IF
            CLOSE (21, STATUS = 'KEEP')
         END DO
      END IF 
   END IF

END SUBROUTINE INITIATE_AVERAGED_SNAPSHOTS

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE ADDITIONAL_PARAMETERS_AVG_SNAPHOTS
!>    @details Read additional parameters related to avg_snapshots
!!    @authors W. Villafana
!!    @date    Apr-17-2024
!-------------------------------------------------------------------------------------------------- 
SUBROUTINE ADDITIONAL_PARAMETERS_AVG_SNAPHOTS

   USE mod_print, ONLY: print_message, print_parser_error
   USE AvgSnapshots, ONLY: plane_x_cuts_location, plane_y_cuts_location,num_plane_x_locations, num_plane_y_locations, save_collision_freq_ee, avg_flux_and_history
   USE CurrentProblemValues, ONLY: string_length, delta_x_m

   IMPLICIT NONE
   INCLUDE 'mpif.h'
   
   ! LOCAL
   INTEGER :: ierr
   CHARACTER (LEN=1000) :: long_buf,line,separator ! long buffer for string
   INTEGER :: i_found ! flag to decide if I found keyword or not.
   REAL(8) :: rval ! buffer for real values
   INTEGER :: ival ! buffer for integer values
   CHARACTER(LEN=string_length) :: caval ! buffer for string values
   CHARACTER(LEN=string_length) :: message, routine  
   LOGICAL :: exists
   INTEGER :: idx
   
   ! Declare routine name and debug level
   routine = 'ADDITIONAL_PARAMETERS_AVG_SNAPHOTS'
  
   separator = '='

   
   WRITE( message, '(A)') "Reading additional parameters for averaged snapshots"//achar(10)
   CALL print_message(message,routine)   

   INQUIRE (FILE = 'init_avgsnapshots.dat', EXIST = exists)
   OPEN (9, file='init_avgsnapshots.dat')

   i_found = 0
   REWIND(9)
   DO
      READ (9,"(A)",iostat=ierr) line ! read line into character variable
      IF ( ierr/=0 ) EXIT
      IF (line == '') CYCLE   ! Skip the rest of the loop if the line is empty. Will cause a crash
      READ (line,*) long_buf ! read first word of line
      IF ( TRIM(long_buf)=="planes_X_locations_to_measure_net_flux" ) THEN ! found search string at beginning of line
         i_found = 1

         CALL count_number_of_inputs(line,long_buf,separator,num_plane_x_locations)
         IF (num_plane_x_locations > 0) THEN
            ALLOCATE(plane_x_cuts_location(num_plane_x_locations))
         ELSE
            WRITE(message, '(A,A)') 'No values found in the line: ',line
            CALL print_parser_error(message)
         END IF            
         
         READ (line,*) long_buf,separator,plane_x_cuts_location(1:num_plane_x_locations)

         WRITE( message,'(A)') "planes_X_locations_to_measure_net_flux keyword present. Will generate fluxes on for cuts at X = ... if averaged snaphots are turned on"
         CALL print_message( message,routine )    
         
         DO idx=1,num_plane_x_locations
            WRITE( message,'(A,ES10.3,A,I10,A,I10)') "Plane cut at X = ",plane_x_cuts_location(idx)," [m] between grid node i = ",FLOOR(plane_x_cuts_location(idx)/delta_x_m)," and i = ",FLOOR(plane_x_cuts_location(idx)/delta_x_m)+1
            plane_x_cuts_location(idx) = plane_x_cuts_location(idx)/delta_x_m
            CALL print_message( message )    
         END DO            

      END IF
   END DO  
   IF ( i_found==0 ) THEN
      WRITE( message,'(A)') "planes_X_locations_to_measure_net_flux keyword not present. I assume we do not want it. "
      CALL print_message( message,routine )       
   END IF                          
   
   i_found = 0
   REWIND(9)
   DO
      READ (9,"(A)",iostat=ierr) line ! read line into character variable
      IF ( ierr/=0 ) EXIT
      IF (line == '') CYCLE   ! Skip the rest of the loop if the line is empty. Will cause a crash
      READ (line,*) long_buf ! read first word of line
      IF ( TRIM(long_buf)=="planes_Y_locations_to_measure_net_flux" ) THEN ! found search string at beginning of line
         i_found = 1
         CALL count_number_of_inputs(line,long_buf,separator,num_plane_y_locations)
         IF (num_plane_y_locations > 0) THEN
            ALLOCATE(plane_y_cuts_location(num_plane_y_locations))
         ELSE
            WRITE(message, '(A,A)') 'No values found in the line: ',line
            CALL print_parser_error(message)
         END IF                   
         READ (line,*) long_buf,separator,plane_y_cuts_location(1:num_plane_y_locations)

         WRITE( message,'(A)') "planes_Y_locations_to_measure_net_flux keyword present. Will generate fluxes on for cuts at Y = ... if averaged snaphots are turned on."
         CALL print_message( message,routine )    
         
         DO idx=1,num_plane_y_locations
            WRITE( message,'(A,ES10.3,A,I10,A,I10)') "Plane cut at Y = ",plane_y_cuts_location(idx)," [m] between grid node j = ",FLOOR(plane_y_cuts_location(idx)/delta_x_m)," and j = ",FLOOR(plane_y_cuts_location(idx)/delta_x_m)+1
            plane_y_cuts_location(idx) = plane_y_cuts_location(idx)/delta_x_m
            CALL print_message( message )    
         END DO            

      END IF
   END DO  
   IF ( i_found==0 ) THEN
      WRITE( message,'(A)') "planes_Y_locations_to_measure_net_flux keyword not present. I assume we do not want it. "
      CALL print_message( message,routine )       
   END IF     
   
   i_found = 0
   REWIND(9)
   DO
      READ (9,"(A)",iostat=ierr) line ! read line into character variable
      IF ( ierr/=0 ) EXIT
      IF (line == '') CYCLE   ! Skip the rest of the loop if the line is empty. Will cause a crash
      READ (line,*) long_buf ! read first word of line
      IF ( TRIM(long_buf)=="coulomb_collision_ee_freq" ) THEN ! found search string at beginning of line
         i_found = 1

         READ (line,*) long_buf,separator,caval            
         
         IF (TRIM(caval)=='yes') THEN
            save_collision_freq_ee = .TRUE.
            WRITE( message, '(A)') "Coulomb collision frequency 2D maps are saved and averaged on the fly."//achar(10)
         ELSE IF (TRIM(caval)=='no') THEN
            save_collision_freq_ee = .FALSE.
            WRITE( message, '(A)') "Coulomb collision frequency 2D maps are NOT saved."//achar(10)
         ELSE 
            WRITE( message, '(A,A,A)') "No sure if you want to save Coulomb collision frequency. For keyword 'coulomb_collision_ee_freq' I received: ",TRIM(caval),". Expected:'yes' or 'no'. Case sensitive"//achar(10)
            CALL print_parser_error(message)
         END IF
         CALL print_message( message )

      END IF
   END DO  
   IF ( i_found==0 ) THEN
      WRITE( message,'(A)') "coulomb_collision_ee_freq keyword not present. I assume we do not want it. "
      CALL print_message( message,routine )       
   END IF         
   
   i_found = 0
   REWIND(9)
   DO
      READ (9,"(A)",iostat=ierr) line ! read line into character variable
      IF ( ierr/=0 ) EXIT
      IF (line == '') CYCLE   ! Skip the rest of the loop if the line is empty. Will cause a crash
      READ (line,*) long_buf ! read first word of line
      IF ( TRIM(long_buf)=="avg_flux_and_history" ) THEN ! found search string at beginning of line
         i_found = 1

         READ (line,*) long_buf,separator,caval            
         
         IF (TRIM(caval)=='yes') THEN
            avg_flux_and_history = .TRUE.
            WRITE( message, '(A)') "Fluxes at the wall and reporting of particles will be averaged on the fly."//achar(10)
         ELSE IF (TRIM(caval)=='no') THEN
            avg_flux_and_history = .FALSE.
            WRITE( message, '(A)') "Fluxes at the wall and reporting of particles will be printed at each iteration."//achar(10)
         ELSE 
            WRITE( message, '(A,A,A)') "No sure if you want to avergae on the fly fluxes and reporting of particles. For keyword 'avg_flux_and_history' I received: ",TRIM(caval),". Expected:'yes' or 'no'. Case sensitive"//achar(10)
            CALL print_parser_error(message)
         END IF
         CALL print_message( message )

      END IF
   END DO  
   IF ( i_found==0 ) THEN
      WRITE( message,'(A)') "avg_flux_and_history keyword not present. I assume we do not want it. "
      CALL print_message( message,routine )       
   END IF                       

   CLOSE (9, STATUS = 'KEEP')   

END SUBROUTINE ADDITIONAL_PARAMETERS_AVG_SNAPHOTS

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE COLLECT_F_EX_EY_FOR_AVERAGED_SNAPSHOT

  USE ParallelOperationValues
  USE AvgSnapshots
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries
  USE IonParticles, ONLY : N_spec
  USE MCCollisions, ONLY : N_neutral_spec, collision_e_neutral
  USE Snapshots, ONLY : diagnostics_neutral


  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER n, p

  INTEGER ALLOC_ERR

  REAL, ALLOCATABLE :: cs_phi(:,:) 

  INTEGER i, j, k
  INTEGER recsize, bufsize
  REAL, ALLOCATABLE :: rbufer(:)

  INTEGER pos1, pos2

! quit if all snapshots were created or if due to any reasons the snapshot counter is 
! larger than the declared number of snapshots (e.g., when no snapshots are requested) 
  IF (current_avgsnap.GT.N_of_all_avgsnaps) RETURN

  IF (T_cntr.LT.avgsnapshot(current_avgsnap)%T_cntr_begin) RETURN

  IF ((T_cntr.EQ.avgsnapshot(current_avgsnap)%T_cntr_begin).AND.(cluster_rank_key.EQ.0)) THEN

     IF (save_avg_data(1)) THEN
        ALLOCATE(cs_avg_phi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_phi = 0.0
     END IF

     IF (save_avg_data(2)) THEN 
        ALLOCATE(cs_avg_EX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_EX = 0.0
     END IF

     IF (save_avg_data(3)) THEN
        ALLOCATE(cs_avg_EY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_EY = 0.0
     END IF

     IF (save_avg_data(7)) THEN 
        ALLOCATE(cs_avg_Ne(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_Ne = 0.0
     END IF

     IF (save_avg_data(4).OR.save_avg_data(8)) THEN
        ALLOCATE(cs_avg_JXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_JXe = 0.0
     END IF

     IF (save_avg_data(5).OR.save_avg_data(9)) THEN
        ALLOCATE(cs_avg_JYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_JYe = 0.0
     END IF

     IF (save_avg_data(6).OR.save_avg_data(10)) THEN
        ALLOCATE(cs_avg_JZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_JZe = 0.0
     END IF

     IF (save_avg_data(11)) THEN
        ALLOCATE(cs_avg_VXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_VXe = 0.0
     END IF

     IF (save_avg_data(12)) THEN
        ALLOCATE(cs_avg_VYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_VYe = 0.0
     END IF

     IF (save_avg_data(13)) THEN
        ALLOCATE(cs_avg_VZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_VZe = 0.0
     END IF

     IF (save_avg_data(14)) THEN
        ALLOCATE(cs_avg_WXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_WXe = 0.0
     END IF

     IF (save_avg_data(15)) THEN
        ALLOCATE(cs_avg_WYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_WYe = 0.0
     END IF

     IF (save_avg_data(16)) THEN
        ALLOCATE(cs_avg_WZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_WZe = 0.0
     END IF

     IF (save_avg_data(17)) THEN
        ALLOCATE(cs_avg_TXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_TXe = 0.0
     END IF

     IF (save_avg_data(18)) THEN
        ALLOCATE(cs_avg_TYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_TYe = 0.0
     END IF

     IF (save_avg_data(19)) THEN
        ALLOCATE(cs_avg_TZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_TZe = 0.0
     END IF

     IF (save_avg_data(20)) THEN
        ALLOCATE(cs_avg_QXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_QXe = 0.0
     END IF

     IF (save_avg_data(21)) THEN
        ALLOCATE(cs_avg_QYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_QYe = 0.0
     END IF

     IF (save_avg_data(22)) THEN
        ALLOCATE(cs_avg_QZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_QZe = 0.0
     END IF

     IF (save_avg_data(23)) THEN
        ALLOCATE(cs_avg_Ni(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_Ni = 0.0
     END IF

     IF (save_avg_data(4).OR.save_avg_data(24)) THEN
        ALLOCATE(cs_avg_JXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_JXi = 0.0
     END IF

     IF (save_avg_data(5).OR.save_avg_data(25)) THEN
        ALLOCATE(cs_avg_JYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_JYi = 0.0
     END IF

     IF (save_avg_data(6).OR.save_avg_data(26)) THEN
        ALLOCATE(cs_avg_JZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_JZi = 0.0
     END IF

     IF (save_avg_data(27)) THEN
        ALLOCATE(cs_avg_VXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_VXi = 0.0
     END IF

     IF (save_avg_data(28)) THEN
        ALLOCATE(cs_avg_VYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_VYi = 0.0
     END IF

     IF (save_avg_data(29)) THEN
        ALLOCATE(cs_avg_VZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_VZi = 0.0
     END IF

     IF (save_avg_data(30)) THEN
        ALLOCATE(cs_avg_WXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_WXi = 0.0
     END IF

     IF (save_avg_data(31)) THEN
        ALLOCATE(cs_avg_WYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_WYi = 0.0
     END IF

     IF (save_avg_data(32)) THEN
        ALLOCATE(cs_avg_WZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_WZi = 0.0
     END IF

     IF (save_avg_data(33)) THEN
        ALLOCATE(cs_avg_TXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_TXi = 0.0
     END IF

     IF (save_avg_data(34)) THEN
        ALLOCATE(cs_avg_TYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_TYi = 0.0
     END IF

     IF (save_avg_data(35)) THEN
        ALLOCATE(cs_avg_TZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_TZi = 0.0
     END IF

     IF (save_avg_data(36)) THEN
        ALLOCATE(cs_avg_QXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_QXi = 0.0
     END IF

     IF (save_avg_data(37)) THEN
        ALLOCATE(cs_avg_QYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_QYi = 0.0
     END IF

     IF (save_avg_data(38)) THEN
        ALLOCATE(cs_avg_QZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_QZi = 0.0
     END IF

     IF (save_avg_data(39)) THEN
        DO n = 1, N_neutral_spec
           DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
              IF (.NOT.collision_e_neutral(n)%colproc_info(p)%save_collfreq_2d) CYCLE
              ALLOCATE(diagnostics_neutral(n)%activated_collision(p)%coll_freq_local(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
              diagnostics_neutral(n)%activated_collision(p)%coll_freq_local = 0.0
           END DO
        END DO
     END IF

      IF ( save_collision_freq_ee ) THEN
         ALLOCATE(cs_avg_coulomb_ee(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
         cs_avg_coulomb_ee = 0.0
      ENDIF

      IF ( num_plane_x_locations+num_plane_y_locations>0 ) THEN
         IF (ALLOCATED(cluster_flux_through_plane_over_one_period)) DEALLOCATE(cluster_flux_through_plane_over_one_period)
         ALLOCATE(cluster_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations))
         ALLOCATE(total_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations))
         total_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations) = 0.0
      ENDIF
  END IF   !### IF ((T_cntr.EQ.avgsnapshot(current_avgsnap)%T_cntr_begin).AND.(cluster_rank_key.EQ.0)) THEN

! potential

  IF (save_avg_data(1)) THEN

     IF (cluster_rank_key.EQ.0) THEN

        ALLOCATE(cs_phi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)

        IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN
! PETSc-based field solver
! account for own contribution of the master as a field calculator
           DO j = indx_y_min, indx_y_max
              DO i = indx_x_min, indx_x_max
                 cs_phi(i,j) = REAL(phi(i,j))
              END DO
           END DO

! account for contributions from other field calculators
           DO k = 2, cluster_N_blocks
              recsize = field_calculator(k)%indx_x_max - field_calculator(k)%indx_x_min + 1
              bufsize = recsize * (field_calculator(k)%indx_y_max - field_calculator(k)%indx_y_min + 1)
              ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)
              CALL MPI_RECV(rbufer, bufsize, MPI_REAL, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
              pos2 = 0
              DO j = field_calculator(k)%indx_y_min, field_calculator(k)%indx_y_max
                 pos1 = pos2 + 1
                 pos2 = pos2 + recsize
                 cs_phi(field_calculator(k)%indx_x_min:field_calculator(k)%indx_x_max,j) = rbufer(pos1:pos2) 
              END DO
              DEALLOCATE(rbufer, STAT = ALLOC_ERR)
           END DO

        ELSE
! FFT-based field solver
! the cluster already knows the potential within its domain, use it -----------------
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_phi(i,j) = REAL(c_phi(i,j))
              END DO
           END DO
        END IF   !### IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

        cs_avg_phi = cs_avg_phi + cs_phi
! cleanup
        DEALLOCATE(cs_phi, STAT = ALLOC_ERR)

     ELSE   !### IF (cluster_rank_key.EQ.0) THEN

        IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

           recsize = indx_x_max - indx_x_min + 1
           bufsize = recsize * (indx_y_max - indx_y_min + 1)

           ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

           pos1 = 1 - indx_x_min
           DO j = indx_y_min, indx_y_max
              DO i = indx_x_min, indx_x_max
                 rbufer(pos1 + i) = REAL(phi(i, j))
              END DO
              pos1 = pos1 + recsize
           END DO

           CALL MPI_SEND(rbufer, bufsize, MPI_REAL, field_master, Rank_of_process, MPI_COMM_WORLD, request, ierr) 

           DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        
        END IF

     END IF   !###   IF (cluster_rank_key.EQ.0) THEN

  END IF   !###   IF (save_avg_data(1)) THEN

  IF (cluster_rank_key.NE.0) THEN
! non-master processes do nothing and can leave
     RETURN
  END IF

! electric field components

  IF (save_avg_data(2)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_EX(i,j) = cs_avg_EX(i,j) + REAL(EX(i,j))
        END DO
     END DO
  END IF

  IF (save_avg_data(3)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_EY(i,j) = cs_avg_EY(i,j) + REAL(EY(i,j))
        END DO
     END DO
  END IF

END SUBROUTINE COLLECT_F_EX_EY_FOR_AVERAGED_SNAPSHOT

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE COLLECT_ELECTRON_DATA_FOR_AVERAGED_SNAPSHOT

  USE AvgSnapshots
  USE IonParticles, ONLY: N_spec
  USE Snapshots, ONLY : cs_N, cs_VX, cs_VY, cs_VZ, cs_WX, cs_WY, cs_WZ, cs_VXVY, cs_VXVZ, cs_VYVZ, cs_QX, cs_QY, cs_QZ
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INTEGER i, j

! checks are performed before calling this procedure
! arrays cs_* are allocated before calling this procedure

! electron VDF moments are already collected in ADVANCE_ELECTRONS_AND_COLLECT_MOMENTS

! electric current along each coordinate direction

  IF (save_avg_data(4).OR.save_avg_data(8)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JXe(i,j) = cs_avg_JXe(i,j) - cs_N(i,j) * cs_VX(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(5).OR.save_avg_data(9)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JYe(i,j) = cs_avg_JYe(i,j) - cs_N(i,j) * cs_VY(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(6).OR.save_avg_data(10)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JZe(i,j) = cs_avg_JZe(i,j) - cs_N(i,j) * cs_VZ(i,j)
        END DO
     END DO
  END IF

! number density

  IF (save_avg_data(7)) THEN
     cs_avg_Ne = cs_avg_Ne + cs_N
  END IF

! average velocity along each coordinate direction

  IF (save_avg_data(11)) THEN
     cs_avg_VXe = cs_avg_VXe + cs_VX
  END IF

  IF (save_avg_data(12)) THEN
     cs_avg_VYe = cs_avg_VYe + cs_VY
  END IF

  IF (save_avg_data(13)) THEN
     cs_avg_VZe = cs_avg_VZe + cs_VZ
  END IF

! average energy of motion along each coordinate direction

  IF (save_avg_data(14)) THEN
     cs_avg_WXe = cs_avg_WXe + cs_WX
  END IF

  IF (save_avg_data(15)) THEN
     cs_avg_WYe = cs_avg_WYe + cs_WY
  END IF

  IF (save_avg_data(16)) THEN
     cs_avg_WZe = cs_avg_WZe + cs_WZ
  END IF

! temperature along each coordinate direction

  IF (save_avg_data(17)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TXe(i,j) = cs_avg_TXe(i,j) + MAX(0.0, cs_WX(i,j) - cs_VX(i,j)**2)
        END DO
     END DO
  END IF

  IF (save_avg_data(18)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TYe(i,j) = cs_avg_TYe(i,j) + MAX(0.0, cs_WY(i,j) - cs_VY(i,j)**2)
        END DO
     END DO
  END IF

  IF (save_avg_data(19)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TZe(i,j) = cs_avg_TZe(i,j) + MAX(0.0, cs_WZ(i,j) - cs_VZ(i,j)**2)
        END DO
     END DO
  END IF

! heat flow along each coordinate direction

  IF (save_avg_data(20)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QXe(i,j) = cs_avg_QXe(i,j) + cs_N(i,j) * cs_QX(i,j)
        END DO
     END DO
  END IF
     
  IF (save_avg_data(21)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QYe(i,j) = cs_avg_QYe(i,j) + cs_N(i,j) * cs_QY(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(22)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QZe(i,j) = cs_avg_QZe(i,j) + cs_N(i,j) * cs_QZ(i,j)
        END DO
     END DO
  END IF

   IF ( num_plane_x_locations+num_plane_y_locations>0 ) THEN 
      total_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations) = total_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations) +&
                                                                                                      cluster_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations)
   END IF

END SUBROUTINE COLLECT_ELECTRON_DATA_FOR_AVERAGED_SNAPSHOT

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE COLLECT_ION_DATA_FOR_AVERAGED_SNAPSHOT(s)

  USE AvgSnapshots
  USE Snapshots, ONLY : cs_N, cs_VX, cs_VY, cs_VZ, cs_WX, cs_WY, cs_WZ, cs_VXVY, cs_VXVZ, cs_VYVZ, cs_QX, cs_QY, cs_QZ
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY: N_spec

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: s

  INTEGER i, j

! checks are performed before calling this procedure
! arrays cs_* are allocated before calling this procedure

! VDF moments for ion species s are already collected in ADVANCE_IONS_AND_COLLECT_MOMENTS

! number density

  IF (save_avg_data(23)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_Ni(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_Ni(c_indx_x_min:c_indx_x_max,j,s) + cs_N(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

! electric current along each coordinate direction

  IF (save_avg_data(4).OR.save_avg_data(24)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JXi(i,j,s) = cs_avg_JXi(i,j,s) + cs_N(i,j) * cs_VX(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(5).OR.save_avg_data(25)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JYi(i,j,s) = cs_avg_JYi(i,j,s) + cs_N(i,j) * cs_VY(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(6).OR.save_avg_data(26)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JZi(i,j,s) = cs_avg_JZi(i,j,s) + cs_N(i,j) * cs_VZ(i,j)
        END DO
     END DO
  END IF

! average velocity along each coordinate direction

  IF (save_avg_data(27)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_VXi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_VXi(c_indx_x_min:c_indx_x_max,j,s) + cs_VX(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

  IF (save_avg_data(28)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_VYi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_VYi(c_indx_x_min:c_indx_x_max,j,s) + cs_VY(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

  IF (save_avg_data(29)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_VZi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_VZi(c_indx_x_min:c_indx_x_max,j,s) + cs_VZ(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

! average energy of motion along each coordinate direction

  IF (save_avg_data(30)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_WXi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_WXi(c_indx_x_min:c_indx_x_max,j,s) + cs_WX(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

  IF (save_avg_data(31)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_WYi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_WYi(c_indx_x_min:c_indx_x_max,j,s) + cs_WY(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

  IF (save_avg_data(32)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_WZi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_WZi(c_indx_x_min:c_indx_x_max,j,s) + cs_WZ(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

! temperature along each coordinate direction

  IF (save_avg_data(33)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TXi(i,j,s) = cs_avg_TXi(i,j,s) + MAX(0.0, cs_WX(i,j) - cs_VX(i,j)**2)
        END DO
     END DO
  END IF

  IF (save_avg_data(34)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TYi(i,j,s) = cs_avg_TYi(i,j,s) + MAX(0.0, cs_WY(i,j) - cs_VY(i,j)**2)
        END DO
     END DO
  END IF

  IF (save_avg_data(35)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TZi(i,j,s) = cs_avg_TZi(i,j,s) + MAX(0.0, cs_WZ(i,j) - cs_VZ(i,j)**2)
        END DO
     END DO
  END IF

! heat flow along each coordinate direction

  IF (save_avg_data(36)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QXi(i,j,s) = cs_avg_QXi(i,j,s) + cs_N(i,j) * cs_QX(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(37)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QYi(i,j,s) = cs_avg_QYi(i,j,s) + cs_N(i,j) * cs_QY(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(38)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QZi(i,j,s) = cs_avg_QZi(i,j,s) + cs_N(i,j) * cs_QZ(i,j)
        END DO
     END DO
  END IF

   IF ( num_plane_x_locations+num_plane_y_locations>0 ) THEN 
      total_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations) = total_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations) +&
                                                                                                   cluster_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations)
   END IF  

END SUBROUTINE COLLECT_ION_DATA_FOR_AVERAGED_SNAPSHOT

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE CREATE_AVERAGED_SNAPSHOT

  USE ParallelOperationValues
  USE AvgSnapshots
  USE CurrentProblemValues, ONLY : N_subcycles, F_scale_V, E_scale_Vm, current_factor_Am2, N_scale_part_m3, V_scale_ms, &
                                 & energy_factor_eV, temperature_factor_eV, heat_flow_factor_Wm2, T_cntr, delta_t_s, string_length, weight_ptcl
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_spec, Qs, Ms
  USE MCCollisions, ONLY : N_neutral_spec, neutral, collision_e_neutral
  USE Snapshots, ONLY : diagnostics_neutral
  USE mod_print, ONLY: print_message

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER N_averaged_timesteps, N_averaged_timesteps_i  ! different for electrons and ions
  REAL avg_factor
  INTEGER ALLOC_ERR

  REAL, ALLOCATABLE :: cs_avg_Jsum(:,:)
  REAL, ALLOCATABLE :: final_flux_through_plane_over_one_period(:,:)

!                                 ! ----x----I----x----I--
  CHARACTER(20) filename_F        ! _NNNN_avg_F_V_2D.bin
  CHARACTER(22) filename_E        ! _NNNN_avg_EX_Vm_2D.bin
!                                 ! ----x----I----x----I----x-
  CHARACTER(26) filename_Jsum     ! _NNNN_avg_JXsum_Am2_2D.bin
!                                 ! ----x----I----x----I----
  CHARACTER(22) filename_Ne       ! _NNNN_avg_Ne_m3_2D.bin
  CHARACTER(24) filename_Je       ! _NNNN_avg_JXe_Am2_2D.bin
  CHARACTER(23) filename_Ve       ! _NNNN_avg_VXe_ms_2D.bin
  CHARACTER(23) filename_We       ! _NNNN_avg_WXe_eV_2D.bin
  CHARACTER(23) filename_Te       ! _NNNN_avg_TXe_eV_2D.bin
  CHARACTER(24) filename_Qe       ! _NNNN_avg_QXe_Wm2_2D.bin
!                                 ! ----x----I----x----I----x-
  CHARACTER(24) filename_Ni       ! _NNNN_avg_Ni_s_m3_2D.bin
  CHARACTER(26) filename_Ji       ! _NNNN_avg_JXi_s_Am2_2D.bin
  CHARACTER(25) filename_Vi       ! _NNNN_avg_VXi_s_ms_2D.bin
  CHARACTER(25) filename_Wi       ! _NNNN_avg_WXi_s_eV_2D.bin
  CHARACTER(25) filename_Ti       ! _NNNN_avg_TXi_s_eV_2D.bin
  CHARACTER(26) filename_Qi       ! _NNNN_avg_QXi_s_Wm2_2D.bin
!                                   ----x----I----x----I----x----I----x----I----x
  CHARACTER(45) filename_encoll   ! _NNNN_avg_frequency_e_n_AAAAAA_coll_id_NN.bin

  CHARACTER(LEN=string_length) :: filename_generic

  INTEGER s, i, j, n, p

  INTEGER :: idx, species
  CHARACTER(LEN=string_length) :: file_name
  CHARACTER(LEN=string_length) :: message

  INTEGER :: avg_output_flag

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  
  CALL DETERMINE_AVG_DATA_CREATION(avg_output_flag) 
  IF (avg_output_flag==0) RETURN
! ! quit if all snapshots were created or if due to any reasons the snapshot counter is 
! ! larger than the declared number of snapshots (e.g., when no snapshots are requested) 
!   IF (current_avgsnap.GT.N_of_all_avgsnaps) RETURN

! ! quit if the current moment of time is not the moment when it is necessary to create the snapshot
!   IF (T_cntr.NE.avgsnapshot(current_avgsnap)%T_cntr_end) RETURN

  IF (cluster_rank_key.NE.0) THEN
! processes which are not cluster masters do not participate in saving data
     current_avgsnap = current_avgsnap + 1           ! increase the snapshots counter 
     RETURN
  END IF

  IF (Rank_of_process.EQ.0) PRINT '("### ^^^^^^^^^^^^^^^ Averaged Snapshot ",i4," will be created now ... ^^^^^^^^^^^^^^^^ ###")', current_avgsnap

  N_averaged_timesteps   =  avgsnapshot(current_avgsnap)%T_cntr_end - avgsnapshot(current_avgsnap)%T_cntr_begin + 1                  ! for fields and electron moments

!  N_averaged_timesteps_i = (avgsnapshot(current_avgsnap)%T_cntr_end - avgsnapshot(current_avgsnap)%T_cntr_begin) / N_subcycles + 1   ! for ion moments
  N_averaged_timesteps_i = (avgsnapshot(current_avgsnap)%T_cntr_end - avgsnapshot(current_avgsnap)%T_cntr_begin + 1) / N_subcycles   ! for ion moments

! potential

  IF (save_avg_data(1)) THEN
     filename_F = '_NNNN_avg_F_V_2D.bin'
     filename_F(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor =  REAL(F_scale_V) / N_averaged_timesteps
     cs_avg_phi = cs_avg_phi * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_phi, filename_F)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_phi, STAT = ALLOC_ERR)
  END IF

! electric field components

  IF (save_avg_data(2)) THEN
     filename_E = '_NNNN_avg_EX_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(E_scale_Vm) / N_averaged_timesteps
     cs_avg_EX = cs_avg_EX * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_EX, filename_E)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_EX, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(3)) THEN
     filename_E = '_NNNN_avg_EY_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(E_scale_Vm) / N_averaged_timesteps
     cs_avg_EY = cs_avg_EY * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_EY, filename_E)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_EY, STAT = ALLOC_ERR)
  END IF

! FULL electric current (sum of electron and ion currents) along each coordinate direction
! restored here from electron and ion currents

  IF (save_avg_data(4)) THEN
     filename_Jsum = '_NNNN_avg_JXsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     ALLOCATE(cs_avg_Jsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps 
     cs_avg_Jsum = cs_avg_JXe * avg_factor
     DO s = 1, N_spec
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_Jsum(i, j) = cs_avg_Jsum(i, j) + cs_avg_JXi(i,j,s) * avg_factor
           END DO
        END DO
     END DO
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Jsum, filename_Jsum)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_Jsum, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(8) ) DEALLOCATE(cs_avg_JXe, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(24)) DEALLOCATE(cs_avg_JXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(5)) THEN
     filename_Jsum = '_NNNN_avg_JYsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     ALLOCATE(cs_avg_Jsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps 
     cs_avg_Jsum = cs_avg_JYe * avg_factor
     DO s = 1, N_spec
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_Jsum(i, j) = cs_avg_Jsum(i, j) + cs_avg_JYi(i,j,s) * avg_factor
           END DO
        END DO
     END DO
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Jsum, filename_Jsum)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_Jsum, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(9) ) DEALLOCATE(cs_avg_JYe, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(25)) DEALLOCATE(cs_avg_JYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(6)) THEN
     filename_Jsum = '_NNNN_avg_JZsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     ALLOCATE(cs_avg_Jsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps 
     cs_avg_Jsum = cs_avg_JZe * avg_factor
     DO s = 1, N_spec
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_Jsum(i, j) = cs_avg_Jsum(i, j) + cs_avg_JZi(i,j,s) * avg_factor
           END DO
        END DO
     END DO
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Jsum, filename_Jsum)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_Jsum, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(10)) DEALLOCATE(cs_avg_JZe, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(26)) DEALLOCATE(cs_avg_JZi, STAT = ALLOC_ERR)
  END IF

! electron number density

  IF (save_avg_data(7)) THEN
     filename_Ne = '_NNNN_avg_Ne_m3_2D.bin'
     filename_Ne(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(N_scale_part_m3) / N_averaged_timesteps
     cs_avg_Ne = cs_avg_Ne * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Ne, filename_Ne)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_Ne, STAT = ALLOC_ERR)
  END IF

! electron electric current along each coordinate direction

  IF (save_avg_data(8)) THEN
     filename_Je = '_NNNN_avg_JXe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps
     cs_avg_JXe = cs_avg_JXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JXe, filename_Je)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_JXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(9)) THEN
     filename_Je = '_NNNN_avg_JYe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps
     cs_avg_JYe = cs_avg_JYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JYe, filename_Je)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_JYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(10)) THEN
     filename_Je = '_NNNN_avg_JZe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps
     cs_avg_JZe = cs_avg_JZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JZe, filename_Je)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_JZe, STAT = ALLOC_ERR)
  END IF

! electron average velocity along each coordinate direction

  IF (save_avg_data(11)) THEN
     filename_Ve = '_NNNN_avg_VXe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps
     cs_avg_VXe = cs_avg_VXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VXe, filename_Ve)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_VXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(12)) THEN
     filename_Ve = '_NNNN_avg_VYe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps
     cs_avg_VYe = cs_avg_VYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VYe, filename_Ve)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_VYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(13)) THEN
     filename_Ve = '_NNNN_avg_VZe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps
     cs_avg_VZe = cs_avg_VZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VZe, filename_Ve)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_VZe, STAT = ALLOC_ERR)
  END IF

! electron average energy of motion along each coordinate direction

  IF (save_avg_data(14)) THEN
     filename_We = '_NNNN_avg_WXe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(energy_factor_eV) / N_averaged_timesteps
     cs_avg_WXe = cs_avg_WXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WXe, filename_We)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_WXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(15)) THEN
     filename_We = '_NNNN_avg_WYe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(energy_factor_eV) / N_averaged_timesteps
     cs_avg_WYe = cs_avg_WYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WYe, filename_We)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_WYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(16)) THEN
     filename_We = '_NNNN_avg_WZe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(energy_factor_eV) / N_averaged_timesteps
     cs_avg_WZe = cs_avg_WZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WZe, filename_We)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_WZe, STAT = ALLOC_ERR)
  END IF

! electron temperature along each coordinate direction

  IF (save_avg_data(17)) THEN
     filename_Te = '_NNNN_avg_TXe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(temperature_factor_eV) / N_averaged_timesteps
     cs_avg_TXe = cs_avg_TXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TXe, filename_Te)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_TXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(18)) THEN
     filename_Te = '_NNNN_avg_TYe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(temperature_factor_eV) / N_averaged_timesteps
     cs_avg_TYe = cs_avg_TYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TYe, filename_Te)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_TYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(19)) THEN
     filename_Te = '_NNNN_avg_TZe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(temperature_factor_eV) / N_averaged_timesteps
     cs_avg_TZe = cs_avg_TZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TZe, filename_Te)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_TZe, STAT = ALLOC_ERR)
  END IF

! heat flow along each coordinate direction

  IF (save_avg_data(20)) THEN
     filename_Qe = '_NNNN_avg_QXe_Wm2_2D.bin'
     filename_Qe(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(heat_flow_factor_Wm2) / N_averaged_timesteps
     cs_avg_QXe = cs_avg_QXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QXe, filename_Qe)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_QXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(21)) THEN
     filename_Qe = '_NNNN_avg_QYe_Wm2_2D.bin'
     filename_Qe(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(heat_flow_factor_Wm2) / N_averaged_timesteps
     cs_avg_QYe = cs_avg_QYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QYe, filename_Qe)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_QYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(22)) THEN
     filename_Qe = '_NNNN_avg_QZe_Wm2_2D.bin'
     filename_Qe(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(heat_flow_factor_Wm2) / N_averaged_timesteps
     cs_avg_QZe = cs_avg_QZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QZe, filename_Qe)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_QZe, STAT = ALLOC_ERR)
  END IF

! ion number density

  IF (save_avg_data(23)) THEN
     avg_factor = REAL(N_scale_part_m3) / N_averaged_timesteps_i
     cs_avg_Ni = cs_avg_Ni * avg_factor
     DO s = 1, N_spec
        filename_Ni = '_NNNN_avg_Ni_s_m3_2D.bin'
        filename_Ni(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ni(14:14) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Ni(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ni)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr) 
     END DO
! cleanup
     DEALLOCATE(cs_avg_Ni, STAT = ALLOC_ERR)
  END IF

! ion electric current along each coordinate direction

  IF (save_avg_data(24)) THEN
     DO s = 1, N_spec
        filename_Ji = '_NNNN_avg_JXi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ji(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_JXi(i,j,s) = cs_avg_JXi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ji)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_JXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(25)) THEN
     DO s = 1, N_spec
        filename_Ji = '_NNNN_avg_JYi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ji(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_JYi(i,j,s) = cs_avg_JYi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ji)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_JYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(26)) THEN
     DO s = 1, N_spec
        filename_Ji = '_NNNN_avg_JZi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ji(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_JZi(i,j,s) = cs_avg_JZi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ji)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_JZi, STAT = ALLOC_ERR)
  END IF

! ion average velocity along each coordinate direction

  IF (save_avg_data(27)) THEN
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps_i
     cs_avg_VXi = cs_avg_VXi * avg_factor
     DO s = 1, N_spec
        filename_Vi = '_NNNN_avg_VXi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Vi(15:15) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Vi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_VXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(28)) THEN
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps_i
     cs_avg_VYi = cs_avg_VYi * avg_factor
     DO s = 1, N_spec
        filename_Vi = '_NNNN_avg_VYi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Vi(15:15) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Vi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_VYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(29)) THEN
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps_i
     cs_avg_VZi = cs_avg_VZi * avg_factor
     DO s = 1, N_spec
        filename_Vi = '_NNNN_avg_VZi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Vi(15:15) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Vi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr) 
     END DO
! cleanup
     DEALLOCATE(cs_avg_VZi, STAT = ALLOC_ERR)
  END IF

! ion average energy of motion along each coordinate direction

  IF (save_avg_data(30)) THEN
     DO s = 1, N_spec
        filename_Wi = '_NNNN_avg_WXi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Wi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * energy_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_WXi(i,j,s) = cs_avg_WXi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Wi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_WXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(31)) THEN
     DO s = 1, N_spec
        filename_Wi = '_NNNN_avg_WYi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Wi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * energy_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_WYi(i,j,s) = cs_avg_WYi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Wi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_WYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(32)) THEN
     DO s = 1, N_spec
        filename_Wi = '_NNNN_avg_WZi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Wi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * energy_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_WZi(i,j,s) = cs_avg_WZi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Wi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_WZi, STAT = ALLOC_ERR)
  END IF

! ion temperature along each coordinate direction

  IF (save_avg_data(33)) THEN
     DO s = 1, N_spec
        filename_Ti = '_NNNN_avg_TXi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ti(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * temperature_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_TXi(i,j,s) = cs_avg_TXi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ti)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_TXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(34)) THEN
     DO s = 1, N_spec
        filename_Ti = '_NNNN_avg_TYi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ti(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * temperature_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_TYi(i,j,s) = cs_avg_TYi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ti)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_TYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(35)) THEN
     DO s = 1, N_spec
        filename_Ti = '_NNNN_avg_TZi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ti(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * temperature_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_TZi(i,j,s) = cs_avg_TZi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ti)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_TZi, STAT = ALLOC_ERR)
  END IF

! heat flow along each coordinate direction

  IF (save_avg_data(36)) THEN
     DO s = 1, N_spec
        filename_Qi = '_NNNN_avg_QXi_s_Wm2_2D.bin'
        filename_Qi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Qi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * heat_flow_factor_Wm2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_QXi(i,j,s) = cs_avg_QXi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Qi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_QXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(37)) THEN
     DO s = 1, N_spec
        filename_Qi = '_NNNN_avg_QYi_s_Wm2_2D.bin'
        filename_Qi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Qi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * heat_flow_factor_Wm2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_QYi(i,j,s) = cs_avg_QYi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Qi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_QYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(38)) THEN
     DO s = 1, N_spec
        filename_Qi = '_NNNN_avg_QZi_s_Wm2_2D.bin'
        filename_Qi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Qi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * heat_flow_factor_Wm2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_QZi(i,j,s) = cs_avg_QZi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Qi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_QZi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(39)) THEN
     DO n = 1, N_neutral_spec
        DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
           IF (.NOT.collision_e_neutral(n)%colproc_info(p)%save_collfreq_2d) CYCLE       
           filename_encoll = '_NNNN_avg_frequency_e_n_AAAAAA_coll_id_NN.bin'
           filename_encoll(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
           filename_encoll(25:30) = neutral(n)%name
           filename_encoll(40:41) = convert_int_to_txt_string(collision_e_neutral(n)%colproc_info(p)%id_number, 2)
           avg_factor = 1.0 / REAL(delta_t_s * N_averaged_timesteps)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 diagnostics_neutral(n)%activated_collision(p)%coll_freq_local(i,j) = diagnostics_neutral(n)%activated_collision(p)%coll_freq_local(i,j) * avg_factor
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(diagnostics_neutral(n)%activated_collision(p)%coll_freq_local, filename_encoll)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
           DEALLOCATE(diagnostics_neutral(n)%activated_collision(p)%coll_freq_local, STAT = ALLOC_ERR)
        END DO
     END DO
  END IF

   IF (save_collision_freq_ee) THEN
      filename_generic = '_NNNN_avg_frequency_e_e_coll.bin'
      filename_generic(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
      avg_factor = 1.0 / REAL(N_averaged_timesteps)
      ! avg_factor = 1.0 / REAL(delta_t_s * N_averaged_timesteps)
      DO j = c_indx_y_min, c_indx_y_max
         DO i = c_indx_x_min, c_indx_x_max
            cs_avg_coulomb_ee(i,j) = cs_avg_coulomb_ee(i,j) * avg_factor
         END DO
      END DO
      CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_coulomb_ee, filename_generic)
      CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
      DEALLOCATE(cs_avg_coulomb_ee, STAT = ALLOC_ERR)
   END IF

   IF (num_plane_x_locations+num_plane_y_locations>0) THEN
      ALLOCATE(final_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations))
      final_flux_through_plane_over_one_period(1:N_spec+1,1:num_plane_x_locations+num_plane_y_locations) = 0.0
      CALL MPI_REDUCE(total_flux_through_plane_over_one_period, final_flux_through_plane_over_one_period, (N_spec+1)*(num_plane_x_locations+num_plane_y_locations), MPI_REAL, MPI_SUM, 0, COMM_HORIZONTAL, ierr)
   ENDIF

   IF (Rank_horizontal==0) THEN 
      avg_factor = 1.0 / REAL(delta_t_s * N_averaged_timesteps)
      DO idx = 1,num_plane_x_locations+num_plane_y_locations
         WRITE( file_name,'(A,I3.3,A)') "Flux_through_plane_",idx,".dat"       
         WRITE( message,'(A)') "Saving files "//TRIM(file_name)
         CALL print_message( message )   
         OPEN (21, file = file_name, position = 'append')
            WRITE (21, '(2x,i10,2x,ES18.10)', ADVANCE='no' ) T_cntr, T_cntr * delta_t_s
            DO species=1,N_spec+1
               WRITE (21, '(2x,ES18.10)', ADVANCE='no' ) final_flux_through_plane_over_one_period(species,idx)*avg_factor*weight_ptcl
            END DO 
         CLOSE (21, status = 'keep')      
      ENDDO
   ENDIF

   IF (ALLOCATED(final_flux_through_plane_over_one_period)) DEALLOCATE(final_flux_through_plane_over_one_period)
   IF (ALLOCATED(total_flux_through_plane_over_one_period)) DEALLOCATE(total_flux_through_plane_over_one_period)

  IF (Rank_of_process.EQ.0) PRINT '(/2x,"### ^^^^^^^^^^^^^^^^^^^^ Averaged Snapshot ",i4," completed :) ^^^^^^^^^^^^^^^^^^^ ###")', current_avgsnap
   
  current_avgsnap = current_avgsnap + 1           ! increase the snapshots counter (cluster masters only, other processes did this already)
  
END SUBROUTINE CREATE_AVERAGED_SNAPSHOT

! --------------------------------------------------------------------------------------------------
!     SUBROUTINE DETERMINE_AVG_DATA_CREATION
! >    @details Determines if averaged snapshots or any quantity using the same time table should be created 
! !    @authors W. Villafana
! !    @date    Aug-25-2024
! -------------------------------------------------------------------------------------------------- 

SUBROUTINE DETERMINE_AVG_DATA_CREATION(avg_output_flag)

   USE AvgSnapshots, ONLY: current_avgsnap, N_of_all_avgsnaps, avgsnapshot
   USE CurrentProblemValues, ONLY: T_cntr, string_length
   USE mod_print, ONLY: print_debug
   IMPLICIT NONE

   !IN/OUT
   INTEGER, INTENT(OUT) :: avg_output_flag

   !LOCAL   
   CHARACTER(LEN=string_length) :: routine
   INTEGER :: local_debug_level   

   routine = "DETERMINE_AVG_DATA_CREATION"
   local_debug_level = 2

   CALL print_debug( routine,local_debug_level)   

   avg_output_flag = 1 ! By default I want an output
   ! quit if all snapshots were created or if due to any reasons the snapshot counter is 
   ! larger than the declared number of snapshots (e.g., when no snapshots are requested) 
   IF (current_avgsnap.GT.N_of_all_avgsnaps) avg_output_flag = 0
   IF (avg_output_flag==0) RETURN 

   ! quit if the current moment of time is not the moment when it is necessary to create the snapshot
   IF (T_cntr.NE.avgsnapshot(current_avgsnap)%T_cntr_end) avg_output_flag = 0   

END SUBROUTINE DETERMINE_AVG_DATA_CREATION

! --------------------------------------------------------------------------------------------------
!     SUBROUTINE DECIDE_COMPUTE_AVG_DATA
! >    @details Determines if averaged data should be computed. Answer can be no after a restart. 
!               In such a case, first iterations will be ignored until I can properly compute the next avg data array (n+1).
!               That implies that avg data array number n can be not computed at all if it could not complete before the restart.
!               Usually that does not occur because I need the code to crash right after checkpoint but before the current avg data array could be computed. 
!               If I compute one more, I am safe. NOT IDEAL
! !    @authors W. Villafana
! !    @date    Aug-29-2024
! -------------------------------------------------------------------------------------------------- 

SUBROUTINE DECIDE_IF_COMPUTE_AVG_DATA_AFTER_RESTART(avg_compute_flag)

   USE AvgSnapshots, ONLY: current_avgsnap, N_of_all_avgsnaps, avgsnapshot
   USE CurrentProblemValues, ONLY: T_cntr, string_length
   USE mod_print, ONLY: print_debug
   IMPLICIT NONE

   !IN/OUT
   INTEGER, INTENT(OUT) :: avg_compute_flag

   !LOCAL   
   CHARACTER(LEN=string_length) :: routine
   INTEGER :: local_debug_level   

   routine = "DECIDE_IF_COMPUTE_AVG_DATA_AFTER_RESTART"
   local_debug_level = 2

   CALL print_debug( routine,local_debug_level)   

   avg_compute_flag = 0
   
   IF ((current_avgsnap.GE.1).AND.(current_avgsnap.LE.N_of_all_avgsnaps)) THEN
      IF ( (T_cntr.GE.avgsnapshot(current_avgsnap)%T_cntr_begin).AND. &
         & (T_cntr.LE.avgsnapshot(current_avgsnap)%T_cntr_end) ) THEN  
            avg_compute_flag = 1 
      END IF
   END IF

END SUBROUTINE DECIDE_IF_COMPUTE_AVG_DATA_AFTER_RESTART


