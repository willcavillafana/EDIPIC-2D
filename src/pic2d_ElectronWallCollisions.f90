!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  particle_not_processed = .TRUE.

  y = MIN(MAX(y,c_Y_area_min),c_Y_area_max)

  DO n = 1, c_N_of_local_object_parts_left

     m = c_index_of_local_object_part_left(n)

     IF ( (y.GE.c_local_object_part(m)%jstart).AND. &
        & (y.LE.c_local_object_part(m)%jend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        CALL ADD_ELECTRON_TO_BO_COLLS_LIST(REAL(y), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) - dqbelow
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) - dqabove
        END SELECT

        IF (whole_object(nwo)%SEE_enabled) THEN
           CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), m, 1)   ! "1" is for a left wall 
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(x, y, vx, vy, vz, tag, x_old,vx_old,vy_old)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC, i_reflection_cyl_electron,string_length
  USE mod_print, ONLY: print_error

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  REAL(8) x, y, vx, vy, vz
  REAL(8), INTENT(IN), OPTIONAL :: x_old,vx_old,vy_old ! original positoin and velocity in r-theta plane (cylindrical only. We need this to properly perform a specular reflection)
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  ! Local
  REAL(8) :: x_reflected,y_reflected,vx_reflected,vy_reflected ! position and veloicty in local Cartesian frame after specular reflection in Cylindrical
  REAL(8) :: x_cart, y_cart, z_cart ! intermediate cartesian coordinates for cylindrical 
  REAL(8) :: alpha_ang ! increment angle for azimuthal coordinate in cylindrical
  REAL(8) :: radius ! radius angle for intermediate calculation in cylindrical system
  REAL(8) :: vx_temp,vy_temp,vz_temp

  CHARACTER(LEN=string_length) :: routine
  CHARACTER(LEN=string_length) :: message

  routine = 'PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT'

  particle_not_processed = .TRUE.

  y = MIN(MAX(y,c_Y_area_min),c_Y_area_max)

  DO n = 1, c_N_of_local_object_parts_right

     m = c_index_of_local_object_part_right(n)

     IF ( (y.GE.c_local_object_part(m)%jstart).AND. &
        & (y.LE.c_local_object_part(m)%jend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        CALL ADD_ELECTRON_TO_BO_COLLS_LIST(REAL(y), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

         ! By default I do not have a specular reflection in cylindrica coordinates
         IF ( i_reflection_cyl_electron==0 ) THEN

             SELECT CASE (whole_object(nwo)%object_type)
                CASE (VACUUM_GAP)
                CASE (METAL_WALL)
                CASE (DIELECTRIC)
                   ! update the surface charge
                   jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
                   jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
                   dqabove = y - jbelow
                   dqbelow = 1.0_8 - dqabove
                   c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) - dqbelow
                   c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) - dqabove
             END SELECT

             IF (whole_object(nwo)%SEE_enabled) THEN
                CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), m, 3)   ! "3" is for a right wall 
             END IF

       !   ! By default I do not have a specular reflection in cylindrica coordinates. Only for domain boundary
         ELSE IF ( i_reflection_cyl_electron==1 ) THEN
            ! If I need reflection I enter here
            IF (whole_object(nwo)%Elast_refl_type==0 .AND. whole_object(nwo)%Emitted_model(1)==1 ) THEN
               IF ( .NOT. PRESENT(x_old) .OR. .NOT. PRESENT(vx_old) .OR. .NOT. PRESENT(vy_old)) THEN
                  message='missing optional paramaters'
                  CALL print_error(message,routine)
               ELSE IF ( m<=0 ) THEN
                  message='specular reflection is not implemented for inner objects'
                  CALL print_error(message,routine)
               END IF

               ! ! Go backward in time to find old radius
               ! x_cart = x - vx ! radius
               ! z_cart = - vz ! in theta direction
               ! x_old =  SQRT( x_cart**2 + z_cart**2 ) ! old radius      
         
               ! ! Deduce old velocity
               ! alpha_ang = DATAN2(z_cart,x_cart) 
         
               ! ! Former velocity in previous coordinate system
               ! vx_old =   COS(alpha_ang)*vx + SIN(alpha_ang)*vz
               ! vy_old = - SIN(alpha_ang)*vx + COS(alpha_ang)*vz            
               ! IF (x_old>456.8361553) print*,'x_old_test',x_old
               CALL REFLECT_CYLINDRICAL ( x_old,vx_old,vy_old,vz,c_X_area_max,x_reflected,y_reflected,vx_reflected,vy_reflected, 1 )

               ! Adjust position in local Cartesian frame after collision
               x_cart = x_reflected
               z_cart = y_reflected
               radius = SQRT( x_cart**2+z_cart**2 )

               ! Velocity in local cartesian frame
               vx = vx_reflected
               vz = vy_reflected  ! theta direction in z  

               ! Then compute increment angle alpha
               alpha_ang = DATAN2(z_cart,x_cart)         

               ! Update radius (X). 
               ! IF ( radius>c_X_area_max) print*,'x_wil',x
               radius = MIN(c_X_area_max,radius) ! SAftey because of error precision
               x = radius
               
               ! Get Final velocities in cylindrical system. (z speed has already been update above)
               vx_temp =   COS(alpha_ang)*vx + SIN(alpha_ang)*vz
               vz_temp = - SIN(alpha_ang)*vx + COS(alpha_ang)*vz

               vx = vx_temp
               vz = vz_temp        
               
               ! IF (radius>456.956 .AND. radius<456.957) print*,'x_old_alone,y_axial',x_old,y
               ! IF (radius>456.956 .AND. radius<456.957) print*,'x_old,x_cart,z_cart,x',x_old,x_cart,z_cart,x
               ! print*,'vx_old,vy_old,vx,vy',vx_old,vy_old,vx,vy
               ! print*,'ec_before,after',vx_old**2+vy_old**2+vy**2,vx**2+vz**2+vy**2

               CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, whole_object(nwo)%object_id_number)
            ELSE ! This obejct has no reflection in cylindrical coordinates
               IF (whole_object(nwo)%SEE_enabled) THEN
                  CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), -1, 3)
               END IF
            ENDIF          
         END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

  particle_not_processed = .TRUE.

  x = MIN(MAX(x,c_X_area_min),c_X_area_max)

  DO n = 1, c_N_of_local_object_parts_below

     m = c_index_of_local_object_part_below(n)

     IF ( (x.GE.c_local_object_part(m)%istart).AND. &
        & (x.LE.c_local_object_part(m)%iend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        CALL ADD_ELECTRON_TO_BO_COLLS_LIST(REAL(x), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft) - dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) - dqright
        END SELECT

        IF (whole_object(nwo)%SEE_enabled) THEN
           CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), m, 4)   ! "4" is for a wall below
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

  particle_not_processed = .TRUE.

  x = MIN(MAX(x,c_X_area_min),c_X_area_max)

  DO n = 1, c_N_of_local_object_parts_above

     m = c_index_of_local_object_part_above(n)

     IF ( (x.GE.c_local_object_part(m)%istart).AND. &
        & (x.LE.c_local_object_part(m)%iend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        CALL ADD_ELECTRON_TO_BO_COLLS_LIST(REAL(x), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft) - dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) - dqright
        END SELECT

        IF (whole_object(nwo)%SEE_enabled) THEN
           CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), m, 2)   ! "2" is for a wall above
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE

!------------------------------------------------------
!
SUBROUTINE COLLECT_ELECTRON_BOUNDARY_HITS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_spec
  USE AvgSnapshots, ONLY: avg_flux_and_history

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER request

  INTEGER, ALLOCATABLE :: ibuf_send(:)
  INTEGER, ALLOCATABLE :: ibuf_receive(:)
  INTEGER ALLOC_ERR

  INTEGER k
  INTEGER :: local_debug_level
  INTEGER :: avg_compute_flag

  local_debug_level = 2

  ALLOCATE(ibuf_send(1:N_of_boundary_and_inner_objects), STAT = ALLOC_ERR)
  ALLOCATE(ibuf_receive(1:N_of_boundary_and_inner_objects), STAT = ALLOC_ERR)

! each cluster adjacent to a boundary assembles electron-boundary hit counters from all cluster members in the master of the cluster

  ibuf_send(1:N_of_boundary_and_inner_objects) = whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count
  ibuf_receive = 0

  CALL MPI_REDUCE(ibuf_send, ibuf_receive, N_of_boundary_and_inner_objects, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)  !??? use Rank_of_bottom_left_cluster_master ???

  IF (Rank_of_process.EQ.0) THEN
! now counters from all processes are assembled in the process with global rank zero
    
     whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count = ibuf_receive(1:N_of_boundary_and_inner_objects)
     IF (debug_level>=local_debug_level) print '("electrons hit boundaries :: ",10(2x,i8))', whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count  

     ! Update avg flux if necessary
      IF (avg_flux_and_history) THEN
         CALL DECIDE_IF_COMPUTE_AVG_DATA_AFTER_RESTART(avg_compute_flag)
         IF (avg_compute_flag==1) THEN 
            whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_flux_avg_per_s = whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_flux_avg_per_s &
                                                                                          + REAL(whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count)     
         END IF
      END IF
     DO k = 1, N_of_boundary_and_inner_objects
        whole_object(k)%ion_hit_count(1:N_spec) = 0
     END DO
     
  END IF

  DEALLOCATE(ibuf_send, STAT = ALLOC_ERR)
  DEALLOCATE(ibuf_receive, STAT = ALLOC_ERR)

END SUBROUTINE COLLECT_ELECTRON_BOUNDARY_HITS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE INITIATE_WALL_DIAGNOSTICS

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : N_of_boundary_and_inner_objects, Start_T_cntr
  USE Checkpoints, ONLY : use_checkpoint
!  USE Diagnostics, ONLY : N_of_saved_records
  USE SetupValues, ONLY : ht_use_e_emission_from_cathode, ht_use_e_emission_from_cathode_zerogradf, ht_emission_constant
  USE AvgSnapshots, ONLY: avg_flux_and_history
  USE CurrentProblemValues, ONLY: string_length

  IMPLICIT NONE

                                    ! ----x----I----x--
  CHARACTER(17) historybo_filename  ! history_bo_NN.dat
  CHARACTER(LEN=string_length) :: file_name
  LOGICAL exists
  INTEGER i, k
  INTEGER i_dummy, ios

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (Rank_of_process.NE.0) RETURN

  IF (ht_use_e_emission_from_cathode.OR.ht_use_e_emission_from_cathode_zerogradf.OR.ht_emission_constant) RETURN

! hardwired for objects #2 (cathode) and #4 (anode)

  IF (use_checkpoint.EQ.1) THEN
   ! start from checkpoint, must trim the time dependences

      IF (avg_flux_and_history) THEN
         DO k = 1, N_of_boundary_and_inner_objects

            file_name = 'history_bo_avg_NN.dat'
            file_name(16:17) = convert_int_to_txt_string(k, 2)

            INQUIRE (FILE = file_name, EXIST = exists)
            IF (exists) THEN                                                       
               OPEN (21, FILE = file_name, STATUS = 'OLD')          
               DO !i = 1, Start_T_cntr   !N_of_saved_records             ! these files are updated at every electron timestep
                  READ (21, '(2x,i8,2x,ES14.7,10(2x,ES14.7))', iostat = ios) i_dummy
                  IF (ios.NE.0) EXIT
                  IF (i_dummy.GE.Start_T_cntr) EXIT
               END DO
               BACKSPACE(21) ! Backspace is restored. I can safely eliminate last poitn because it should be recomputed as checkpoint automatically restarts at start of average window
               ENDFILE 21       
               CLOSE (21, STATUS = 'KEEP')        
            ELSE! Start a new one. Back compatibility
         
                  OPEN  (21, FILE = file_name, STATUS = 'REPLACE')          
                  CLOSE (21, STATUS = 'KEEP')
         
            END IF

         END DO

      ELSE
         DO k = 1, N_of_boundary_and_inner_objects

            historybo_filename = 'history_bo_NN.dat'
            historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

            INQUIRE (FILE = historybo_filename, EXIST = exists)
            IF (exists) THEN                                                       
               OPEN (21, FILE = historybo_filename, STATUS = 'OLD')          
               DO !i = 1, Start_T_cntr   !N_of_saved_records             ! these files are updated at every electron timestep
                  READ (21, '(2x,i8,10(2x,i8))', iostat = ios) i_dummy
                  IF (ios.NE.0) EXIT
                  IF (i_dummy.GE.Start_T_cntr) EXIT
               END DO
               BACKSPACE(21)
               ENDFILE 21       
               CLOSE (21, STATUS = 'KEEP')        
            ELSE ! Start a new one. Back compatibility
      
                  OPEN  (21, FILE = historybo_filename, STATUS = 'OLD')          
                  CLOSE (21, STATUS = 'KEEP')
      
            END IF

         END DO
      END IF

  ELSE
   ! fresh start, empty files, clean up whatever garbage there might be

      IF (avg_flux_and_history) THEN
         DO k = 1, N_of_boundary_and_inner_objects

            file_name = 'history_bo_avg_NN.dat'
            file_name(16:17) = convert_int_to_txt_string(k, 2)

            OPEN  (21, FILE = file_name, STATUS = 'REPLACE')          
            CLOSE (21, STATUS = 'KEEP')

         END DO
      ELSE

         DO k = 1, N_of_boundary_and_inner_objects

            historybo_filename = 'history_bo_NN.dat'
            historybo_filename(12:13) = convert_int_to_txt_string(k, 2)
   
            OPEN  (21, FILE = historybo_filename, STATUS = 'REPLACE')          
            CLOSE (21, STATUS = 'KEEP')
   
         END DO
      END IF   
  END IF

END SUBROUTINE INITIATE_WALL_DIAGNOSTICS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE SetupValues, ONLY : ht_use_e_emission_from_cathode, ht_use_e_emission_from_cathode_zerogradf, ht_emission_constant
  USE IonParticles, ONLY : N_spec, Qs
  USE ExternalCircuit
  USE AvgSnapshots, ONLY: avg_flux_and_history, avgsnapshot, current_avgsnap
  USE mod_print, ONLY: print_message

  IMPLICIT NONE

  INTEGER nn, noi, s
  INTEGER k
                                    ! ----x----I----x--
  CHARACTER(17) historybo_filename  ! history_bo_NN.dat
  INTEGER :: avg_output_flag
  INTEGER :: N_averaged_timesteps
  REAL(8) :: time_window
  CHARACTER(LEN=string_length) :: message
  CHARACTER(LEN=string_length) :: file_name

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  avg_output_flag = 0
  IF (Rank_of_process.NE.0) RETURN

  IF (ht_use_e_emission_from_cathode.OR.ht_use_e_emission_from_cathode_zerogradf.OR.ht_emission_constant) RETURN

  DO nn = 1, N_of_object_potentials_to_solve
     noi = object_charge_calculation(1)%noi
     dQ_plasma_of_object(nn) = -whole_object(noi)%electron_hit_count + &
                             &  whole_object(noi)%electron_emit_count                      ! include electron emission
     DO s = 1, N_spec
        dQ_plasma_of_object(nn) = dQ_plasma_of_object(nn) + Qs(s) * whole_object(noi)%ion_hit_count(s)
     END DO
  END DO

  IF (avg_flux_and_history) THEN
      CALL DETERMINE_AVG_DATA_CREATION(avg_output_flag)
      IF (avg_output_flag==1) THEN
         N_averaged_timesteps   =  avgsnapshot(current_avgsnap)%T_cntr_end - avgsnapshot(current_avgsnap)%T_cntr_begin + 1
         time_window = REAL(N_averaged_timesteps*delta_t_s)
         WRITE( message,'(A,I4,A)') "### ^^^^^^^^^^^^^^^^^^^^ Averaged flux ",current_avgsnap," will be created ^^^^^^^^^^^^^^^^^^^ ###"
         CALL print_message( message )   
         DO k = 1, N_of_boundary_and_inner_objects

            file_name = 'history_bo_avg_NN.dat'
            file_name(16:17) = convert_int_to_txt_string(k, 2)
            
            WRITE( message,'(A)') "Saving file "//TRIM(file_name)
            CALL print_message( message )   
            
            OPEN (21, FILE = file_name, POSITION = 'APPEND')
            WRITE (21, '(2x,i8,2x,ES14.7,10(2x,ES14.7))') &
                 & T_cntr, &
                 & REAL(T_cntr)*delta_t_s, &
                 & whole_object(k)%electron_hit_flux_avg_per_s*weight_ptcl/time_window , &
                 & whole_object(k)%ion_hit_flux_avg_per_s(1:N_spec)*weight_ptcl/time_window, &
                 & whole_object(k)%electron_emission_flux_avg_per_s*weight_ptcl/time_window
       
            CLOSE (21, STATUS = 'KEEP')

            ! Reset everything
            whole_object(k)%electron_hit_flux_avg_per_s = zero
            whole_object(k)%ion_hit_flux_avg_per_s(1:N_spec) = zero
            whole_object(k)%electron_emission_flux_avg_per_s = zero
         END DO         
         WRITE( message,'(A,I4,A)') "### ^^^^^^^^^^^^^^^^^^^^ Averaged flux ",current_avgsnap," saved ^^^^^^^^^^^^^^^^^^^ ###"
         CALL print_message( message )            
      END IF
      RETURN ! If I need average on the fly data, I am done
  END IF  

  DO k = 1, N_of_boundary_and_inner_objects

     historybo_filename = 'history_bo_NN.dat'
     historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

     OPEN (21, FILE = historybo_filename, POSITION = 'APPEND')
     WRITE (21, '(2x,i8,10(2x,i8))') &
          & T_cntr, &
          & whole_object(k)%electron_hit_count , &
          & whole_object(k)%ion_hit_count(1:N_spec), &
          & whole_object(k)%electron_emit_count

     CLOSE (21, STATUS = 'KEEP')
  
  END DO

END SUBROUTINE SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS

!-----------------------------------------
!
SUBROUTINE TRY_ELECTRON_COLL_WITH_INNER_OBJECT(x, y, vx, vy, vz, tag, n_obj_collision) !, myobject)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues !, ONLY : inner_object, METAL_WALL, DIELECTRICJ_ext

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

!  INTEGER nio  ! number of the inner object
  REAL(8) x, y, vx, vy, vz
  INTEGER tag
  REAL(8) :: r_old, vx_old, vy_old, x_cart, z_cart ! old radius, readial and axial velocity in r-z
  REAL(8) :: alpha_ang
  REAL(8) :: vx_new, vy_new, vz_new
  INTEGER, INTENT(OUT) :: n_obj_collision
  !  TYPE(boundary_object) myobject

  REAL(8) xorg, yorg

  INTEGER n_do     ! number of inner object that the particle collided with
  INTEGER mcross   ! number of the segment of inner object number n_do that the particle collided with
  REAL(8) xcross, ycross  ! coordinates of the crossing
  REAL(8) distorg         ! distance from the origin to the crossing (we keep the crossing with the smallest distance from the origin)

  INTEGER n_try
  INTEGER mcross_try
  REAL(8) xcross_try, ycross_try, distorg_try

  INTEGER coll_direction_flag
  REAL(8) :: t_star, R_max

  REAL coll_coord   ! coordinate of collision point, y/x for collisions with vertical/horizontal segments, respectively

  n_obj_collision = 0
  ! Find previous position
  IF (i_cylindrical==0) THEN
      xorg = x - vx
      yorg = y - vy

  ELSE IF ( i_cylindrical==2 ) THEN
      ! Go backward in time to find old radius
      x_cart = x - vx ! radius
      z_cart = - vz ! in theta direction
      xorg =  SQRT( x_cart**2 + z_cart**2 ) ! old radius      
      yorg = y - vy

      ! Check if old particle position is not outside simulation domain. In such a case I need to correct initial origin point. This isutation can happen when I have cylindrical reflections
      R_max = DBLE(c_indx_x_max_total)
      IF ( i_cylindrical==2 .AND. xorg> R_max ) THEN

         xorg = R_max-1.0D-6

         
         ! Compute actual time you spent from external boundary. Positive solution should be OK
         t_star = (x*vx+SQRT( (R_max*vx)**2 + vz**2*(R_max**2 - x**2) ))/(vx**2+vz**2)

         ! Correct axial position 
         yorg = y - vy*t_star

         ! Correct init positions in local frame
         x_cart = x - vx*t_star
         z_cart = -vz*t_star

         ! print*,'xorg,t_star,t_star_2,x_cart,z_cart',xorg,t_star,(x*vx-SQRT( (xorg*vx)**2 + vz**2*(xorg**2 - x**2) ))/(vx**2+vz**2),x_cart,z_cart
      END IF

      ! Deduce old velocity
      alpha_ang = DATAN2(z_cart,x_cart) 

      ! Former velocity in previous coordinate system
      vx_old =   COS(alpha_ang)*vx + SIN(alpha_ang)*vz
      vy_old = - SIN(alpha_ang)*vx + COS(alpha_ang)*vz
   ENDIF

   ! By default the velocity at impact is the one I have right now (will be different in cylindrical)
   vx_new = vx ! r in cylindrical rz
   vy_new = vy ! z in cylindrical rz
   vz_new = vz ! theta in cylindrical rz

  n_do = -1
  mcross = -1

  DO n_try = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     CALL FIND_CLOSEST_INTERSECTION_WITH_OBJECT(xorg, yorg, x, y, n_try, mcross_try, xcross_try, ycross_try, distorg_try,vx_old,vy_old, vy, vx_new, vz_new)

     IF (mcross_try.LT.0) CYCLE  ! no crossing found

     IF (mcross.EQ.-1) THEN
! the very first crossing was found
        n_do = n_try
        mcross = mcross_try
        xcross = xcross_try
        ycross = ycross_try
        distorg = distorg_try
     ELSE
        IF (distorg_try.GE.distorg) CYCLE
! the new crossing is closer to the origin than the previously found one
        n_do = n_try
        mcross = mcross_try
        xcross = xcross_try
        ycross = ycross_try
        distorg = distorg_try
     END IF

  END DO   !### DO n_try = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

  IF (mcross.EQ.-1) THEN
     print*,'xorg,yorg,xnew,ynew,vx.vy.vz,vx_old,vy_old',xorg,yorg,x,y,vx,vy,vz,vx_old,vy_old
     PRINT '("Error-1 in TRY_ELECTRON_COLL_WITH_INNER_OBJECT. If this is a problem with cylindrical and specular reflection, remove specular walls behind inner object. Top and bottom right corner can sometimes present problems (including leaks). Might try to add some pace between inner objects and these corners as well. ",4(2x,f10.4))', xorg, yorg, x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  SELECT CASE (mcross)
     CASE (1)
        coll_direction_flag = 3
        coll_coord = REAL(ycross)
     CASE (2)
        coll_direction_flag = 4
        coll_coord = REAL(xcross)
     CASE (3)
        coll_direction_flag = 1
        coll_coord = REAL(ycross)
     CASE (4)
        coll_direction_flag = 2
        coll_coord = REAL(xcross)
  END SELECT

  CALL ADD_ELECTRON_TO_BO_COLLS_LIST(coll_coord, REAL(vx_new), REAL(vy_new), REAL(vz_new), tag, n_do, mcross)
  n_obj_collision = n_do

  CALL DO_ELECTRON_COLL_WITH_INNER_OBJECT(xcross, ycross, vx_new, vy_new, vz_new, tag, whole_object(n_do), coll_direction_flag,xorg,vx_old,vy_old)

END SUBROUTINE TRY_ELECTRON_COLL_WITH_INNER_OBJECT

!-----------------------------------------
!
SUBROUTINE DO_ELECTRON_COLL_WITH_INNER_OBJECT(x, y, vx, vy, vz, tag, myobject, coll_direction_flag,x_old,vx_old,vy_old)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues !, ONLY : inner_object, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

!  INTEGER nio  ! number of the inner object
  REAL(8) :: x_old, vx_old, vy_old ! old radius, readial and azimuthal velocity in r-z
  REAL(8) x, y, vx, vy, vz
  INTEGER tag
  TYPE(boundary_object) myobject
  INTEGER coll_direction_flag

  REAL(8) xmin, xmax, ymin, ymax

  INTEGER i_left_top, i_right_top, i_right_bottom, i_left_bottom_bis, i, ip1
  REAL(8) :: x_reflected,y_reflected,vx_reflected,vy_reflected ! position and veloicty in local Cartesian frame after specular reflection in Cylindrical
  REAL(8) :: R_object, alpha_ang, vx_temp, vz_temp, x_cart, z_cart, radius
  REAL(8) dqip1, dqi

! identify side of the inner object hit by the particle

  xmin = myobject%Xmin
  xmax = myobject%Xmax
  ymin = myobject%Ymin
  ymax = myobject%Ymax

  myobject%electron_hit_count = myobject%electron_hit_count + 1

  IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge

! ilt   -  ---- irt
!  |             |
!  |             |
!  1 ilbb ---- irb
!
     i_left_top        =                  myobject%jtop   - myobject%jbottom + 1
     i_right_top       = i_left_top     + myobject%iright - myobject%ileft
     i_right_bottom    = i_right_top    + myobject%jtop   - myobject%jbottom
     i_left_bottom_bis = i_right_bottom + myobject%iright - myobject%ileft - 1

     SELECT CASE (coll_direction_flag)
        CASE (3)
! left wall of the object
           i = MIN(INT(y - ymin) + 1, i_left_top - 1)
           dqip1 = y - INT(y)
           dqi = 1.0_8 - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   - dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) - dqip1
           
        CASE (4)
! top wall of the object
           i = MIN(INT(x - xmin) + i_left_top, i_right_top - 1)
           dqip1 = x - INT(x)
           dqi = 1.0_8 - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   - dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) - dqip1

        CASE (1)
! right wall of the object
           i = MIN(INT(ymax - y) + i_right_top, i_right_bottom - 1)
           dqi = y - INT(y)
           dqip1 = 1.0_8 - dqi
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   - dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) - dqip1

        CASE (2)
! bottom wall of the object
           i = MIN(INT(xmax - x) + i_right_bottom, i_left_bottom_bis)
           dqi = x - INT(x)
           dqip1 = 1.0_8 - dqi
           ip1 = i+1
           IF (i.EQ.i_left_bottom_bis) ip1 = 1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   - dqi
           myobject%surface_charge_variation(ip1) = myobject%surface_charge_variation(ip1) - dqip1

     END SELECT
  END IF   !### IF (myobject%object_type.EQ.DIELECTRIC) THEN
! print*,'hello,i_reflection_cyl_electron,coll_direction_flag',i_reflection_cyl_electron,coll_direction_flag
   IF (i_reflection_cyl_electron==0 ) THEN

      IF (myobject%SEE_enabled) THEN
         CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, myobject, -1, coll_direction_flag)
      END IF
   ! Cylindrical coordinates, I need a special treatment
   ELSE IF ( i_cylindrical==2 ) THEN
      ! If this object has cylindrical reflection, I need to be carfeul
      IF ( myobject%Elast_refl_type==0 .AND. myobject%Emitted_model(1)==1) THEN
         ! ELSE IF ( i_reflection_cyl_electron==1 ) THEN
         ! Left or right
         IF ( coll_direction_flag==1 .OR. coll_direction_flag==3 ) THEN

            ! right
            IF ( coll_direction_flag==1 ) THEN
               R_object = xmax
            ! left
            ELSE IF ( coll_direction_flag==3 ) THEN
               R_object = xmin
            ENDIF
            CALL REFLECT_CYLINDRICAL( x_old,vx_old,vy_old,vz,R_object,x_reflected,y_reflected,vx_reflected,vy_reflected, 1 )

            ! Adjust position in local Cartesian frame after collision
            x_cart = x_reflected
            z_cart = y_reflected
            radius = SQRT( x_cart**2+z_cart**2 )

            ! Velocity in local cartesian frame
            vx = vx_reflected
            vz = vy_reflected  ! theta direction in z  

            ! Then compute increment angle alpha
            alpha_ang = DATAN2(z_cart,x_cart)         

            ! Update radius (X). 
            ! IF ( radius>c_X_area_max) print*,'x_wil',x
            radius = MIN(c_X_area_max,radius) ! SAftey because of error precision
            x = radius
            
            ! Get Final velocities in cylindrical system. (z speed has already been update above)
            vx_temp =   COS(alpha_ang)*vx + SIN(alpha_ang)*vz
            vz_temp = - SIN(alpha_ang)*vx + COS(alpha_ang)*vz

            vx = vx_temp
            vz = vz_temp        
            ! print*,'x_old,x',x_old,x
            ! IF (x_old > 457.01 .AND. x_old<457.15) print*,'x_old,x,wil_1',x_old,x
            ! IF (x > 456.97 .AND. x<456.98) print*,'x_old,x,wil_2',x_old,x
            ! print*,'x_old,x_cart,z_cart,x',x_old,x_cart,z_cart,x
            ! print*,'vx_old,vy_old,vx,vy',vx_old,vy_old,vx,vy
            ! print*,'ec_before,after',vx_old**2+vy_old**2+vy**2,vx**2+vz**2+vy**2

            CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, myobject%object_id_number)
         ! For other directions, I need to check
         ELSE
            IF (myobject%SEE_enabled) THEN
               CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, myobject, -1, coll_direction_flag)
               ! print*,'x_old,x_2',x_old,x
            END IF         
         ENDIF
      ELSE ! This obejct has no reflection in cylindrical coordinates
         IF (myobject%SEE_enabled) THEN
            CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, myobject, -1, coll_direction_flag)
         END IF
      ENDIF         
   ENDIF
         

END SUBROUTINE DO_ELECTRON_COLL_WITH_INNER_OBJECT

!--------------------------
!
SUBROUTINE FIND_CLOSEST_INTERSECTION_WITH_OBJECT(xorg, yorg, x, y, n, mcross, xcross, ycross, distorg, vx_old, vy_old, vz_axial, vx_new, vz_new)

  USE CurrentProblemValues !, ONLY : inner_object, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: xorg, yorg, x, y
  REAL(8), INTENT(IN) :: vx_old, vy_old, vz_axial
  INTEGER, INTENT(IN) :: n

  INTEGER, INTENT(OUT) :: mcross
  REAL(8), INTENT(OUT) :: xcross, ycross
  REAL(8), INTENT(OUT) :: distorg
  REAL(8), INTENT(INOUT) :: vx_new, vz_new ! New velocities at impact with the wall (in cylindrical, it will be different). Will not be updated if Cartesian

  INTEGER m

  INTEGER jbot, jtop, jcross
  INTEGER ileft, iright, icross
  REAL(8) myxcross, myycross
  INTEGER mystatus
  REAL(8) mydistorg

  mcross = -1

  DO m = 1, whole_object(n)%number_of_segments

     IF (whole_object(n)%segment(m)%istart.EQ.whole_object(n)%segment(m)%iend) THEN
! vertical segment
         
        jbot = MIN(whole_object(n)%segment(m)%jstart, whole_object(n)%segment(m)%jend)
        jtop = MAX(whole_object(n)%segment(m)%jstart, whole_object(n)%segment(m)%jend)
        myxcross = DBLE(whole_object(n)%segment(m)%istart)
        CALL CHECK_INTERSECTION_WITH_VERTICAL_SEGMENT( xorg, yorg, x, y, myxcross, DBLE(jbot), DBLE(jtop), mystatus, myycross, vx_old, vy_old, vz_axial, vx_new, vz_new )
        IF (mystatus.NE.0) CYCLE

        jcross = MAX(jbot, MIN(jtop-1, INT(myycross)))

!print *,"aa ", jcross, n, m, whole_object(n)%segment(m)%cell_is_covered(jcross)

! check that the intersection point is not in the prohibited (covered) part of the segment
        IF (whole_object(n)%segment(m)%cell_is_covered(jcross)) CYCLE

        mydistorg = (myxcross-xorg)**2 + (myycross-yorg)**2

        IF (mcross.EQ.-1) THEN
           mcross = m
           xcross = myxcross
           ycross = myycross
           distorg = mydistorg
        ELSE
           IF (mydistorg.GE.distorg) CYCLE
           mcross = m
           xcross = myxcross
           ycross = myycross
           distorg = mydistorg
        END IF

     ELSE IF (whole_object(n)%segment(m)%jstart.EQ.whole_object(n)%segment(m)%jend) THEN
! horizontal segment
      
        ileft  = MIN(whole_object(n)%segment(m)%istart, whole_object(n)%segment(m)%iend)
        iright = MAX(whole_object(n)%segment(m)%istart, whole_object(n)%segment(m)%iend)
        myycross = DBLE(whole_object(n)%segment(m)%jstart)

        CALL CHECK_INTERSECTION_WITH_HORIZONTAL_SEGMENT( xorg, yorg, x, y, myycross, DBLE(ileft), DBLE(iright), mystatus, myxcross, vx_old, vy_old, vz_axial, vx_new, vz_new )
        IF (mystatus.NE.0) CYCLE

        icross = MAX(ileft, MIN(iright-1, INT(myxcross)))

!print *,"bb ", icross, n, m, whole_object(n)%segment(m)%cell_is_covered(icross)

! check that the intersection point is not in the prohibited (covered) part of the segment
        IF (whole_object(n)%segment(m)%cell_is_covered(icross)) CYCLE

        mydistorg = (myxcross-xorg)**2 + (myycross-yorg)**2

        IF (mcross.EQ.-1) THEN
           mcross = m
           xcross = myxcross
           ycross = myycross
           distorg = mydistorg
        ELSE
           IF (mydistorg.GE.distorg) CYCLE
           mcross = m
           xcross = myxcross
           ycross = myycross
           distorg = mydistorg
        END IF

     END IF  !### ELSE IF (whole_object(n)%segment(m)%jstart.EQ.whole_object(n)%segment(m)%jend) THEN

  END DO   !### DO m = 1, whole_object(n)%number_of_segments

END SUBROUTINE FIND_CLOSEST_INTERSECTION_WITH_OBJECT

!----------------------------------------------------------
!### assume that yminseg < ymaxseg
!
SUBROUTINE CHECK_INTERSECTION_WITH_VERTICAL_SEGMENT( xorg, yorg, x, y, xseg, yminseg, ymaxseg, mystatus, ycross, vx_old, vy_old, vz_axial, vx_new, vz_new )

  use, intrinsic :: ieee_arithmetic
  USE CurrentProblemValues, ONLY: i_cylindrical, four, one, two, zero

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: xorg, yorg, x, y          ! coordinates of the ends of particle trajectory segment
  REAL(8), INTENT(IN) :: xseg, yminseg, ymaxseg    ! coordinates of the ends of vertical boundary segment
  REAL(8), INTENT(IN) :: vx_old, vy_old, vz_axial        ! old vr and vtheta speed with axial speed vz (not rotated so I keep that)
  INTEGER, INTENT(OUT) :: mystatus                 ! zero if crossing found, nonzero otherwise
  REAL(8), INTENT(OUT) :: ycross                   ! y-coordinate of the crossing
  REAL(8), INTENT(INOUT) :: vx_new, vz_new   ! Velocities at time of impact in cylindrical

  ! IN/OUT
  REAL(8) :: a,b,delta ! trajectory parameters and intermediate computation variables
  REAL(8) :: x_star,y_star ! Intersection point with outer radius
  REAl(8) :: t_star, x_start
  REAL(8) :: alpha_ang
  REAL(8) :: sol_sign

  mystatus = -1

! check obvious things first
  IF (MAX(xorg, x).LT.xseg) RETURN
  IF (MIN(xorg, x).GT.xseg) RETURN
  IF (MAX(yorg, y).LT.yminseg) RETURN
  IF (MIN(yorg, y).GT.ymaxseg) RETURN
  !IF (xorg>456.999) print*,'xorg_vert,yorg',xorg,yorg
! extremely unlikely situation, particle goes exactly along the surface of the object
  IF (xorg.EQ.x) RETURN

! since we are here, ends of segment {xorg,yorg}-{x,y} are on different sides of segment {xseg,yminseg}-{xseg,ymaxseg}
   IF ( i_cylindrical==0 ) THEN
      ycross = yorg + (y-yorg) * (xseg-xorg) / (x-xorg)
   ELSE IF ( i_cylindrical==2 ) THEN
      ! Step 1: compute coefficients of straight line if I had no impact, y= ax+b
      ! Compute initial trajectory
      x_start = xorg

      !!! Need to differentiate some cases 
      ! Vertical line (improbable case)
      IF ( vx_old==zero ) THEN 
          x_star = x_start
          y_star = SIGN(one,vy_old)*SQRT(xseg**2-x_start**2)
          t_star = y_star/vy_old
      ELSE ! Line with finite slope  
         sol_sign = SIGN(one,vx_old) ! By default we assume a wall on the right. The largest solution should be picked if we go to thr right in the circle. Otherwise it should be the smallest one.
         IF (xseg<x_start) sol_sign = one ! if this is a wall on the left, then the largest solution of quadratic equation should be picked 

          a = vy_old/vx_old
          b = -a*x_start
       
          ! Step 2: get coordinates of interection point with radius
          delta = four*(a**2*(xseg**2-x_start**2)+xseg**2)
          !delta = a**2*(xseg-x_start)*(xseg+x_start)+xseg**2
          x_star = (-two*a*b+sol_sign*SQRT(delta))/(two*(a**2+one)) ! Correct solution is autmoatically chosen 
          !x_star = (-a*b+SQRT(delta))/(a**2+one) ! Take positive solution (should be OK)
          y_star = a*x_star+b   
          
          ! Step 3: deduce spent time
          t_star = (x_star-x_start)/vx_old
          IF (t_star<zero) print*,'Time is negative in collision with vertical wall: t_star,xseg,x_star,y_star,x_start,vx_old,vy_old',t_star,xseg,x_star,y_star,x_start,vx_old,vy_old

      END IF
      ! Step 4: compute axial location of intersection
      ycross = vz_axial*t_star+yorg

      ! Step 5: Deduce final velocity at time of impact
      alpha_ang = DATAN2(y_star,x_star)
      vx_new =  COS(alpha_ang)*vx_old + SIN(alpha_ang)*vy_old
      vz_new = -SIN(alpha_ang)*vx_old + COS(alpha_ang)*vy_old ! theta
      !  IF (xorg > 1024.1077 .AND. xorg<1024.1079) print*,'xorg,y_star,tnew,vy_old,ycross,t_star,yorg,vz_axial,x_star,x_start',xorg,y_star,(y_star)/vy_old,vy_old,ycross,t_star,yorg,vz_axial,x_star,x_start
   END IF

  IF (.NOT.ieee_is_finite(ycross)) THEN
     mystatus = 1
     RETURN
  END IF
!   IF (xorg > 1024.1077 .AND. xorg<1024.1079) print*,'xorg,ycross,yminseg,ymaxseg,yorg,y',xorg,ycross,yminseg,ymaxseg,yorg,y
  IF (ycross.LT.yminseg) RETURN
  IF (ycross.GT.ymaxseg) RETURN
  IF (ycross.LT.MIN(y,yorg)) RETURN  ! paranoidal failsafe check?
  IF (ycross.GT.MAX(y,yorg)) RETURN  !

  mystatus = 0
  RETURN

END SUBROUTINE CHECK_INTERSECTION_WITH_VERTICAL_SEGMENT

!----------------------------------------------------------
!### assume that xminseg < xmagseg
!
SUBROUTINE CHECK_INTERSECTION_WITH_HORIZONTAL_SEGMENT( xorg, yorg, x, y, yseg, xminseg, xmaxseg, mystatus, xcross, vx_old, vy_old, vz_axial, vx_new, vz_new )

  use, intrinsic :: ieee_arithmetic
  USE CurrentProblemValues, ONLY: i_cylindrical, four, one, two

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: xorg, yorg, x, y          ! coordinates of the ends of particle trajectory segment
  REAL(8), INTENT(IN) :: yseg, xminseg, xmaxseg    ! coordinates of the ends of horizontal boundary segment
  REAL(8), INTENT(IN) :: vx_old, vy_old, vz_axial        ! old vr and vtheta speed with axial speed vz (not rotated so I keep that)
  INTEGER, INTENT(OUT) :: mystatus                 ! zero if crossing found, nonzero otherwise
  REAL(8), INTENT(OUT) :: xcross                   ! x-coordinate of the crossing
  REAL(8), INTENT(INOUT) :: vx_new, vz_new   ! Velocities at time of impact in cylindrical

  ! IN/OUT
  REAL(8) :: a,b,delta ! trajectory parameters and intermediate computation variables
  REAL(8) :: x_star,y_star ! Intersection point with outer radius
  REAl(8) :: t_star, x_start
  REAL(8) :: alpha_ang

  mystatus = -1

! check obvious things first
  IF (MAX(yorg, y).LT.yseg) RETURN
  IF (MIN(yorg, y).GT.yseg) RETURN
  IF (MAX(xorg, x).LT.xminseg) RETURN
  IF (MIN(xorg, x).GT.xmaxseg) RETURN
!   IF (xorg>456.999) print*,'xorg_horizontal,yorg',xorg,yorg
! extremely unlikely situation, particle goes exactly along the surface of the object
  IF (yorg.EQ.y) RETURN

! since we are here, ends of segment {xorg,yorg}-{x,y} are on different sides of segment {xseg,yminseg}-{xseg,ymaxseg}
  IF ( i_cylindrical==0 ) THEN
      xcross = xorg + (x-xorg) * (yseg-yorg) / (y-yorg)
   ELSE IF ( i_cylindrical==2 ) THEN
      ! Step 1: deduce how long I need to reach the horizontal wall
      t_star = (yseg-yorg)/vz_axial ! vz cannot be zero here (or we could not reach the wall 
      ! IF (xorg > 329.7661 .AND. xorg<329.7663) t_star = one
      ! Step 2: Compute new x and y coordinates in the r-theta plane
      x_start = xorg
      x_star = x_start+vx_old*t_star
      y_star = vy_old*t_star

      ! Step 3: Compute new radius when the particle hits the wall
      xcross = SQRT(x_star**2+y_star**2)

      ! Step 4 Compute velocities vr and vtheta at time of imnpact
      alpha_ang = DATAN2(y_star,x_star)
      vx_new =  COS(alpha_ang)*vx_old + SIN(alpha_ang)*vy_old
      vz_new = -SIN(alpha_ang)*vx_old + COS(alpha_ang)*vy_old ! theta
      ! IF (xorg > 456 .AND. xorg<457.1) print*,'xorg,xcross,x,xminseg,xmaxseg',xorg,xcross,x,xminseg,xmaxseg
   END IF

  IF (.NOT.ieee_is_finite(xcross)) THEN
     mystatus = 1
     RETURN
  END IF

  IF (xcross.LT.xminseg) RETURN
  IF (xcross.GT.xmaxseg) RETURN
  IF (xcross.LT.MIN(x,xorg) .AND.i_cylindrical==0 ) RETURN  ! paranoidal failsafe check? In cylindrical I can have this situation because the radius is not a linear function.
  IF (xcross.GT.MAX(x,xorg).AND.i_cylindrical==0 ) RETURN  !

  mystatus = 0
  RETURN

END SUBROUTINE CHECK_INTERSECTION_WITH_HORIZONTAL_SEGMENT

!-------------------------------------------------------------------------------------------
! Prepares the tabulated values of integral of the maxwell distribution function
!  
SUBROUTINE PrepareMaxwellDistribIntegral

  USE ParallelOperationValues
  USE MaxwellVelocity
!  USE CurrentProblemValues, ONLY : N_box_vel
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER i
  INTEGER N_pnts
  INTEGER count
  REAL(8) V_min, V_max
  REAL(8) F(0:180003)      ! to be sure that we overcome V_max
  REAL(8) temp
  REAL(8) dV

  LOGICAL check1, check2

  check1 = .FALSE.
  check2 = .FALSE.

! ------- for symmetrical maxwellian
  N_pnts = 180000  !30000
  V_min  = -U_max
  V_max  =  U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     F(i) = F(i-1) + EXP( - (V_min + (DBLE(i)-0.5_8) * dV)**2 )
  END DO

  temp = F(N_pnts)
  F = F * R_max / temp   ! normalize integral such that F(N_pnts) = R_max

  v(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v(count) = V_min + i * dV
        IF (count.EQ.R_max) THEN
           check1 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

!  v = v * N_box_vel

!--------- for asymmetrical maxwellian * v (used for injection, v > 0)

  N_pnts = 90000   !15000
  V_min  = 0.0_8
  V_max  = U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     temp = V_min + (REAL(i)-0.5_8) * dV
     F(i) = F(i-1) + EXP( - temp**2 ) * temp
  END DO

  temp = F(N_pnts)
  F(1:(N_pnts+3)) = F(1:(N_pnts+3)) * R_max_inj / temp   ! normalize integral such that F(N_pnts) = R_max_inj

  v_inj(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v_inj(count) = V_min + i * dV
        IF (count.EQ.R_max_inj) THEN
           check2 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

!  v_inj = v_inj * N_box_vel

  IF (check1.AND.check2) THEN
!     PRINT '(2x,"Process ",i3," : Integrals for producing maxwell distributions are successfully obtained ...")', &
!                                                                                                  & Rank_of_process
  ELSE
     PRINT '(2x,"Process ",i3," : ERROR in PrepareMaxwellDistribIntegral !!!")', Rank_of_process
     PRINT '(2x,"The initialization in PrepareMaxwellDistribIntegral is not performed !!!")'
     PRINT '(2x,"The program will be terminated now :(")'
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE PrepareMaxwellDistribIntegral

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetInjMaxwellVelocity(U) 

  USE MaxwellVelocity
 
  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) U

  REAL(8) R
  INTEGER indx
  
  R = R_max_inj * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max_inj) THEN
     U = v_inj(indx) + (R - indx) * (v_inj(indx+1) - v_inj(indx))
  ELSE
     U = v_inj(R_max_inj)
  END IF
  RETURN
  
END SUBROUTINE GetInjMaxwellVelocity

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetMaxwellVelocity(U) 

  USE MaxwellVelocity

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) U

  REAL(8) R
  INTEGER indx
  
  R = R_max * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max) THEN
     U = v(indx) + (R - indx) * (v(indx+1) - v(indx))
  ELSE
     U = v(R_max)
  END IF
  RETURN
  
END SUBROUTINE GetMaxwellVelocity

!-------------------------
!
REAL(8) FUNCTION vector_product_z(ax, ay, bx, by)

  IMPLICIT NONE
  REAL(8) ax, ay, bx, by

  vector_product_z = ax * by - ay * bx

END FUNCTION vector_product_z

!---------------------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_BO_COLLS_LIST(coll_coord, vx, vy, vz, tag, nwo, nseg)

  USE CurrentProblemValues, ONLY : e_colls_with_bo
  USE Snapshots

  IMPLICIT NONE

  REAL coll_coord, vx, vy, vz
  INTEGER tag
  
  INTEGER nwo   ! number of the boundary object
  INTEGER nseg  ! number of the segment of the boundary object
  
  TYPE collided_particle
     INTEGER token
     REAL coll_coord
     REAL VX
     REAL VY
     REAL VZ
  END TYPE collided_particle

  TYPE(collided_particle), ALLOCATABLE :: bufer(:)
  INTEGER ALLOC_ERR

  INTEGER k
  INTEGER current_N

  IF (.NOT.e_colls_with_bo(nwo)%must_be_saved) RETURN

  IF (current_snap.GT.N_of_all_snaps) RETURN

  IF (.NOT.save_e_collided_with_bo(current_snap)) RETURN

  e_colls_with_bo(nwo)%N_of_saved_parts = e_colls_with_bo(nwo)%N_of_saved_parts+1

  IF (e_colls_with_bo(nwo)%N_of_saved_parts.GT.e_colls_with_bo(nwo)%max_N_of_saved_parts) THEN
! increase the size of the array
     current_N = e_colls_with_bo(nwo)%max_N_of_saved_parts
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%token      = e_colls_with_bo(nwo)%part(k)%token
        bufer(k)%coll_coord = e_colls_with_bo(nwo)%part(k)%coll_coord
        bufer(k)%VX         = e_colls_with_bo(nwo)%part(k)%VX
        bufer(k)%VY         = e_colls_with_bo(nwo)%part(k)%VY
        bufer(k)%VZ         = e_colls_with_bo(nwo)%part(k)%VZ
     END DO
     IF (ALLOCATED(e_colls_with_bo(nwo)%part)) DEALLOCATE(e_colls_with_bo(nwo)%part, STAT=ALLOC_ERR)
     e_colls_with_bo(nwo)%max_N_of_saved_parts = e_colls_with_bo(nwo)%max_N_of_saved_parts + MAX(50, e_colls_with_bo(nwo)%max_N_of_saved_parts/10)
     ALLOCATE(e_colls_with_bo(nwo)%part(1:e_colls_with_bo(nwo)%max_N_of_saved_parts), STAT=ALLOC_ERR)
     DO k = 1, current_N
        e_colls_with_bo(nwo)%part(k)%token      = bufer(k)%token
        e_colls_with_bo(nwo)%part(k)%coll_coord = bufer(k)%coll_coord
        e_colls_with_bo(nwo)%part(k)%VX         = bufer(k)%VX
        e_colls_with_bo(nwo)%part(k)%VY         = bufer(k)%VY
        e_colls_with_bo(nwo)%part(k)%VZ         = bufer(k)%VZ
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=ALLOC_ERR)
  END IF

  k = e_colls_with_bo(nwo)%N_of_saved_parts

  tag = tag - (tag/10000)*10000   ! to clear tracing increment if there is any

  e_colls_with_bo(nwo)%part(k)%token = tag + 100 * nseg

  e_colls_with_bo(nwo)%part(k)%coll_coord = coll_coord
  e_colls_with_bo(nwo)%part(k)%VX = vx
  e_colls_with_bo(nwo)%part(k)%VY = vy
  e_colls_with_bo(nwo)%part(k)%VZ = vz

END SUBROUTINE ADD_ELECTRON_TO_BO_COLLS_LIST

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE REFLECT_CYLINDRICAL
!>    @details Compute new position and speed of specularly reflected particles
!!    @authors W. Villafana
!!    @date    Dec-19-2022
!-------------------------------------------------------------------------------------------------- 

SUBROUTINE REFLECT_CYLINDRICAL ( x_start,vx,vy,vz_axial,R_max ,xf,yf,dot_prod_i,dot_prod_j,n_subcycles)
   
   USE mod_print, ONLY: print_debug
   USE CurrentProblemValues, ONLY: string_length,two,four,one, zero
   USE ParallelOperationValues, ONLY: Rank_of_process
   IMPLICIT NONE

   !IN/OUT
   REAL(8), INTENT(IN) :: x_start ! x start in current iteration
   REAL(8), INTENT(IN) :: vx ! radial normalized velocity before collision. It is the increment in radial direction
   REAL(8), INTENT(IN) :: vy ! azimuthal normalized before collision. It is the increment in azimuthal direction
   REAL(8), INTENT(IN) :: vz_axial ! axial velocity normalized. It is the increment in axial direction
   REAL(8), INTENT(IN) :: R_max ! normalized radius at which we have a reflection
   INTEGER, INTENT(IN) :: n_subcycles ! subcycle frequency (for ions only)
   REAL(8), INTENT(OUT) :: dot_prod_i,dot_prod_j,xf,yf ! final position in local Cartesian frame  

   !LOCAL
   REAL(8) :: a,b,delta ! trajectory parameters and intermediate computation variables
   REAL(8) :: x_star,y_star ! Intersection point with outer radius
   REAL(8) :: d_star,d_remaining ! travelled and remaining distance after collision
   REAL(8) :: beta,v_perp,v_par ! polar angle at intersection point and velocity at intersection point
   REAL(8) :: v_perp_new ! new velocity after collision
   REAL(8) :: module_velocity,dt_remaining ! velocity module and remaining time 
   REAL(8) :: n_sub ! subcycle frequency for ions (will be one for electrons)
   REAL(8) :: t_star! amount of time I need to hit the wall
   REAl(8) :: sol_sign ! Sign to pick correct solution

   CHARACTER(LEN=string_length) :: routine
   INTEGER :: local_debug_level

   routine = "REFLECT_CYLINDRICAL"
   local_debug_level = 2

   CALL print_debug( routine,local_debug_level)
   
   ! Subcycle freq if needed
   n_sub = REAL(n_subcycles) ! this is for ions normally 

   ! Step 1: compute coefficients of straight line if I had no refelction, y= ax+b
   ! Compute initial trajectory
   !!! Need to differentiate some cases
   ! Vertical line (improbable case)
   IF ( vx==zero ) THEN
       x_star = x_start
       y_star = SIGN(one,vy)*SQRT(R_max**2-x_start**2)
       t_star = y_star/vy  
   ELSE ! Line with finite slope
      sol_sign = SIGN(one,vx) ! By default we assume a wall on the right. The largest solution should be picked if we go to thr right in the circle. Otherwise it should be the smallest one. 
      IF (R_max<x_start) sol_sign = one ! if this is a wall on the left, then the largest solution of quadratic equation should be picked 
      ! sol_sign = one !SIGN(one,vx) ! This will determine which solution of the quadratice equation we retain    
       a = vy/vx
       b = -a*x_start

       ! Step 2: get coordinates of interection point with outer radius
       delta = four*(a**2*(R_max**2-x_start**2)+R_max**2)
       x_star = (-two*a*b+sol_sign*SQRT(delta))/(two*(a**2+one)) ! Take correcvt solution depending on sign
       y_star = a*x_star+b  
   

       ! Step 3: compute remaining distance after collision
       t_star = (x_star-x_start)/(vx*n_sub) ! Portion of time I have used
   END IF

   ! IF (x_start>456.8361553 ) print*,'x_start,x_star,y_star,xneg,R,t_star',x_start,x_star,y_star,(-two*a*b-SQRT(delta))/(two*(a**2+one)),SQRT(x_star**2+y_star**2),t_star
   d_star = SQRT((x_star-x_start)**2+y_star**2+(vz_axial*n_sub*t_star)**2)
   d_remaining = SQRT((vx*n_sub)**2+(vy*n_sub)**2+(vz_axial*n_sub)**2) - d_star  
   
   ! Step 4: compute orthogonal and parallel velocity at intersection point, with respect to tangent  
   beta = DATAN2(y_star,x_star) ! polar angle at intersection point
   v_perp = vx*COS(beta)+vy*SIN(beta) !radial
   v_par = -vx*SIN(beta)+vy*COS(beta) ! azimuthal   

   ! Step 5: Define new velocity vector after collision. In local polar coordinate from intersection point
   v_perp_new = - v_perp ! v_par is unchanged   

   ! Step 6: Compute equivalent time to spend of the remaining distance
   module_velocity = SQRT(v_perp**2+v_par**2+(vz_axial)**2)
   dt_remaining = d_remaining/module_velocity   

   ! Step 7: compute back velocity in current Cartesian frame
   dot_prod_i = v_perp_new*COS(beta)-v_par*SIN(beta) ! dot product v_prime with i vector (x)
   dot_prod_j = v_perp_new*SIN(beta)+v_par*COS(beta) ! dot product v_prime with j vector (y)
   xf = x_star+dot_prod_i*dt_remaining*n_sub
   yf = y_star + dot_prod_j*dt_remaining*n_sub   


   IF ( (SQRT(xf**2+yf**2)>R_max .AND. x_start<R_max) .OR. (SQRT(xf**2+yf**2)<R_max .AND. x_start>R_max) ) THEN
      print*,'R_max',R_max
      print*,'t_star',t_star
      print*,'x_start,',x_start
      print*,'vx,vy',vx,vy
      print*,'xf,yf,',xf,yf
      print*,'dot_prod_i,dot_prod_j',dot_prod_i,dot_prod_j
      print*,'beta,',beta
      print*,'v_par',v_par
      print*,'v_perp',v_perp
      print*,'vz_axial',vz_axial
      print*,'d_remaining',d_remaining
      print*,'module_velocity',module_velocity
      print*,'d_star',d_star
      print*,'xstar,dt_remaining',x_star,dt_remaining
      print*,'n_sub',n_sub


   ENDIF

END SUBROUTINE        
