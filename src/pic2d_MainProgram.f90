!===========================================
PROGRAM MainProg

  USE CurrentProblemValues
  USE ParallelOperationValues
  USE LoadBalancing
  USE ClusterAndItsBoundaries
  USE Checkpoints
  USE mod_def_timers
  USE mod_timers, ONLY: start_timer, end_timer, print_timer, pic_loop_timer
  USE mod_print, ONLY: print_message, print_git_info
  USE IonParticles, ONLY: i_freeze_ions

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  REAL(8) start, finish
  INTEGER n_sub
  LOGICAL ions_moved
  INTEGER :: local_debug_level

  ! Timer
  REAL(8)                      :: loop_time,time_measured
  CHARACTER(LEN=string_length) :: message

  CHARACTER(54) rmandmkdir_command    ! rm -rfv checkdir_TTTTTTTT ; mkdir -v checkdir_TTTTTTTT
                                      ! ----x----I----x----I----x----I----x----I----x----I----

  !REAL(8) t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, Rank_of_process, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_of_processes, ierr)

  ! Print current version of the code 
  CALL print_git_info

  CALL SET_PHYSICAL_CONSTANTS

  CALL PrepareMaxwellDistribIntegral

  local_debug_level = 3
  Start_T_cntr = 0
  T_cntr_global_load_balance = Start_T_cntr
  T_cntr_cluster_load_balance = Start_T_cntr

  CALL INITIATE_PARAMETERS
!print *, "did INITIATE_PARAMETERS"

  CALL ANALYZE_BOUNDARY_OBJECTS

  CALL CALCULATE_OBJECT_POTENTIALS_2D

  CALL CALCULATE_OBJECT_POTENTIAL_CHARGE_COEFFS

  CALL INITIATE_ELECTRON_NEUTRAL_COLLISIONS

  CALL INITIATE_ION_NEUTRAL_COLLISIONS

!print *, "did INITIATE_ELECTRON_NEUTRAL_COLLISIONS"
!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
!CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

  CALL INITIATE_PROBE_DIAGNOSTICS

!###  CALL INITIATE_WALL_DIAGNOSTICS_HT_SETUP   ! only one of the two actually works
  CALL INITIATE_WALL_DIAGNOSTICS            !

  CALL INITIATE_GENERAL_DIAGNOSTICS

  CALL INITIATE_EXT_CIRCUIT_DIAGNOSTICS

  CALL INITIATE_en_COLL_DIAGNOSTICS

  CALL INITIATE_in_COLL_DIAGNOSTICS

  CALL INITIATE_AVERAGED_SNAPSHOTS

  CALL INITIATE_SNAPSHOTS

  CALL ADJUST_T_CNTR_SAVE_CHECKPOINT

  IF ((.NOT.use_mpiio_checkpoint).AND.(Rank_of_process.EQ.0).AND.(T_cntr_save_checkpoint.GE.Start_T_cntr)) THEN
! create a directory for a future checkpoint to make sure that the directory exists before the checkpoint files are saved
     rmandmkdir_command = 'rm -rfv checkdir_TTTTTTTT ; mkdir -v checkdir_TTTTTTTT'
                         ! ----x----I----x----I----x----I----x----I----x----I----
     rmandmkdir_command(18:25) = convert_int_to_txt_string(T_cntr_save_checkpoint, 8)
     rmandmkdir_command(47:54) = convert_int_to_txt_string(T_cntr_save_checkpoint, 8)
     CALL SYSTEM(rmandmkdir_command)
  END IF

  start = MPI_WTIME()

  n_sub = 0
  ions_moved = .TRUE.

  CALL start_timer( total_timer )
  CALL ADJUST_HARMONIC_OSCILLATIONS_PHASE

  CALL ADJUST_ECPS_HARMONIC_OSCILLATIONS_PHASE

  DO T_cntr = Start_T_cntr, Max_T_cntr

      CALL start_timer( single_pic_loop_timer )
   !   if (Rank_of_process.eq.0) print *, "Iteration, time (ns) # ",T_cntr,T_cntr*delta_t_s*1e9

     CALL start_timer( save_checkpoint_timer )
     !t0 = MPI_WTIME()
     IF (T_cntr.EQ.T_cntr_save_checkpoint) THEN
        IF (use_mpiio_checkpoint) THEN
           CALL SAVE_CHECKPOINT_MPIIO_2(n_sub)
        ELSE
           CALL SAVE_CHECKPOINT_POSIX(n_sub)
        END IF
        T_cntr_save_checkpoint = T_cntr_save_checkpoint + dT_save_checkpoint
        CALL ADJUST_T_CNTR_SAVE_CHECKPOINT
        IF ((.NOT.use_mpiio_checkpoint).AND.(Rank_of_process.EQ.0)) THEN
! create a directory for a future checkpoint to make sure that the directory exists before the checkpoint files are saved
           rmandmkdir_command = 'rm -rfv checkdir_TTTTTTTT ; mkdir -v checkdir_TTTTTTTT'
                               ! ----x----I----x----I----x----I----x----I----x----I----
           rmandmkdir_command(18:25) = convert_int_to_txt_string(T_cntr_save_checkpoint, 8)
           rmandmkdir_command(47:54) = convert_int_to_txt_string(T_cntr_save_checkpoint, 8)
           CALL SYSTEM(rmandmkdir_command)
        END IF
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     CALL end_timer( save_checkpoint_timer )
     CALL start_timer( global_load_balance_timer )
     !t1 = MPI_WTIME()
     call report_total_number_of_particles
     IF (T_cntr.EQ.T_cntr_global_load_balance) THEN
        IF (n_sub.NE.0 .AND. i_freeze_ions==0) THEN
           PRINT '("Process ",i5," :: ERROR-1 in MainProg :: GLOBAL_LOAD_BALANCE is about to be called at wrong time :: T_cntr = ",i8," n_sub = ",i8)', Rank_of_process, T_cntr, n_sub
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        END IF
      !   IF (i_cylindrical/=0 .OR. Delta_r/=one .OR. Delta_z/=one ) THEN!IF (T_cntr==Start_T_cntr .AND. i_cylindrical/=0 ) THEN
      !    ! WRITE( message,'(A)') achar(10)//"No global load balancing in cylindrical coordinates for first iteration. Avoids bad redistributions if arrays have not been initialized"
      !    WRITE( message,'(A)') achar(10)//"No global load balancing in cylindrical coordinates or partial cartesian. I need to fix this at some point"
      !    CALL print_message( message )
      !   ELSE
         CALL GLOBAL_LOAD_BALANCE  ! includes calls to SET_COMMUNICATIONS 
                                  !                   DISTRIBUTE_CLUSTER_PARAMETERS
      !   ENDIF
        T_cntr_global_load_balance = T_cntr_global_load_balance + dT_global_load_balance
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL end_timer( global_load_balance_timer )
     CALL start_timer( internal_load_balance_timer )
     !t2 = MPI_WTIME()

     IF (T_cntr.EQ.T_cntr_cluster_load_balance) THEN
        CALL BALANCE_LOAD_WITHIN_CLUSTER
        T_cntr_cluster_load_balance = T_cntr_cluster_load_balance + dT_cluster_load_balance
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL end_timer( internal_load_balance_timer )
     CALL start_timer( gather_ion_charge_density_timer )    
     !t3 = MPI_WTIME()

     IF (n_sub.EQ.0) CALL GATHER_ION_CHARGE_DENSITY 

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL end_timer( gather_ion_charge_density_timer )
     CALL start_timer( gather_electron_charge_density_timer )    
     !t4 = MPI_WTIME()

     CALL GATHER_ELECTRON_CHARGE_DENSITY      ! here surface charge density on inner dielectric objects is subtracted from the electron volume charge density

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!###     CALL PERFORM_ELECTRON_EMISSION_HT_SETUP_ZERO_GRAD_F   ! this procedure is used only when axial-azimuthal periodic model of a Hall thruster is simulated
!###     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL UPDATE_WALL_POTENTIALS(T_cntr)

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL end_timer( gather_electron_charge_density_timer )
     CALL start_timer( poisson_solver_timer )    
     !t5 = MPI_WTIME()

     IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

        CALL SOLVE_POTENTIAL_WITH_PETSC
        
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL end_timer( poisson_solver_timer )
        CALL start_timer( calculate_electric_field_timer )    
        !t6 = MPI_WTIME()

        CALL start_timer( externa_circuit_timer )    
        CALL SOLVE_EXTERNAL_CONTOUR
        CALL end_timer( externa_circuit_timer )    

        CALL CALCULATE_ELECTRIC_FIELD               ! n

     ELSE IF (periodicity_flag.EQ.PERIODICITY_X) THEN

        CALL SOLVE_POISSON_FFTX_LINSYSY

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL end_timer( poisson_solver_timer )
        CALL start_timer( calculate_electric_field_timer )    
        !t6 = MPI_WTIME()

        CALL CALCULATE_ELECTRIC_FIELD_FFTX_LINSYSY
     
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     
     CALL end_timer( calculate_electric_field_timer )
     CALL start_timer( compute_averaged_snapshot_timer )    
     !t7 = MPI_WTIME()

     CALL COLLECT_F_EX_EY_FOR_AVERAGED_SNAPSHOT

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL end_timer( compute_averaged_snapshot_timer )
     CALL start_timer( create_instantaneous_snapshot_timer )    
     !t8 = MPI_WTIME()

     CALL CREATE_SNAPSHOT

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL end_timer( create_instantaneous_snapshot_timer )
     CALL start_timer( electron_pusher_with_collisions_inner_object_timer )    
     !t9 = MPI_WTIME()

     CALL ADVANCE_ELECTRONS_PLUS                      !   velocity: n-1/2 ---> n+1/2
                                                      ! coordinate: n     ---> n+1 
                                                      ! may calculate electron moments for averaged snapshots

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL SAVE_ELECTRONS_COLLIDED_WITH_BOUNDARY_OBJECTS

     IF ((n_sub+1).NE.N_subcycles) CALL FIND_ALIENS_IN_ELECTRON_ADD_LIST        ! when n_sub+1==N_subcycles, the ions will be advanced below
                                                                                ! there may be more electrons in the electron_to_add array due to ion-induced SEE
                                                                                ! so at this timestep we call this procedure later, inside the ion IF clause

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL end_timer( electron_pusher_with_collisions_inner_object_timer )
     CALL start_timer( ions_pusher_with_collisions_inner_object_timer )    
     !t10 = MPI_WTIME()

     n_sub = n_sub + 1
     IF (n_sub.EQ.N_subcycles .AND. i_freeze_ions==0) THEN            ! N_subcycles is odd

        if (Rank_of_process.eq.0 .AND. debug_level>=local_debug_level) print '("----- doing ions at step ",i6," ------")', T_cntr

        CALL ADVANCE_IONS_PLUS                      !   velocity: n-N_e_subcycles+1/2 ---> n+1/2
                                                    ! coordinate: n-int(N_e_subcycles/2) ---> n-int(N_e_subcycles/2)+N_e_subcycles

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL SAVE_IONS_COLLIDED_WITH_BOUNDARY_OBJECTS

        ions_moved = .TRUE.
        CALL FIND_ALIENS_IN_ION_ADD_LIST
        CALL FIND_ALIENS_IN_ELECTRON_ADD_LIST

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL end_timer( ions_pusher_with_collisions_inner_object_timer )
        CALL start_timer( transfer_particle_after_pusher_timer )   
        !t11 = MPI_WTIME()

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN  ! send/receive BOTH electrons and ions crossing the borders
           CALL EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS            ! need only one X-pass for self-connected X-periodic clusters
        ELSE                                                              !
           CALL EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS             ! in general, three passes X-Y-X are needed
           CALL EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS            !
           CALL EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS             !
        END IF

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL FIND_INNER_OBJECT_COLL_IN_ELECTRON_ADD_LIST
        CALL FIND_INNER_OBJECT_COLL_IN_ION_ADD_LIST

        CALL PROCESS_ADDED_ELECTRONS                ! add the new electrons to the main array   !### NEW   
    
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        CALL end_timer( transfer_particle_after_pusher_timer )
        CALL start_timer( collect_particles_hitting_with_bo_timer )   
        !t12 = MPI_WTIME()

        CALL COLLECT_PARTICLE_BOUNDARY_HITS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      !   IF ( T_cntr>1000)STOP

        CALL end_timer( collect_particles_hitting_with_bo_timer )
        CALL start_timer( compute_mcc_timer ) 
        !t13 = MPI_WTIME()

        CALL COLLECT_ELECTRON_DENSITY_FOR_COLL_FREQS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL PERFORM_ELECTRON_NEUTRAL_COLLISIONS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL SAVE_en_COLLISIONS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL SAVE_en_COLLISIONS_2D

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL PERFORM_RESONANT_CHARGE_EXCHANGE    ! the only ion-neutral collision kind for now, does not produce new ions

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL SAVE_in_COLLISIONS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL PERFORM_PLASMA_EMISSION_ECR_SETUP ! This is for the ECR configuration
!###        CALL PERFORM_IONIZATION_HT_SETUP
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        CALL end_timer( compute_mcc_timer )
        CALL start_timer( add_ions_after_collisions_timer ) 
        !t14 = MPI_WTIME()

        CALL PROCESS_ADDED_IONS                  ! add the new ions to the main array

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL end_timer( add_ions_after_collisions_timer )
        CALL start_timer( clear_accumulated_fields_timer ) 
        !t15 = MPI_WTIME()

        CALL CLEAR_ACCUMULATED_FIELDS
        n_sub = 0                                 !### n_sub reset to zero here

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL end_timer( clear_accumulated_fields_timer )
        CALL start_timer( create_averaged_snapshot_and_compute_ptcl_emission_timer ) 
        !t16 = MPI_WTIME()

     ELSE

        CALL end_timer( ions_pusher_with_collisions_inner_object_timer )
        CALL start_timer( transfer_particle_after_pusher_timer ) 
        !t11 = t10

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN  ! send/receive ONLY electrons crossing the borders
           CALL EXCHANGE_ELECTRONS_WITH_ABOVE_BELOW_NEIGHBOURS            ! need only one X-pass for self-connected X-periodic clusters
        ELSE                                                              !
           CALL EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS             ! in general, three passes X-Y-X are needed
           CALL EXCHANGE_ELECTRONS_WITH_ABOVE_BELOW_NEIGHBOURS            !
           CALL EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS             !
        END IF

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        
        CALL end_timer( transfer_particle_after_pusher_timer )
        CALL start_timer( collect_particles_hitting_with_bo_timer )   
        !t12 = MPI_WTIME()

        CALL FIND_INNER_OBJECT_COLL_IN_ELECTRON_ADD_LIST

        CALL COLLECT_ELECTRON_BOUNDARY_HITS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL end_timer( collect_particles_hitting_with_bo_timer )
        CALL start_timer( compute_mcc_timer )    
        CALL end_timer( compute_mcc_timer )
        CALL start_timer( add_ions_after_collisions_timer )     
        CALL end_timer( add_ions_after_collisions_timer )
        CALL start_timer( clear_accumulated_fields_timer ) 
        CALL end_timer( clear_accumulated_fields_timer )
        CALL start_timer( create_averaged_snapshot_and_compute_ptcl_emission_timer )                          
      !   t13 = MPI_WTIME()
      !   t14 = t13
      !   t15 = t13
      !   t16 = t13

     END IF

     CALL DO_PROBE_DIAGNOSTICS

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL CREATE_AVERAGED_SNAPSHOT

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!###     CALL PERFORM_ELECTRON_EMISSION_HT_SETUP        ! either this or
                                                    ! PERFORM_ELECTRON_EMISSION_HT_SETUP_ZERO_GRAD_F (called above) works
                                                    ! not both
     CALL PERFORM_ELECTRON_EMISSION_SETUP           ! this works for non-HT setup

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL PERFORM_ELECTRON_EMISSION_SETUP_INNER_OBJECTS

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
     
     CALL end_timer( create_averaged_snapshot_and_compute_ptcl_emission_timer )
     CALL start_timer( add_electrons_after_emission_timer )   
     !t17 = MPI_WTIME()

     CALL PROCESS_ADDED_ELECTRONS                ! add the new electrons to the main array

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL end_timer( add_electrons_after_emission_timer )
     CALL start_timer( save_bo_particle_hits_emissions_timer )   
     !t18 = MPI_WTIME()

!###     CALL SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS_HT_SETUP  ! only one of the two will work
     CALL SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS           ! 

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL end_timer( save_bo_particle_hits_emissions_timer )
     CALL start_timer( gather_surface_charge_density_timer )   
     !t19 = MPI_WTIME()

     CALL GATHER_SURFACE_CHARGE_DENSITY      ! 

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL GATHER_SURFACE_CHARGE_DENSITY_INNER_OBJECTS   ! whole_object%surface_charge_variation=0 is done here

     CALL end_timer( gather_surface_charge_density_timer )
 

     !t20 = MPI_WTIME()
     ! Output time info
   !   loop_time = t20 - t0

   !   IF ( Rank_of_process==0 ) THEN
   !    message = "Save checkpoint"
   !    time_measured = t1  - t0
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"
   !    message = "Global load balance"
   !    time_measured = t2  - t1
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"
   !    message = "Internal load balance"
   !    time_measured = t3  - t2
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"                                                                                              
   !    message = "Gather ion charge density"
   !    time_measured = t4  - t3
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Gather electron charge density"
   !    time_measured = t5  - t4
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Solve Poisson equation"
   !    time_measured = t6  - t5
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Calculate electric field"
   !    time_measured = t7  - t6
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Compute averaged snapshot"
   !    time_measured = t8  - t7
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Create instantaneous snapshot"
   !    time_measured = t9  - t8
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Electron pusher + collisions with internal object"
   !    time_measured = t10  - t9
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Ion pusher + collisions with internal object"
   !    time_measured = t11  - t10
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"                                                    
   !    message = "Transfer particles after pusher"
   !    time_measured = t12  - t11
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Collect particles hitting with boundary"
   !    time_measured = t13  - t12
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Compute MCC"
   !    time_measured = t14  - t13
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Add ions after collisions"
   !    time_measured = t15  - t14
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Clear accumulated fields"
   !    time_measured = t16  - t15
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Create average snapshot + particle emission"
   !    time_measured = t17  - t16
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"     
   !    message = "Add electrons after emission"
   !    time_measured = t18  - t17
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"  
   !    message = "SAVE BOUNDARY PARTICLE HITS EMISSIONS"
   !    time_measured = t19  - t18
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"  
   !    message = "Gather surface charge density"
   !    time_measured = t20  - t19
   !    WRITE(*,'(A,T50,A,T120,E12.4,A,F6.2,A)') "CPU INFO ",">>> "//TRIM(message)//" (s) : ",&
   !                                              time_measured, " ==> ", time_measured/loop_time*100.0_8, "%"  
   !  ENDIF
   !   IF (Rank_of_process.EQ.0) PRINT '(2x,i9,2x,f8.3,2x,20(1x,f5.1))', &
   !        & T_cntr, &
   !        & REAL(t20 - t0), &
   !        & 100.0 * REAL((t1  - t0 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t2  - t1 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t3  - t2 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t4  - t3 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t5  - t4 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t6  - t5 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t7  - t6 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t8  - t7 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t9  - t8 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t10 - t9 ) / (t20 - t0)), &
   !        & 100.0 * REAL((t11 - t10) / (t20 - t0)), &
   !        & 100.0 * REAL((t12 - t11) / (t20 - t0)), &
   !        & 100.0 * REAL((t13 - t12) / (t20 - t0)), &
   !        & 100.0 * REAL((t14 - t13) / (t20 - t0)), &
   !        & 100.0 * REAL((t15 - t14) / (t20 - t0)), &
   !        & 100.0 * REAL((t16 - t15) / (t20 - t0)), &
   !        & 100.0 * REAL((t17 - t16) / (t20 - t0)), &
   !        & 100.0 * REAL((t18 - t17) / (t20 - t0)), &
   !        & 100.0 * REAL((t19 - t18) / (t20 - t0)), &
   !        & 100.0 * REAL((t20 - t19) / (t20 - t0))

     CALL end_timer( single_pic_loop_timer )
     CALL print_iteration_info(T_cntr)
  END DO
  CALL end_timer( total_timer )     
  pic_loop_timer = total_timer  

  finish = MPI_WTIME()

  !PRINT '(2x,"**** Process ",i3" : Simulation time is  : ", f12.3," sec")', Rank_of_process, finish - start

  CALL FINISH_SNAPSHOTS

!---------------------------------------------------------------------
!  Print CPU info
!---------------------------------------------------------------------
CALL print_cpu_time
!   message = "Total PIC loop"
!   CALL print_timer( pic_loop_timer,message )
!   message = "Save checkpoint"
!   CALL print_timer( save_checkpoint_timer,message )
!   message = "Global load balance"        
!   CALL print_timer( global_load_balance_timer,message )
!   message = "Internal load balance"        
!   CALL print_timer( internal_load_balance_timer,message )
!   message = "Gather ion charge density"
!   CALL print_timer( gather_ion_charge_density_timer,message )
!   message = "Gather electron charge density"
!   CALL print_timer( gather_electron_charge_density_timer,message )
!   message = "Solve Poisson equation"
!   CALL print_timer( poisson_solver_timer,message )
!   message = "Calculate electric field"
!   CALL print_timer( calculate_electric_field_timer,message )
!   message = "Compute averaged snapshot"
!   CALL print_timer( compute_averaged_snapshot_timer,message )
!   message = "Create instantaneous snapshot"
!   CALL print_timer( create_instantaneous_snapshot_timer,message )
!   message = "Electron pusher + collisions with internal object"
!   CALL print_timer( electron_pusher_with_collisions_inner_object_timer,message )
!   message = "Ion pusher + collisions with internal object"
!   CALL print_timer( ions_pusher_with_collisions_inner_object_timer,message )
!   message = "Transfer particles after pusher"
!   CALL print_timer( transfer_particle_after_pusher_timer,message )
!   message = "Collect particles hitting with boundary"
!   CALL print_timer( collect_particles_hitting_with_bo_timer,message )
!   message = "Compute MCC"
!   CALL print_timer( compute_mcc_timer,message )
!   message = "Add ions after collisions"
!   CALL print_timer( add_ions_after_collisions_timer,message )
!   message = "Clear accumulated fields"
!   CALL print_timer( clear_accumulated_fields_timer,message )
!   message = "Create average snapshot + particle emission"
!   CALL print_timer( create_averaged_snapshot_and_compute_ptcl_emission,message )
!   message = "Add electrons after emission"
!   CALL print_timer( add_electrons_after_emission_timer,message )
!   message = "Save bo particle hits emissions"
!   CALL print_timer( save_bo_particle_hits_emissions_timer,message )
!   message = "Gather surface charge density"
!   CALL print_timer( gather_surface_charge_density_timer,message )

!---------------------------------------------------------------------
!  End program
!---------------------------------------------------------------------

  CALL MPI_FINALIZE(ierr)

END PROGRAM MainProg

