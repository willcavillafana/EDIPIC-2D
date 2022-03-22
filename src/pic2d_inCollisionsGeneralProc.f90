
!===================================================================================================
SUBROUTINE INITIATE_ION_NEUTRAL_COLLISIONS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : kB_JK, e_Cl, N_max_vel, T_e_eV  !, amu_kg, m_e_kg, V_scale_ms
  USE IonParticles, ONLY : N_spec, Ms
!  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER n

  character(33) init_rcx_param_filename      ! init_neutral_AAAAAA_rcx_param.dat
                                             ! ----x----I----x----I----x----I---
 
  LOGICAL exists

  CHARACTER(1) buf
  INTEGER s

  INTEGER ALLOC_ERR

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  no_rcx_collisions = .TRUE.

! N_neutral_spec and other general neutral parameters are acquired from init_neutrals.dat [if it exists] in INITIATE_ELECTRON_NEUTRAL_COLLISIONS
  IF (N_neutral_spec.EQ.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### Neutrals deactivated, ion-neutral collisions are turned off ###")'
     RETURN
  END IF

  DO n = 1, N_neutral_spec
! reading the parameters for the resonant charge exchange model:
     init_rcx_param_filename = 'init_neutral_AAAAAA_rcx_param.dat'
     init_rcx_param_filename(14:19) = neutral(n)%name
     INQUIRE (FILE = init_rcx_param_filename, EXIST = exists)
     if (.not.exists) then
        if (Rank_of_process.eq.0) print '("### file ",A33," not found, RCX with neutral species ", i2," (",A6,") is turned off")', init_rcx_param_filename, n, neutral(n)%name
        neutral(n)%rcx_on = .false.
        neutral(n)%sigma_rcx_m2_1eV = 0.0_8 !just in case
        neutral(n)%alpha_rcx = 0.0_8        !just in case
        cycle
     end if
     open(9, file = init_rcx_param_filename)
     read(9, *) buf
     read(9, *) s, neutral(n)%sigma_rcx_m2_1eV, neutral(n)%alpha_rcx
     close(9, status = 'keep')
     if ((s.gt.0).and.(s.le.N_spec)) then
        no_rcx_collisions = .false.
        neutral(n)%rcx_on = .true.     
        neutral(n)%rcx_ion_species_index = s
     else
        neutral(n)%rcx_on = .false.
        neutral(n)%rcx_ion_species_index = 0
     end if     
  END DO

  IF (no_rcx_collisions) RETURN
    
  allocate (collision_rcx(1:N_spec), stat = alloc_err)
  do s = 1, N_spec
     collision_rcx(s)%rcx_on = .false.
  end do

  do n = 1, N_neutral_spec
     if (neutral(n)%rcx_on) then
        s = neutral(n)%rcx_ion_species_index
        collision_rcx(s)%rcx_on = .true.
        collision_rcx(s)%neutral_species_index = n
        collision_rcx(s)%vfactor = SQRT(neutral(n)%T_K * kB_JK / (T_e_eV * e_Cl * Ms(s))) / DBLE(N_max_vel)
!??        collision_rcx(s)%factor_eV = Ms(s) * energy_factor_eV !0.5_8 * m_e_kg * Ms(s) * V_scale_ms**2 / e_Cl
     end if
  end do

  CALL calculate_thermal_cx_probab

END SUBROUTINE INITIATE_ION_NEUTRAL_COLLISIONS

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE INITIATE_IONS_NEUTRAL_COLLISIONS_ELASTIC
!>    @details Reads parameters for elastic ions-neutral collisions. Should be merged with INITIATE_ION_NEUTRAL_COLLISIONS at some point 
!!    @authors W. Villafana
!!    @date    03-21-2022
!--------------------------------------------------------------------------------------------------

SUBROUTINE INITIATE_IONS_NEUTRAL_COLLISIONS_ELASTIC

   USE ParallelOperationValues
   USE MCCollisions, ONLY: in_collisions_turned_off, no_ionization_collisions, neutral, collision_i_neutral, N_neutral_spec
   USE CurrentProblemValues, ONLY : kB_JK, e_Cl, N_max_vel, T_e_eV, N_max_vel
   USE IonParticles, ONLY : N_spec, Ms
 !  USE ClusterAndItsBoundaries
 
   IMPLICIT NONE
 
   INCLUDE 'mpif.h'
 
   INTEGER ierr
  
   LOGICAL exists
 
   CHARACTER(1) buf
   INTEGER n
 
   CHARACTER(23) initneutral_filename   ! init_neutral_AAAAAA.dat
                                        ! ----x----I----x----I---
   INTEGER ALLOC_ERR
   INTEGER p, i, s
   INTEGER colflag
 
   CHARACTER(49) initneutral_crsect_filename  ! init_neutral_AAAAAA_crsect_coll_id_NN_type_NN.dat
                                              ! ----x----I----x----I----x----I----x----I----x----
 
   INTEGER count
 
   INTEGER j
 
   REAL(8) energy_eV
   INTEGER indx_energy
   INTEGER indx_energy_max_prob
   REAL(8) temp
 
 ! functions
   REAL(8) frequency_of_in_collision
 
   INTERFACE
      FUNCTION convert_int_to_txt_string(int_number, length_of_string)
        CHARACTER*(length_of_string) convert_int_to_txt_string
        INTEGER int_number
        INTEGER length_of_string
      END FUNCTION convert_int_to_txt_string
   END INTERFACE
 
   
   in_collisions_turned_off = .TRUE.
 
   INQUIRE (FILE = 'init_neutrals.dat', EXIST = exists)
 
   IF (.NOT.exists) THEN
      IF (Rank_of_process.EQ.0) PRINT '("### file init_neutrals.dat not found, ion-neutral collisions are turned off ###")'
      RETURN
   END IF
 
 ! specify neutral species included in simulation
 
   ! OPEN (9, FILE = 'init_neutrals.dat')
   ! READ (9, '(A1)') buf    !----------dd--- number of neutral species
   ! READ (9, '(10x,i2)') N_neutral_spec
 
   ! IF (N_neutral_spec.LT.0) N_neutral_spec = 0  ! will be used for checks later
 
   IF (N_neutral_spec.EQ.0) THEN
      IF (Rank_of_process.EQ.0) PRINT '("### file init_neutrals.dat found, neutrals deactivated, ion-neutral collisions are turned off ###")'
      RETURN
   END IF
 
   ! ALLOCATE(neutral(1:N_neutral_spec), STAT=ALLOC_ERR)
 
   ! DO n = 1, N_neutral_spec
   !    READ (9, '(A1)') buf !------AAAAAA--- code/abbreviation of the neutral gas, character string 
   !    READ (9, '(6x,A6)') neutral(n)%name
   !    READ (9, '(A1)') buf !--#d.dddE#dd--- Density [m^-3]
   !    READ (9, '(2x,e10.3)') neutral(n)%N_m3
   !    READ (9, '(A1)') buf !-----ddddd.d--- Temperature [K]
   !    READ (9, '(5x,f7.1)') neutral(n)%T_K
   ! END DO
   ! CLOSE (9, STATUS = 'KEEP')
 
 !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
 !print *, "read init_neutrals.dat"
 
 ! for each neutral species included, specify activated collisions
 
   DO n = 1, N_neutral_spec
      initneutral_filename = 'init_neutral_AAAAAA_ions.dat'
      initneutral_filename(14:19) = neutral(n)%name
      INQUIRE (FILE = initneutral_filename, EXIST = exists)
      IF (.NOT.exists) THEN
         IF (Rank_of_process.EQ.0) PRINT '("### file ",A23," NOT found, ion-neutral collisions for neutral species ",i2," (",A6,") are turned OFF ###")', &
              & initneutral_filename, n, neutral(n)%name
         neutral(n)%N_in_colproc = 0
         CYCLE
      END IF
 
      OPEN (9, FILE = initneutral_filename)
      READ (9, '(A1)') buf !---ddd.ddd--- mass [a.m.u.]
      READ (9, '(3x,f7.3)') neutral(n)%M_amu
      READ (9, '(A1)') buf !--------dd--- number of all possible collisional processes
      READ (9, '(8x,i2)') neutral(n)%N_in_colproc
      IF (neutral(n)%N_in_colproc.LE.0) THEN
         IF (Rank_of_process.EQ.0) PRINT '("### file ",A23," does NOT specify ion-neutral collisions for neutral species ",i2," (",A6,"), the collisions are turned OFF ###")', &
              & initneutral_filename, n, neutral(n)%name
         CLOSE (9, STATUS = 'KEEP')
         CYCLE
      END IF
 
      ALLOCATE(neutral(n)%in_colproc(1:neutral(n)%N_in_colproc), STAT=ALLOC_ERR)
      DO p = 1, neutral(n)%N_in_colproc
         READ (9, '(A1)') buf !---dd--d--dd--- collision #NN :: type / activated (1/0 = Yes/No) 
         READ (9, '(3x,i2,2x,i1,2x)') neutral(n)%in_colproc(p)%type, colflag
         neutral(n)%in_colproc(p)%activated = .FALSE.
         IF (colflag.NE.0) neutral(n)%in_colproc(p)%activated = .TRUE.
      END DO
 
      READ (9, '(A1)') buf !--------dd--- number of energy segments for collision probabilities (>0)
      READ (9, '(8x,i2)') neutral(n)%N_of_energy_segments_ions
      ALLOCATE(neutral(n)%energy_segment_boundary_value_ions(0:neutral(n)%N_of_energy_segments_ions), STAT=ALLOC_ERR)
      ALLOCATE(neutral(n)%energy_segment_step_ions(1:neutral(n)%N_of_energy_segments_ions), STAT=ALLOC_ERR)
      READ (9, '(A1)') buf !--ddddd.ddd------------- minimal energy [eV]
      READ (9, '(2x,f9.3)') neutral(n)%energy_segment_boundary_value_ions(0)
      DO i = 1, neutral(n)%N_of_energy_segments_ions
         READ (9, '(A1)') buf !--ddddd.ddd---ddd.ddd--- energy segment NN :: upper boundary [eV] / resolution [eV]
         READ (9, '(2x,f9.3,3x,f7.3)') neutral(n)%energy_segment_boundary_value_ions(i), neutral(n)%energy_segment_step_ions(i)
      END DO
 
      CLOSE (9, STATUS = 'KEEP')
   END DO
 
 !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
 !print *, "read ", initneutral_filename
 
 ! collision types
 ! 10 = elastic
 ! 20 = inelastic
 !   21, 22, 23, etc.
 ! 30 = ionization, ion+ and e-
 ! 40 = ionization, ion++ and e- e-
 ! etc, we do the 3 first types only, though allow multiple inelastic collisions
 
 ! init_neutral_Xenon-.dat
 ! cross sections
 ! init_neutral_Xenon-_crsect_colltype_10.dat
 ! init_neutral_Xenon-_crsect_colltype_20.dat
 ! init_neutral_Xenon-_crsect_colltype_30.dat
 ! init_neutral_Xenon-_crsect_colltype_40.dat etc
 ! 
 
 ! for each neutral species included, for each activated collisional process, read cross sections
 
   DO n = 1, N_neutral_spec
      DO p = 1, neutral(n)%N_in_colproc
 
         IF (.NOT.neutral(n)%in_colproc(p)%activated) CYCLE
 
         in_collisions_turned_off = .FALSE.                                           !### flip the general collision switch
  
         initneutral_crsect_filename = 'init_neutral_AAAAAA_crsect_coll_id_NN_type_NN_ions.dat'
         initneutral_crsect_filename(14:19) = neutral(n)%name
         initneutral_crsect_filename(36:37) = convert_int_to_txt_string(p, 2)
         initneutral_crsect_filename(44:45) = convert_int_to_txt_string(neutral(n)%in_colproc(p)%type, 2)
         INQUIRE( FILE = initneutral_crsect_filename, EXIST = exists)
         IF (.NOT.exists) THEN
            IF (Rank_of_process.EQ.0) PRINT '("###ERROR :: file ",A49," not found, for neutrals ",A6," while collisions ",i2," of type ",i2," are activated, program terminated")', &
                 & initneutral_crsect_filename, neutral(n)%name, p, neutral(n)%in_colproc(p)%type
            CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
         END IF
 
         OPEN(9, FILE = initneutral_crsect_filename)
         READ (9, '(2x,i4)') neutral(n)%in_colproc(p)%N_crsect_points
         ALLOCATE(neutral(n)%in_colproc(p)%energy_eV(1:neutral(n)%in_colproc(p)%N_crsect_points), STAT = ALLOC_ERR)
         ALLOCATE(neutral(n)%in_colproc(p)%crsect_m2(1:neutral(n)%in_colproc(p)%N_crsect_points), STAT = ALLOC_ERR)
         DO i = 1, neutral(n)%in_colproc(p)%N_crsect_points
            READ (9, '(4x,f9.3,2x,e10.3)') neutral(n)%in_colproc(p)%energy_eV(i), neutral(n)%in_colproc(p)%crsect_m2(i)
         END DO
         CLOSE (9, STATUS = 'KEEP')
 
 ! set thresholds for inelastic/ionization collisions
         neutral(n)%in_colproc(p)%threshold_energy_eV = 0.0_8   ! for elastic collisions
         ! For now discard other collisions than elastic
         !IF (neutral(n)%en_colproc(p)%type.GE.20) neutral(n)%en_colproc(p)%threshold_energy_eV = neutral(n)%en_colproc(p)%energy_eV(1)
 
      END DO
   END DO
 
   ALLOCATE(collision_i_neutral(1:N_neutral_spec), STAT = ALLOC_ERR)
   DO n = 1, N_neutral_spec
 
 ! count activated collision procedures
      collision_i_neutral(n)%N_of_activated_colproc = 0
      DO p = 1, neutral(n)%N_in_colproc
         IF (neutral(n)%in_colproc(p)%activated) collision_i_neutral(n)%N_of_activated_colproc = collision_i_neutral(n)%N_of_activated_colproc + 1
      END DO
 
      IF (collision_i_neutral(n)%N_of_activated_colproc.EQ.0) CYCLE
 
 ! setup information about the activated collisions
      ALLOCATE(collision_i_neutral(n)%colproc_info(1:collision_i_neutral(n)%N_of_activated_colproc), STAT = ALLOC_ERR)     
      count=0
      DO p = 1, neutral(n)%N_in_colproc
         IF (.NOT.neutral(n)%in_colproc(p)%activated) CYCLE
         count=count+1 
         collision_i_neutral(n)%colproc_info(count)%id_number = p
         collision_i_neutral(n)%colproc_info(count)%type = neutral(n)%in_colproc(p)%type
         collision_i_neutral(n)%colproc_info(count)%threshold_energy_eV = neutral(n)%in_colproc(p)%threshold_energy_eV 
 
      END DO
 
 ! setup general counters of collision events for minimal diagnostics
      ALLOCATE(collision_i_neutral(n)%counter(1:collision_i_neutral(n)%N_of_activated_colproc), STAT = ALLOC_ERR)     
      DO p = 1, neutral(n)%N_in_colproc
         collision_i_neutral(n)%counter(p) = 0
      END DO
 
 ! copy energy segments and steps
      collision_i_neutral(n)%N_of_energy_segments = neutral(n)%N_of_energy_segments_ions
      ALLOCATE(collision_i_neutral(n)%energy_segment_boundary_value(0:collision_i_neutral(n)%N_of_energy_segments), STAT=ALLOC_ERR)
      collision_i_neutral(n)%energy_segment_boundary_value(0:collision_i_neutral(n)%N_of_energy_segments) = &
                & neutral(n)%energy_segment_boundary_value_ions(0:            neutral(n)%N_of_energy_segments_ions)
      ALLOCATE(collision_i_neutral(n)%energy_segment_step(1:collision_i_neutral(n)%N_of_energy_segments), STAT=ALLOC_ERR)
      collision_i_neutral(n)%energy_segment_step(1:collision_i_neutral(n)%N_of_energy_segments) = &
                & neutral(n)%energy_segment_step_ions(1:            neutral(n)%N_of_energy_segments_ions)
 
 ! identify total number of probability points / energy values
      collision_i_neutral(n)%N_of_energy_values = 0
      DO i = 1, collision_i_neutral(n)%N_of_energy_segments
         collision_i_neutral(n)%N_of_energy_values = collision_i_neutral(n)%N_of_energy_values + &
                                                   & INT( ( collision_i_neutral(n)%energy_segment_boundary_value(i) - &
                                                   &        collision_i_neutral(n)%energy_segment_boundary_value(i-1) ) / &
                                                   &      collision_i_neutral(n)%energy_segment_step(i) + 0.001_8 )
      END DO
 
 ! set index values for energy segment boundaries
      ALLOCATE(collision_i_neutral(n)%energy_segment_boundary_index(0:collision_i_neutral(n)%N_of_energy_segments), STAT=ALLOC_ERR)
      collision_i_neutral(n)%energy_segment_boundary_index(0) = 0
      DO i = 1, collision_i_neutral(n)%N_of_energy_segments
         collision_i_neutral(n)%energy_segment_boundary_index(i) = collision_i_neutral(n)%energy_segment_boundary_index(i-1) + &
                                                                 & INT( ( collision_i_neutral(n)%energy_segment_boundary_value(i) - &
                                                                 &        collision_i_neutral(n)%energy_segment_boundary_value(i-1) ) / &
                                                                 &      collision_i_neutral(n)%energy_segment_step(i) + 0.001_8 )
      END DO
 
 ! fool proof (paranoidal)
      IF (collision_i_neutral(n)%N_of_energy_values.NE.collision_i_neutral(n)%energy_segment_boundary_index(collision_i_neutral(n)%N_of_energy_segments)) THEN
 ! error
         PRINT '("error")'
         CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
      END IF
 
 ! create array of energy values corresponding to the probability array
      ALLOCATE(collision_i_neutral(n)%energy_eV(0:collision_i_neutral(n)%N_of_energy_values), STAT=ALLOC_ERR)
      collision_i_neutral(n)%energy_eV(0) = collision_i_neutral(n)%energy_segment_boundary_value(0)
      DO i = 1, collision_i_neutral(n)%N_of_energy_segments
         DO j = collision_i_neutral(n)%energy_segment_boundary_index(i-1)+1, collision_i_neutral(n)%energy_segment_boundary_index(i)-1
            collision_i_neutral(n)%energy_eV(j) = collision_i_neutral(n)%energy_segment_boundary_value(i-1) + &
                                                & (j-collision_i_neutral(n)%energy_segment_boundary_index(i-1)) * &
                                                & collision_i_neutral(n)%energy_segment_step(i)
         END DO
         collision_i_neutral(n)%energy_eV(collision_i_neutral(n)%energy_segment_boundary_index(i)) = collision_i_neutral(n)%energy_segment_boundary_value(i)
      END DO
 
 ! set the probability arrays --------------------------------------------------------->>>>>>>>>>>>>>>
 
      ALLOCATE(collision_i_neutral(n)%prob_colproc_energy(1:collision_i_neutral(n)%N_of_activated_colproc, 0:collision_i_neutral(n)%N_of_energy_values), STAT = ALLOC_ERR)
 
      DO indx_energy = 0, collision_i_neutral(n)%N_of_energy_values
 
         energy_eV = collision_i_neutral(n)%energy_eV(indx_energy)
 
         collision_i_neutral(n)%prob_colproc_energy(1, indx_energy) = frequency_of_in_collision(energy_eV, n, collision_i_neutral(n)%colproc_info(1)%id_number)
         DO p = 2, collision_i_neutral(n)%N_of_activated_colproc
            collision_i_neutral(n)%prob_colproc_energy(p, indx_energy) = collision_i_neutral(n)%prob_colproc_energy(p-1, indx_energy) + &
                                                                    & frequency_of_in_collision(energy_eV, n, collision_i_neutral(n)%colproc_info(p)%id_number)
         END DO
      END DO
 
 ! find the maximum
      indx_energy_max_prob = 0
      DO indx_energy = 1, collision_i_neutral(n)%N_of_energy_values
         IF ( collision_i_neutral(n)%prob_colproc_energy(collision_i_neutral(n)%N_of_activated_colproc, indx_energy).GT. &
            & collision_i_neutral(n)%prob_colproc_energy(collision_i_neutral(n)%N_of_activated_colproc, indx_energy_max_prob) ) indx_energy_max_prob = indx_energy
      END DO
 
 ! find maximal fraction of colliding particles for the neutral density neutral(n)%N_m3
      collision_i_neutral(n)%max_colliding_fraction = 1.0_8 - EXP(-collision_i_neutral(n)%prob_colproc_energy(collision_i_neutral(n)%N_of_activated_colproc, indx_energy_max_prob))
 
 ! convert accumulated collision frequencies into boundaries of probability ranges for collision processes
      temp = collision_i_neutral(n)%prob_colproc_energy(collision_i_neutral(n)%N_of_activated_colproc, indx_energy_max_prob)
      DO indx_energy = 0, collision_i_neutral(n)%N_of_energy_values
         DO p = 1, collision_i_neutral(n)%N_of_activated_colproc
            collision_i_neutral(n)%prob_colproc_energy(p, indx_energy) = MAX(0.0_8, MIN(1.0_8, collision_i_neutral(n)%prob_colproc_energy(p, indx_energy)/temp))
         END DO
      END DO
      collision_i_neutral(n)%prob_colproc_energy(collision_i_neutral(n)%N_of_activated_colproc, indx_energy_max_prob) = 1.0_8
 
 ! fool proof (paranoidal again) - make sure that boundaries of probability ranges do not decrease
      DO indx_energy = 0, collision_i_neutral(n)%N_of_energy_values
         DO p = 1, collision_i_neutral(n)%N_of_activated_colproc-1
            collision_i_neutral(n)%prob_colproc_energy(p+1, indx_energy) = MAX( collision_i_neutral(n)%prob_colproc_energy(p,   indx_energy), &
                                                                              & collision_i_neutral(n)%prob_colproc_energy(p+1, indx_energy) )
         END DO
      END DO
 
   END DO   !###  DO n = 1, N_neutral_spec
 
 END SUBROUTINE INITIATE_IONS_NEUTRAL_COLLISIONS_ELASTIC

!-------------------------------------------------------------------------------------------
!
SUBROUTINE INITIATE_in_COLL_DIAGNOSTICS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : Start_T_cntr, N_subcycles, delta_t_s
  USE Checkpoints, ONLY : use_checkpoint
  USE IonParticles, ONLY : N_spec

  IMPLICIT NONE

  INTEGER s, n, i
                                      ! ----x----I----x----I----x----
  CHARACTER(29) historycoll_filename  ! history_coll_i_S_n_AAAAAA.dat

  LOGICAL exists
  CHARACTER(1) buf
  INTEGER i_dummy

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (Rank_of_process.NE.0) RETURN

  IF (no_rcx_collisions) RETURN

  IF (use_checkpoint.EQ.1) THEN
! start from checkpoint, must trim the time dependences

     DO s = 1, N_spec
        IF (.NOT.collision_rcx(s)%rcx_on) CYCLE

        n = collision_rcx(s)%neutral_species_index

        historycoll_filename = 'history_coll_i_S_n_AAAAAA.dat'
        historycoll_filename(16:16) = convert_int_to_txt_string(s, 1)
        historycoll_filename(20:25) = neutral(n)%name

        INQUIRE (FILE = historycoll_filename, EXIST = exists)
        IF (exists) THEN                                                       
           OPEN (21, FILE = historycoll_filename, STATUS = 'OLD')          
! skip header
           READ (21, '(A1)') buf ! WRITE (21, '("# electron time step is ",e14.7," s")'), delta_t_s
           READ (21, '(A1)') buf ! WRITE (21, '("#      ion time step is ",e14.7," s")'), delta_t_s * N_subcycles
           READ (21, '(A1)') buf ! WRITE (21, '("# column  1 is the electron step counter")')
           READ (21, '(A1)') buf ! WRITE (21, '("# column  2 is the total number of ion macroparticles of species SS in the whole system")')
           READ (21, '(A1)') buf ! WRITE (21, '("# column  3 is the number of rezonance charge exchange collision events during past ion time step")')

           DO i = 1, Start_T_cntr / N_subcycles            ! these files are updated at every ion timestep
              READ (21, '(2x,i9,2x,i9,2x,i8)') i_dummy
           END DO
           ENDFILE 21       
           CLOSE (21, STATUS = 'KEEP')
        ELSE

           OPEN  (21, FILE = historycoll_filename, STATUS = 'REPLACE')
! save header, for now resonance charge exchange only
           WRITE (21, '("# electron time step is ",e14.7," s")') delta_t_s
           WRITE (21, '("#      ion time step is ",e14.7," s")') delta_t_s * N_subcycles
           WRITE (21, '("# column  1 is the electron step counter")')
           WRITE (21, '("# column  2 is the total number of ion macroparticles of species ",i2," in the whole system")') s
           WRITE (21, '("# column  3 is the number of rezonance charge exchange collision events during past ion time step")')
           CLOSE (21, STATUS = 'KEEP')

        END IF

     END DO

  ELSE
! fresh start - create empty files with a header

     DO s = 1, N_spec
        IF (.NOT.collision_rcx(s)%rcx_on) CYCLE
    
        n = collision_rcx(s)%neutral_species_index

        historycoll_filename = 'history_coll_i_S_n_AAAAAA.dat'
        historycoll_filename(16:16) = convert_int_to_txt_string(s, 1)
        historycoll_filename(20:25) = neutral(n)%name

        OPEN  (21, FILE = historycoll_filename, STATUS = 'REPLACE')
! save header, for now resonance charge exchange only
        WRITE (21, '("# electron time step is ",e14.7," s")') delta_t_s
        WRITE (21, '("#      ion time step is ",e14.7," s")') delta_t_s * N_subcycles
        WRITE (21, '("# column  1 is the electron step counter")')
        WRITE (21, '("# column  2 is the total number of ion macroparticles of species ",i2," in the whole system")') s
        WRITE (21, '("# column  3 is the number of rezonance charge exchange collision events during past ion time step")')
        CLOSE (21, STATUS = 'KEEP')

     END DO

  END IF

END SUBROUTINE INITIATE_in_COLL_DIAGNOSTICS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_in_COLLISIONS

  USE ParallelOperationValues
  USE MCCollisions
  USE IonParticles, ONLY : N_spec, N_ions
  USE CurrentProblemValues, ONLY : T_cntr

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER request

  INTEGER buflen, s, n
  INTEGER, ALLOCATABLE :: ibufer_send(:)
  INTEGER, ALLOCATABLE :: ibufer_receive(:)
  INTEGER ALLOC_ERR

                                      ! ----x----I----x----I----x----
  CHARACTER(29) historycoll_filename  ! history_coll_i_S_n_AAAAAA.dat

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (no_rcx_collisions) RETURN

! report all collision counters to the process with zero global rank

  buflen=2*N_spec             ! include N_ions(1:N_spec) and collision_rcx(1:N_spec)%counter
                              ! note that presently rcx is built so that there is only one typr of rcx collisions and only with the same neutral kind
                              ! that is rcx does not create new ion species, only replaces velocity components of the collided ion

  ALLOCATE (ibufer_send(1:buflen), STAT = ALLOC_ERR)
  ALLOCATE (ibufer_receive(1:buflen), STAT = ALLOC_ERR)

  ibufer_send(1:N_spec) = N_ions(1:N_spec)
  DO s = 1, N_spec
     ibufer_send(N_spec+s) = collision_rcx(s)%counter
  END DO

  ibufer_receive = 0

  CALL MPI_REDUCE(ibufer_send, ibufer_receive, buflen, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) THEN

     PRINT '("Total charge exchange collisions for all ions species :: ",10(2x,i6))', ibufer_receive(N_spec+1:N_spec+N_spec)

     DO s = 1, N_spec

        IF (.NOT.collision_rcx(s)%rcx_on) CYCLE  

        n = collision_rcx(s)%neutral_species_index

        historycoll_filename = 'history_coll_i_S_n_AAAAAA.dat'
        historycoll_filename(16:16) = convert_int_to_txt_string(s, 1)
        historycoll_filename(20:25) = neutral(n)%name

        OPEN (21, FILE = historycoll_filename, POSITION = 'APPEND')
        WRITE (21, '(2x,i9,2x,i9,2x,i8)') &
             & T_cntr, &
             & ibufer_receive(s), &
             & ibufer_receive(N_spec+s)
        CLOSE (21, STATUS = 'KEEP')
  
     END DO
  END IF

  DEALLOCATE(ibufer_send, STAT = ALLOC_ERR)
  DEALLOCATE(ibufer_receive, STAT = ALLOC_ERR)

END SUBROUTINE SAVE_in_COLLISIONS
