
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
!     SUBROUTINE INITIATE_IONS_NEUTRAL_COLLISIONS_GENERAL
!>    @details Reads parameters for elastic and charge exchange ions-neutral collisions. Collisions are based on cross sections tables. Should be merged with INITIATE_ION_NEUTRAL_COLLISIONS at some point
!!    @authors W. Villafana
!!    @date    03-21-2022
!--------------------------------------------------------------------------------------------------

SUBROUTINE INITIATE_IONS_NEUTRAL_COLLISIONS_GENERAL

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
 
   CHARACTER(28) initneutral_filename   ! init_neutral_AAAAAA_ions.dat
                                        ! ----x----I----x----I---
   INTEGER ALLOC_ERR
   INTEGER p, i, s
   INTEGER colflag
 
   CHARACTER(54) initneutral_crsect_filename  ! init_neutral_AAAAAA_crsect_coll_id_NN_type_NN_ions.dat
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
         READ (9, '(A1)') buf !---dd--d--- collision #NN :: type / activated (1/0 = Yes/No) 
         READ (9, '(3x,i2,2x,i1)') neutral(n)%in_colproc(p)%type, colflag
         neutral(n)%in_colproc(p)%activated = .FALSE.
         IF (colflag.NE.0) neutral(n)%in_colproc(p)%activated = .TRUE.

         ! Check if charge exchange has been requested for analytical model 
         IF ( neutral(n)%rcx_on .AND. Rank_of_process==0 .AND. neutral(n)%in_colproc(p)%type==20 .AND. colflag==1 ) THEN
            PRINT*,'You are requesting ion-neutral charge exchange from cross section table. However you have already requested charge exchange as an analytical model. Choose one.'
            STOP
         END IF
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
 ! 20 = charge exchange
 
 
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
            PRINT '("###ERROR :: file ",A54," not found, for neutrals ",A6," while collisions ",i2," of type ",i2," are activated, program terminated")', &
                 & initneutral_crsect_filename, neutral(n)%name, p, neutral(n)%in_colproc(p)%type
            CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
         END IF
         OPEN(9, FILE = initneutral_crsect_filename)
         READ (9, '(2x,i5)') neutral(n)%in_colproc(p)%N_crsect_points
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
 
END SUBROUTINE INITIATE_IONS_NEUTRAL_COLLISIONS_GENERAL

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE PERFORM_ION_NEUTRAL_COLLISIONS_ELASTIC
!>    @details Performs elastic and charge exchange ions-neutral collisions from cross section tables
!!    @authors W. Villafana
!!    @date    08-24-2022
!--------------------------------------------------------------------------------------------------

SUBROUTINE PERFORM_ION_NEUTRAL_COLLISIONS

   USE ParallelOperationValues
   USE MCCollisions
   USE IonParticles
   USE CurrentProblemValues, ONLY : V_scale_ms, m_e_kg, e_Cl
   USE ClusterAndItsBoundaries
   USE rng_wrapper
   USE Snapshots
 
 !  USE ParallelOperationValues
 
   IMPLICIT NONE
 
   INCLUDE 'mpif.h'
 
   INTEGER ierr
   INTEGER stattus(MPI_STATUS_SIZE)
   INTEGER request
 
   INTEGER ALLOC_ERR
 
   INTEGER n, p
 
   REAL(8) R_collided, F_collided
   INTEGER I_collided
 
   INTEGER j
 
   REAL(8) random_r, random_n
   INTEGER random_j
 
   REAL(8) energy_eV
   INTEGER indx_energy
   REAL(8) a1, a0
   INTEGER i, indx_segment
 
   INTEGER k, indx_coll
 
   INTEGER n1, n2, n3, bufsize, pos
   REAL, ALLOCATABLE :: rbufer(:), rbufer2(:)

   INTEGER :: ion_index ! index of ion species that will undergo elactis collision. For now only one ion specie can have elastic collisions
   REAL(8) :: Vx_n,Vy_n,Vz_n ! speed of neutral particle

   ! functions
   LOGICAL Find_in_stored_list_ions
   REAL(8) neutral_density_normalized

   INTERFACE
      RECURSIVE SUBROUTINE Node_Killer_ion(node)
        USE MCCollisions
        TYPE (binary_tree), POINTER :: node
      END SUBROUTINE Node_Killer_ion
 
      RECURSIVE SUBROUTINE Transfer_collisions_from_stored_list_ions(node, n_neutral, indx_coll, bufsize, n1, n2, n3, rbufer)
        USE MCCollisions
        USE IonParticles
        USE ClusterAndItsBoundaries
        TYPE (binary_tree), POINTER :: node
        INTEGER, INTENT(IN) :: n_neutral, indx_coll, bufsize, n1, n2, n3
        REAL, DIMENSION(bufsize), INTENT(INOUT) :: rbufer
      END SUBROUTINE Transfer_collisions_from_stored_list_ions
   END INTERFACE

   IF (in_collisions_turned_off) RETURN

   ion_index = 1 ! By default first ions species will have collision.

 ! Allocate binary tree to store the numbers of particles which have collided already
   NULLIFY(Collided_particle_ions)
   IF (.NOT.ASSOCIATED(Collided_particle_ions)) THEN
      ALLOCATE(Collided_particle_ions, STAT=ALLOC_ERR)
   END IF
   NULLIFY(Collided_particle_ions%Larger)
   NULLIFY(Collided_particle_ions%Smaller)
 
 ! clear all collision counters
   DO n = 1, N_neutral_spec
      DO p = 1, collision_i_neutral(n)%N_of_activated_colproc
         collision_i_neutral(n)%counter(p) = 0
      END DO
   END DO
   
   DO n = 1, N_neutral_spec
 
      IF (collision_i_neutral(n)%N_of_activated_colproc.EQ.0) CYCLE
      R_collided = collision_i_neutral(n)%max_colliding_fraction * N_ions(ion_index)
      I_collided = INT(R_collided)
      F_collided = R_collided - I_collided

      DO j = 0, I_collided
 
         IF (j.EQ.0) THEN
 ! here we process the "fractional" collisional event
            IF (well_random_number().GT.F_collided) CYCLE
         END IF

 !------------- Determine the kind of collision for the selected particle
         random_r = well_random_number()
         random_n = well_random_number()
 
         DO                        ! search will be repeated until a number will be successfully obtained
            random_j = INT(well_random_number() * N_ions(ion_index))
            random_j = MIN(MAX(random_j, 1), N_ions(ion_index))
            ! print*,'rand,hello',random_j
            IF (.NOT.Find_in_stored_list_ions(random_j)) EXIT    !#### needs some safety mechanism to avoid endless cycling
         END DO
 ! account for reduced neutral density
         IF (random_n.GT.neutral_density_normalized(n, ion(ion_index)%part(random_j)%X, ion(ion_index)%part(random_j)%Y)) CYCLE      
 
         !!!! Compute incident energy of ions with respect to neutral speed
         ! Compute neutral speed in dimensionless form
         CALL GetMaxwellVelocity(Vx_n)
         CALL GetMaxwellVelocity(Vy_n)
         CALL GetMaxwellVelocity(Vz_n)
         ! Compute energy of incident ion
         energy_eV = (  (ion(ion_index)%part(random_j)%VX-Vx_n)**2 + &
                      & (ion(ion_index)%part(random_j)%VY-Vy_n)**2 + &
                      & (ion(ion_index)%part(random_j)%VZ-Vz_n)**2) * V_scale_ms**2 * Ms(ion_index) * m_e_kg * 0.5_8 / e_Cl
 
         IF (energy_eV.GE.collision_i_neutral(n)%energy_segment_boundary_value(collision_i_neutral(n)%N_of_energy_segments)) THEN
            indx_energy = collision_i_neutral(n)%N_of_energy_values-1
            a1 = 1.0_8
            a0 = 0.0_8
         ELSE IF (energy_eV.LT.collision_e_neutral(n)%energy_segment_boundary_value(0)) THEN
            indx_energy = 0
            a1 = 0.0_8
            a0 = 1.0_8
         ELSE
            DO i = collision_i_neutral(n)%N_of_energy_segments-1, 0, -1
               IF (energy_eV.GE.collision_i_neutral(n)%energy_segment_boundary_value(i)) THEN
                  indx_segment = i+1
                  EXIT
               END IF
            END DO
            indx_energy =              collision_i_neutral(n)%energy_segment_boundary_index(indx_segment-1) + &
                        & (energy_eV - collision_i_neutral(n)%energy_segment_boundary_value(indx_segment-1)) / collision_i_neutral(n)%energy_segment_step(indx_segment)
            indx_energy = MAX(0, MIN(indx_energy, collision_i_neutral(n)%N_of_energy_values-1))
 ! double check
            IF ((energy_eV.LT.collision_i_neutral(n)%energy_eV(indx_energy)).OR.(energy_eV.GT.collision_i_neutral(n)%energy_eV(indx_energy+1))) THEN
 ! error
               PRINT '("Proc ",i4," error in PERFORM_ION_NEUTRAL_COLLISIONS")', Rank_of_process
               CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
            END IF
            a1 = (energy_eV - collision_i_neutral(n)%energy_eV(indx_energy)) / collision_i_neutral(n)%energy_segment_step(indx_segment)
            a0 = 1.0_8 - a1
         END IF
 
         DO k = collision_i_neutral(n)%N_of_activated_colproc, 1, -1 
            IF (random_r.GT.(a0 * collision_i_neutral(n)%prob_colproc_energy(k, indx_energy) + &
                          &  a1 * collision_i_neutral(n)%prob_colproc_energy(k, indx_energy+1)) ) EXIT
         END DO
         indx_coll = k + 1        
 
         IF (indx_coll.GT.collision_i_neutral(n)%N_of_activated_colproc) CYCLE   ! the null collision
         CALL Add_to_stored_list_ions(random_j, n, indx_coll)

         SELECT CASE (collision_i_neutral(n)%colproc_info(indx_coll)%type)

         ! scattering with the angle Ksi as in EDIPIC 1D
         CASE (10)
            CALL in_Collision_Elastic( ion_index, random_j, collision_i_neutral(n)%counter(indx_coll),Vx_n,Vy_n,Vz_n)
            
         ! Charge Exchange from given table
         CASE (20)
            CALL in_Collision_Charge_Exchange( ion_index, random_j, collision_i_neutral(n)%counter(indx_coll),Vx_n,Vy_n,Vz_n)            

         END SELECT
      END DO
   END DO
 
   ! print*,'ici'
   CALL Node_Killer_ion(Collided_particle_ions)
   ! print*,'ici1'

END SUBROUTINE PERFORM_ION_NEUTRAL_COLLISIONS

LOGICAL FUNCTION Find_in_stored_list_ions(number)

  USE MCCollisions
  IMPLICIT NONE

  INTEGER number
  TYPE (binary_tree), POINTER :: current

  Find_in_stored_list_ions = .FALSE.

  current => Collided_particle_ions

  DO 
      ! print*,'find number, current%number',number, current%number
     IF (number.GT.current%number) THEN
        IF (ASSOCIATED(current%Larger)) THEN
           current => current%Larger               ! 
           CYCLE                                   ! go to the next node, with larger "number"
        ELSE
           EXIT
        END IF
     END IF

     IF (number.LT.current%number) THEN
        IF (ASSOCIATED(current%Smaller)) THEN
           current => current%Smaller              ! 
           CYCLE                                   ! go to the next node, with smaller "number"
        ELSE
           EXIT
        END IF
     END IF

     Find_in_stored_list_ions = .TRUE.                  ! number.EQ.current%number
   !   print*,'find number, current%number',number, current%number
     EXIT                                          ! if we are here, then we found the match
     
  END DO

END FUNCTION Find_in_stored_list_ions

!-----------------------------------------------------------------
! subroutine adds number to the binary tree
! ! we assume that there are no nodes in the tree with the same value yet
SUBROUTINE Add_to_stored_list_ions(number, n_neutral, indx_coll)

   USE MCCollisions
 
   IMPLICIT NONE
 
   INTEGER number
   INTEGER n_neutral
   INTEGER indx_coll
 
   TYPE (binary_tree), POINTER :: current
   INTEGER ALLOC_ERR
 
   current => Collided_particle_ions                  ! start from the head node of the binary tree
 
   DO                                            ! go through the allocated nodes to the end of the branch
      IF (number.GT.current%number) THEN         
         IF (ASSOCIATED(current%Larger)) THEN        
            current => current%Larger               ! 
            CYCLE                                   ! go to the next node, with larger "number"
         ELSE
            ALLOCATE(current%Larger, STAT=ALLOC_ERR)
            current => current%Larger
            EXIT
         END IF
      END IF
 
      IF (number.LT.current%number) THEN
         IF (ASSOCIATED(current%Smaller)) THEN        
            current => current%Smaller              ! 
            CYCLE                                   ! go to the next node, with smaller "number"
         ELSE
            ALLOCATE(current%Smaller, STAT=ALLOC_ERR)
            current => current%Smaller
            EXIT
         END IF
      END IF
      
   END DO
 
   current%number    = number                       ! collided electron particle number
   current%neutral   = n_neutral                    ! neutral species number
   current%indx_coll = indx_coll                    ! index of activated collision [not collision id]
 
   NULLIFY(current%Larger)
   NULLIFY(current%Smaller)
 
END SUBROUTINE Add_to_stored_list_ions

RECURSIVE SUBROUTINE Transfer_collisions_from_stored_list_ions(node, n_neutral, indx_coll, bufsize, n1, n2, n3, rbufer)

 USE MCCollisions
 USE IonParticles
 USE ClusterAndItsBoundaries

 IMPLICIT NONE

 TYPE (binary_tree), POINTER :: node

 INTEGER, INTENT(IN) :: n_neutral, indx_coll, bufsize, n1, n2, n3
 REAL, DIMENSION(bufsize), INTENT(INOUT) :: rbufer

 INTEGER i, j, k
 INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
 REAL ax_ip1, ax_i, ay_jp1, ay_j
 REAL vij, vip1j, vijp1

 INTEGER :: ion_index

 ion_index = 1 ! By default, only first ion species is considered

 IF ((node%neutral.EQ.n_neutral).AND.(node%indx_coll.EQ.indx_coll)) THEN

    k = node%number

    i = INT(ion(ion_index)%part(k)%X)
    j = INT(ion(ion_index)%part(k)%Y)
    IF (ion(ion_index)%part(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
    IF (ion(ion_index)%part(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1
    
    pos_i_j     = i + j * n3 + n2
    pos_ip1_j   = pos_i_j + 1
    pos_i_jp1   = pos_i_j + n3
    pos_ip1_jp1 = pos_i_jp1 + 1

    ax_ip1 = REAL(ion(ion_index)%part(k)%X) - REAL(i)
    ax_i   = 1.0 - ax_ip1

    ay_jp1 = REAL(ion(ion_index)%part(k)%Y) - REAL(j)
    ay_j   = 1.0 - ay_jp1

    vij   = ax_i   * ay_j
    vip1j = ax_ip1 * ay_j
    vijp1 = ax_i   * ay_jp1

    rbufer(pos_i_j)     = rbufer(pos_i_j)     + vij                         !ax_i   * ay_j
    rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   + vip1j                       !ax_ip1 * ay_j
    rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   + vijp1                       !ax_i   * ay_jp1
    rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) + 1.0 - vij - vip1j - vijp1   !ax_ip1 * ay_jp1

 END IF

 IF (ASSOCIATED(node%Larger)) CALL Transfer_collisions_from_stored_list_ions(node%Larger, n_neutral, indx_coll, bufsize, n1, n2, n3, rbufer)

 IF (ASSOCIATED(node%Smaller)) CALL Transfer_collisions_from_stored_list_ions(node%Smaller, n_neutral, indx_coll, bufsize, n1, n2, n3, rbufer)

 RETURN

END SUBROUTINE Transfer_collisions_from_stored_list_ions

!---------------------------------------------------
! this subroutine kills the nodes of the binary tree
RECURSIVE SUBROUTINE Node_Killer_ion(node)

 USE MCCollisions
 IMPLICIT NONE
 
 TYPE (binary_tree), POINTER :: node
 INTEGER DEALLOC_ERR

 IF (ASSOCIATED(node%Larger))  CALL Node_Killer_ion(node%Larger)
 IF (ASSOCIATED(node%Smaller)) CALL Node_Killer_ion(node%Smaller)

 DEALLOCATE(node, STAT=DEALLOC_ERR)

 RETURN

END SUBROUTINE Node_Killer_ion 

!-------------------------------------------------------------------------------------------

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
