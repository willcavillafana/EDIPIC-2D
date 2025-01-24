
!===================================================================================================
SUBROUTINE INITIATE_ION_NEUTRAL_COLLISIONS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : kB_JK, e_Cl, N_max_vel, T_e_eV, pi  !, amu_kg, m_e_kg, V_scale_ms
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
  INTEGER s, cnt
  real(8) Fel, dtheta, theta, xi, dxi


  INTEGER ALLOC_ERR

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  no_in_collisions = .TRUE.

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
   else
      open(9, file = init_rcx_param_filename)
      read(9, *) buf
      read(9, *) s, neutral(n)%sigma_rcx_m2_1eV, neutral(n)%alpha_rcx
      close(9, status = 'keep')
      if ((s.gt.0).and.(s.le.N_spec)) then
         no_in_collisions = .false.
         neutral(n)%rcx_on = .true.     
         neutral(n)%rcx_ion_species_index = s
      else
         neutral(n)%rcx_on = .false.
         neutral(n)%rcx_ion_species_index = 0
      end if   
   end if

   if ((neutral(n)%rcx_on).AND.(neutral(n)%in_cols_on).AND.(Rank_of_process.eq.0)) PRINT '("### Both the old model for RCX as well as the new one for polarization scattering + RCX are switched on. Only the new model will be used. ###")'
   if (neutral(n)%in_cols_on) no_in_collisions = .false. ! if we are using the new in-col model, we still need the preparation below for diagnostics

  END DO
  

  IF (no_in_collisions) RETURN
    
  allocate (collision_rcx(1:N_spec), stat = alloc_err)
  do s = 1, N_spec
     collision_rcx(s)%rcx_on = .false.
  end do

  do n = 1, N_neutral_spec
     if ((neutral(n)%rcx_on).OR.(neutral(n)%in_cols_on)) then
        s = neutral(n)%rcx_ion_species_index
        collision_rcx(s)%rcx_on = .true.
        collision_rcx(s)%neutral_species_index = n
        collision_rcx(s)%vfactor = SQRT(neutral(n)%T_K * kB_JK / (T_e_eV * e_Cl * Ms(s))) / DBLE(N_max_vel)
!??        collision_rcx(s)%factor_eV = Ms(s) * energy_factor_eV !0.5_8 * m_e_kg * Ms(s) * V_scale_ms**2 / e_Cl
     end if
  end do

  do n = 1, N_neutral_spec
   if (neutral(n)%in_cols_on) then ! Precalculate the collision integral
      ALLOCATE(neutral(n)%Fel_precalc(0:(50000-1)), STAT=ALLOC_ERR)
      dtheta = 0.0001
      dxi = 1.0/50000.0
      cnt = 0
      do xi = 0.0_8, 1.0-dxi, dxi
         Fel = 0.0
         do theta = 0.0_8, pi/2.0, dtheta ! incomplete elliptic integral of first kind 
            Fel = Fel + 1.0/sqrt((1.0 - (xi**2)*sin(theta)**2))*dtheta 
         end do
         !PRINT *, "Cnt: ", cnt  
         neutral(n)%Fel_precalc(cnt) = Fel
         cnt = cnt + 1
      end do

   end if
  end do

  CALL calculate_thermal_cx_probab

END SUBROUTINE INITIATE_ION_NEUTRAL_COLLISIONS

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
  INTEGER i_dummy, ios

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (Rank_of_process.NE.0) RETURN

  IF (no_in_collisions) RETURN

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

            DO
               READ (21, '(2x,i9,2x,i9,2x,i8)', iostat = ios) i_dummy
               IF (ios.NE.0) EXIT
               IF (i_dummy.GE.Start_T_cntr) EXIT
            ENDDO
            BACKSPACE(21)            
         !   DO i = 1, Start_T_cntr / N_subcycles            ! these files are updated at every ion timestep
         !      READ (21, '(2x,i9,2x,i9,2x,i8)') i_dummy
         !   END DO
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
  USE CurrentProblemValues, ONLY : T_cntr, debug_level
  USE AvgSnapshots, ONLY: avg_flux_and_history, current_avgsnap

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
  INTEGER :: local_debug_level
  INTEGER :: avg_output_flag

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  local_debug_level = 3

  IF (no_in_collisions) RETURN

   IF (avg_flux_and_history) THEN
      CALL DETERMINE_AVG_DATA_CREATION(avg_output_flag)
      IF (avg_output_flag==0) RETURN
   END IF  

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

   IF (debug_level>=local_debug_level) PRINT '("Total charge exchange collisions for all ions species :: ",10(2x,i6))', ibufer_receive(N_spec+1:N_spec+N_spec)

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
