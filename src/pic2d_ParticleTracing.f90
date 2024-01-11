!-----------------------------------------------
!
SUBROUTINE SET_PARTICLE_TRACING_e

  USE ParticleTracing
  USE Checkpoints
  USE ParallelOperationValues
  USE CurrentProblemValues,  ONLY : delta_t_s, delta_x_m, global_maximal_i, global_maximal_j

!   USE MPI

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'

  LOGICAL exists

  CHARACTER(1) buf
  INTEGER content_flag(6:18)
  INTEGER n
  REAL(8) dtemp
  INTEGER ios
  REAL tempx, tempy
  INTEGER tempnpart, flag1, flag2

! default values
  Tcntr_e_tracing_start = -1
  T_cntr_save_traced_electrons = -1

  INQUIRE (FILE = 'init_trace_electrons.dat', EXIST = exists)

  IF (.NOT.exists) THEN
     IF (Rank_of_process.EQ.0) PRINT '("###### File init_trace_electrons.dat not found, NO TRACING of electrons in this simulation ######")'
     RETURN
  END IF
     
  OPEN (9, FILE = 'init_trace_electrons.dat')

  READ (9, '(A1)') buf ! # default record content: tag, x, y, vx, and vy are saved always
  READ (9, '(A1)') buf ! # additional record content: to save/omit values below enter 1/0 under each value name
  READ (9, '(A1)') buf ! # vz    irrEx   irrEy   extEz   solEx   solEy   solEz   solBx   solBy   solBz   extBx   extBy   extBz
  content_flag = 0

  READ (9, *, iostat = ios) &
       &  content_flag(6), &   ! vz
       &  content_flag(7), &   ! irrEx
       &  content_flag(8), &   ! irrEy
       &  content_flag(9), &   ! extEz
       &  content_flag(10), &  ! solEx
       &  content_flag(11), &  ! solEy
       &  content_flag(12), &  ! solEz
       &  content_flag(13), &  ! solBx
       &  content_flag(14), &  ! solBy
       &  content_flag(15), &  ! solBz
       &  content_flag(16), &  ! extBx
       &  content_flag(17), &  ! extBy
       &  content_flag(18)     ! extBz

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("###### Error in line 4 of init_trace_electrons.dat, NO TRACING of electrons in this simulation ######")'
     RETURN
  END IF

  save_e_vz = .FALSE.
  save_e_irrEx = .FALSE.
  save_e_irrEy = .FALSE.
  save_e_extEz = .FALSE.
  save_e_solEx = .FALSE.
  save_e_solEy = .FALSE.
  save_e_solEz = .FALSE.
  save_e_solBx = .FALSE.
  save_e_solBy = .FALSE.
  save_e_solBz = .FALSE.
  save_e_extBx = .FALSE.
  save_e_extBy = .FALSE.
  save_e_extBz = .FALSE.

  IF (content_flag(6) .NE.0) save_e_vz    = .TRUE.
  IF (content_flag(7) .NE.0) save_e_irrEx = .TRUE.
  IF (content_flag(8) .NE.0) save_e_irrEy = .TRUE.
  IF (content_flag(9) .NE.0) save_e_extEz = .TRUE.
  IF (content_flag(10).NE.0) save_e_solEx = .TRUE.
  IF (content_flag(11).NE.0) save_e_solEy = .TRUE.
  IF (content_flag(12).NE.0) save_e_solEz = .TRUE.
  IF (content_flag(13).NE.0) save_e_solBx = .TRUE.
  IF (content_flag(14).NE.0) save_e_solBy = .TRUE.
  IF (content_flag(15).NE.0) save_e_solBz = .TRUE.
  IF (content_flag(16).NE.0) save_e_extBx = .TRUE.
  IF (content_flag(17).NE.0) save_e_extBy = .TRUE.
  IF (content_flag(18).NE.0) save_e_extBz = .TRUE.

  e_record_length = 5
  DO n = 6, 18
     IF (content_flag(n).NE.0) e_record_length = e_record_length+1
  END DO

  READ (9, '(A1)') buf ! # time when tracing begins, in nanoseconds if >=0, timesteps if <0
  READ (9, *, iostat = ios)  dtemp

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("###### Error in line 6 of init_trace_electrons.dat, NO TRACING of electrons in this simulation ######")'
     RETURN
  END IF

  IF (dtemp.ge.0.0) THEN
     Tcntr_e_tracing_start = INT(dtemp * 1.0d-9 / delta_t_s)
  ELSE
     Tcntr_e_tracing_start = INT(ABS(dtemp))
  END IF

  READ (9, '(A1)') buf ! # time when tracing ends, in nanoseconds if >=0, timesteps if <0
  READ (9, *, iostat = ios)  dtemp

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("###### Error in line 8 of init_trace_electrons.dat, NO TRACING of electrons in this simulation ######")'
     RETURN
  END IF

  IF (dtemp.ge.0.0) THEN
     Tcntr_e_tracing_end = INT(dtemp * 1.0d-9 / delta_t_s)
  ELSE
     Tcntr_e_tracing_end = INT(ABS(dtemp))
  END IF

  READ (9, '(A1)') buf ! # interval between saving particles to file, in nanoseconds if >=0, timesteps if <0 (0 to save each timestep)
  READ (9, *, iostat = ios)  dtemp

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("###### Error in line 10 of init_trace_electrons.dat, NO TRACING of electrons in this simulation ######")'
     RETURN
  END IF

  IF (dtemp.ge.0.0) THEN
     Tcntr_e_tracing_step = INT(dtemp * 1.0d-9 / delta_t_s)
  ELSE
     Tcntr_e_tracing_step = INT(ABS(dtemp))
  END IF
  Tcntr_e_tracing_step = MAX(1, Tcntr_e_tracing_step)

  READ (9, '(A1)') buf ! # definitions for particle selection:
  READ (9, '(A1)') buf ! # one can specify up to 50 locations (cells) where particles to be traced will be selected
  READ (9, '(A1)') buf ! # in each location one can select up to 20 particles
  READ (9, '(A1)') buf ! # location coordinates x/y represent the left bottom corner of the cell containing the particles to be traced
  READ (9, '(A1)') buf ! # the coordinate is in meters if >=0, the node number if <0
  READ (9, '(A1)') buf ! # rule 1 for particle selection: to choose particles with VX<0, any VX, or VX>0 set flag 1 to -1/0/1 respectively
  READ (9, '(A1)') buf ! # rule 2 for particle selection: to choose particles with VY<0, any VY, or VY>0 set flag 2 to -1/0/1 respectively

  READ (9, '(A1)') buf ! # below, for each location type location x/y, number of particles to be traced, flag 1, flag 2

  n = 0
  DO
     READ (9, *, iostat = ios) tempx, tempy, tempnpart, flag1, flag2

     IF (ios.NE.0) EXIT

     n = n + 1

     IF (tempx.GE.0.0) THEN
        cell_e_tracing(n)%icell = INT(tempx / delta_x_m)
     ELSE
        cell_e_tracing(n)%icell = INT(ABS(tempx))
     END IF

     IF (cell_e_tracing(n)%icell.GE.global_maximal_i) THEN 
        n = n-1
        CYCLE
     END IF

     IF (tempy.GE.0.0) THEN
        cell_e_tracing(n)%jcell = INT(tempy / delta_x_m)
     ELSE
        cell_e_tracing(n)%jcell = INT(ABS(tempy))
     END IF

     IF (cell_e_tracing(n)%jcell.GE.global_maximal_j) THEN 
        n = n-1
        CYCLE
     END IF

     cell_e_tracing(n)%N_to_trace = MAX(0, MIN(maximal_N_part_to_trace_in_cell, tempnpart))
     IF (cell_e_tracing(n)%N_to_trace.EQ.0) THEN
        n = n-1
        CYCLE
     END IF

     cell_e_tracing(n)%choose_negative_VX = .FALSE.
     cell_e_tracing(n)%choose_positive_VX = .FALSE.
     cell_e_tracing(n)%choose_negative_VY = .FALSE.
     cell_e_tracing(n)%choose_positive_VY = .FALSE.

     IF (flag1.LT.0) cell_e_tracing(n)%choose_negative_VX = .TRUE.
     IF (flag1.GT.0) cell_e_tracing(n)%choose_positive_VX = .TRUE.

     IF (flag2.LT.0) cell_e_tracing(n)%choose_negative_VY = .TRUE.
     IF (flag2.GT.0) cell_e_tracing(n)%choose_positive_VY = .TRUE.

     IF (n.EQ.maximal_N_cells_tracing) EXIT
     
  END DO

  N_of_cells_e_tracing = n

  CLOSE (9, STATUS = 'KEEP')

  IF (N_of_cells_e_tracing.EQ.0) THEN
! no valid cells provided
     Tcntr_e_tracing_start = -1   ! restore the default value which cancels tracing
     IF (Rank_of_process.EQ.0) THEN
        PRINT '("### Electron tracing disabled because no valid locations for electron tracing was provided")'
     END IF
     RETURN
  END IF

  tempnpart = 0
  DO n = 1, N_of_cells_e_tracing
     tempnpart = tempnpart + cell_e_tracing(n)%N_to_trace
  END DO
  IF (tempnpart.GT.0) THEN
     IF (Rank_of_process.EQ.0) THEN
        PRINT '("### Electron tracing enabled :: requested total ",i4," particles in ",i2," locations, start/end/step = ",3(2x,i8))', &
             & tempnpart, N_of_cells_e_tracing, Tcntr_e_tracing_start, Tcntr_e_tracing_end, Tcntr_e_tracing_step
        PRINT '("### Expected size of file electron_particles_vst.dat is ",f10.1," kB or less")', REAL(e_record_length * 4 * tempnpart) * REAL((Tcntr_e_tracing_end - Tcntr_e_tracing_start)/Tcntr_e_tracing_step) / 1024.0
     END IF
  ELSE
! no particles in valid cells requested
     Tcntr_e_tracing_start = -1   ! restore the default value which cancels tracing
     IF (Rank_of_process.EQ.0) THEN
        PRINT '("### Electron tracing disabled because no particles in valid locations were requested")'
     END IF
  END IF

END SUBROUTINE SET_PARTICLE_TRACING_e

!-----------------------------------------------
!
SUBROUTINE START_PARTICLE_TRACING_e

  USE ParticleTracing
  USE Checkpoints
  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : T_cntr
  USE ElectronParticles
  USE ClusterAndItsBoundaries

!   USE MPI

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr, errcode
  INTEGER stattus(MPI_STATUS_SIZE)

  INTEGER k
  INTEGER i_limit_right, j_limit_above
  INTEGER ibufer(2)

  INTEGER ALLOC_ERR
  INTEGER, ALLOCATABLE :: N_particles_cell_proc(:)
  INTEGER, ALLOCATABLE :: N_distribute_per_process(:)
  INTEGER, ALLOCATABLE :: start_tag(:)

  INTEGER n
  INTEGER my_N_particles_cell

  INTEGER m
  INTEGER N_total_in_cell
  INTEGER N_to_distribute

  INTEGER my_N_to_distribute
  INTEGER my_start_tag

  INTEGER count
 
  IF (T_cntr.NE.Tcntr_e_tracing_start) RETURN  ! ensures that the procedure is either applied once or not applied at all

print '("proc ",i4," entered START_PARTICLE_TRACING_e")', Rank_of_process

  T_cntr_save_traced_electrons = T_cntr

! clear the tracing increment (if simulation starts with a checkpoint created while another tracing was applied)
  DO k = 1, N_electrons
     IF (electron(k)%tag.GT.9999) THEN
        electron(k)%tag = electron(k)%tag - (electron(k)%tag / 10000) * 10000   ! 10,000 -> 10,000 - 10,000 = 0
                                                                                ! 10,001 -> 10,001 - 10,000 = 1, etc.
     END IF
  END DO

  IF (Rank_cluster.EQ.0) THEN
     i_limit_right = c_indx_x_max-2
     IF (Rank_of_master_right.LT.0) i_limit_right = c_indx_x_max-1   ! wall on the right

     j_limit_above = c_indx_y_max-2
     IF (Rank_of_master_above.LT.0) j_limit_above = c_indx_y_max-1   ! wall above

     ibufer(1) = i_limit_right
     ibufer(2) = j_limit_above
     
     CALL MPI_BCAST(ibufer, 2, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
  ELSE
     CALL MPI_BCAST(ibufer, 2, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
     i_limit_right = ibufer(1)
     j_limit_above = ibufer(2)
  END IF

  IF (Rank_cluster.EQ.0) THEN
     ALLOCATE(N_distribute_per_process(0:N_processes_cluster-1), STAT = ALLOC_ERR)
     ALLOCATE(   N_particles_cell_proc(0:N_processes_cluster-1), STAT = ALLOC_ERR)
     ALLOCATE(               start_tag(0:N_processes_cluster-1), STAT = ALLOC_ERR)
  END IF

  DO n = 1, N_of_cells_e_tracing

! print '(2x,i4,2x,i2,4x,3(2x,i4),4x,3(2x,i4))', Rank_of_process, n, c_indx_x_min, cell_e_tracing(n)%icell, i_limit_right, c_indx_y_min, cell_e_tracing(n)%jcell, j_limit_above

! skip locations outside the cluster's domain
     IF (cell_e_tracing(n)%icell.LT.c_indx_x_min) CYCLE
     IF (cell_e_tracing(n)%icell.GT.i_limit_right) CYCLE
     IF (cell_e_tracing(n)%jcell.LT.c_indx_y_min) CYCLE
     IF (cell_e_tracing(n)%jcell.GT.j_limit_above) CYCLE

! find the number of own particles in the given cell
! note that the particles must satisfy all the criteria
! and they may not be already chosen
     my_N_particles_cell = 0
     DO k = 1, N_electrons

        IF (INT(electron(k)%X).NE.cell_e_tracing(n)%icell) CYCLE
        IF (INT(electron(k)%Y).NE.cell_e_tracing(n)%jcell) CYCLE
        IF (electron(k)%tag.GT.9999) CYCLE
        IF (cell_e_tracing(n)%choose_negative_VX) THEN
           IF (electron(k)%VX.GT.0.0_8) CYCLE
        ELSE IF (cell_e_tracing(n)%choose_positive_VX) THEN
           IF (electron(k)%VX.LT.0.0_8) CYCLE
        END IF
        IF (cell_e_tracing(n)%choose_negative_VY) THEN
           IF (electron(k)%VY.GT.0.0_8) CYCLE
        ELSE IF (cell_e_tracing(n)%choose_positive_VY) THEN
           IF (electron(k)%VY.LT.0.0_8) CYCLE
        END IF

        my_N_particles_cell = my_N_particles_cell + 1

     END DO

print '("proc ",i4," location ",i2," my_N_particles_cell ",i4)', Rank_of_process, n, my_N_particles_cell

     IF (Rank_cluster.EQ.0) THEN      

        N_particles_cell_proc(0) = my_N_particles_cell

! get the number of particles in the given cell in each process       
        DO m = 1, N_processes_cluster-1
           CALL MPI_RECV(N_particles_cell_proc(m), 1, MPI_INTEGER, m, m, COMM_CLUSTER, stattus, ierr)
        END DO     

! prepare the number of particles to be traced for each process belonging to the cluster
        N_total_in_cell = N_particles_cell_proc(0)
        DO m = 1, N_processes_cluster-1
           N_total_in_cell = N_total_in_cell + N_particles_cell_proc(m)
        END DO

        N_to_distribute = MIN(N_total_in_cell, cell_e_tracing(n)%N_to_trace)

        N_distribute_per_process = 0   ! array
        DO WHILE (N_to_distribute.GT.0)
           DO m = 0, N_processes_cluster-1
              IF (N_particles_cell_proc(m).GT.0) THEN
                 N_distribute_per_process(m) = N_distribute_per_process(m) + 1
                 N_particles_cell_proc(m) = N_particles_cell_proc(m) - 1
                 N_to_distribute = N_to_distribute-1
              END IF
              IF (N_to_distribute.EQ.0) EXIT
           END DO
        END DO

! consistency check    
        N_to_distribute = 0
        DO m = 0, N_processes_cluster-1
           N_to_distribute = N_to_distribute + N_distribute_per_process(m)
        END DO
        
        IF (N_to_distribute.NE.MIN(N_total_in_cell, cell_e_tracing(n)%N_to_trace)) THEN
! ERROR, report
           errcode = 1020
           CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        END IF

        start_tag(0) = 1 + (n-1) * maximal_N_part_to_trace_in_cell   ! n is the number of the cell where tracing particles are selected
        DO m = 1, N_processes_cluster-1
           start_tag(m) = start_tag(m-1) + N_distribute_per_process(m-1)
        END DO

        DO m = 1, N_processes_cluster-1
           ibufer(1) = N_distribute_per_process(m)
           ibufer(2) = start_tag(m)
           CALL MPI_SEND(ibufer, 2, MPI_INTEGER, m, 0, COMM_CLUSTER, ierr)
        END DO

        my_N_to_distribute = N_distribute_per_process(0)
        my_start_tag = start_tag(0)

     ELSE   !###   IF (Rank_cluster.EQ.0) THEN

        CALL MPI_SEND(my_N_particles_cell, 1, MPI_INTEGER, 0, Rank_cluster, COMM_CLUSTER, ierr)

        CALL MPI_RECV(ibufer, 2, MPI_INTEGER, 0, 0, COMM_CLUSTER, stattus, ierr)
        my_N_to_distribute = ibufer(1)
        my_start_tag = ibufer(2)

        IF (my_N_to_distribute.GT.my_N_particles_cell) THEN
! ERROR, report
           errcode = 1030
           CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        END IF

     END IF   !###   IF (Rank_cluster.EQ.0) THEN

print '("proc ",i4," my_N_to_distribute = ",i4)', Rank_of_process, my_N_to_distribute

! label particles to be traced
     count = 0
     DO k = 1, N_electrons

        IF (INT(electron(k)%X).NE.cell_e_tracing(n)%icell) CYCLE
        IF (INT(electron(k)%Y).NE.cell_e_tracing(n)%jcell) CYCLE
        IF (electron(k)%tag.GT.9999) CYCLE
        IF (cell_e_tracing(n)%choose_negative_VX) THEN
           IF (electron(k)%VX.GT.0.0_8) CYCLE
        ELSE IF (cell_e_tracing(n)%choose_positive_VX) THEN
           IF (electron(k)%VX.LT.0.0_8) CYCLE
        END IF
        IF (cell_e_tracing(n)%choose_negative_VY) THEN
           IF (electron(k)%VY.GT.0.0_8) CYCLE
        ELSE IF (cell_e_tracing(n)%choose_positive_VY) THEN
           IF (electron(k)%VY.LT.0.0_8) CYCLE
        END IF
! since we are here, the particle will be selected
        electron(k)%tag = electron(k)%tag + (my_start_tag + count) * 10000
        count = count + 1
        IF ( count.EQ.my_N_to_distribute ) EXIT

     END DO

  END DO   !###   DO n = 1, N_of_cells_e_tracing

! cleanup
  IF (ALLOCATED(N_distribute_per_process)) DEALLOCATE(N_distribute_per_process, STAT = ALLOC_ERR)
  IF (ALLOCATED(N_particles_cell_proc)) DEALLOCATE(N_particles_cell_proc, STAT = ALLOC_ERR)
  IF (ALLOCATED(start_tag)) DEALLOCATE(start_tag, STAT = ALLOC_ERR)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

END SUBROUTINE START_PARTICLE_TRACING_e

!----------------------------------------------
!
SUBROUTINE SAVE_TRACED_PARTICLES_e

  USE ParticleTracing
  USE ParallelOperationValues
  USE ElectronParticles
  USE CurrentProblemValues, ONLY : delta_x_m, delta_t_s, T_cntr, V_scale_ms, E_scale_Vm, B_scale_T, EX, EY, zero 
  USE SolenoidalFields
  USE ClusterAndItsBoundaries

!   USE MPI

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr, errcode
  INTEGER stattus(MPI_STATUS_SIZE)

  INTEGER bufsize
  REAL, ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  INTEGER my_N_to_save
  INTEGER pos, shift
  INTEGER k
  INTEGER i, j, ii, jj
  REAL(8) X_ii, Y_jj
  REAL(8) aii_ip1, aii_i
  REAL(8) ajj_jp1, ajj_j
  REAL(8) a_ip1, a_i
  REAL(8) a_jp1, a_j
  REAL(8) temp

  INTEGER, ALLOCATABLE :: N_to_save(:)

  INTEGER m
  INTEGER bigbufsize
  REAL, ALLOCATABLE :: bigbufer(:)

  INTEGER kk

! functions
  REAL(8) Ez, Bx, By, Bz

  IF (T_cntr.NE.T_cntr_save_traced_electrons) RETURN

  bufsize = e_record_length * maximal_N_cells_tracing * maximal_N_part_to_trace_in_cell  ! assume maximal size
  ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

  my_N_to_save = 0
  pos=1

  DO k = 1, N_electrons

     IF (electron(k)%tag.LT.10000) CYCLE

! default content
     rbufer(pos)   = REAL(electron(k)%tag) + 0.0001
     rbufer(pos+1) = REAL(electron(k)%X * delta_x_m)
     rbufer(pos+2) = REAL(electron(k)%Y * delta_x_m)
     rbufer(pos+3) = REAL(electron(k)%VX * V_scale_ms)
     rbufer(pos+4) = REAL(electron(k)%VY * V_scale_ms)

! additional content
     shift=5
     IF (save_e_vz) THEN
        rbufer(pos+shift) = REAL(electron(k)%VZ * V_scale_ms)
        shift = shift+1
     END IF

     i = FLOOR(electron(k)%X)
     j = FLOOR(electron(k)%Y)

     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

     ii = FLOOR(electron(k)%X - 0.5_8)   ! sol_BY/EX x-indices limits: c_indx_x_min - 1 : c_indx_x_max
     jj = FLOOR(electron(k)%Y - 0.5_8)   ! sol_BX/EY y-indices limits: c_indx_y_min - 1 : c_indx_y_max

     X_ii = DBLE(ii) + 0.5_8             ! x-coordinate of sol_BY node ii
     Y_jj = DBLE(jj) + 0.5_8             ! y-coordinate of sol_BX node jj

     aii_ip1 = electron(k)%X - X_ii
     aii_i   = 1.0_8 - aii_ip1

     ajj_jp1 = electron(k)%Y - Y_jj
     ajj_j   = 1.0_8 - ajj_jp1

     a_ip1 =  electron(k)%X - DBLE(i)
     a_i   = 1.0_8 - a_ip1

     a_jp1 =  electron(k)%Y - DBLE(j)
     a_j   = 1.0_8 - a_jp1

     IF (save_e_irrEx) THEN
        temp = EX(ii,j) * aii_i * a_j + EX(ii+1,j) * aii_ip1 * a_j + EX(ii,j+1) * aii_i * a_jp1 + EX(ii+1,j+1) * aii_ip1 * a_jp1
        rbufer(pos+shift) = REAL(temp * E_scale_Vm)
        shift = shift+1
     END IF

     IF (save_e_irrEy) THEN
        temp = EY(i,jj) * a_i * ajj_j + EY(i+1, jj) * a_ip1 * ajj_j + EY(i,jj+1) * a_i * ajj_jp1 + EY(i+1,jj+1) * a_ip1 * ajj_jp1
        rbufer(pos+shift) = REAL(temp * E_scale_Vm)
        shift = shift+1
     END IF

     IF (save_e_extEz) THEN
        rbufer(pos+shift) = REAL(Ez(electron(k)%X, electron(k)%Y) * E_scale_Vm)
        shift = shift+1
     END IF

     IF (save_e_solEx) THEN
        temp = zero !sol_EX(ii,j) * aii_i * a_j + sol_EX(ii+1,j) * aii_ip1 * a_j + sol_EX(ii,j+1) * aii_i * a_jp1 + sol_EX(ii+1,j+1) * aii_ip1 * a_jp1
        rbufer(pos+shift) = REAL(temp * E_scale_Vm)
        shift = shift+1
     END IF

     IF (save_e_solEy) THEN
        temp = zero !sol_EY(i,jj) * a_i * ajj_j + sol_EY(i+1,jj) * a_ip1 * ajj_j + sol_EY(i,jj+1) * a_i * ajj_jp1 + sol_EY(i+1,jj+1) * a_ip1 * ajj_jp1
        rbufer(pos+shift) = REAL(temp * E_scale_Vm)
        shift = shift+1
     END IF

     IF (save_e_solEz) THEN
        temp = zero !sol_EZ(i,j) * a_i * a_j + sol_EZ(i+1,j) * a_ip1 * a_j + sol_EZ(i,j+1) * a_i * a_jp1 + sol_EZ(i+1,j+1) * a_ip1 * a_jp1
        rbufer(pos+shift) = REAL(temp * E_scale_Vm)
        shift = shift+1
     END IF

     IF (save_e_solBx) THEN
        temp = zero !sol_BX(i,jj) * a_i * ajj_j + sol_BX(i+1,jj) * a_ip1 * ajj_j + sol_BX(i,jj+1) * a_i * ajj_jp1 + sol_BX(i+1,jj+1) * a_ip1 * ajj_jp1
        rbufer(pos+shift) = REAL(temp * B_scale_T)
        shift = shift+1
     END IF

     IF (save_e_solBy) THEN
        temp = zero !sol_BY(ii,j) * aii_i * a_j + sol_BY(ii+1,j) * aii_ip1 * a_j + sol_BY(ii,j+1) * aii_i * a_jp1 + sol_BY(ii+1,j+1) * aii_ip1 * a_jp1
        rbufer(pos+shift) = REAL(temp * B_scale_T)
        shift = shift+1
     END IF

     IF (save_e_solBz) THEN
        temp = zero !sol_BZ(ii,jj) * aii_i * ajj_j + sol_BZ(ii+1,jj) * aii_ip1 * ajj_j + sol_BZ(ii,jj+1) * aii_i * ajj_jp1 + sol_BZ(ii+1,jj+1) * aii_ip1 * ajj_jp1
        rbufer(pos+shift) = REAL(temp * B_scale_T)
        shift = shift+1
     END IF

     IF (save_e_extBx) THEN
        rbufer(pos+shift) = REAL(Bx(electron(k)%X, electron(k)%Y) * B_scale_T)
        shift = shift+1
     END IF

     IF (save_e_extBy) THEN
        rbufer(pos+shift) = REAL(By(electron(k)%X, electron(k)%Y) * B_scale_T)
        shift = shift+1
     END IF

     IF (save_e_extBz) THEN
        rbufer(pos+shift) = REAL(Bz(electron(k)%X, electron(k)%Y) * B_scale_T)
        shift = shift+1
     END IF

! paranoidal check
     IF (shift.NE.e_record_length) THEN
        errcode=1010
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
     END IF

     pos = pos + e_record_length

     my_N_to_save = my_N_to_save + 1   ! count particles

  END DO   !###  DO k = 1, N_electrons

  IF (Rank_of_process.EQ.0) THEN

     ALLOCATE(N_to_save(1:N_of_processes-1), STAT = ALLOC_ERR)

     bigbufsize = my_N_to_save
     DO m = 1, N_of_processes-1
        CALL MPI_RECV(N_to_save(m), 1, MPI_INTEGER, m, 0, MPI_COMM_WORLD, stattus, ierr)
        bigbufsize = bigbufsize + N_to_save(m)
     END DO

! print '(16(2x,i4))', my_N_to_save, N_to_save

     bigbufsize = bigbufsize * e_record_length
     ALLOCATE(bigbufer(1:bigbufsize), STAT = ALLOC_ERR)
     bufsize = e_record_length * my_N_to_save
     IF (bufsize.GT.0) bigbufer(1:bufsize) = rbufer(1:bufsize)
     
     pos = bufsize + 1
     DO m = 1, N_of_processes-1
        IF (N_to_save(m).EQ.0) CYCLE
        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        bufsize = e_record_length * N_to_save(m)
        ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, m, m, MPI_COMM_WORLD, stattus, ierr)
        bigbufer(pos:pos+bufsize-1) = rbufer(1:bufsize)
        pos = pos + bufsize
     END DO

     IF (T_cntr.EQ.Tcntr_e_tracing_start) THEN
        OPEN (19, FILE = 'electron_particles_vst.dat', STATUS = 'REPLACE', ACCESS = 'STREAM')
        WRITE (19) save_e_vz
        WRITE (19) save_e_irrEx
        WRITE (19) save_e_irrEy
        WRITE (19) save_e_extEz
        WRITE (19) save_e_solEx
        WRITE (19) save_e_solEy
        WRITE (19) save_e_solEz
        WRITE (19) save_e_solBx
        WRITE (19) save_e_solBy
        WRITE (19) save_e_solBz
        WRITE (19) save_e_extBx
        WRITE (19) save_e_extBy
        WRITE (19) save_e_extBz
     ELSE
        OPEN (19, FILE = 'electron_particles_vst.dat', POSITION = 'APPEND', ACCESS = 'STREAM')
     END IF

     bigbufsize = bigbufsize / e_record_length

     WRITE (19) T_cntr
     WRITE (19) REAL(T_cntr * delta_t_s * 1.0d9)   ! time in nanoseconds
     WRITE (19) bigbufsize       ! number of particles saved at this step
     pos=1
     DO k = 1, bigbufsize
        WRITE (19) INT(bigbufer(pos))   ! tag
        DO kk = 1, e_record_length-1
           WRITE (19) bigbufer(pos+kk)  ! x, y, vx, vy, &&& vz, Ex, Ey, extEz, solEx, solEy, solEz, solBx, solBy, solBz, extBx, extBy, extBz
        END DO
        pos = pos + e_record_length
     END DO

     CLOSE (19, STATUS = 'KEEP')

! cleanup
     IF (ALLOCATED(bigbufer)) DEALLOCATE(bigbufer, STAT = ALLOC_ERR)

  ELSE   !###   IF (Rank_of_process.EQ.0) THEN

     CALL MPI_SEND(my_N_to_save, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)

     IF (my_N_to_save.GT.0) THEN
        bufsize = e_record_length * my_N_to_save  ! actual size of saved particles data
        CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_REAL, 0, Rank_of_process, MPI_COMM_WORLD, ierr)
     END IF

  END IF   !###   IF (Rank_of_process.EQ.0) THEN

! cleanup
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)

! update next time to save particles
  T_cntr_save_traced_electrons = T_cntr_save_traced_electrons + Tcntr_e_tracing_step

  IF (T_cntr_save_traced_electrons.GT.Tcntr_e_tracing_end) THEN
! tracing is over
     Tcntr_e_tracing_start = -1
     T_cntr_save_traced_electrons = -1
! clear the tracing increment from tags 
     DO k = 1, N_electrons
        IF (electron(k)%tag.GT.9999) THEN
           electron(k)%tag = electron(k)%tag - (electron(k)%tag / 10000) * 10000   ! 10,000 -> 10,000 - 10,000 = 0
                                                                                   ! 10,001 -> 10,001 - 10,000 = 1, etc.
        END IF
     END DO
  END IF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

END SUBROUTINE SAVE_TRACED_PARTICLES_e
