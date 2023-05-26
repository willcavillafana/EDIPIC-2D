!---------------------------------------
!
SUBROUTINE PERFORM_ELECTRON_COULOMB_SCATTERING

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  USE rng_wrapper
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER k
  INTEGER i, j
  REAL(8) ax_ip1, ax_i, ay_jp1, ay_j, a00, a01, a10, a11
  REAL(8) n_loc, vx_loc, vy_loc, vz_loc, v2_loc, T_loc_eV
  REAL(8) vel_factor, dens_ratio, temp_ratio
  REAL(8) Vx, Vy, Vz, Vx_sample, Vy_sample, Vz_sample, Vframe_x, Vframe_y, Vframe_z
  REAL(8) Ux, Uy, Uz, Usq, Usq_scaled, U_scaled ! relative velocity
  REAL(8) Wx, Wy, Wz, Wsq, Wsq_scaled
  REAL(8) V, V_xy, a, b, Vx_s, Vy_s, Vz_s
  REAL(8) L_ee, A_norm, s_coll, CosKsi, SinKsi, Fi, CosFi, SinFi, R

  INTEGER m, bufsize

!  write (*, *) "enter Coulomb", Rank_of_process
  IF (.NOT.Coulomb_flag) RETURN

  IF (cluster_rank_key.EQ.0) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO m = 0, 4 
           acc_rho_e(m, c_indx_x_min:c_indx_x_max, j) = acc_rho_e(m, c_indx_x_min:c_indx_x_max, j) / DBLE(N_subcycles)  
        END DO
     END DO   

!     bufsize = 5 * (c_indx_x_max - c_indx_x_min + 1) * (c_indx_y_max - c_indx_y_min + 1)
!     CALL MPI_BCAST(acc_rho_e, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)
!  ELSE 
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
  bufsize = 5 * (c_indx_x_max - c_indx_x_min + 1) * (c_indx_y_max - c_indx_y_min + 1)
  CALL MPI_BCAST(acc_rho_e, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)     
!  END IF   
! all processes:
  k = 0

  DO WHILE (k.LT.N_electrons)

     k = k + 1

     if ( (electron(k)%X.lt.c_X_area_min).or. &
        & (electron(k)%X.gt.c_X_area_max).or. &
        & (electron(k)%Y.lt.c_Y_area_min).or. &
        & (electron(k)%Y.gt.c_Y_area_max) ) then
        print '("Process ",i4," : Error-1 in PERFORM_COULOMB_SCATTERING : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
        print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
        print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%Vx, electron(k)%Vy, electron(k)%Vz, electron(k)%tag
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if

! interpolate the moments

     i = INT(electron(k)%X)
     j = INT(electron(k)%Y)

     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

if ((i.lt.c_indx_x_min).or.(i.gt.(c_indx_x_max-1)).or.(j.lt.c_indx_y_min).or.(j.gt.(c_indx_y_max-1))) then
   print '("Process ",i4," : Error-2 in PERFORM_COULOMB_SCATTERING : index out of bounds")', Rank_of_process
   print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
   print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%Vx, electron(k)%Vy, electron(k)%Vz, electron(k)%tag
   print '("Process ",i4," : minx/maxx/miny/maxy : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
   CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
end if

     ax_ip1 = electron(k)%X - DBLE(i) ! consider setting all ax, ay to 0.5 to smooth out the moments?
     ax_i   = 1.0_8 - ax_ip1

     ay_jp1 = electron(k)%Y - DBLE(j)
     ay_j = 1.0_8 - ay_jp1

!     a00 =   ax_i * ay_j
!     a01 =   ax_i * ay_jp1
!     a10 = ax_ip1 * ay_j
!     a11 = ax_ip1 * ay_jp1
     a00 = 0.25_8
     a01 = 0.25_8
     a10 = 0.25_8
     a11 = 0.25_8
! these moments have already been time-averaged over the ion time step:
     n_loc = acc_rho_e(0,i,j) * a00 + acc_rho_e(0,i+1,j) * a10 + acc_rho_e(0,i,j+1) * a01 + acc_rho_e(0,i+1,j+1) * a11
     vx_loc= acc_rho_e(1,i,j) * a00 + acc_rho_e(1,i+1,j) * a10 + acc_rho_e(1,i,j+1) * a01 + acc_rho_e(1,i+1,j+1) * a11
     vy_loc= acc_rho_e(2,i,j) * a00 + acc_rho_e(2,i+1,j) * a10 + acc_rho_e(2,i,j+1) * a01 + acc_rho_e(2,i+1,j+1) * a11
     vz_loc= acc_rho_e(3,i,j) * a00 + acc_rho_e(3,i+1,j) * a10 + acc_rho_e(3,i,j+1) * a01 + acc_rho_e(3,i+1,j+1) * a11
     v2_loc= acc_rho_e(4,i,j) * a00 + acc_rho_e(4,i+1,j) * a10 + acc_rho_e(4,i,j+1) * a01 + acc_rho_e(4,i+1,j+1) * a11
      
     IF (n_loc .EQ. 0.0_8) CYCLE

        vx_loc = vx_loc / n_loc
        vy_loc = vy_loc / n_loc
        vz_loc = vz_loc / n_loc     

        v2_loc = v2_loc / n_loc - (vx_loc**2 + vy_loc**2 + vz_loc**2)      
        T_loc_eV  = 0.6666666667_8 * v2_loc * (N_max_vel ** 2) * T_e_eV 
! normalized density and temperature:
        dens_ratio = n_loc / DBLE(N_of_particles_cell)
        temp_ratio = T_loc_eV / T_e_eV

     Vx = electron(k)%VX
     Vy = electron(k)%VY
     Vz = electron(k)%VZ
! sample a collision partner from the field particle distribution:     
     call GetMaxwellVelocity(Vx_sample)
     call GetMaxwellVelocity(Vy_sample)
     call GetMaxwellVelocity(Vz_sample)
     vel_factor = SQRT(temp_ratio) / N_max_vel ! convert to internal units
     Vx_sample = Vx_sample * vel_factor + vx_loc
     Vy_sample = Vy_sample * vel_factor + vy_loc
     Vz_sample = Vz_sample * vel_factor + vz_loc

! relative velocity in internal units of the code:
     Ux = Vx - Vx_sample
     Uy = Vy - Vy_sample
     Uz = Vz - Vz_sample
     Usq = Ux**2 + Uy**2 + Uz**2
     Usq_scaled = Usq * N_max_vel**2 * 2.0_8 !normalized by the standard value of v_th, sqrt(T/m)
     U_scaled = SQRT(Usq_scaled)

     Wx = Vx - vx_loc
     Wy = Vy - vy_loc
     Wz = Vz - vz_loc
     Wsq = Wx**2 + Wy**2 + Wz**2
     Wsq_scaled = Wsq * N_max_vel**2 * 2.0_8 

     Vframe_x = 0.5_8 * (Vx + Vx_sample)
     Vframe_y = 0.5_8 * (Vy + Vy_sample)
     Vframe_z = 0.5_8 * (Vz + Vz_sample)

! Nanbu [1997] is the reference work applied here:
!     L_ee = L_ee_0 * (temp_ratio / dens_ratio)**0.5 * (3.0_8 * temp_ratio + Wsq_scaled) ! account for possible beam electron
     L_ee = L_ee_0 * (temp_ratio / dens_ratio)**0.5 * (6.0_8 * temp_ratio)
     IF (U_scaled.GE.1.e-5) THEN
        s_coll = LOG(L_ee) / pi * base_Coulomb_probab * dens_ratio / U_scaled**3
     ELSE
        s_coll = 10.0_8 ! results in isotropic scattering
     END IF  

     A_norm = 3.0_8 * EXP(-s_coll) / (1.0_8 - EXP(-3.0_8 * s_coll)) ! new approximation for Nanbu normalization factor

! sample the scattering angle from the random walk process:      
     R = well_random_number()
     IF (R.LE.1.e-20) THEN 
        CosKsi = -1.0_8     
     ELSEIF (s_coll.LE.2.e-2) THEN
             CosKsi = 1.0_8 + s_coll * LOG(R) !this will be the prevailing case 
     ELSEIF (s_coll.ge.6.0_8) THEN 
        CosKsi = 2.0_8 * R - 1.0_8
     ELSE   
     CosKsi = LOG( EXP(-A_norm) + 2.* R * SINH(A_norm)) / A_norm 
     END IF

     CosKsi = MAX(-1.0_8, MIN(CosKsi, 1.0_8))
     SinKsi = SQRT(1.0_8 - CosKsi ** 2)
     Fi = 2.0_8 * pi * well_random_number()
     CosFi = COS(Fi)
     SinFi = SIN(Fi)
! velocity in the center-of-mass frame:
     Vx = 0.5_8 * Ux
     Vy = 0.5_8 * Uy
     Vz = 0.5_8 * Uz
     
! Scattering the electron in the CM frame.
! Turn the velocity in that frame, using Ksi and Fi. After that transform the velocity to the laboratory frame.

    V =    0.5_8 * SQRT(Usq) ! (U/2)
    V_xy = 0.5_8 * SQRT(Ux*Ux + Uy*Uy)

  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_s = Vx * CosKsi + (SinFi * V * b + CosFi * Vz * a) * SinKsi
     Vy_s = Vy * CosKsi - (SinFi * V * a - CosFi * Vz * b) * SinKsi
     Vz_s = Vz * CosKsi - V_xy * CosFi * SinKsi
  ELSE
     Vx_s = ABS(Vz) * SinKsi * CosFi
     Vy_s = ABS(Vz) * SinKsi * SinFi
     Vz_s = Vz * CosKsi
  END IF
! back to "lab" frame after scattering:  
  electron(k)%VX = Vframe_x + Vx_s 
  electron(k)%VY = Vframe_y + Vy_s
  electron(k)%VZ = Vframe_z + Vz_s

  END DO ! particle loop

!  print '("Coulomb scattering done ",i4)', Rank_of_process
  
END SUBROUTINE PERFORM_ELECTRON_COULOMB_SCATTERING

SUBROUTINE CLEAR_COULOMB_ARRAYS

  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INTEGER i, j, m
   
  IF(.NOT.Coulomb_flag) RETURN

  DO j = c_indx_y_min, c_indx_y_max
     DO m = 0, 4
        acc_rho_e(m, c_indx_x_min : c_indx_x_max, j) = 0.0_8
     END DO
  END DO

END SUBROUTINE CLEAR_COULOMB_ARRAYS


SUBROUTINE INITIATE_COULOMB_SCATTERING
  
  USE ParallelOperationValues
  USE CurrentProblemValues
!  USE IonParticles, ONLY : N_spec, Ms
!  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER flag_scatter

!  REAL(8) L_ee_0

  LOGICAL exists

! Coulomb_flag is read in INITIATE_PARAMETERS to set the required arrays of velocity moments
  IF (.NOT.Coulomb_flag) RETURN

  L_debye_0 = L_debye_m / sqrt(2.0_8)
  base_plasma_param = N_plasma_m3 * L_debye_0**3 ! L_debye_m has extra square root of 2 as defined. 
  base_Coulomb_probab = W_plasma_s1 * delta_t_s / base_plasma_param * DBLE(N_subcycles)
  L_ee_0 = pi * base_plasma_param   
  
END SUBROUTINE INITIATE_COULOMB_SCATTERING         

!--------------------------------------------------
!
subroutine get_global_electron_moments

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles, ONLY : N_electrons, electron, kin_energy_e_global

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  REAL(8), allocatable :: rbufer(:), pwbufer(:), totpwbufer(:)
  INTEGER ALLOC_ERR
  INTEGER k

  ALLOCATE(    rbufer(1:4), STAT = ALLOC_ERR)
  ALLOCATE(   pwbufer(1:4), STAT = ALLOC_ERR)
  ALLOCATE(totpwbufer(1:4), STAT = ALLOC_ERR)

!print '("get_global_electron_moments:: Rank_of_process ",i4)', Rank_of_process

  rbufer = 0.0_8
  pwbufer = 0.0_8
  totpwbufer = 0.0_8

  DO k = 1, N_electrons ! each process
     rbufer(1) = rbufer(1) + electron(k)%VX
     rbufer(2) = rbufer(2) + electron(k)%VY
     rbufer(3) = rbufer(3) + electron(k)%VZ
     rbufer(4) = rbufer(4) + electron(k)%VX**2 + electron(k)%VY**2 + electron(k)%VZ**2
  END DO

! collect electron moments from all processes in a cluster
  CALL MPI_REDUCE(rbufer, pwbufer, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (cluster_rank_key.EQ.0) THEN

     CALL MPI_REDUCE(pwbufer, totpwbufer, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_HORIZONTAL, ierr)

     IF (Rank_horizontal.EQ.0) THEN
!        PRINT '("Total :  momentum X/Y/Z = ", 3(2x,e16.9)," energy = ", e16.9)', totpwbufer(1:4)
         kin_energy_e_global = totpwbufer(4)
     END IF    
     CALL MPI_BCAST(kin_energy_e_global, 1, MPI_DOUBLE_PRECISION, 0, COMM_HORIZONTAL, ierr)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(kin_energy_e_global, 1, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)

!print '("get_global_electron_moments:: Rank_of_process ",i4," cluster_rank_key ",i4," done")', Rank_of_process, cluster_rank_key

end subroutine get_global_electron_moments