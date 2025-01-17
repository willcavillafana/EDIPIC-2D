
!---------------------------------------------------------------------------------------------------------------------
!
subroutine calculate_thermal_cx_probab

  USE CurrentProblemValues, ONLY: delta_t_s, N_subcycles, e_Cl, kB_JK, m_e_kg
  USE MCCollisions
  USE IonParticles, ONLY : N_spec, Ms

  implicit none

  integer s, n
  real(8) sigma_m2_1eV, alpha, Tgas_eV, ngas_m3, sigma_m2_therm, Vmean_ms
!function
  real(8) sigma_rcx_m2

  if (no_in_collisions) return

  DO s = 1, N_spec
     if (.not.collision_rcx(s)%rcx_on) cycle
     n = collision_rcx(s)%neutral_species_index
     if (neutral(n)%in_cols_on) cycle ! skip if we use the new in-col model
     sigma_m2_1eV = neutral(n)%sigma_rcx_m2_1eV
     alpha = neutral(n)%alpha_rcx
     Tgas_eV = neutral(n)%T_K * kB_JK / e_Cl
     ngas_m3 = neutral(n)%N_m3
     sigma_m2_therm = sigma_rcx_m2(1.5_8 * Tgas_eV, sigma_m2_1eV, alpha)
     Vmean_ms = sqrt(2.54647908947033_8 * Tgas_eV * e_Cl / (Ms(s) * m_e_kg))   ! 2.5464... = 8/pi, note that Ms(s) is (amu_kg * M_i_amu(s)) / m_e_kg
     collision_rcx(s)%probab_thermal = ngas_m3 * sigma_m2_therm * Vmean_ms * delta_t_s * N_subcycles
  END DO

end subroutine calculate_thermal_cx_probab

!------------------------------------------------------------------------- 
!
real(8) function sigma_rcx_m2(energy_eV, sigma_m2_1eV, alpha)

  implicit none

  real(8) energy_eV, sigma_m2_1eV, alpha, Tgas_eV

  sigma_rcx_m2 = sigma_m2_1eV * (1.0_8 - alpha * log(MAX(energy_eV,0.001_8)))**2   ! the limiter avoids singularity
                                                                                   ! for comparison, 300 K = 0.026 eV

end function sigma_rcx_m2

!---------------------------------------------------------------------------------------------------------------------
!
subroutine PERFORM_RESONANT_CHARGE_EXCHANGE

!  USE ParallelOperationValues

  USE MCCollisions
  USE IonParticles
  USE CurrentProblemValues, ONLY : energy_factor_eV, delta_t_s, N_subcycles, V_scale_ms
  USE rng_wrapper

  IMPLICIT NONE

!  INCLUDE 'mpif.h'
!  INTEGER ierr
!  INTEGER stattus(MPI_STATUS_SIZE)

  INTEGER s, n, i
  real(8) ngas_m3, sigma_m2_1eV, alpha, probab_rcx_therm_2
  real(8) factor_eV, vfactor, prob_factor

  real(8) vx, vy, vz, vsq, energy_eV, vabs_ms
  real(8) probab_rcx

! functions
  real(8) neutral_density_normalized, sigma_rcx_m2
  
  if (no_in_collisions) return

! clear collision counters
  DO s = 1, N_spec
     collision_rcx(s)%counter = 0
  END DO

  DO s = 1, N_spec

     if (.not.collision_rcx(s)%rcx_on) cycle

     n = collision_rcx(s)%neutral_species_index
     if (neutral(n)%in_cols_on) cycle ! skip if we use the new in-col model

     ngas_m3 = neutral(n)%N_m3
     sigma_m2_1eV = neutral(n)%sigma_rcx_m2_1eV
     alpha =        neutral(n)%alpha_rcx
     probab_rcx_therm_2  = (collision_rcx(s)%probab_thermal)**2

     factor_eV = Ms(s) * energy_factor_eV         ! instead of collision_rcx(s)%factor_eV
     vfactor = collision_rcx(s)%vfactor           ! to convert Maxwellian sample
     prob_factor = ngas_m3 * delta_t_s * N_subcycles

     DO i = 1, N_ions(s)


        vx = ion(s)%part(i)%VX
        vy = ion(s)%part(i)%VY
        vz = ion(s)%part(i)%VZ
        vsq = vx**2 + vy**2 +vz**2
        energy_eV = vsq * factor_eV
        vabs_ms = sqrt(vsq) * V_scale_ms 

        probab_rcx = prob_factor * vabs_ms * sigma_rcx_m2(energy_eV, sigma_m2_1eV, alpha)
        probab_rcx = neutral_density_normalized(n, ion(s)%part(i)%x, ion(s)%part(i)%y) * sqrt(probab_rcx**2 + probab_rcx_therm_2)  ! account for the nonuniform density and the low-energy correction

        if (well_random_number().le.probab_rcx) then
           call GetMaxwellVelocity(VX)
           call GetMaxwellVelocity(VY)
           call GetMaxwellVelocity(VZ)
           ion(s)%part(i)%VX = VX * vfactor
           ion(s)%part(i)%VY = VY * vfactor
           ion(s)%part(i)%VZ = VZ * vfactor
!           ion(s)%part(k)%tag = CXtag
           collision_rcx(s)%counter = collision_rcx(s)%counter + 1
!        else
        end if
     END DO

  END DO   

END SUBROUTINE PERFORM_RESONANT_CHARGE_EXCHANGE



SUBROUTINE PERFORM_ION_NEUTRAL_COLLISION
  ! Julian Held, j.held@tue.nl, 2025
  ! After: Nanbu and Kitatani: J. Phys. D: Appl. Phys. 28 (1995) 324-330
  
    USE MCCollisions
    USE IonParticles
    USE CurrentProblemValues, ONLY : energy_factor_eV, delta_t_s, N_subcycles, V_scale_ms, T_cntr, pi, e_Cl, true_eps_0_Fm, amu_kg, kB_JK, T_e_eV, N_max_vel
    USE rng_wrapper
    !use stdlib_specialfunctions_gamma, only: gamma
  
    IMPLICIT NONE
    INCLUDE 'mpif.h'
  
  
    INTEGER s, n, i
    INTEGER ierr
    LOGICAL CX ! did charge exchange occur?
    real(8) vfactor, ngas_m3
  
    real(8) vx, vy, vz, vx_, vy_, vz_, Ekin, vr_ms
    real(8) vxn, vyn, vzn, vxn_, vyn_, vzn_ ! neutral velocities before and after
    real(8) Rx, Ry, Rz
    real(8) gx, gy, gz, g, g_perp, hx, hy, hz
    real(8) Mi, Mn
    real(8) E_ratio, E_ratio_ion
  
    real(8) neutral_density_normalized ! function
  
    real(8) beta_inf, alpha0, q, mr, beta_cx
    real(8) sigmaL, sigmaP, sigmaT, sigma_cx, d0, a, Acx
    real(8) p_col, bmax_col, bmax_cx, bmax, b
    real(8) xi0, xi1, xi, chi, beta, beta0, theta, theta0, Fel, dtheta
    real(8) cos_chi, sin_chi, phi
  
    beta0 = 1.001 ! beta > beta0 -> spiraling 

    if (no_in_collisions) return

    ! clear collision counters
    DO s = 1, N_spec
        collision_rcx(s)%counter = 0
    END DO

    DO s = 1, N_spec
      DO n = 1, N_neutral_spec
        if (.not.neutral(n)%in_cols_on) cycle

        ngas_m3 = neutral(n)%N_m3
        beta_inf = neutral(n)%beta_inf ! beta_inf may be lowered for slow ions or small E/N to allow for larger timesteps
        alpha0 = neutral(n)%alpha * 1D-30 ! polarizability of neutral
        Acx = neutral(n)%Axc ! resonant charge exchange parameter (check for same species later)
        Mn = neutral(n)%M_amu * amu_kg
        vfactor = SQRT(neutral(n)%T_K * kB_JK / (T_e_eV * e_Cl * Ms(s))) / DBLE(N_max_vel) 

        q = Qs(s) * e_Cl 
        Mi = M_i_amu(s) * amu_kg
        mr = Mi*Mn/(Mi+Mn) 
    
        DO i = 1, N_ions(s)
          CX = .False.
          a = alpha0 * (q**2) / (2.0*(4.0*pi*true_eps_0_Fm))
          p_col = ngas_m3 * sqrt(8.0*a/mr) * pi * (beta_inf**2) * delta_t_s * N_subcycles * neutral_density_normalized(n, ion(s)%part(i)%x, ion(s)%part(i)%y)
    
          if (well_random_number().le.p_col) then ! perform collision
            collision_rcx(s)%counter = collision_rcx(s)%counter + 1
            beta = beta_inf * sqrt(well_random_number())
            beta_cx = Acx * (Ekin/e_Cl)**(0.25)
    
            ! create a virtual neutral particle
            call GetMaxwellVelocity(vxn)
            call GetMaxwellVelocity(vyn)
            call GetMaxwellVelocity(vzn)
            vxn = vxn * vfactor * V_scale_ms
            vyn = vyn * vfactor * V_scale_ms
            vzn = vzn * vfactor * V_scale_ms
      
            vx = ion(s)%part(i)%VX * V_scale_ms
            vy = ion(s)%part(i)%VY * V_scale_ms
            vz = ion(s)%part(i)%VZ * V_scale_ms
      
            ! relative velocity between colliding particles
            gx = vxn-vx
            gy = vyn-vy
            gz = vzn-vz
            g = sqrt((gx)**2 + (gy)**2 + (gz)**2)
            g_perp = sqrt(gy**2 + gz**2)
    
            Ekin = 0.5 * mr * g**2
    
            phi = 2.0*pi*well_random_number()
            hx = g_perp * cos(phi)
            hy = -(gy*gx*cos(phi) + g*gz*sin(phi)) / g_perp                                                
            hz = -(gz*gx*cos(phi) - g*gy*sin(phi)) / g_perp
    
            if (beta.GE.beta0) then ! beta > 1 ->  polarization scattering and maybe CX
              xi0 = sqrt(beta**2 - sqrt(beta**4 - 1))
              xi1 = sqrt(beta**2 + sqrt(beta**4 - 1))
              xi = xi0/xi1
    
              ! find scattering angle chi
              if (beta.LE.3) then
                Fel = neutral(n)%Fel_precalc( nint(xi*(50000.0-1.0)) )
                theta0 = sqrt(2.0)*beta/xi1 * Fel
                chi = pi - 2.0 * theta0
              else
                chi = -(3.0*pi/16.0)*beta**(-4)
              end if
    
              cos_chi = cos(chi)
              sin_chi = abs(sin(chi))
              
              ! post-collision velocities
              vx_ =  vx +  Mn/(Mi + Mn) * (gx*(1-cos_chi) + hx*sin_chi)
              vy_ =  vy +  Mn/(Mi + Mn) * (gy*(1-cos_chi) + hy*sin_chi)
              vz_ =  vz +  Mn/(Mi + Mn) * (gz*(1-cos_chi) + hz*sin_chi)
              vxn_ = vxn - Mi/(Mi + Mn) * (gx*(1-cos_chi) + hx*sin_chi)
              vyn_ = vyn - Mi/(Mi + Mn) * (gy*(1-cos_chi) + hy*sin_chi)
              vzn_ = vzn - Mi/(Mi + Mn) * (gz*(1-cos_chi) + hz*sin_chi)
      
              ! CX: identity switch (50% chance if beta < beta_cx)
              if ((well_random_number().le.0.5).AND.(beta.LE.beta_cx)) then  
                CX = .True.
              end if
    
            end if  
    
    
            if (beta.LT.beta0) then ! spiraling motion, VHS (variable hard sphere) model
              
              ! VHS model -> random direction
              theta = acos(1.0 - 2.0*well_random_number())
              phi = 2.0*pi*well_random_number()
              Rx = cos(theta)
              Ry = sin(theta) * cos(phi)
              Rz = sin(theta) * sin(phi)
              
              ! post-collision velocities
              vx_ =  (1/(Mi + Mn)) * (Mi*vx + Mn*vxn - Mn*g*Rx)
              vy_ =  (1/(Mi + Mn)) * (Mi*vy + Mn*vyn - Mn*g*Ry)
              vz_ =  (1/(Mi + Mn)) * (Mi*vz + Mn*vzn - Mn*g*Rz)
              vxn_ = (1/(Mi + Mn)) * (Mi*vx + Mn*vxn + Mi*g*Rx)
              vyn_ = (1/(Mi + Mn)) * (Mi*vy + Mn*vyn + Mi*g*Ry)
              vzn_ = (1/(Mi + Mn)) * (Mi*vz + Mn*vzn + Mi*g*Rz)
    
              if (well_random_number().le.0.5) then ! CX: (50% chance when spiraling)
                CX = .TRUE.
              end if 

            end if
            
            ! check if ion and neutral species are the same (via mass)
            ! if not, disable charge exchange after all (ignoring the non-resonant case)
            if (ABS(Mn-Mi).GT.0.1) then
              CX = .FALSE.
            end if

            if (CX) then ! assign neutral post-collision velcoity to the ion (col + identity switch)
              ion(s)%part(i)%VX = vxn_/V_scale_ms
              ion(s)%part(i)%VY = vyn_/V_scale_ms
              ion(s)%part(i)%VZ = vzn_/V_scale_ms
            else ! assign ion post-collision velocity (no identity switch)
              ion(s)%part(i)%VX = vx_/V_scale_ms
              ion(s)%part(i)%VY = vy_/V_scale_ms
              ion(s)%part(i)%VZ = vz_/V_scale_ms
            end if
    
            ! Check for energy conservation. 
            E_ratio = (0.5*Mi*(vx_**2+vy_**2+vz_**2) + 0.5*Mn*(vxn_**2 + vyn_**2 + vzn_**2))/(0.5*Mi*(vx**2+vy**2+vz**2) + 0.5*Mn*(vxn**2 + vyn**2 + vzn**2))
            E_ratio_ion = (0.5*Mi*(vx_**2+vy_**2+vz_**2))/(0.5*Mi*(vx**2+vy**2+vz**2))
            if ((E_ratio-1.0).GT.1E-6) then
              PRINT *, "Error in PERFORM_ION_NEUTRAL_COLLISION. Energy not conserved. This should never happen. Energy lost (1-E_before/E_after): ", (1-E_ratio)
              PRINT *, "Beta: ", beta  
              PRINT *, "Charge exchange occured?: ", CX  
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
            end if
    
            ! check if collision frequency is okay
            if (p_col.GT.0.5) then
              PRINT *, "Error in PERFORM_ION_NEUTRAL_COLLISION. Collision rate too high. Need to use smaller time steps. Collision probability was p_col =  ", p_col
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
            end if
    
          end if ! if col
        END DO ! ion loop
      END DO ! neutral species loop
    END DO ! ion species loop
  
  END SUBROUTINE PERFORM_ION_NEUTRAL_COLLISION