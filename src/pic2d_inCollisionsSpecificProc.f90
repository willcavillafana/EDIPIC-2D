
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

  if (no_rcx_collisions) return

  DO s = 1, N_spec
     if (.not.collision_rcx(s)%rcx_on) cycle
     n = collision_rcx(s)%neutral_species_index
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
  
  if (no_rcx_collisions) return

! clear collision counters
  DO s = 1, N_spec
     collision_rcx(s)%counter = 0
  END DO

  DO s = 1, N_spec

     if (.not.collision_rcx(s)%rcx_on) cycle

     n = collision_rcx(s)%neutral_species_index

     ngas_m3 = neutral(n)%N_m3
     sigma_m2_1eV = neutral(n)%sigma_rcx_m2_1eV
     alpha =        neutral(n)%alpha_rcx
     probab_rcx_therm_2  = (collision_rcx(s)%probab_thermal)**2

     factor_eV = Ms(s) * energy_factor_eV         ! instead of collision_rcx(s)%factor_eV
     vfactor = collision_rcx(s)%vfactor           ! to convert Maxwellian sample
     prob_factor = ngas_m3 * delta_t_s * N_subcycles

     DO i = 1, N_ions(s)

!        if (well_random_number().GT.neutral_density_normalized(n, ion(s)%part(i)%x, ion(s)%part(i)%y)) cycle   ! for uniform density profile this is not necessary

        vx = ion(s)%part(i)%VX
        vy = ion(s)%part(i)%VY
        vz = ion(s)%part(i)%VZ
        vsq = vx**2 + vy**2 +vz**2
        energy_eV = vsq * factor_eV
        vabs_ms = sqrt(vsq) * V_scale_ms 

!        probab_rcx = ngas_m3 * vabs_ms * sigma_rcx_m2(energy_eV, sigma_m2_1eV, alpha) * delta_t_s * N_subcycles
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

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE in_Collision_Elastic
!>    @details Computes scattered velocity for elastic ions-neutral collisions. 
!              Same formula as in EDIPIPC 1D. Maybe same formula as in Vahedi and Surendra, “A Monte Carlo Collision Model for the Particle-in-Cell Method.”
!!    @authors W. Villafana
!!    @date    08-25-2022
!--------------------------------------------------------------------------------------------------

SUBROUTINE in_Collision_Elastic(ion_index, indx_particle, counter,Vx_n,Vy_n,Vz_n)

   USE MCCollisions
   USE IonParticles
 
   USE rng_wrapper
 
   IMPLICIT NONE
 
   ! INTEGER, INTENT(IN)   ::indx_neutral   ! ordering numeber of neutral species
   INTEGER, INTENT(IN)   ::ion_index      ! species  ion index
   INTEGER, INTENT(IN)   ::indx_particle  ! ordering number of particle 
   INTEGER, INTENT(INOUT)::counter        ! counter of collisions
   REAL(8), INTENT(IN)   ::Vx_n,Vy_n,Vz_n ! vecloity component of neutral particle
 
   REAL(8) R                ! random number             
   REAL(8) Ksi              ! scattering angle (relative to the initial direction)
   REAL(8) CosKsi, SinKsi     
   REAL(8) Fi               ! azimuthal scattering angle 
   REAL(8) CosFi, SinFi     
   REAL(8) Vx, Vy, Vz       ! velocity components, before scattering
   REAL(8) Vx_s, Vy_s, Vz_s ! velocity components, after scattering
   REAL(8) V, V_xy, a, b            
   !REAL(8) delta_energy     ! electron energy drop due to collision
   REAL(8) alpha            ! coefficient, accounting the ion energy drop
 
 
   ! We assume that the masses of ion and neutral atom are equal, 
   ! in this case the scattering angle in laboratory frame is half of scattering angle in the center of mass (COM) frame
   ! In C.O.M. frame scattering is uniform...

   ! Calculate the scattering angle Ksi for the ion
   R      = well_random_number()
   Ksi    = R * 1.5707963268_8
   CosKsi = COS(Ksi)
   SinKsi = SIN(Ksi)   
 
   ! Calculate the azimuthal scattering angle for the ion
   R     = well_random_number()
   Fi    = R * 6.28318530718_8
   CosFi = COS(Fi)
   SinFi = SIN(Fi)

   ! Take the initial ion velocity in the FRAME where the NEUTRAL atom is initially at REST. This frame moves with the initial velocity of neutral atom.
   Vx = ion(ion_index)%part(indx_particle)%VX - Vx_n
   Vy = ion(ion_index)%part(indx_particle)%VY - Vy_n
   Vz = ion(ion_index)%part(indx_particle)%VZ - Vz_n

   ! We have scattered the ion, i.e. calculated the angles Ksi and Fi, in the frame where the neutral atom is initially at rest.
   ! We must turn the velocity in that frame, using Ksi and Fi. After that we must transform the velocity to the initial laboratory frame.   
 
   ! Precalculate   
   V    = SQRT(Vx * Vx + Vy * Vy + Vz * Vz) 
   V_xy = SQRT(Vx * Vx + Vy * Vy)
 
   ! Calculate the velocity of ion after scattering (turn the velocity) in the frame, where the neutral atom was initially at rest.
   IF (V_xy.GT.0.0_8) THEN
      a    = SinKsi * SinFi * V / V_xy
      b    = SinKsi * CosFi * Vz / V_xy 
      Vx_s = Vx * CosKsi + Vy * a + Vx * b
      Vy_s = Vy * CosKsi - Vx * a + Vy * b
      Vz_s = Vz * CosKsi          - V_xy * SinKsi * CosFi
   ELSE
      Vx_s = ABS(Vz) * SinKsi * CosFi
      Vy_s = ABS(Vz) * SinKsi * SinFi
      Vz_s = Vz * CosKsi
   END IF
   
   ! Energy loss. Probably under the assumption of mass(neutral)=mass(ion). See 3.2.1 Vahedi and Surendra, “A Monte Carlo Collision Model for the Particle-in-Cell Method.”
   alpha = CosKsi

   ! Renormalize the scattered ion velocity in order to account the energy drop and transform the velocity back to the laboratory frame
   Vx_s = Vx_s * alpha
   Vy_s = Vy_s * alpha
   Vz_s = Vz_s * alpha
   ion(ion_index)%part(indx_particle)%VX = Vx_s + Vx_n !VX_of_spec(2)%part(num) = Vx_s + Vx_n
   ion(ion_index)%part(indx_particle)%VY = Vy_s + Vy_n !VY_of_spec(2)%part(num) = Vy_s + Vy_n
   ion(ion_index)%part(indx_particle)%VZ = Vz_s + Vz_n !VZ_of_spec(2)%part(num) = Vz_s + Vz_n
 
   counter = counter + 1

END SUBROUTINE in_Collision_Elastic

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE in_Collision_Charge_Exchange
!>    @details Computes the new velocity of ion after charge exchange. 
!              Same formula as in EDIPIPC 1D. Based on cross section tables
!!    @authors W. Villafana
!!    @date    08-25-2022
!--------------------------------------------------------------------------------------------------

SUBROUTINE in_Collision_Charge_Exchange(ion_index, indx_particle, counter,Vx_n,Vy_n,Vz_n)!num, Vx_n, Vy_n, Vz_n)

   USE MCCollisions
   USE IonParticles
   USE rng_wrapper

   IMPLICIT NONE
   ! INTEGER, INTENT(IN)   ::indx_neutral   ! ordering numeber of neutral species
   INTEGER, INTENT(IN)   ::ion_index      ! species  ion index
   INTEGER, INTENT(IN)   ::indx_particle  ! ordering number of particle 
   INTEGER, INTENT(INOUT)::counter        ! counter of collisions
   REAL(8), INTENT(IN)   ::Vx_n,Vy_n,Vz_n ! vecloity component of neutral particle
 
   ! Set the ion velocity equal to the velocity of atom (we assume that neutral and ion just exchange one electron)
   ion(ion_index)%part(indx_particle)%VX = Vx_n !VX_of_spec(2)%part(num) = Vx_n
   ion(ion_index)%part(indx_particle)%VY = Vy_n !VY_of_spec(2)%part(num) = Vy_n
   ion(ion_index)%part(indx_particle)%VZ = Vz_n !VZ_of_spec(2)%part(num) = Vz_n
 
   ! Update counter value
   counter = counter +1 
 
END SUBROUTINE in_Collision_Charge_Exchange

!--------------------------------------------------------------------------------------------------
!     FUNCTION frequency_of_in_collision
!>    @details Computes ions-neutral collision frequency. Similar to frequency_of_en_collision
!!    @authors W. Villafana
!!    @date    03-22-2022
!--------------------------------------------------------------------------------------------------

REAL(8) FUNCTION frequency_of_in_collision(energy_eV, indx_neutral, colproc_id)

  USE MCCollisions, ONLY: neutral
  USE CurrentProblemValues, ONLY : V_scale_ms, e_Cl, N_subcycles, delta_t_s, amu_kg
  USE IonParticles, ONLY: M_i_amu

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: energy_eV
  INTEGER, INTENT(IN) :: indx_neutral
  INTEGER, INTENT(IN) :: colproc_id

  INTEGER N_crsect_points
  REAL(8) f_temp
  INTEGER j
  REAL(8) energy_j_eV, energy_jp1_eV, f_j, f_jp1
  INTEGER :: index_ion ! selected ion species

  index_ion = 1

  N_crsect_points = neutral(indx_neutral)%in_colproc(colproc_id)%N_crsect_points


  IF (energy_eV.GE.neutral(indx_neutral)%in_colproc(colproc_id)%energy_eV(N_crsect_points)) THEN

     f_temp = neutral(indx_neutral)%in_colproc(colproc_id)%crsect_m2(N_crsect_points) * SQRT(2.0_8 * energy_eV * e_Cl / (M_i_amu(index_ion)*amu_kg))

  ELSE IF (energy_eV.LT.neutral(indx_neutral)%in_colproc(colproc_id)%energy_eV(1)) THEN

     IF (neutral(indx_neutral)%in_colproc(colproc_id)%type.LT.30) THEN ! For ions: 10 = elastic, 20 = charge exchange
! for elastic collisions only
        f_temp = neutral(indx_neutral)%in_colproc(colproc_id)%crsect_m2(1) * SQRT(2.0_8 * energy_eV * e_Cl / (M_i_amu(index_ion)*amu_kg))
     ELSE
! no inelastic and ionization collisions if energy below threshold. Case that should be implemented for ions
        f_temp = 0.0_8
     END IF

  ELSE
     
     DO j = 1, N_crsect_points-1
        IF ( (energy_eV.GE.neutral(indx_neutral)%in_colproc(colproc_id)%energy_eV(j)).AND. &
           & (energy_eV.LT.neutral(indx_neutral)%in_colproc(colproc_id)%energy_eV(j+1)) ) EXIT
     END DO

     j = MIN(j, N_crsect_points-1)

     energy_j_eV   = neutral(indx_neutral)%in_colproc(colproc_id)%energy_eV(j)
     energy_jp1_eV = neutral(indx_neutral)%in_colproc(colproc_id)%energy_eV(j+1)

     f_j   = neutral(indx_neutral)%in_colproc(colproc_id)%crsect_m2(j)   * SQRT(2.0_8 * energy_j_eV   * e_Cl / (M_i_amu(index_ion)*amu_kg)) ! not sure about this formula
     f_jp1 = neutral(indx_neutral)%in_colproc(colproc_id)%crsect_m2(j+1) * SQRT(2.0_8 * energy_jp1_eV * e_Cl / (M_i_amu(index_ion)*amu_kg))

     f_temp = f_j + (f_jp1 - f_j) * (energy_eV - energy_j_eV) / (energy_jp1_eV - energy_j_eV)

  END IF

  frequency_of_in_collision = f_temp * neutral(indx_neutral)%N_m3 * N_subcycles * delta_t_s 
  
END FUNCTION frequency_of_in_collision