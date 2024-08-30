MODULE mod_def_timers
    USE CurrentProblemValues, ONLY: zero,string_length
    USE mod_timers, ONLY: T_TIMER,print_timer,pic_loop_timer
!--------------------------------------------------------------------------------------------------
!     Timer definition and initialization
!>    @details Define and initialize timers
!!    @authors W. Villafana
!!    @date    Apr-25-2022
!-------------------------------------------------------------------------------------------------- 
    
    TYPE(T_TIMER) :: total_timer = t_timer ( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: single_pic_loop_timer = t_timer ( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: save_checkpoint_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: global_load_balance_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: internal_load_balance_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: gather_ion_charge_density_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: gather_electron_charge_density_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: poisson_solver_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: externa_circuit_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: calculate_electric_field_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: compute_averaged_snapshot_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: create_instantaneous_snapshot_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: electron_pusher_with_collisions_inner_object_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: ions_pusher_with_collisions_inner_object_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: transfer_particle_after_pusher_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: collect_particles_hitting_with_bo_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: compute_mcc_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: add_ions_after_collisions_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: clear_accumulated_fields_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: compute_ptcl_emission_timer = t_timer( start=zero,end=zero,total=zero )    
    TYPE(T_TIMER) :: add_electrons_after_emission_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: save_bo_particle_hits_emissions_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: gather_surface_charge_density_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: ptcl_tracing_timer = t_timer( start=zero,end=zero,total=zero )
    TYPE(T_TIMER) :: create_averaged_snapshot_timer = t_timer( start=zero,end=zero,total=zero )

    CONTAINS
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_cpu_time
!>    @details Print all CPU info at the end
!!    @authors W. Villafana
!!    @date    Apr-25-2022
!-------------------------------------------------------------------------------------------------- 
    SUBROUTINE print_cpu_time 

        IMPLICIT NONE

        !LOCAL
        CHARACTER(LEN=string_length) :: message

        message = "Total PIC loop"
        CALL print_timer( pic_loop_timer,message )
        message = "Save checkpoint"
        CALL print_timer( save_checkpoint_timer,message )
        message = "Global load balance"        
        CALL print_timer( global_load_balance_timer,message )
        message = "Internal load balance"        
        CALL print_timer( internal_load_balance_timer,message )
        message = "Gather ion charge density"
        CALL print_timer( gather_ion_charge_density_timer,message )
        message = "Gather electron charge density"
        CALL print_timer( gather_electron_charge_density_timer,message )
        message = "Solve Poisson equation"
        CALL print_timer( poisson_solver_timer,message )
        message = "Adjust potential accounting for external circuit"
        CALL print_timer( externa_circuit_timer,message )        
        message = "Calculate electric field"
        CALL print_timer( calculate_electric_field_timer,message )
        message = "Compute averaged snapshot"
        CALL print_timer( compute_averaged_snapshot_timer,message )
        message = "Create instantaneous snapshot"
        CALL print_timer( create_instantaneous_snapshot_timer,message )
        message = "Electron pusher + collisions with internal object"
        CALL print_timer( electron_pusher_with_collisions_inner_object_timer,message )
        message = "Ion pusher + collisions with internal object"
        CALL print_timer( ions_pusher_with_collisions_inner_object_timer,message )
        message = "Transfer particles after pusher"
        CALL print_timer( transfer_particle_after_pusher_timer,message )
        message = "Collect particles hitting with boundary"
        CALL print_timer( collect_particles_hitting_with_bo_timer,message )
        message = "Compute MCC"
        CALL print_timer( compute_mcc_timer,message )
        message = "Add ions after collisions"
        CALL print_timer( add_ions_after_collisions_timer,message )
        message = "Clear accumulated fields"
        CALL print_timer( clear_accumulated_fields_timer,message )
        message = "Compute particle emission"
        CALL print_timer( compute_ptcl_emission_timer,message )
        message = "Add electrons after emission"
        CALL print_timer( add_electrons_after_emission_timer,message )
        message = "Save bo particle hits emissions"
        CALL print_timer( save_bo_particle_hits_emissions_timer,message )
        message = "Gather surface charge density"
        CALL print_timer( gather_surface_charge_density_timer,message )
        message = "Particle tracing"
        CALL print_timer( ptcl_tracing_timer,message )
        message = "Create average snapshot"
        CALL print_timer( create_averaged_snapshot_timer,message )        
        
    END SUBROUTINE 
    
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_iteration_info
!>    @details Print iteration information
!!    @authors W. Villafana
!!    @date    Mar-02-2023
!-------------------------------------------------------------------------------------------------- 
    SUBROUTINE print_iteration_info ( time_counter )
        
        USE CurrentProblemValues, ONLY: delta_t_s,one, Start_T_cntr
        USE mod_print, ONLY: print_message
    
        IMPLICIT NONE
        INCLUDE 'mpif.h'
    
        !IN/OUT
        INTEGER, INTENT(IN) :: time_counter
        
        ! LOCAL
        REAL(8) :: new_time, old_time, speed_loop
        CHARACTER(LEN=string_length) :: message
    
        ! First iteration
        IF ( time_counter==Start_T_cntr ) THEN
          WRITE(message,'(A)') "Iteration #        Total_time [s]        Iteration/second [s-1]"//achar(10)
          CALL print_message(message)
        END IF
    
        ! Compute time step ion the current iteration
        new_time = single_pic_loop_timer%end
        old_time = single_pic_loop_timer%start
        speed_loop = one/(new_time-old_time)
    
        WRITE(message,'(I10,9x,ES18.7,9x,ES12.4)') time_counter,DBLE(time_counter)*delta_t_s,speed_loop
        CALL print_message(message)
    
      END SUBROUTINE        

END MODULE