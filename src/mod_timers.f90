!--------------------------------------------------------------------------------------------------
!     MODULE timers
!>    @details Measures time spent in EDIPIC kernels
!!    @authors W. Villafana
!!    @date    Apr-25-2022
!--------------------------------------------------------------------------------------------------

MODULE mod_timers

    USE CurrentProblemValues, ONLY: zero,string_length
    IMPLICIT NONE

    !INCLUDE 'mpif.h'

    ! Define timer structure
    TYPE :: T_TIMER
        REAL(8) :: start ! starting time
        REAL(8) :: end   ! ending time
        REAL(8) :: total ! total time
    END TYPE   

    TYPE(T_TIMER) :: pic_loop_timer

!--------------------------------------------------------------------------------------------------
!     Timer definition and initialization
!>    @details Define and initialize timers
!!    @authors W. Villafana
!!    @date    Apr-25-2022
!-------------------------------------------------------------------------------------------------- 
    
    ! TYPE(T_TIMER) :: pic_loop_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: save_checkpoint_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: global_load_balance_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: internal_load_balance_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: gather_ion_charge_density_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: gather_electron_charge_density_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: poisson_solver_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: calculate_electric_field_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: compute_averaged_snapshot_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: create_instantaneous_snapshot_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: electron_pusher_with_collisions_inner_object_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: ions_pusher_with_collisions_inner_object_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: transfer_particle_after_pusher_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: collect_particles_hitting_with_bo_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: compute_mcc_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: add_ions_after_collisions_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: clear_accumulated_fields_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: create_averaged_snapshot_and_compute_ptcl_emission = t_timer( start=zero,end=zero,total=zero )    
    ! TYPE(T_TIMER) :: add_electrons_after_emission_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: save_bo_particle_hits_emissions_timer = t_timer( start=zero,end=zero,total=zero )
    ! TYPE(T_TIMER) :: gather_surface_charge_density_timer = t_timer( start=zero,end=zero,total=zero )

    CONTAINS

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE start_timer
!>    @details Set starting time
!!    @authors W. Villafana
!!    @date    Apr-25-2022
!-------------------------------------------------------------------------------------------------- 
    SUBROUTINE start_timer ( timer_t )


        IMPLICIT NONE
        INCLUDE 'mpif.h'

        !IN/OUT
        TYPE(T_TIMER), INTENT(INOUT) :: timer_t

        timer_t%start = MPI_WTIME()

    END SUBROUTINE

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE end_timer
!>    @details Stop timer increment total time
!!    @authors W. Villafana
!!    @date    Apr-25-2022
!-------------------------------------------------------------------------------------------------- 
    SUBROUTINE end_timer ( timer_t )

        IMPLICIT NONE
        INCLUDE 'mpif.h'
        
        !IN/OUT
        TYPE(T_TIMER), INTENT(INOUT) :: timer_t

        timer_t%end = MPI_WTIME()
        timer_t%total = timer_t%total + timer_t%end - timer_t%start
        
    END SUBROUTINE    

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_timer
!>    @details Print timer at the end
!!    @authors W. Villafana
!!    @date    Apr-25-2022
!-------------------------------------------------------------------------------------------------- 
    SUBROUTINE print_timer ( timer_t,timer_name )
        
        USE ParallelOperationValues, ONLY: Rank_of_process
        IMPLICIT NONE

        !IN/OUT
        TYPE(T_TIMER), INTENT(IN)                :: timer_t
        CHARACTER(LEN=string_length), INTENT(IN) :: timer_name

        IF ( Rank_of_process==0 ) THEN
            WRITE(*,'(A,T50,A,T120,ES12.4,A,F6.2,A)') "CPU TIME ",">>> "//TRIM(timer_name)//": ",&
                    timer_t%total, " (s) ==> ", timer_t%total/pic_loop_timer%total*100.0_8, "%"  
        ENDIF 
    END SUBROUTINE        

END MODULE
