!------------------------------------
!
MODULE ExternalFields

  REAL(8) Bx_ext
  REAL(8) By_ext
  REAL(8) Bz_ext

  REAL(8) Ez_ext

! parameters magnetic field as in the paper of Boeuf and Garrigues
  REAL(8) y_Bmax
  REAL(8) Bz_0
  REAL(8) Bz_max
  REAL(8) Bz_Lsys
  REAL(8) half_over_sigma2_1
  REAL(8) half_over_sigma2_2
  REAL(8) a1, a2
  REAL(8) b1, b2

! parameters of set of wires with JZ current to get magnetic field with Bx, By
  INTEGER N_JZ_wires
  REAL(8), ALLOCATABLE :: JZwire_X(:)    ! x-coordinate of the wire
  REAL(8), ALLOCATABLE :: JZwire_Y(:)    ! y-coordinate of the wire
  REAL(8), ALLOCATABLE :: JZwire_JZ(:)    ! JZ current in the wire

  ! Parameter for ECR cathode
  INTEGER :: i_mag_profile ! 0: off, 1 gaussian in x direction, 1 double direction in y direction
  REAL(8) :: top_B_val_1, bottom_B_val_1, y_discontinuity_1 ! for i_mag=2 (double unifrom profile in y). y_discontinuity_1 is where the field goes from top_B_val_1 to bottom_B_val_1

END MODULE ExternalFields

!------------------------------------
!
MODULE ParallelOperationValues

  INTEGER, PARAMETER :: SHIFT1 = 0 !1000000
  INTEGER, PARAMETER :: SHIFT2 = 0 !2000000
  INTEGER, PARAMETER :: SHIFT3 = 0 !2000000

  INTEGER Rank_of_process
  INTEGER N_of_processes

  INTEGER Rank_of_process_left
  INTEGER Rank_of_process_right
  INTEGER Rank_of_process_above
  INTEGER Rank_of_process_below

  LOGICAL WHITE           ! for field solver, block-related

  LOGICAL WHITE_CLUSTER   ! for particle mover, cluster-related
  
  INTEGER field_master
  INTEGER particle_master
  INTEGER cluster_rank_key

  INTEGER Rank_of_master_left
  INTEGER Rank_of_master_right
  INTEGER Rank_of_master_above
  INTEGER Rank_of_master_below

  INTEGER COMM_BOUNDARY
  INTEGER Rank_boundary
  INTEGER N_processes_boundary

  INTEGER COMM_CLUSTER
  INTEGER Rank_cluster
  INTEGER N_processes_cluster

  INTEGER COMM_HORIZONTAL
  INTEGER Rank_horizontal
  INTEGER N_processes_horizontal

  INTEGER N_processes_cluster_left
  INTEGER N_processes_cluster_right
  INTEGER N_processes_cluster_above
  INTEGER N_processes_cluster_below

  INTEGER Rank_horizontal_left
  INTEGER Rank_horizontal_right
  INTEGER Rank_horizontal_above
  INTEGER Rank_horizontal_below

  INTEGER Rank_of_bottom_left_cluster_master

END MODULE  ParallelOperationValues

!------------------------------------
!
MODULE ParallelFFTX

  INTEGER maxbufsize

  INTEGER N_of_comm_steps_X
  INTEGER N_of_comm_steps_Y

  INTEGER fftx_band_jmin       ! are used in forward fft
  INTEGER fftx_band_jmax       ! no overlapping in j

  INTEGER invfftx_band_jmin    ! are used in reverse fft
  INTEGER invfftx_band_jmax    ! overlapping in j permitted, include endpoints c_indx_x_min/max where applicable

  INTEGER fftx_strip_jmin      ! are used in forward fft
  INTEGER fftx_strip_jmax      ! no overlapping in j

  INTEGER invfftx_strip_jmin   ! are used in inverse fft
  INTEGER invfftx_strip_jmax   ! overlapping in j permitted, include endpoints c_indx_x_min/max where applicable

  INTEGER sysy_band_nmin
  INTEGER sysy_band_nmax

  INTEGER sysy_strip_nmin
  INTEGER sysy_strip_nmax

  TYPE comm_proc_fftx
     INTEGER rank
     INTEGER fftx_band_jmin      ! these limits are used for calculation of forward fft
     INTEGER fftx_band_jmax      ! they avoid overlapping [between neighbor cluster lines]
     INTEGER invfftx_band_jmin   ! these limits are use for calculation of inverse fft
     INTEGER invfftx_band_jmax   ! they allow overlapping [between neighbor cluster lines]
     INTEGER c_indx_x_min
     INTEGER c_indx_x_max
  END TYPE comm_proc_fftx

  TYPE(comm_proc_fftx), ALLOCATABLE :: comm_FFTX_step(:)

  TYPE comm_proc_sysy
     INTEGER rank
     INTEGER sysy_band_nmin
     INTEGER sysy_band_nmax
     INTEGER c_indx_y_max
     INTEGER c_indx_y_min
  END TYPE comm_proc_sysy

  TYPE(comm_proc_sysy), ALLOCATABLE :: comm_SYSY_step(:)

  REAL(8), ALLOCATABLE :: a_eq(:,:)    ! pre-calculated upper-diagonal coefficients of the linear y-system of equations
                                       ! [for several fft-n, second index]
  LOGICAL two_dielectric_walls

END MODULE ParallelFFTX

!------------------------------------
!
MODULE CurrentProblemValues

  INTEGER, PARAMETER :: VACUUM_GAP = 0
  INTEGER, PARAMETER :: METAL_WALL = 1
  INTEGER, PARAMETER :: PERIODIC_PIPELINE_X = 2
  INTEGER, PARAMETER :: PERIODIC_PIPELINE_Y = 3
  INTEGER, PARAMETER :: DIELECTRIC = 4
  INTEGER, PARAMETER :: SYMMETRY_PLANE = 5
  INTEGER, PARAMETER :: NEUMANN = 6

  INTEGER, PARAMETER :: PERIODICITY_NONE    = 0
  INTEGER, PARAMETER :: PERIODICITY_X       = 1
  INTEGER, PARAMETER :: PERIODICITY_X_PETSC = 10
  INTEGER, PARAMETER :: PERIODICITY_X_Y     = 2

  INTEGER, PARAMETER :: END_FLAT = 0
  INTEGER, PARAMETER :: END_CORNER_CONCAVE = 1
  INTEGER, PARAMETER :: END_CORNER_CONVEX = 2
  INTEGER, PARAMETER :: CONCAVE_CORNER = 3
  INTEGER, PARAMETER :: CONVEX_CORNER = 4

  REAL(8), PARAMETER :: e_Cl     = 1.602189d-19      ! Charge of single electron [Cl]
  REAL(8), PARAMETER :: m_e_kg   = 9.109534d-31      ! Mass of single electron [kg]
  REAL(8), PARAMETER :: true_eps_0_Fm = 8.854188d-12      ! The dielectric constant [F/m]
  REAL(8), PARAMETER :: mu_0_Hm  = 1.25663706d-6     ! vacuum permeability [H/m] or [m kg s^-2 A^-2]
  REAL(8), PARAMETER :: amu_kg   = 1.660565d-27      ! atomic mass unit [kg]
  REAL(8), PARAMETER :: kB_JK    = 1.38064852d-23    ! Boltzmann constant [J/K]

  REAL(8), PARAMETER :: pi = 3.141592653589793_8
  REAL(8), PARAMETER :: zero = 0.0_8
  REAL(8), PARAMETER :: half = 1.0_8/2.0_8
  REAL(8), PARAMETER :: third = 1.0_8/3.0_8
  REAL(8), PARAMETER :: fourth = 1.0_8/4.0_8
  REAL(8), PARAMETER :: one = 1.0_8
  REAL(8), PARAMETER :: two = 2.0_8
  REAL(8), PARAMETER :: three = 3.0_8
  REAL(8), PARAMETER :: four = 4.0_8
  REAL(8), PARAMETER :: five = 5.0_8
  REAL(8), PARAMETER :: bessel_first_zero_order_zero = 2.4048255_8
  REAL(8) eps_0_Fm

  INTEGER, PARAMETER :: string_length = 300
  INTEGER, PARAMETER :: dble_precision = selected_real_kind(15, 307)

  CHARACTER(len=string_length), PARAMETER :: GIT_BRANCH="GIT_BRANCH"
  CHARACTER(len=string_length), PARAMETER :: GIT_HASH="GIT_HASH"
  CHARACTER(len=string_length), PARAMETER :: GIT_DATE="GIT_DATE"

  INTEGER ::  i_cylindrical ! Choose if this is a Cartesian (=0) or a Cylindrical (>0) case. r_theta plane=> 1,  r_z plane=>2, 
                            ! For the z_theta plane we can take a Cartesian geometry for now
  INTEGER :: debug_level

  REAL(8) :: weight_ptcl ! Statistical weight of particles

  INTEGER :: i_no_poisson                                    ! Deactivate poisson   (no  e_field)
  INTEGER :: i_reflection_cyl_electron, i_reflection_cyl_ion ! Indicates if   we asked for a specular reflection for electrons and ions in cylindrical geometry. Inelastic collisions are not implemented for now (Dec 20, 2022)
  INTEGER :: i_empty_domain                                  ! Initialize domain with no particles (=1, 0 otherwise by default)
  INTEGER :: i_one_particle                                  ! Initialize domain with one particle (=1, 0 otherwise by default)
  REAL(8) :: vthx_factor,vthy_factor,vthz_factor             ! If one particle in domain, these are factor of of the thermal speed for each direction (you can choose the velocity of the imposed particle)
  INTEGER :: i_velocity_one_particle                         ! I impose velocity of the only particle in domain (=1) or I do nothing(=0)
  INTEGER :: i_customed_profile_initialization                 ! Customed initial profiles imposed. =1: n=N_uniform*J0(alpha*X)*(ch(alpha*(Y-Y_shift_exp_Y))-ch(alpha*Y0_exp)/sh(alpha*Y0_exp)*sh(alpha*(Y-Y_shift_exp_Y))
  REAL(8) :: Rw_bessel                                       ! i_customed_profile_initialization=1: side wall. Profile valid for 0<=X<=Rw_bessel
  REAL(8) :: Y_shift_exp_Y                                   ! i_customed_profile_initialization=1: Y shift. Profile valid for Y_shift_exp_Y<Y<Y0_exp+Y_shift_exp_Y
  REAL(8) :: Y0_exp                                          ! i_customed_profile_initialization=1: total distance over which profile is valid. Profile valid for Y_shift_exp_Y<Y<Y0_exp+Y_shift_exp_Y

  REAL(8) :: LX, LY ! maximal dimensions of domain (dimensions)

  INTEGER :: c_indx_x_min_total, c_indx_x_max_total, c_indx_y_min_total, c_indx_y_max_total ! Maximal and minal index of simulation doain in X and Y directions

  ! Coulomb collisions stuff
  REAL(8) L_ee_0 ! base argument of Coulomb logarithm
  REAL(8) base_Coulomb_probab
  REAL(8) base_plasma_param

  INTEGER i_given_F_double_period_sys    ! in a system which is periodic in both X and Y directions, if there is no metal objects with given potential
  INTEGER j_given_F_double_period_sys    ! it is necessary to specify a point with some given potential, otherwise the field solver converges only if no dielectric objects are inside (pure plasma)
  REAL(8) given_F_double_period_sys      ! so here we specify the node (i,j) and the potential in the node

  INTEGER periodicity_flag   ! shows the presence of periodicity, is used to switch between periodic and non-periodic field solvers

  REAL(8) T_e_eV        ! scale electron temperature [eV]
  REAL(8) N_plasma_m3   ! scale electron density [m^-3]

  INTEGER N_of_cells_debye   ! number of cells per scale electron Debye length
  INTEGER N_max_vel          ! maximal expected velocity [units of scale thermal electron velocity]
  
  INTEGER N_blocks_x    ! number of blocks (processes) along the X (horizontal) direction (numbering starts at 1)
  INTEGER N_blocks_y    ! number of blocks (processes) along the Y (vertical) direction (numbering starts at 1)

  INTEGER N_grid_block_x  ! number of cells along the X-direction in a block
  INTEGER N_grid_block_y  ! number of cells along the Y-direction in a block

  INTEGER(8) :: N_of_particles_cell  ! number of macroparticles per cell for the scale density
  REAL(KIND=dble_precision) :: N_of_particles_cell_dble

  INTEGER cluster_N_blocks_x  ! number of blocks along the X-direction in a cluster
  INTEGER cluster_N_blocks_y  ! number of blocks along the Y-direction in a cluster
  INTEGER cluster_N_blocks    ! total number of processes (blocks) in a cluster

  INTEGER global_maximal_i    ! maximal values of indices in the X- and Y-directions
  INTEGER global_maximal_j    ! note that the origin has i/j=0/0

  INTEGER N_clusters_x        ! for rectangular domains, number of clusters along X
  INTEGER N_clusters_y        ! for rectangular domains, number of clusters along Y

  REAL(8) init_Te_eV   ! initial electron temperature [eV]. Assumed to be equal to init_TXe_eV in case of anisotropy 
  REAL(8) :: init_TXe_eV, init_TYe_eV, init_TZe_eV ! Possible anisotropy now
  REAL(8) init_Ne_m3   ! initial electron density [m^-3]

  REAL(8) v_Te_ms
  REAL(8) W_plasma_s1
  REAL(8) L_debye_m
  REAL(8) delta_x_m
  REAL(8) delta_t_s
  REAL(8) :: Delta_r,Delta_z
  REAL(8) :: shift_center_z

  REAL(8) E_scale_Vm
  REAL(8) B_scale_T
  REAL(8) F_scale_V
  REAL(8) V_scale_ms

  REAL(8) N_scale_part_m3
  REAL(8) current_factor_Am2
  REAL(8) energy_factor_eV
  REAL(8) temperature_factor_eV
  REAL(8) heat_flow_factor_Wm2

  INTEGER T_cntr
  INTEGER Start_T_cntr
  INTEGER Max_T_cntr
  INTEGER N_subcycles              ! number of electron cycles per ion cycle (must be odd)

  LOGICAL :: Coulomb_flag

  TYPE segment_of_object
     INTEGER istart
     INTEGER jstart
     INTEGER iend
     INTEGER jend
!     REAL(8), ALLOCATABLE :: surface_charge(:)

     INTEGER start_type
     INTEGER end_type

     LOGICAL, ALLOCATABLE :: cell_is_covered(:)

! surface_n_flux_e(:)
! surface_n_flux_i(:,:)
! surface_energy_flux_e(:)
! surface_energy_flux_i(:)

  END TYPE segment_of_object

  TYPE boundary_object
     INTEGER object_id_number
     INTEGER object_type
     INTEGER electron_hit_count
     INTEGER electron_emit_count
     INTEGER, ALLOCATABLE :: ion_hit_count(:)
     REAL(8) :: electron_hit_flux_avg_per_s
     REAL(8) :: electron_emission_flux_avg_per_s
     REAL(8), ALLOCATABLE :: ion_hit_flux_avg_per_s(:)

     REAL(8) total_charge
     REAL(8) phi                  ! electrostatic potential, used with metal walls and to calculate external field with dielectric walls
     REAL(8) phi_const            ! constant part of the electrostatic potential
     REAL(8) phi_var              ! time varying part of the electrostatic potential
     REAL(8) omega                ! frequency of time varying part
     REAL(8) phase                ! phase of time varying part
     REAL(8) phase_adjusted       ! phase of time varying part adjusted to the beginning of the amplitude profile period

     ! Customized waveform
     INTEGER :: i_customized_waveform ! 0 = no, 1 = yes
     INTEGER :: i_customized_waveform_name ! 1 = cosinus series as defined in DOI 10.1088/1361-6463/abf229 
     INTEGER :: nb_harmonics ! number of harmonics to be used ( see DOI 10.1088/1361-6463/abf229 )
     REAL(8) :: customized_waveform_phi ! phi_lf as in DOI 10.1088/1361-6463/abf229
     REAL(8) :: customized_waveform_freq ! freq as in DOI 10.1088/1361-6463/abf229

     ! linear model for ion induced model
     INTEGER :: i_linear_model_ion_induced ! 0 = no, 1 = yes
     REAL(8), ALLOCATABLE :: min_yield_ion_induced_ee_linear_model(:)
     REAL(8), ALLOCATABLE :: max_yield_ion_induced_ee_linear_model(:)

! waveform defines periodic non-harmonic variation of potential, the shape is defined by a user via data file
     LOGICAL use_waveform
     INTEGER N_wf_points                    ! number of waveform data points, must be no less than 2
     REAL,    ALLOCATABLE :: wf_phi(:)      ! array of potential values of waveform data points
     INTEGER, ALLOCATABLE :: wf_T_cntr(:)   ! array of times (in units of timesteps) of waveform data points

! amplitude profile for the oscillatory potential (includes the harmonic potential and the waveform)
     LOGICAL use_amplitude_profile
     INTEGER N_ap_points                    ! number of oscillation amplitude profile data points, must be no less than 2
     REAL,    ALLOCATABLE :: ap_factor(:)   ! array of factor values which will be multiplied by the oscillatory potential
     INTEGER, ALLOCATABLE :: ap_T_cntr(:)   ! array of times (in units of timesteps) of amplitude profile data points

     LOGICAL potential_must_be_solved       ! .TRUE. for [metal] electrodes connected to external circuits, .FALSE. otherwise
     INTEGER :: N_objects_connected_to_it   ! Number of objects connected to the preset=nt object (useful if I want to connect electrically object A with object B. Object B can be connected to external circuit)
     INTEGER, ALLOCATABLE :: connected_object_no(:)

     REAL    N_electron_constant_emit      ! constant number of electron macroparticles to be injected each time step (for example due to emission from a thermocathode)
     INTEGER model_constant_emit
     REAL(8) factor_convert_vinj_normal_constant_emit    ! factor to be used to convert values provided by Get*Velocity procedures to desired temperature
     REAL(8) factor_convert_vinj_parallel_constant_emit  ! factor to be used to convert values provided by Get*Velocity procedures to desired temperature
     REAL(8) v_ebeam_constant_emit                       ! velocity of the electron beam

     REAL(8) eps_diel                 ! relative dielectric permittivity, used with dielectric walls only

     INTEGER number_of_segments
     INTEGER L                    ! total length of the boundary object
     INTEGER n_connected_to_start ! number of the boundary object connected to the start point of this boundary object
     INTEGER n_connected_to_end   ! number of the boundary object connected to the start point of this boundary object
     TYPE(segment_of_object), ALLOCATABLE :: segment(:)
     REAL(8), ALLOCATABLE :: phi_profile(:)   ! ######## must be redone ??? #######, just put it here to compile now, must be in the segment

     INTEGER ileft               ! for a rectangular inner object we specify only left, right, bottom, and top index limits
     INTEGER iright
     INTEGER jbottom
     INTEGER jtop
     INTEGER N_boundary_nodes    ! number of nodes 

     REAL(8) Xmin
     REAL(8) Xmax
     REAL(8) Ymin
     REAL(8) Ymax

! are -1 by default
! for inner objects placed across periodic boundaries
     INTEGER object_copy_periodic_X_right  ! number of a copy inner object placed across opposite X periodic boundary (on the right)
     INTEGER object_copy_periodic_Y_above  ! number of a copy inner object placed across opposite Y periodic boundary (above)

     LOGICAL object_does_NOT_cross_symmetry_plane_X  ! default is .TRUE.

     REAL(8), ALLOCATABLE :: surface_charge(:)
     REAL(8), ALLOCATABLE :: surface_charge_variation(:)

     CHARACTER(6) material

! variables below define electron-induced emission of secondary electrons
! based on the 1D EDIPIC
     
! we have three possible processes when an electron hits the wall, producing a secondary electron:
! 1 - elastic backscattering
! 2 - inelastic backscattering
! 3 - true secondary emission

     LOGICAL SEE_enabled              ! switches on/off secondary electron emission

     INTEGER Emitted_model(1:3)          ! for each possible process determines the way of processing
     
     INTEGER Elast_refl_type             ! type of reflection for elastically reflected electron: specular / random

! ELASTIC, MODEL 1
     REAL(8) setD_see_elastic            ! constant coefficient (ratio of emitted and incident electron fluxess)
     REAL(8) minE_see_elastic            ! LOWER energy boundary, [dimensionless] 
     REAL(8) maxE_see_elastic            ! UPPER energy boundary, [dimensionless]

! ELASTIC, MODEL 2
     REAL(8) E_elast_0                   ! threshold energy, [dimensionless]
     REAL(8) E_elast_max                 ! maximum emission energy, [dimensionless]
     REAL(8) maxD_elast                  ! maximum emission yield
     REAL(8) dE_elast                    ! rate of decaying (like half-width) at energies > E_elast_max_eV, [dimensionless]
     REAL(8) Frac_elast_highenergy       ! fraction of total SEE yield at high energies, >=0, <<1
     
! INELASTIC, MODEL 1
     REAL(8) setD_see_inelastic          ! constant coefficient (ratio of emitted and incident electron fluxes)
     REAL(8) minE_see_inelastic          ! LOWER energy boundary, [dimensionless] 
     REAL(8) maxE_see_inelastic          ! UPPER energy boundary, [dimensionless]

! INELASTIC, MODEL 2
     REAL(8) Frac_inelastic              ! fraction of total SEE yield, >=0, <<1

! TRUE SECONDARY, MODEL 1
     REAL(8) setD_see_true               ! the coefficient (ratio of emitted to incident electrons)
     REAL(8) minE_see_true               ! LOWER energy boundary, [dimensionless] 
     REAL(8) maxE_see_true               ! UPPER energy boundary, [dimensionless]

! TRUE SECONDARY, MODEL 2 / CLASSIC SEE COEFFICIENT
     REAL(8) E_see_0                     ! threshold energy, [dimensionless]
     REAL(8) E_see_max                   ! maximal emission energy, [dimensionless]
     REAL(8) maxD_see_classic            ! maximal emission coefficient (normal to the surface)
     REAL(8) k_smooth                    ! smoothness factor (from 0 for very rough to 2 for polished surfaces)

     REAL(8) T_see_true_eV               ! Temperature of injected true secondary electrons, [eV]
     REAL(8) factor_convert_seetrue_vinj ! factor to be used to convert values provided by Get*Velocity procedures to desired temperature
     REAL(8) lowest_energy_for_see       ! threshold energy below which no emission occurs

! variables below define emission of electrons caused by ion impact on the material surface
! based on the 1D EDIPIC
! here ii_ee stands for ion-induced-electron-emission

     LOGICAL reflects_all_ions        ! switches on/off reflection of all ions (presently specular)
     LOGICAL ion_induced_EE_enabled   ! switches on/off electron emission caused by ions hitting the wall
     LOGICAL :: ion_thermalization    ! switches on/off thermalization of ions based on initial temperature

     REAL(8), ALLOCATABLE :: setD_ii_ee_true(:)                    ! the coefficient (ratio of emitted to incident electrons)
     REAL(8), ALLOCATABLE :: minE_ii_ee_true(:)                    ! LOWER energy boundary, [dimensionless] 
     REAL(8), ALLOCATABLE :: maxE_ii_ee_true(:)                    ! UPPER energy boundary, [dimensionless]
     REAL(8), ALLOCATABLE :: T_ii_ee_true_eV(:)                    ! Temperature of injected true secondary electrons, [eV]
     REAL(8), ALLOCATABLE :: factor_convert_ii_ee_true_vinj(:)     ! factor to be used to convert values provided by Get*Velocity procedures to desired temperature

  END TYPE boundary_object

  INTEGER N_of_boundary_objects
  INTEGER N_of_inner_objects
! inner objects can only be of type METAL_WALL or DIELECTRIC
! inner objects may overlap
! inner objects may touch the boundary of the domain
  INTEGER N_of_boundary_and_inner_objects  ! = N_of_boundary_objects + N_of_inner_objects

  TYPE(boundary_object), ALLOCATABLE :: whole_object(:)

  TYPE collided_particle
     INTEGER token
     REAL coll_coord
     REAL VX
     REAL VY
     REAL VZ
  END TYPE collided_particle

  TYPE boundary_object_statistics
     LOGICAL must_be_saved
     INTEGER max_N_of_saved_parts
     INTEGER N_of_saved_parts
     TYPE(collided_particle), ALLOCATABLE :: part(:)
  END TYPE boundary_object_statistics

  TYPE(boundary_object_statistics), ALLOCATABLE :: ion_colls_with_bo(:)
  TYPE(boundary_object_statistics), ALLOCATABLE :: e_colls_with_bo(:)

  REAL(8), ALLOCATABLE :: EX(:,:)        ! these arrays cover the whole cluster
  REAL(8), ALLOCATABLE :: EY(:,:)        !
  REAL(8), ALLOCATABLE :: BX_grid(:,:)        ! these arrays cover the whole cluster
  REAL(8), ALLOCATABLE :: BY_grid(:,:)        !
 
  REAL(8), ALLOCATABLE :: acc_EX(:,:)    ! these arrays cover the whole cluster
  REAL(8), ALLOCATABLE :: acc_EY(:,:)    ! 

  REAL(8), ALLOCATABLE :: c_rho(:,:)     ! these arrays cover the whole cluster
  REAL(8), ALLOCATABLE :: c_rho_i(:,:)   ! they are used in [semi]-periodic systems
  REAL(8), ALLOCATABLE :: c_phi(:,:)     !

  REAL(8), ALLOCATABLE :: acc_rho_e(:,:,:) ! to accumulate higher moments with Coulomb scattering enabled 
  REAL(8), ALLOCATABLE :: c_rho_ext(:,:,:) ! to accumulate higher moments with Coulomb scattering enabled 

  REAL(8), ALLOCATABLE :: phi(:,:)       !
  REAL(8), ALLOCATABLE :: rho_i(:,:)     ! these arrays cover a single field-calculating block, they are used when SOR is involved
  REAL(8), ALLOCATABLE :: rho_e(:,:)     ! 

!  REAL(8), ALLOCATABLE :: ext_phi(:)     ! linear potential profile corresponding to the external voltage applied between the two electrodes 
!                                         ! in case when the system is periodic in X and non-periodic in Y

END MODULE CurrentProblemValues

!-----------------------------   ????????????????????
!
MODULE ArraysOfGridValues

END MODULE ArraysOfGridValues

!-----------------------------
!
MODULE ElectronParticles

  TYPE particle
     REAL(8) X
     REAL(8) Y
     REAL(8) VX
     REAL(8) VY
     REAL(8) VZ
     INTEGER tag
  END TYPE particle

  INTEGER     N_electrons ! number of electron macroparticles
  INTEGER max_N_electrons ! current size (number of particles) of array of electrons macroparticles

  INTEGER     N_e_to_add  ! number of electron macroparticles to be added
  INTEGER max_N_e_to_add  ! current size (number of particles) of array electron_to_add

  REAL(8) :: kin_energy_e_global
  REAL(8) :: moments_global(0:4) ! for forcing energy and momentum conservation in Coulomb collisions  

  TYPE(particle), ALLOCATABLE :: electron(:)
  TYPE(particle), ALLOCATABLE :: electron_to_add(:)

! electron counters
  INTEGER N_e_to_send_left
  INTEGER N_e_to_send_right
  INTEGER N_e_to_send_above
  INTEGER N_e_to_send_below

! sizes of arrays of sent particles
  INTEGER max_N_e_to_send_left
  INTEGER max_N_e_to_send_right
  INTEGER max_N_e_to_send_above
  INTEGER max_N_e_to_send_below
  
! electron arrays of particles to be sent to neighbors
  TYPE(particle), ALLOCATABLE :: electron_to_send_left(:)
  TYPE(particle), ALLOCATABLE :: electron_to_send_right(:)
  TYPE(particle), ALLOCATABLE :: electron_to_send_above(:)
  TYPE(particle), ALLOCATABLE :: electron_to_send_below(:)

END MODULE ElectronParticles

!-----------------------------
!
MODULE IonParticles

  LOGICAL ions_sense_magnetic_field
  LOGICAL ions_sense_EZ
  INTEGER :: i_freeze_ions ! indicates if ions are moving or frozen in current simulation

!  INTEGER, PARAMETER :: N_spec = 1         ! number of ion species
  INTEGER N_spec        ! number of ion species

  TYPE particle
     REAL(8) X
     REAL(8) Y
     REAL(8) VX
     REAL(8) VY
     REAL(8) VZ
     INTEGER tag
  END TYPE particle

  TYPE particle_array
!     TYPE(particle), POINTER :: part(:)
     TYPE(particle), ALLOCATABLE :: part(:)
  END TYPE particle_array

  INTEGER, ALLOCATABLE :: Qs(:) !(1:N_spec)               ! q/e, [relative] charge of each ion species (e positive)
                                                          ! note, in the present version the ion charge can be either positive or negative 
                                                          ! but the absolute value should not exceed 3e
  REAL(8), ALLOCATABLE :: M_i_amu(:) !(1:N_spec)          ! ion mass, atomic mass unit

  REAL(8), ALLOCATABLE :: init_Ti_eV(:)                   ! initial ion temperature [eV]
  REAL(8), ALLOCATABLE :: init_NiNe(:)                    ! initial relative density

  REAL(8), ALLOCATABLE :: Ms(:) !(1:N_spec)               ! Mi/me
  REAL(8), ALLOCATABLE :: QM2s(:) !(1:N_spec)             ! (q/e)(me/Mi)/2
  REAL(8), ALLOCATABLE :: QM2sNsub(:) !(1:N_spec)         ! (q/e)(me/Mi)(N_subcycles/2)

  INTEGER, ALLOCATABLE ::      N_ions(:) !(1:N_spec)       ! number of ion macroparticles for each ion species
  INTEGER, ALLOCATABLE ::  max_N_ions(:) !(1:N_spec)       ! current size (number of particles) of array of macroparticles for each ion species

  INTEGER, ALLOCATABLE ::      N_ions_to_add(:) !(1:N_spec)  ! number of ion macroparticles to be added for each ion species
  INTEGER, ALLOCATABLE ::  max_N_ions_to_add(:) !(1:N_spec)  ! current size (number of particles) of array ion_to_add

  TYPE(particle_array), ALLOCATABLE :: ion(:)          ! Array of pointers on arrays of ions of different species
  TYPE(particle_array), ALLOCATABLE :: ion_to_add(:)   ! Array of pointers on arrays of ions of different species to be added

! ion species counters for particles sent to neighbors
  INTEGER, ALLOCATABLE ::  N_ions_to_send_left(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  N_ions_to_send_right(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  N_ions_to_send_above(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  N_ions_to_send_below(:) !(1:N_spec)

! sizes of arrays of sent particles
  INTEGER, ALLOCATABLE ::  max_N_ions_to_send_left(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  max_N_ions_to_send_right(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  max_N_ions_to_send_above(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  max_N_ions_to_send_below(:) !(1:N_spec)
  
! ion species arrays of particles to be sent to neighbors
  TYPE(particle_array), ALLOCATABLE :: ion_to_send_left(:)
  TYPE(particle_array), ALLOCATABLE :: ion_to_send_right(:)
  TYPE(particle_array), ALLOCATABLE :: ion_to_send_above(:)
  TYPE(particle_array), ALLOCATABLE :: ion_to_send_below(:)

END MODULE IonParticles

!----------------------------
!
MODULE ClusterAndItsBoundaries

  INTEGER, PARAMETER :: HAS_TWO_NEIGHBORS       = 10
  INTEGER, PARAMETER :: SURROUNDED_BY_WALL      = 20
  INTEGER, PARAMETER :: FLAT_WALL_LEFT          = 31
  INTEGER, PARAMETER :: FLAT_WALL_RIGHT         = 32
  INTEGER, PARAMETER :: FLAT_WALL_ABOVE         = 33
  INTEGER, PARAMETER :: FLAT_WALL_BELOW         = 34
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_LEFT  = 41
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_RIGHT = 42
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_BELOW = 43
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_ABOVE = 44

  INTEGER c_row
  INTEGER c_column

  REAL(8) c_X_area_min
  REAL(8) c_X_area_max
  REAL(8) c_Y_area_min
  REAL(8) c_Y_area_max

  REAL(8) :: index_maxi_r, index_maxi_z

  INTEGER c_indx_x_min
  INTEGER c_indx_x_max
  INTEGER c_indx_y_min
  INTEGER c_indx_y_max

  INTEGER c_left_bottom_corner_type
  INTEGER c_left_top_corner_type
  INTEGER c_right_bottom_corner_type
  INTEGER c_right_top_corner_type

  LOGICAL periodic_boundary_X_left
  LOGICAL periodic_boundary_X_right
  LOGICAL periodic_boundary_Y_below
  LOGICAL periodic_boundary_Y_above

  INTEGER i_period_x   ! integer values are easier to transmit to particle calculators
  INTEGER i_period_y   ! 
  REAL(8) L_period_X
  REAL(8) L_period_Y

  TYPE object_link
     INTEGER object_number
     INTEGER segment_number
     INTEGER istart
     INTEGER jstart
     INTEGER iend
     INTEGER jend
!     TYPE(segment_of_object) segment     
     REAL(8), ALLOCATABLE :: surface_charge(:) 

  END TYPE object_link

! if the object is bigger than the cluster area these flags show how to connect with neighbor clusters to get their contributions to the surface charge

  LOGICAL connect_left      ! if a cluster has connect_left=TRUE,  it has a neighbor on the left  with connect_right = TRUE
  LOGICAL connect_right     ! if a cluster has connect_right=TRUE, it has a neighbor on the right with connect_left  = TRUE
  LOGICAL connect_below     ! if a cluster has connect_below=TRUE, it has a neighbor below with connect_above = TRUE
  LOGICAL connect_above     ! if a cluster has connect_above=TRUE, it has a neighbor above with connect_below = TRUE

  LOGICAL symmetry_plane_X_left

  ! Neumann boundary imposed. Cluster variable
  LOGICAL neumann_X_cluster_left
  LOGICAL neumann_X_cluster_right
  LOGICAL neumann_Y_cluster_bottom
  LOGICAL neumann_Y_cluster_top

  INTEGER n_left(1:2)       ! contain indices of objects in array c_local_object_part
  INTEGER n_right(1:2)      ! which endpoints will participate in surface charge exchange
  INTEGER n_below(1:2)      ! with left/right/below/above neighbor cluster
  INTEGER n_above(1:2)

  INTEGER, PARAMETER :: c_max_N_of_local_object_parts = 20
  INTEGER c_N_of_local_object_parts
  TYPE(object_link) c_local_object_part(1:c_max_N_of_local_object_parts)

  INTEGER c_N_of_local_object_parts_left
  INTEGER c_N_of_local_object_parts_right
  INTEGER c_N_of_local_object_parts_above
  INTEGER c_N_of_local_object_parts_below

  INTEGER c_index_of_local_object_part_left( 1:c_max_N_of_local_object_parts)
  INTEGER c_index_of_local_object_part_right(1:c_max_N_of_local_object_parts)
  INTEGER c_index_of_local_object_part_above(1:c_max_N_of_local_object_parts)
  INTEGER c_index_of_local_object_part_below(1:c_max_N_of_local_object_parts)

  TYPE field_calc_proc
     INTEGER rank

     INTEGER indx_x_min
     INTEGER indx_x_max
     INTEGER indx_y_min
     INTEGER indx_y_max

     INTEGER fftx_strip_jmin
     INTEGER fftx_strip_jmax
     INTEGER invfftx_strip_jmin
     INTEGER invfftx_strip_jmax

     INTEGER sysy_strip_nmin
     INTEGER sysy_strip_nmax
  END TYPE field_calc_proc

  TYPE(field_calc_proc), ALLOCATABLE :: field_calculator(:)

  REAL(8), ALLOCATABLE :: vol_r_m3(:) ! node volume in cylindrical 
  REAL(8), ALLOCATABLE :: vol_cart(:) ! node volume in cartesian. By default it is dx**2. It is locally corrected after on the fly with IF conditions (not great) when there are wall boundaries 
  REAL(8), ALLOCATABLE :: factor_cyl_vol(:) ! corrective factory in claculation of density in cylindrical coordites: V_cart/V_cyl

END MODULE ClusterAndItsBoundaries

!-----------------------------------
!
MODULE LoadBalancing

  INTEGER T_cntr_global_load_balance
  INTEGER T_cntr_cluster_load_balance
  INTEGER dT_global_load_balance            ! time step interval for global load balancing, must be an integer number of N_subcycles
  INTEGER dT_cluster_load_balance           ! time step interval for balancing load within a cluster, must be an integer number of dT_cluster_load_balance

  TYPE cluster_params
     INTEGER particle_master
     INTEGER N_processes
     INTEGER N_particles

     INTEGER N_processes_balanced
     REAL    avg_load_balanced
  END TYPE cluster_params

  TYPE(cluster_params), ALLOCATABLE :: cluster(:)

END MODULE LoadBalancing


!----------------------------
!
MODULE BlockAndItsBoundaries

  USE CurrentProblemValues, ONLY: string_length
!  INTEGER, PARAMETER :: VACUUM_GAP = 0
!  INTEGER, PARAMETER :: METAL_WALL = 1

  INTEGER, PARAMETER :: HAS_TWO_NEIGHBORS       = 10
  INTEGER, PARAMETER :: SURROUNDED_BY_WALL      = 20
  INTEGER, PARAMETER :: FLAT_WALL_LEFT          = 31
  INTEGER, PARAMETER :: FLAT_WALL_RIGHT         = 32
  INTEGER, PARAMETER :: FLAT_WALL_ABOVE         = 33
  INTEGER, PARAMETER :: FLAT_WALL_BELOW         = 34
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_LEFT  = 41
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_RIGHT = 42
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_BELOW = 43
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_ABOVE = 44

  INTEGER block_row         ! from 1 to N_blocks_y, number of row in the matrix of processes
  INTEGER block_column     ! from 1 to N_blocks_x, number of column in the matrix of processes

  REAL(8) X_area_min
  REAL(8) X_area_max
  REAL(8) Y_area_min
  REAL(8) Y_area_max

  INTEGER indx_x_min
  INTEGER indx_x_max
  INTEGER indx_y_min
  INTEGER indx_y_max

  INTEGER left_bottom_corner_type
  INTEGER left_top_corner_type
  INTEGER right_bottom_corner_type
  INTEGER right_top_corner_type

  TYPE object_link
     INTEGER object_number
     INTEGER istart
     INTEGER jstart
     INTEGER iend
     INTEGER jend
!     TYPE(segment_of_object) segment     
!     REAL(8), ALLOCATABLE :: surface_charge(:) 
  END TYPE object_link

  INTEGER, PARAMETER :: max_N_of_local_object_parts = 20
  INTEGER N_of_local_object_parts
  TYPE(object_link) local_object_part(1:max_N_of_local_object_parts)

  LOGICAL block_has_symmetry_plane_X_left

  LOGICAL block_has_neumann_bc_X_left
  LOGICAL block_has_neumann_bc_X_right
  LOGICAL block_has_neumann_bc_Y_bottom
  LOGICAL block_has_neumann_bc_Y_top

  INTEGER N_of_local_object_parts_left
  INTEGER N_of_local_object_parts_right
  INTEGER N_of_local_object_parts_above
  INTEGER N_of_local_object_parts_below

  INTEGER index_of_local_object_part_left(1:max_N_of_local_object_parts)
  INTEGER index_of_local_object_part_right(1:max_N_of_local_object_parts)
  INTEGER index_of_local_object_part_above(1:max_N_of_local_object_parts)
  INTEGER index_of_local_object_part_below(1:max_N_of_local_object_parts)

  LOGICAL block_periodic_boundary_X_left
  LOGICAL block_periodic_boundary_X_right
  LOGICAL block_periodic_boundary_Y_below
  LOGICAL block_periodic_boundary_Y_above

  INTEGER N_to_solve_total
  INTEGER block_N_of_nodes_to_solve
  INTEGER global_offset

  INTEGER process_left_bottom_right_inner_node
  INTEGER process_left_solved_nodes_row_length

  INTEGER process_right_bottom_left_inner_node
  INTEGER process_right_solved_nodes_row_length

  INTEGER process_above_left_bottom_inner_node
  INTEGER process_below_left_top_inner_node

  CHARACTER(LEN=string_length) :: work_dir_partition_and_fields_files ! working directory for files used for partitionning and external fields
  
END MODULE BlockAndItsBoundaries


!----------------------------
!
MODULE Diagnostics

  USE CurrentProblemValues, ONLY: string_length

  CHARACTER(LEN=string_length) :: work_dir_probes  ! wkdir for probes

! diagnostic control
  INTEGER, PARAMETER :: Max_N_of_probes = 100

  INTEGER Save_probes_e_data_step
  INTEGER Save_probes_i_data_step

  INTEGER Save_probes_e_data_T_cntr 
  INTEGER Save_probes_i_data_T_cntr 

!  INTEGER Save_probes_data_T_cntr
  INTEGER Save_probes_data_step        ! WriteOut_step   ! Time interval (in steps) for writing into the file
  INTEGER Save_probes_data_T_cntr_rff  ! the value of Save_probes_data_T_cntr read from file (rff) init_probes.dat
                                       ! we need to save this value for consistent initialization of snapshots 
                                       ! when checkpoints are used to continue simulation

  INTEGER TextOut_skip    ! Periods of writing to be skipped between text outputs, 0 = write each, 1 = skip one, etc.

  INTEGER text_output_counter    ! this counter is used to skip extra diagnostics periods between two text outputs

  INTEGER N_of_saved_records     ! number of records saved in each time dependence data file, is used to trim protocol data files
                                 ! when checkpoint is used to initialize the system

  INTEGER N_of_probes         ! Total number of probes
  INTEGER N_of_probes_cluster ! number of probes in a particular cluster
  INTEGER N_of_probes_block   ! number of probes in a field calculator

  INTEGER, ALLOCATABLE :: Probe_position(:,:)           ! 1:2,1:N_of_probes_all ; (1,n)=x_n, (2,n)=y_n
  INTEGER, ALLOCATABLE :: List_of_probes_cluster(:)     ! 1:N_of_probes_cluster
  INTEGER, ALLOCATABLE :: Probe_params_block_list(:,:)  ! 1:3,1:N_of_probes_block

  REAL, ALLOCATABLE :: probe_F_cluster(:)         ! used by cluster masters in FFT-based field solver
                                                  ! and to assemble potential from blocks when PETSc-based solver is involved

  REAL, ALLOCATABLE :: probe_F_block(:)           ! used by a field calculator in PETSc-based field solver

  REAL, ALLOCATABLE :: probe_Ne_cluster(:)        ! these arrays keep diagnostics values obtained in different subroutines

  REAL, ALLOCATABLE :: probe_VXe_cluster(:)       ! till they are saved in the main probe diagnostics routine
  REAL, ALLOCATABLE :: probe_VYe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_VZe_cluster(:)       !

  REAL, ALLOCATABLE :: probe_VXVYe_cluster(:)     !
  REAL, ALLOCATABLE :: probe_VXVZe_cluster(:)     !
  REAL, ALLOCATABLE :: probe_VYVZe_cluster(:)     !

  REAL, ALLOCATABLE :: probe_WXe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_WYe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_WZe_cluster(:)       !

  REAL, ALLOCATABLE :: probe_QXe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_QYe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_QZe_cluster(:)       !

  REAL, ALLOCATABLE :: probe_Ni_cluster(:,:)      !

  REAL, ALLOCATABLE :: probe_VXi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_VYi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_VZi_cluster(:,:)     !

  REAL, ALLOCATABLE :: probe_VXVYi_cluster(:,:)   !
  REAL, ALLOCATABLE :: probe_VXVZi_cluster(:,:)   !
  REAL, ALLOCATABLE :: probe_VYVZi_cluster(:,:)   !

  REAL, ALLOCATABLE :: probe_WXi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_WYi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_WZi_cluster(:,:)     !

  REAL, ALLOCATABLE :: probe_QXi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_QYi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_QZi_cluster(:,:)     !

  REAL, ALLOCATABLE :: probe_JXsum(:)       !
  REAL, ALLOCATABLE :: probe_JYsum(:)       !
  REAL, ALLOCATABLE :: probe_JZsum(:)       !

  TYPE diagnostic_process
     INTEGER rank
     INTEGER N_probes
     INTEGER, ALLOCATABLE :: probe_number(:)
  END TYPE diagnostic_process

  TYPE(diagnostic_process), ALLOCATABLE :: diag_cluster(:)  ! this array wil be used by process with rank 0 which writes to the data files

  TYPE(diagnostic_process), ALLOCATABLE :: diag_block(:)    ! these arrays wil be used by master processes

END MODULE Diagnostics

!--------------------------
!
MODULE Snapshots

  USE CurrentProblemValues, ONLY: string_length

  INTEGER N_of_all_snaps                        ! number of all snapshots  
  INTEGER, ALLOCATABLE ::     Tcntr_snapshot(:)     ! timesteps when the snapshot files are written
  INTEGER, ALLOCATABLE :: save_evdf_snapshot(:)     ! flags controlling how evdf is saved
  INTEGER, ALLOCATABLE ::   save_pp_snapshot(:)     ! flags controlling saving of phase planes

  LOGICAL, ALLOCATABLE :: save_ionization_rates_2d(:)   ! turns on/off accumulation and saving of ioniization rates

  LOGICAL, ALLOCATABLE :: save_ions_collided_with_bo(:)    ! turns on/off saving of ions collided with boundary objects
                                                           ! (this must be confirmed by individual request for each boundary object)

  LOGICAL, ALLOCATABLE :: save_e_collided_with_bo(:)    ! turns on/off saving of electrons collided with boundary objects
                                                        ! (this must be confirmed by individual request for each boundary object)

  INTEGER current_snap                          ! index of current snapshot (which must be created)

! below, prefix cs stands for c-luster s-napshot
 
! these arrays are made global to simplify synchronization of overlapping nodes
  REAL, ALLOCATABLE :: cs_N(:,:)

  REAL, ALLOCATABLE :: cs_VX(:,:)
  REAL, ALLOCATABLE :: cs_VY(:,:)
  REAL, ALLOCATABLE :: cs_VZ(:,:)
 
  REAL, ALLOCATABLE :: cs_WX(:,:)
  REAL, ALLOCATABLE :: cs_WY(:,:)
  REAL, ALLOCATABLE :: cs_WZ(:,:)

  REAL, ALLOCATABLE :: cs_VXVY(:,:)  ! arrays required to calculate heat flow vector
  REAL, ALLOCATABLE :: cs_VXVZ(:,:)
  REAL, ALLOCATABLE :: cs_VYVZ(:,:)

  REAL, ALLOCATABLE :: cs_QX(:,:)    ! arrays required to calculate heat flow vector
  REAL, ALLOCATABLE :: cs_QY(:,:)
  REAL, ALLOCATABLE :: cs_QZ(:,:)

! flags for turning output of various parameters on/off

  LOGICAL save_data(1:38)  ! 1+2+3+16+16

! variables for calculation of 1D and 2D velocity distribution functions

  INTEGER, PARAMETER :: NOANYVDF = 0
  INTEGER, PARAMETER :: ONLY1D = 1
  INTEGER, PARAMETER :: ONLY2D = 2
  INTEGER, PARAMETER :: BOTH1DAND2D = 3

  INTEGER N_vdfbox_x       ! number of spatial boxes in the X-direction in a cluster
  INTEGER N_vdfbox_y       !           same as above in the Y-direction
  INTEGER N_vdfbox_all     ! total number of spatial boxes in a cluster

  INTEGER N_max_vel_e      ! maximal electron velocity in units of v_Te_ms
  INTEGER N_vbins_e        ! number of velocity bins per v_Te_ms for electrons

  INTEGER N_max_vel_i      ! maximal electron velocity in units of v_Te_ms * sqrt(me/Ms)
  INTEGER N_vbins_i        ! number of velocity bins per v_Te_ms * sqrt(me/Ms) for ions

  INTEGER indx_v_min_e     ! limits of velocity bin index for electrons
  INTEGER indx_v_max_e     !

  INTEGER indx_v_min_i     ! limits of velocity bin index for ions
  INTEGER indx_v_max_i     !

  INTEGER, ALLOCATABLE :: evxdf(:,:)
  INTEGER, ALLOCATABLE :: evydf(:,:)
  INTEGER, ALLOCATABLE :: evzdf(:,:)

  INTEGER, ALLOCATABLE :: evxvydf(:,:,:)

  INTEGER, ALLOCATABLE :: isvxdf(:,:,:)
  INTEGER, ALLOCATABLE :: isvydf(:,:,:)
  INTEGER, ALLOCATABLE :: isvzdf(:,:,:)

! variables for saving phase planes

  INTEGER, PARAMETER :: NOANYPP = 0
  INTEGER, PARAMETER :: ONLYelectronPP = 1
  INTEGER, PARAMETER :: ONLYionPP = 2
  INTEGER, PARAMETER :: BOTHelectronANDionPP = 3

  CHARACTER(LEN=string_length) :: work_dir_2d_map  ! working directory for 2D maps (instantaneous and averaged)

  INTEGER N_pp_boxes

  TYPE index_limits
     INTEGER imin
     INTEGER jmin
     INTEGER imax
     INTEGER jmax
  END TYPE index_limits

  TYPE(index_limits), ALLOCATABLE :: pp_box(:)

  TYPE specific_coll_diag
     REAL, ALLOCATABLE :: ionization_rate_local(:,:)
     REAL, ALLOCATABLE :: coll_freq_local(:,:)
  END TYPE specific_coll_diag

  TYPE coll_diag
!     INTEGER N_of_activated_colproc
     TYPE(specific_coll_diag), ALLOCATABLE :: activated_collision(:)
  END TYPE coll_diag

  TYPE(coll_diag), ALLOCATABLE :: diagnostics_neutral(:)

END MODULE Snapshots

!*******************************************************************************************
! This module is used by two subroutines, calculating the arbitrary velocity
! according to the isotropic maxwell distribution
MODULE MaxwellVelocity
  REAL(8), PARAMETER :: U_max = 3.0_8
  INTEGER, PARAMETER :: R_max     = 9000 !300
  INTEGER, PARAMETER :: R_max_inj = 4500 !150
  REAL(8) v(0:R_max)
  REAL(8) v_inj(0:R_max_inj)
END MODULE MaxwellVelocity

!--------------------------------------------------------------
!
MODULE SetupValues

!  USE IonParticles, ONLY : N_spec

  LOGICAL ht_use_e_emission_from_cathode_zerogradf
  LOGICAL ht_use_e_emission_from_cathode
  LOGICAL ht_emission_constant
  LOGICAL ht_grid_requested
  LOGICAL ht_soft_grid_requested
  LOGICAL ht_injection_inside
  LOGICAL :: ecr_injection_inside
  LOGICAL ht_use_ionization_source

  INTEGER N_macro_constant_injection
  REAL(8) injection_y
  REAL(8) :: injection_y_ecr
  INTEGER grid_j
  REAL(8) F_grid

  INTEGER total_cathode_N_e_to_inject   ! this variable is used to calculate total value of negatively charged particles
                                        ! that left the system without being balanced by an escaping ion
                                        ! this value becomes negative if the number of escaping ions exceeds 
                                        ! that of escaping electrons

  INTEGER total_cathode_N_e_injected    ! actual number of electron macroparticles injected, used for diagnostics only

  INTEGER j_ion_source_1   ! full source boundaries
  INTEGER j_ion_source_2   !

  INTEGER c_j_ion_source_1   ! source boundaries within a cluster 
  INTEGER c_j_ion_source_2   !

  INTEGER, ALLOCATABLE :: N_to_ionize_total(:)      !#########(1:N_spec)
  INTEGER, ALLOCATABLE :: N_to_ionize_cluster(:)    !#########(1:N_spec)
  INTEGER :: N_to_ionize_total_ecr      !#########(1:N_spec). For ECR setup 
  INTEGER :: N_to_ionize_cluster_ecr    !#########(1:N_spec). For ECR setup 

  INTEGER, PARAMETER :: c_R_max = 1000              ! maximal value of integral of ionization rate distribution within the range allocated to a cluster
  REAL(8) yi(0:c_R_max)
  REAL(8) yi_ecr(0:c_R_max)                         ! y position of particles that are injected with the ionization thoughts for ECR. 
                                                    !As a reminder, in the ECR setup you can control the size in the x direction (doesn't have to be the whole width of the domain)

  REAL(8) factor_convert_vinj                       ! velocity conversion factors, injected electrons
  REAL(8) :: factor_convert_vinj_i                       ! velocity conversion factors, injected ions. One species for now
  REAL(8) factor_convert_vion_e                     ! ionization electrons
  REAL(8), ALLOCATABLE ::  factor_convert_vion_i(:) !########(1:N_spec)  ! ionization ions
  INTEGER :: use_ecr_injection ! 0=off, 1=ON (e-i injection pair)


  ! ECR ionization moduel
  INTEGER :: j_ion_source_start_ecr, j_ion_source_end_ecr, i_ion_source_start_ecr, i_ion_source_end_ecr
  INTEGER :: c_j_ion_source_start_ecr, c_j_ion_source_end_ecr, c_i_ion_source_start_ecr, c_i_ion_source_end_ecr ! boundaries for cluster

  INTEGER :: i_neutral_profile ! For ECR project. Selects predefined neutral profiles
  REAL(8) :: ys_neutral_1,ye_neutral_1,nn_neutral_1 ! For i_neutral_profile=1 (graident in Y direction). This is starting and ending points of gradient. nn_neutral_1 is the final neutral density profile. Below that it remains the profiles remains at this value 
  INTEGER :: i_ionize_source_ecr ! For Ecr, volumic ionization source term 
  REAL(8) :: ioniz_ecr_vol_I_injected ! For ECR volumic ionization source: injected current in A
  REAL(8) :: xs_ioniz,xe_ioniz,ys_ioniz,ye_ioniz ! For ECR volumic ionization source: position of ionization volumic source

  ! Ionization source term that is flux dependent
  INTEGER :: i_ionize_source_flux_dependent ! if equal to 1, the number of pairs of e-i to be injected is given by a measured ion flux in the domain
  REAL(8) :: xs_ioniz_flux_dependent,xe_ioniz_flux_dependent,ys_ioniz_flux_dependent,ye_ioniz_flux_dependent ! Volume of ionization source term 
  ! INTEGER :: j_ion_source_start_flux_dependent, j_ion_source_end_flux_dependent, i_ion_source_start_flux_dependent, i_ion_source_end_flux_dependent
  INTEGER :: c_j_ion_source_start_flux_dependent, c_j_ion_source_end_flux_dependent, c_i_ion_source_start_flux_dependent, c_i_ion_source_end_flux_dependent ! boundaries for cluster
  INTEGER :: i_ionize_source_flux_dependent_plane_X, i_ionize_source_flux_dependent_plane_Y ! only one of them can be equal to 1. Will mean that plane is at location X= or Y= something
  REAL(8) :: cut_location_flux_plane_for_ionization ! position of flux plane in m 
  REAL(8) :: number_of_ions_crossing_plane, cluster_number_of_ions_crossing_plane, total_number_of_ions_crossing_plane
  REAL(8) :: portion_injections_per_cluster_flux_dependent_ionization, leftovers_injection!, total_number_injected
  ! INTEGER :: total_pairs_injected

END MODULE SetupValues

!--------------------------------------------------------------
!
MODULE Checkpoints

  LOGICAL use_mpiio_checkpoint

  INTEGER use_checkpoint             ! 0/1/2 = don't use / use to continue older run / use to start a new run
  INTEGER T_cntr_to_continue

  INTEGER dT_save_checkpoint
  INTEGER T_cntr_save_checkpoint

  INTEGER Save_probes_data_T_cntr_check
  INTEGER N_of_saved_records_check
  INTEGER current_snap_check

END MODULE Checkpoints

!-------------------------
!
MODULE MCCollisions

  LOGICAL en_collisions_turned_off
  LOGICAL no_ionization_collisions
  LOGICAL no_in_collisions

  INTEGER N_neutral_spec


  TYPE collision_type
     LOGICAL activated
     INTEGER type
     INTEGER ion_species_produced
     LOGICAL save_collfreq_2d
     REAL(8) threshold_energy_eV
     INTEGER N_crsect_points
     REAL(8), ALLOCATABLE :: energy_eV(:)
     REAL(8), ALLOCATABLE :: crsect_m2(:)
  END TYPE collision_type

  TYPE neutral_type
     CHARACTER(6) name
     REAL(8) M_amu
     REAL(8) N_m3
     REAL(8) T_K
     INTEGER N_en_colproc
     INTEGER N_of_energy_segments
     TYPE(collision_type), ALLOCATABLE :: en_colproc(:)
     REAL(8), ALLOCATABLE :: energy_segment_boundary_value(:)
     REAL(8), ALLOCATABLE :: energy_segment_step(:)

     logical rcx_on
     integer rcx_ion_species_index
     real(8) sigma_rcx_m2_1eV
     real(8) alpha_rcx
     logical in_cols_on ! switch ion-neutral collisions on/off
     real(8) alpha ! polarizability of neutral (in units of Angstrom^3)
     real(8) beta_inf ! cutoff for impact parameter of ion-neutral collisions, lower for speed, raise for accuracy
     real(8) Axc ! resonant charge exchange parameter (determine by comparison with swarm experiments)
     real(8), ALLOCATABLE :: Fel_precalc(:) ! careful when changing size!
  END TYPE neutral_type

  TYPE(neutral_type), ALLOCATABLE :: neutral(:)

  TYPE brief_collision_type
     INTEGER id_number
     INTEGER type
     INTEGER ion_species_produced
     REAL(8) threshold_energy_eV
     REAL(8) ion_velocity_factor
     LOGICAL save_collfreq_2d
  END TYPE brief_collision_type

  TYPE selected_collision_probability
     INTEGER N_of_energy_segments
     INTEGER N_of_energy_values
     INTEGER N_of_activated_colproc
     REAL(8) max_colliding_fraction
     TYPE(brief_collision_type), ALLOCATABLE :: colproc_info(:)
     INTEGER, ALLOCATABLE :: counter(:)
     REAL(8), ALLOCATABLE :: energy_segment_boundary_value(:)
     REAL(8), ALLOCATABLE :: energy_segment_step(:)
     INTEGER, ALLOCATABLE :: energy_segment_boundary_index(:)
     REAL(8), ALLOCATABLE :: prob_colproc_energy(:,:)
     REAL(8), ALLOCATABLE :: energy_eV(:)
  END TYPE selected_collision_probability

  TYPE(selected_collision_probability), ALLOCATABLE :: collision_e_neutral(:)

  TYPE rcx_collision_data
     logical rcx_on
     integer neutral_species_index
     real(8) probab_thermal
!     real(8) factor_eV
     real(8) vfactor !to scale the Maxwellian distribution
     INTEGER counter  ! to be used for diagnostics
  END TYPE rcx_collision_data

  TYPE(rcx_collision_data), allocatable :: collision_rcx(:)

! LINKED LIST, which store the numbers of particles participating in collisions

  TYPE binary_tree
     INTEGER number

     INTEGER neutral
     INTEGER indx_coll

     TYPE (binary_tree), POINTER :: Larger
     TYPE (binary_tree), POINTER :: Smaller     
  END TYPE binary_tree

  TYPE(binary_tree), POINTER :: Collided_particle

END MODULE MCCollisions

!--------------------------------
!
MODULE ExternalCircuit

  INTEGER N_of_object_potentials_to_solve

  REAL(8), ALLOCATABLE :: phi_due_object(:,:,:)   ! potential shape function
                                                  ! for the potential created by object nn when it's potential is 1 while all other walls are grounded
                                                  ! indices (i,j,nn) with i,j the spatial indices, nn the number of the object
                                                  ! note: object numbering is different from the one of whole_object
                                                  ! because it only accounts for objects connected to the circuit

  REAL(8), ALLOCATABLE :: potential_of_object(:)  ! actual value of the potential of the object
                                                  ! numbering only accounts for objects connected to the circuit
                                                  ! need this to have potential values from preivious time step

  REAL(8), ALLOCATABLE :: charge_of_object(:)     ! full charge of the object
                                                  ! need this to have charge values from previous step

  REAL(8), ALLOCATABLE :: object_charge_coeff(:,:)   ! dimension (0:N_of_object_potentials_to_solve, 1:N_of_object_potentials_to_solve)
                                                     ! object_charge(nn) = object_charge_coeff(0,nn) +
                                                     !                     object_charge_coeff(1,nn) * potential_of_object(1) + 
                                                     !                     object_charge_coeff(2,nn) * potential_of_object(2) + etc
                                                     ! object_charge_coeff(0,nn) depends on charge density and related potential profile,
                                                     ! it must be updated each timestep
                                                     ! object_charge_coeff(1:N_of_object_potentials_to_solve,nn) are precalculated

  REAL(8), ALLOCATABLE :: dQ_plasma_of_object(:)  ! change of charge of object due to plasma electrons and ions collided with the object
                                                  ! and emission of electrons only

  REAL(8) :: potential_of_object_avg ! for diagnostics, average on the fly potential of first object
  REAL(8) :: charge_of_object_avg ! for diagnostics, average on the fly charge of first object
  REAL(8) :: dQ_full_avg ! for diagnostics, average on the fly cumulated charge of first object
  REAL(8) :: dQ_plasma_of_object_avg ! for diagnostics, average on the fly change of charge due to plasma of first object

! for a point on the surface of an object we need a template to define which nodes to use
! and which coefficients to use
! in the structure below, there are arrays of size 3
! element 1 is the node on the surface
! elements 2 and 3 are additional nodes to be used
! if only one additional node is used, element 3 is not used
  TYPE node_control
     INTEGER N_of_nodes_to_use     ! from 0 to 5
     INTEGER, ALLOCATABLE :: use_i(:)              ! index i of node to use
     INTEGER, ALLOCATABLE :: use_j(:)              ! index j of node to use
     REAL(8) use_alpha_rho                         ! coefficient for charge density (the surface node only)
     REAL(8), ALLOCATABLE :: use_alpha_phi(:)      ! coefficient for potential in the node
  END TYPE node_control

  TYPE calculation_control
     INTEGER noi                           ! actual index of boundary object (index in whole_object array)
     INTEGER N_of_points_to_process
     TYPE(node_control), ALLOCATABLE :: control(:)
  END TYPE calculation_control

  TYPE(calculation_control), ALLOCATABLE :: object_charge_calculation(:)

! circuit parameters (for the test case only, must be edited for different systems)
!  REAL(8) source_U
!  REAL(8) source_omega
!  REAL(8) source_phase

  INTEGER :: circuit_type
  INTEGER N_of_power_supplies   ! number of power supplies in the exernal circuit

  TYPE PSU_type
! similar to the corresponding section of the boundary_object type

     REAL(8) phi_const            ! constant part of the electrostatic potential
     REAL(8) phi_var              ! time varying part of the electrostatic potential
     REAL(8) omega                ! frequency of time varying part
     REAL(8) phase                ! phase of time varying part
     REAL(8) phase_adjusted       ! phase of time varying part adjusted to the beginning of the amplitude profile period

! waveform defines periodic non-harmonic variation of potential, the shape is defined by a user via data file
     LOGICAL use_waveform
     INTEGER N_wf_points                    ! number of waveform data points, must be no less than 2
     REAL,    ALLOCATABLE :: wf_phi(:)      ! array of potential values of waveform data points
     INTEGER, ALLOCATABLE :: wf_T_cntr(:)   ! array of times (in units of timesteps) of waveform data points

! amplitude profile for the oscillatory potential (includes the harmonic potential and the waveform)
     LOGICAL use_amplitude_profile
     INTEGER N_ap_points                    ! number of oscillation amplitude profile data points, must be no less than 2
     REAL,    ALLOCATABLE :: ap_factor(:)   ! array of factor values which will be multiplied by the oscillatory potential
     INTEGER, ALLOCATABLE :: ap_T_cntr(:)   ! array of times (in units of timesteps) of amplitude profile data points
  END TYPE PSU_type

  TYPE(PSU_type), ALLOCATABLE :: EC_power_supply(:)
  
  INTEGER N_of_resistors    ! number of resistors
  INTEGER N_of_capacitors   ! number of capacitors
  INTEGER N_of_inductors    ! number of inductors

  REAL(8), ALLOCATABLE :: resistor_R_Ohm(:)
  REAL(8), ALLOCATABLE :: capacitor_C_F(:)
  REAL(8), ALLOCATABLE :: inductor_L_H(:)

  REAL(8), ALLOCATABLE :: J_ext(:) ! for current-driven circuit

END MODULE ExternalCircuit

!--------------------------
!
MODULE AvgSnapshots

  INTEGER N_of_all_avgsnaps                        ! number of all snapshots

  INTEGER current_avgsnap                          ! index of current snapshot (which must be created)

  TYPE average_snapshot_timing 
     INTEGER T_cntr_begin          ! time step when accumulation of average data for this snapshot begins
     INTEGER T_cntr_end            ! time step when it ends and the data are saved
  END TYPE average_snapshot_timing

  TYPE(average_snapshot_timing), ALLOCATABLE :: avgsnapshot(:)

!  INTEGER avg_data_collection_offset               ! offset relative to timestep following the ion move timestep
                                                   ! 0 means that data will be collected at the timestep immediately after the timestep when ions moved
                                                   ! N_subcycles-1 means that data will be collected at the timestep when the ions move
                                                   ! [note that average data are collected after electric field is calculated but before particles are advanced]

! arrays for accumulation of data

  REAL, ALLOCATABLE :: cs_avg_phi(:,:)    !  1

  REAL, ALLOCATABLE :: cs_avg_EX(:,:)     !  2
  REAL, ALLOCATABLE :: cs_avg_EY(:,:)     !  3

! total current is calculated immediately before writing to file, no need for the global array 

  REAL, ALLOCATABLE :: cs_avg_Ne(:,:)     !  7

  REAL, ALLOCATABLE :: cs_avg_JXe(:,:)    !  8
  REAL, ALLOCATABLE :: cs_avg_JYe(:,:)    !  9
  REAL, ALLOCATABLE :: cs_avg_JZe(:,:)    ! 10

  REAL, ALLOCATABLE :: cs_avg_VXe(:,:)    ! 11
  REAL, ALLOCATABLE :: cs_avg_VYe(:,:)    ! 12
  REAL, ALLOCATABLE :: cs_avg_VZe(:,:)    ! 13
 
  REAL, ALLOCATABLE :: cs_avg_WXe(:,:)    ! 14
  REAL, ALLOCATABLE :: cs_avg_WYe(:,:)    ! 15
  REAL, ALLOCATABLE :: cs_avg_WZe(:,:)    ! 16

  REAL, ALLOCATABLE :: cs_avg_TXe(:,:)    ! 17
  REAL, ALLOCATABLE :: cs_avg_TYe(:,:)    ! 18
  REAL, ALLOCATABLE :: cs_avg_TZe(:,:)    ! 19

  REAL, ALLOCATABLE :: cs_avg_QXe(:,:)    ! 20
  REAL, ALLOCATABLE :: cs_avg_QYe(:,:)    ! 21
  REAL, ALLOCATABLE :: cs_avg_QZe(:,:)    ! 22

  REAL, ALLOCATABLE :: cs_avg_Ni(:,:,:)   ! 23

  REAL, ALLOCATABLE :: cs_avg_JXi(:,:,:)  ! 24
  REAL, ALLOCATABLE :: cs_avg_JYi(:,:,:)  ! 25
  REAL, ALLOCATABLE :: cs_avg_JZi(:,:,:)  ! 26

  REAL, ALLOCATABLE :: cs_avg_VXi(:,:,:)  ! 27
  REAL, ALLOCATABLE :: cs_avg_VYi(:,:,:)  ! 28
  REAL, ALLOCATABLE :: cs_avg_VZi(:,:,:)  ! 29
 
  REAL, ALLOCATABLE :: cs_avg_WXi(:,:,:)  ! 30
  REAL, ALLOCATABLE :: cs_avg_WYi(:,:,:)  ! 31
  REAL, ALLOCATABLE :: cs_avg_WZi(:,:,:)  ! 32

  REAL, ALLOCATABLE :: cs_avg_TXi(:,:,:)  ! 33
  REAL, ALLOCATABLE :: cs_avg_TYi(:,:,:)  ! 34
  REAL, ALLOCATABLE :: cs_avg_TZi(:,:,:)  ! 35

  REAL, ALLOCATABLE :: cs_avg_QXi(:,:,:)  ! 36
  REAL, ALLOCATABLE :: cs_avg_QYi(:,:,:)  ! 37
  REAL, ALLOCATABLE :: cs_avg_QZi(:,:,:)  ! 38

  REAL, ALLOCATABLE :: cs_avg_coulomb_ee(:,:)  
  LOGICAL :: save_collision_freq_ee 

! flags for turning output of various parameters on/off (see above)

  LOGICAL save_avg_data(1:39)  ! 1+2+3+16+16+1

  REAL, ALLOCATABLE :: cs_Npart_coll(:,:)    ! electron densities (in units of macroparticles) immediately before the e-neutral collisions are applied
                                             ! is used to calculate e-neutral collision frequencies

  INTEGER :: num_plane_x_locations, num_plane_y_locations
  REAL, ALLOCATABLE :: flux_through_plane_over_one_period(:,:), cluster_flux_through_plane_over_one_period(:,:), total_flux_through_plane_over_one_period(:,:) ! For electrons and s species of ions (first dim) and N planes X or Y (second dim)
  REAL(8), ALLOCATABLE :: plane_x_cuts_location(:), plane_y_cuts_location(:)

  LOGICAL :: avg_flux_and_history

END MODULE AvgSnapshots

!--------------------------------------------------------------------------------------------------
!     MODULE print
!>    @details Prints readable messages for user and debug
!!    @authors W. Villafana
!!    @date    Dec-5-2022
!--------------------------------------------------------------------------------------------------

MODULE mod_print
    
    USE CurrentProblemValues, ONLY: string_length
    IMPLICIT NONE

    CONTAINS
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_message_cluster
!>    @details Print message by cluster master. 
!!    @authors W. Villafana
!!    @date    Nov-25-2022
!-------------------------------------------------------------------------------------------------- 
    SUBROUTINE print_message_cluster ( message, debug_level )
      
      USE ParallelOperationValues, ONLY: Rank_of_process
      IMPLICIT NONE

      !IN/OUT
      CHARACTER(LEN=string_length), INTENT(IN) :: message
      INTEGER, INTENT(IN) :: debug_level

      IF ( debug_level>0) WRITE(*,'(T8,A,I3,A)') "Cluster master ",Rank_of_process,": "//TRIM(message) 

  END SUBROUTINE    
  
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_message
!>    @details Print message by master proc and indicates current subroutine 
!!    @authors W. Villafana
!!    @date    Nov-25-2022
!-------------------------------------------------------------------------------------------------- 
  SUBROUTINE print_message ( message,routine )
      
    USE ParallelOperationValues, ONLY: Rank_of_process
    IMPLICIT NONE

    !IN/OUT
    CHARACTER(LEN=string_length), INTENT(IN) :: message
    CHARACTER(LEN=string_length), INTENT(IN), OPTIONAL :: routine
    
    IF ( Rank_of_process==0) THEN
        IF (PRESENT( routine )) THEN 
            WRITE(*,'(T8,A,T18,A)') "In subroutine "//TRIM(routine)//": "//TRIM(message) 
        ELSE
            WRITE(*,'(T8,A)') TRIM(message) 
        END IF
    END IF

  END SUBROUTINE  

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_debug
!>    @details Returns name of subroutine depending on debug_level
!!    @authors W. Villafana
!!    @date    Nov-25-2022
!-------------------------------------------------------------------------------------------------- 
  SUBROUTINE print_debug ( routine,local_debug_level )
      
    USE ParallelOperationValues, ONLY: Rank_of_process
    USE CurrentProblemValues, ONLY: debug_level
    IMPLICIT NONE

    !IN/OUT
    CHARACTER(LEN=string_length), INTENT(IN) :: routine
    INTEGER, INTENT(IN) :: local_debug_level
    
    IF ( Rank_of_process==0) THEN
        IF ( local_debug_level<debug_level) THEN 
            WRITE(*,'(T8,A)') "In subroutine "//TRIM(routine)
        END IF
    END IF

  END SUBROUTINE    

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_output
!>    @details Print generic output by master proc 
!!    @authors W. Villafana
!!    @date    Nov-25-2022
!-------------------------------------------------------------------------------------------------- 
  SUBROUTINE print_output ( message )
        
    USE ParallelOperationValues, ONLY: Rank_of_process
    IMPLICIT NONE

    !IN/OUT
    CHARACTER(LEN=string_length), INTENT(IN) :: message
    
    IF ( Rank_of_process==0) WRITE(*,'(A)') TRIM(message)

  END SUBROUTINE      

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_output_all
!>    @details Print generic output by all called proc 
!!    @authors W. Villafana
!!    @date    Nov-25-2022
!-------------------------------------------------------------------------------------------------- 
  SUBROUTINE print_output_all ( message,debug_level )
        
    IMPLICIT NONE

    !IN/OUT
    CHARACTER(LEN=string_length), INTENT(IN) :: message
    INTEGER, INTENT(IN) :: debug_level
    
    IF ( debug_level>0 ) WRITE(*,'(A)') TRIM(message)

  END SUBROUTINE      
  
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_error
!>    @details Print error message encountered by proc 
!!    @authors W. Villafana
!!    @date    Dec-8-2022
!-------------------------------------------------------------------------------------------------- 
  SUBROUTINE print_error ( message,routine )
        
    USE ParallelOperationValues, ONLY: Rank_of_process
    
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    !IN/OUT
    CHARACTER(LEN=string_length), INTENT(IN) :: message
    CHARACTER(LEN=string_length), INTENT(IN), OPTIONAL :: routine
    
    ! LOCAL
    INTEGER :: ierr
    WRITE(*,'(T8,A,I4,A)') ">>> In subroutine "//TRIM(routine)//", error PROC ",Rank_of_process," "//TRIM(message)
    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

  END SUBROUTINE       
  
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_parser_error
!>    @details Print error message by root proc when reading input files
!!    @authors W. Villafana
!!    @date    Dec-23-2022
!-------------------------------------------------------------------------------------------------- 
  SUBROUTINE print_parser_error ( message )
        
    USE ParallelOperationValues, ONLY: Rank_of_process
    
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    !IN/OUT
    CHARACTER(LEN=string_length), INTENT(IN) :: message
    
    ! LOCAL
    INTEGER :: ierr

    IF ( Rank_of_process==0 ) THEN
      WRITE(*,'(T8,A,I4,A)') ">>> Parser error: "//TRIM(message)
      CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
    END IF

  END SUBROUTINE       
  
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_warning
!>    @details Print error message by root proc when reading input files
!!    @authors W. Villafana
!!    @date    Mar-01-2023
!-------------------------------------------------------------------------------------------------- 
  SUBROUTINE print_warning ( message )
        
    USE ParallelOperationValues, ONLY: Rank_of_process
    
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    !IN/OUT
    CHARACTER(LEN=string_length), INTENT(IN) :: message
    
    ! LOCAL
    INTEGER :: ierr

    IF ( Rank_of_process==0 ) THEN
      WRITE(*,'(T8,A,I4,A)') ">>> WARNING: "//TRIM(message)
      ! CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
    END IF

  END SUBROUTINE     

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE print_git_info
!>    @details Print last commit ID, branch and compilation time
!!    @authors W. Villafana
!!    @date    Mar-03-2023
!-------------------------------------------------------------------------------------------------- 
  SUBROUTINE print_git_info
        
    USE ParallelOperationValues, ONLY: Rank_of_process
    USE CurrentProblemValues, ONLY: string_length, GIT_BRANCH, GIT_HASH, GIT_DATE
    
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    
    ! LOCAL
    CHARACTER(LEN=string_length) :: message
    
    IF ( Rank_of_process==0 ) THEN

      ! Print Branch
      WRITE(message,'(A,A)') "GIT BRANCH: ",TRIM(GIT_BRANCH)
      CALL print_message( message )

      ! Print Commit ID
      WRITE(message,'(A,A)') "GIT COMMIT: ",TRIM(GIT_HASH)
      CALL print_message( message )

      ! Print Commit date
      WRITE(message,'(A,A,A)') "GIT DATE: ",TRIM(GIT_DATE),achar(10)
      CALL print_message( message )              
      
    END IF

  END SUBROUTINE     

END MODULE mod_print

!*******************************************************************************
!> author: Jacob Williams
!  license: BSD
!  date: 2/7/2016
!
!  Carlson symmetric forms of elliptic integrals.
!
!  These routines are refactored versions of the ones from [SLATEC](http://www.netlib.org/slatec/).
!  They have been converted into modern Fortran, and the documentation has been
!  converted to FORD syntax.
! From https://github.com/jacobwilliams/carlson-elliptic-integrals
module carlson_elliptic_module

  use iso_fortran_env, only: error_unit, wp => real64

  implicit none

  private

  integer,parameter,public :: carlson_elliptic_module_wp = wp

  !**************************************************************
  !>
  !  Machine constants (replaces the old SLATEC [D1MACH](http://www.netlib.org/slatec/src/d1mach.f) function)
  !
  !  The traditional D1MACH constants are:
  !  * `D1MACH( 1) = B**(EMIN-1)`,           the smallest positive magnitude.
  !  * `D1MACH( 2) = B**EMAX*(1 - B**(-T))`, the largest magnitude.
  !  * `D1MACH( 3) = B**(-T)`,               the smallest relative spacing.
  !  * `D1MACH( 4) = B**(1-T)`,              the largest relative spacing.
  !  * `D1MACH( 5) = LOG10(B)`

  real(wp),dimension(5),parameter :: d1mach = &
      [  tiny(1.0_wp), &
         huge(1.0_wp), &
         real(radix(1.0_wp),wp)**(-digits(1.0_wp)), &
         epsilon(1.0_wp), &
         log10(real(radix(1.0_wp),wp)) ]

  !**************************************************************

  public :: drf,drd

  contains
!*******************************************************************************



  !*******************************************************************************
!>
!  Compute an approximation for the incomplete or
!  complete elliptic integral of the 2nd kind:
!  $$ R_D(x,y,z) = \frac{3}{2} \int_{0}^{\infty}
!                       (t+x)^{-1/2}
!                       (t+y)^{-1/2}
!                       (t+z)^{-3/2} dt $$
!  Where \(x\ge0\), \(y\ge0\), \(x+y>0\), and \(z>0\).
!
!  If \(x=0\) or \(y=0\), the integral is complete.
!
!  The duplication theorem is iterated until the variables are
!  nearly equal, and the function is then expanded in Taylor
!  series to fifth order.
!
!### DRD Special Comments
!
!  $$
!    \begin{array}{rl}
!      R_D(x,y,z) + R_D(y,z,x) + R_D(z,x,y) = \frac{3}{\sqrt{x y z}}  &  x>0, y>0, z>0
!    \end{array}
!  $$
!
!### Special functions via DRD and DRF
!
!  * Legendre form of ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      E(\phi,k) = \sin \phi  R_F(\cos^2 \phi,1-k^2 \sin^2 \phi,1)
!      -\frac{k^2}{3} \sin^3 \phi R_D(\cos^2 \phi,1-k^2 \sin^2 \phi,1)
!    $$
!    When \( \phi = \pi /2 \) the integral is complete:
!    $$
!      \begin{array}{rcl}
!       E(k) &=& R_F(0,1-k^2 ,1) - \frac{k^2}{3} R_D(0,1-k^2 ,1) \\
!            &=& \int_{0}^{\pi/2} \sqrt{1-k^2 \sin^2 \phi}  d \phi
!      \end{array}
!    $$
!
!  * Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      \mathrm{EL2}(x,k_c,a,b) = ax R_F(1,1+k_c^2 x^2 ,1+x^2 )
!        + \frac{1}{3}(b-a) x^3 R_D(1,1+k_c^2 x^2 ,1+x^2 )
!    $$
!
!  * Legendre form of alternative ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      \begin{array}{rcl}
!      D(q,k) &=& \int_{0}^{q} \sin^2 p (1-k^2 \sin^2 p)^{-1/2} dp \\
!             &=& \frac{1}{3} (\sin^3 q) R_D(\cos^2 q,1-k^2 \sin^2 q,1)
!      \end{array}
!    $$
!
!  * Lemniscate constant B:
!
!    $$
!      \begin{array}{rcl}
!      B &=& \int_{0}^{1} s^2 (1-s^4)^{-1/2} ds \\
!        &=& \frac{1}{3} R_D (0,2,1)
!      \end{array}
!    $$
!
!  * Heuman's LAMBDA function:
!
!    $$
!      \begin{array}{rcl}
!      \frac{\pi}{2} \Lambda_0(a,b) &=&
!      \sin b \left(R_F(0,\cos^2 a,1)-\frac{1}{3} \sin^2 a
!      R_D(0,\cos^2 a,1) \right) R_F(\cos^2 b,1-\cos^2 a \sin^2 b,1) \\
!      & &-\frac{1}{3} \cos^2 a \sin^3 b R_F(0,\cos^2 a,1)
!      R_D(\cos^2 b,1-\cos^2 a \sin^2 b,1)
!      \end{array}
!    $$
!
!  * Jacobi ZETA function:
!
!    $$
!      \begin{array}{rcl}
!      Z(b,k) &=& \frac{k^2}{3} \sin b R_F(\cos^2 b,1-k^2 \sin^2 b,1)
!               R_D(0,1-k^2 ,1)/R_F(0,1-k^2 ,1) \\
!                & & -\frac{k^2}{3} \sin^3 b R_D(\cos^2 b,1-k^2 \sin^2 b,1)
!      \end{array}
!    $$
!
!### Authors
!  * Carlson, B. C. Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Notis, E. M., Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Pexton, R. L., Lawrence Livermore National Laboratory, Livermore, CA  94550
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Modify calls to XERMSG to put in standard form.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drd.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

  real(wp) function drd(x,y,z,ier)

  implicit none

  real(wp),intent(in) :: x    !! nonnegative variable (\(x+y>0\))
  real(wp),intent(in) :: y    !! nonnegative variable (\(x+y>0\))
  real(wp),intent(in) :: z    !! positive variable
  integer,intent(out) :: ier  !! indicates normal or abnormal termination:
                              !!
                              !! * `IER = 0`: Normal and reliable termination of the
                              !!   routine. It is assumed that the requested
                              !!   accuracy has been achieved.
                              !! * `IER > 0`: Abnormal termination of the routine:
                              !! * `IER = 1`: `min(x,y) < 0`
                              !! * `IER = 2`: `min(x + y, z ) < LOLIM`
                              !! * `IER = 3`: `max(x,y,z) > UPLIM`

  character(len=16) :: xern3 , xern4 , xern5 , xern6
  real(wp) :: epslon, ea , eb , ec , ed , ef , lamda
  real(wp) :: mu , power4 , sigma , s1 , s2 , xn , xndev
  real(wp) :: xnroot , yn , yndev , ynroot , zn , zndev , znroot

  real(wp),parameter :: errtol = (d1mach(3)/3.0_wp)**(1.0_wp/6.0_wp)
      !! Determines the accuracy of the answer.
      !! The value assigned by the routine will result
      !! in solution precision within 1-2 decimals of
      !! machine precision.
      !!
      !! Relative error due to truncation is less than
      !! `3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2`.
      !!
      !! The accuracy of the computed approximation to the integral
      !! can be controlled by choosing the value of ERRTOL.
      !! Truncation of a Taylor series after terms of fifth order
      !! introduces an error less than the amount shown in the
      !! second column of the following table for each value of
      !! ERRTOL in the first column.  In addition to the truncation
      !! error there will be round-off error, but in practice the
      !! total error from both sources is usually less than the
      !! amount given in the table.
      !!
      !! Sample choices:
      !! (ERRTOL, Relative truncation error less than):
      !! (1.0e-3, 4.0e-18),
      !! (3.0e-3, 3.0e-15),
      !! (1.0e-2, 4.0e-12),
      !! (3.0e-2, 3.0e-9),
      !! (1.0e-1, 4.0e-6)
      !!
      !! Decreasing ERRTOL by a factor of 10 yields six more
      !! decimal digits of accuracy at the expense of one or
      !! two more iterations of the duplication theorem.

  real(wp),parameter :: lolim  = 2.0_wp/(d1mach(2))**(2.0_wp/3.0_wp) !! Lower limit of valid arguments
  real(wp),parameter :: tuplim = (0.10_wp*errtol)**(1.0_wp/3.0_wp)/&
                                  d1mach(1)**(1.0_wp/3.0_wp)
  real(wp),parameter :: uplim  = tuplim**2  !! Upper limit of valid arguments
  real(wp),parameter :: c1     = 3.0_wp/14.0_wp
  real(wp),parameter :: c2     = 1.0_wp/6.0_wp
  real(wp),parameter :: c3     = 9.0_wp/22.0_wp
  real(wp),parameter :: c4     = 3.0_wp/26.0_wp

  ! initialize:
  drd = 0.0_wp

  ! check for errors:
  if ( min(x,y)<0.0_wp ) then
      ier = 1
      write (xern3,'(1pe15.6)') x
      write (xern4,'(1pe15.6)') y
      write(error_unit,'(a)') 'drd: min(x,y)<0 where x = '//xern3// &
                              ' and y = '//xern4
      return
  endif

  if ( max(x,y,z)>uplim ) then
      ier = 3
      write (xern3,'(1pe15.6)') x
      write (xern4,'(1pe15.6)') y
      write (xern5,'(1pe15.6)') z
      write (xern6,'(1pe15.6)') uplim
      write(error_unit,'(a)') 'drd: max(x,y,z)>uplim where x = '// &
                              xern3//' y = '//xern4//' z = '//xern5// &
                              ' and uplim = '//xern6
      return
  endif

  if ( min(x+y,z)<lolim ) then
      ier = 2
      write (xern3,'(1pe15.6)') x
      write (xern4,'(1pe15.6)') y
      write (xern5,'(1pe15.6)') z
      write (xern6,'(1pe15.6)') lolim
      write(error_unit,'(a)') 'drd: min(x+y,z)<lolim where x = '// &
                              xern3//' y = '//xern4//' z = '//xern5// &
                              ' and lolim = '//xern6
      return
  endif

  ier    = 0
  xn     = x
  yn     = y
  zn     = z
  sigma  = 0.0_wp
  power4 = 1.0_wp

  do
      mu     = (xn+yn+3.0_wp*zn)*0.20_wp
      xndev  = (mu-xn)/mu
      yndev  = (mu-yn)/mu
      zndev  = (mu-zn)/mu
      epslon = max(abs(xndev),abs(yndev),abs(zndev))
      if ( epslon<errtol ) exit
      xnroot = sqrt(xn)
      ynroot = sqrt(yn)
      znroot = sqrt(zn)
      lamda  = xnroot*(ynroot+znroot) + ynroot*znroot
      sigma  = sigma + power4/(znroot*(zn+lamda))
      power4 = power4*0.250_wp
      xn     = (xn+lamda)*0.250_wp
      yn     = (yn+lamda)*0.250_wp
      zn     = (zn+lamda)*0.250_wp
  end do

  ea  = xndev*yndev
  eb  = zndev*zndev
  ec  = ea - eb
  ed  = ea - 6.0_wp*eb
  ef  = ed + ec + ec
  s1  = ed*(-c1+0.250_wp*c3*ed-1.50_wp*c4*zndev*ef)
  s2  = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea))
  drd = 3.0_wp*sigma + power4*(1.0_wp+s1+s2)/(mu*sqrt(mu))

  end function drd
!*******************************************************************************

!*******************************************************************************
!>
!  Compute an approximation for the incomplete or
!  complete elliptic integral of the 1st kind:
!  $$ R_F(x,y,z) = \frac{1}{2} \int_{0}^{\infty}
!                       (t+x)^{-1/2}
!                       (t+y)^{-1/2}
!                       (t+z)^{-1/2} dt $$
!  Where \(x\ge0\), \(y\ge0\), \(z\ge0\), and at most one of
!  them is \(=0\).
!
!  If \(x=0\), \(y=0\), or \(z=0\), the integral is complete.
!
!  The duplication theorem is iterated until the variables are
!  nearly equal, and the function is then expanded in Taylor
!  series to fifth order.
!
!### DRF Special Comments
!
!  $$
!    \begin{array}{rl}
!    R_F(x,x+z,x+w) + R_F(y,y+z,y+w) = R_F(0,z,w)
!    & x>0, y>0, z>0, x y = z w
!    \end{array}
!  $$
!
!### Special functions via DRF
!
!  * Legendre form of ELLIPTIC INTEGRAL of 1st kind:
!
!    $$
!    \begin{array}{rl}
!      F(\phi,k) &= \sin \phi R_F( \cos^2 \phi,1-k^2 \sin^2 \phi,1) \\
!           K(k) &= R_F(0,1-k^2 ,1) = \int_{0}^{\pi/2} (1-k^2 sin^2 \phi )^{-1/2} d \phi
!    \end{array}
!    $$
!
!  * Bulirsch form of ELLIPTIC INTEGRAL of 1st kind:
!
!    $$
!      \mathrm{EL1}(x,k_c) = x R_F(1,1+k_c^2 x^2 ,1+x^2 )
!    $$
!
!  * Lemniscate constant A:
!
!    $$
!      A = \int_{0}^{1} (1-s^4 )^{-1/2}    ds = R_F(0,1,2) = R_F(0,2,1)
!    $$
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891009  Removed unreferenced statement labels.  (WRB)
!  * 891009  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Changed calls to XERMSG to standard form, and some editorial changes.  (RWC))
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drf.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

  real(wp) function drf(x,y,z,ier)

  implicit none

  real(wp),intent(in) :: x    !! nonnegative variable
  real(wp),intent(in) :: y    !! nonnegative variable
  real(wp),intent(in) :: z    !! nonnegative variable
  integer,intent(out) :: ier  !! indicates normal or abnormal termination:
                              !!
                              !! * `IER = 0`: Normal and reliable termination of the
                              !!   routine. It is assumed that the requested
                              !!   accuracy has been achieved.
                              !! * `IER > 0`: Abnormal termination of the routine:
                              !! * `IER = 1`: `min(x,y,z) < 0`
                              !! * `IER = 2`:` min(x+y,x+z,y+z) < LOLIM`
                              !! * `IER = 3`: `max(x,y,z) > UPLIM`

  character(len=16) :: xern3 , xern4 , xern5 , xern6
  real(wp) :: epslon, e2 , e3 , lamda
  real(wp) :: mu , s , xn , xndev
  real(wp) :: xnroot , yn , yndev , ynroot , zn , zndev , znroot

  real(wp),parameter :: errtol = (4.0_wp*d1mach(3))**(1.0_wp/6.0_wp)
      !! Determines the accuracy of the answer.
      !! The value assigned by the routine will result
      !! in solution precision within 1-2 decimals of
      !! machine precision.
      !!
      !! Relative error due to truncation is less than
      !! `ERRTOL ** 6 / (4 * (1-ERRTOL)`.
      !!
      !! The accuracy of the computed approximation to the integral
      !! can be controlled by choosing the value of ERRTOL.
      !! Truncation of a Taylor series after terms of fifth order
      !! introduces an error less than the amount shown in the
      !! second column of the following table for each value of
      !! ERRTOL in the first column.  In addition to the truncation
      !! error there will be round-off error, but in practice the
      !! total error from both sources is usually less than the
      !! amount given in the table.
      !!
      !! Sample choices:
      !! (ERRTOL, Relative truncation error less than):
      !! (1.0e-3, 3.0e-19),
      !! (3.0e-3, 2.0e-16),
      !! (1.0e-2, 3.0e-13),
      !! (3.0e-2, 2.0e-10),
      !! (1.0e-1, 3.0e-7)
      !!
      !! Decreasing ERRTOL by a factor of 10 yields six more
      !! decimal digits of accuracy at the expense of one or
      !! two more iterations of the duplication theorem.

  real(wp),parameter :: lolim  = 5.0_wp*d1mach(1) !! Lower limit of valid arguments
  real(wp),parameter :: uplim  = d1mach(2)/5.0_wp !! Upper limit of valid arguments
  real(wp),parameter :: c1     = 1.0_wp/24.0_wp
  real(wp),parameter :: c2     = 3.0_wp/44.0_wp
  real(wp),parameter :: c3     = 1.0_wp/14.0_wp

  ! initialize:
  drf = 0.0_wp

  ! check for errors:
  if ( min(x,y,z)<0.0_wp ) then
      ier = 1
      write (xern3,'(1pe15.6)') x
      write (xern4,'(1pe15.6)') y
      write (xern5,'(1pe15.6)') z
      write(error_unit,'(a)') 'drf: min(x,y,z)<0 where x = '// &
              xern3//' y = '//xern4//' and z = '//xern5
      return
  endif

  if ( max(x,y,z)>uplim ) then
      ier = 3
      write (xern3,'(1pe15.6)') x
      write (xern4,'(1pe15.6)') y
      write (xern5,'(1pe15.6)') z
      write (xern6,'(1pe15.6)') uplim
      write(error_unit,'(a)') 'drf: max(x,y,z)>uplim where x = '// &
              xern3//' y = '//xern4//' z = '//xern5// &
              ' and uplim = '//xern6
      return
  endif

  if ( min(x+y,x+z,y+z)<lolim ) then
      ier = 2
      write (xern3,'(1pe15.6)') x
      write (xern4,'(1pe15.6)') y
      write (xern5,'(1pe15.6)') z
      write (xern6,'(1pe15.6)') lolim
      write(error_unit,'(a)') 'drf: min(x+y,x+z,y+z)<lolim where x = '//xern3// &
              ' y = '//xern4//' z = '//xern5//' and lolim = '//xern6
      return
  endif

  ier = 0
  xn  = x
  yn  = y
  zn  = z

  do
      mu     = (xn+yn+zn)/3.0_wp
      xndev  = 2.0_wp - (mu+xn)/mu
      yndev  = 2.0_wp - (mu+yn)/mu
      zndev  = 2.0_wp - (mu+zn)/mu
      epslon = max(abs(xndev),abs(yndev),abs(zndev))
      if ( epslon<errtol ) exit
      xnroot = sqrt(xn)
      ynroot = sqrt(yn)
      znroot = sqrt(zn)
      lamda  = xnroot*(ynroot+znroot) + ynroot*znroot
      xn     = (xn+lamda)*0.250_wp
      yn     = (yn+lamda)*0.250_wp
      zn     = (zn+lamda)*0.250_wp
  end do
  e2  = xndev*yndev - zndev*zndev
  e3  = xndev*yndev*zndev
  s   = 1.0_wp + (c1*e2-0.10_wp-c2*e3)*e2 + c3*e3
  drf = s/sqrt(mu)

  end function drf
!*******************************************************************************  

!*******************************************************************************
end module carlson_elliptic_module
!*******************************************************************************

!--------------------------------------------------------------------------------------------------
!     MODULE print
!>    @details Set up everything you need to perform particle tracing (electrons as of Jan 8 2024)
!!    @authors W. Villafana
!!    @date    Jan-8-2024
!--------------------------------------------------------------------------------------------------

MODULE ParticleTracing

  INTEGER, PARAMETER :: maximal_N_cells_tracing = 50
  INTEGER, PARAMETER :: maximal_N_part_to_trace_in_cell = 20

  LOGICAL save_e_vz
  LOGICAL save_e_irrEx
  LOGICAL save_e_irrEy
  LOGICAL save_e_extEz
  LOGICAL save_e_solEx
  LOGICAL save_e_solEy
  LOGICAL save_e_solEz
  LOGICAL save_e_solBx
  LOGICAL save_e_solBy
  LOGICAL save_e_solBz
  LOGICAL save_e_extBx
  LOGICAL save_e_extBy
  LOGICAL save_e_extBz

  INTEGER Tcntr_e_tracing_start
  INTEGER Tcntr_e_tracing_end
  INTEGER Tcntr_e_tracing_step
  INTEGER T_cntr_save_traced_electrons

  INTEGER e_record_length

  INTEGER N_of_cells_e_tracing

  TYPE ParticleTracingSetup
     INTEGER icell
     INTEGER jcell
     INTEGER N_to_trace
     LOGICAL choose_negative_VX
     LOGICAL choose_positive_VX
     LOGICAL choose_negative_VY
     LOGICAL choose_positive_VY
  END TYPE ParticleTracingSetup

  TYPE(ParticleTracingSetup) cell_e_tracing(1:maximal_N_cells_tracing)

END MODULE ParticleTracing

!--------------------------------------------------------------------------------------------------
!     MODULE print
!>    @details Stuff needed for Darwin (some variables present by default in particle tracing)
!!    @authors W. Villafana
!!    @date    Jan-8-2024
!--------------------------------------------------------------------------------------------------

MODULE SolenoidalFields

  REAL(8), ALLOCATABLE :: sol_EX(:,:)
  REAL(8), ALLOCATABLE :: sol_EY(:,:)
  REAL(8), ALLOCATABLE :: sol_EZ(:,:)

  REAL(8), ALLOCATABLE :: sol_BX(:,:)
  REAL(8), ALLOCATABLE :: sol_BY(:,:)
  REAL(8), ALLOCATABLE :: sol_BZ(:,:)

  REAL(8), ALLOCATABLE :: old_sol_EX(:,:)  ! these fields are used by master processes only
  REAL(8), ALLOCATABLE :: old_sol_EY(:,:)  ! they are updated in CALCULATE_BC_FOR_BZ

  REAL(8), ALLOCATABLE :: old_sol_BZ(:,:)

  REAL(8), ALLOCATABLE :: delta_JX(:,:)
  REAL(8), ALLOCATABLE :: delta_JY(:,:)
  REAL(8), ALLOCATABLE :: delta_JZ(:,:)
  
  TYPE electric_current_patch
     REAL(8) amplitude
     REAL(8) omega
     REAL(8) phase
     INTEGER imin
     INTEGER imax
     INTEGER jmin
     INTEGER jmax
  END TYPE electric_current_patch

  INTEGER N_of_JZ_patches
  TYPE(electric_current_patch), ALLOCATABLE :: JZpatch(:)

  TYPE contour_box
     INTEGER imin
     INTEGER jmin
     INTEGER imax
     INTEGER jmax
  END TYPE contour_box

  TYPE contour
     INTEGER material_object_number
!     INTEGER N_boxes ! presently 2, box 1 is the main contour, box 2 is the primed contour which has common side with the prime one but is smaller
     TYPE(contour_box) box(1:2)
  END TYPE contour

  INTEGER N_Faraday_contours
  TYPE(contour), ALLOCATABLE :: Faraday_contour(:)

  REAL(8), ALLOCATABLE :: EXEY_circulation_contour_field(:,:)
  REAL(8), ALLOCATABLE :: EXEY_circulation_primed_contour_field(:,:)

  REAL(8), ALLOCATABLE :: sol_EX_1(:,:,:)
  REAL(8), ALLOCATABLE :: sol_EY_1(:,:,:)

  REAL(8), ALLOCATABLE :: A(:)
  REAL(8), ALLOCATABLE :: B(:,:)

  REAL(8), ALLOCATABLE :: sol_BX_partial(:,:,:)
  REAL(8), ALLOCATABLE :: sol_BY_partial(:,:,:)

  REAL(8), ALLOCATABLE :: integral_solBXBY_loop_field(:,:)

  REAL(8), ALLOCATABLE :: sol_EZ_partial(:,:,:)

END MODULE SolenoidalFields

!--------------------------------------------------------------------------------------------------
!     MODULE IOs
!>    @details Deal with IOs subroutines
!!    @authors W. Villafana
!!    @date    Jan-29-2024
!--------------------------------------------------------------------------------------------------

MODULE mod_input_output
    
    USE CurrentProblemValues, ONLY: string_length
    IMPLICIT NONE

    CONTAINS
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE check_file_existence
!>    @details check if file exists and deletes it if it does. 
!!    @authors W. Villafana
!!    @date    Dec-29-2024
!-------------------------------------------------------------------------------------------------- 

  SUBROUTINE check_file_existence  ( filename )
        
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    !IN/OUT
    CHARACTER(LEN=string_length), INTENT(IN) :: filename
    
    ! LOCAL
    LOGICAL :: exists, ierr

    INQUIRE (FILE = TRIM(filename), EXIST = exists)
    IF (exists) CALL MPI_FILE_DELETE( TRIM(filename), MPI_INFO_NULL, ierr ) 

  END SUBROUTINE check_file_existence
 END MODULE mod_input_output

!--------------------------------------------------------------------------------------------------
!     MODULE numerics
!>    @details Useful numerical techniques, such as integrals
!!    @authors W. Villafana
!!    @date    Oct-28-2024
!--------------------------------------------------------------------------------------------------

MODULE mod_numerics

  USE CurrentProblemValues, ONLY: zero, half
  IMPLICIT NONE
  CONTAINS

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE mid_point_integration
!>    @details This subroutine calculates the integral of a function using the midpoint rule.
!!    @authors W. Villafana
!!    @date    Oct-28-2024
!     @input 
!       f       - Function to integrate (passed as an argument).
!       a, b    - Limits of integration.
!       N       - Number of intervals (must be positive).
!     @output
!       result  - Result of the integration (output).
!--------------------------------------------------------------------------------------------------   
  
  SUBROUTINE mid_point_integration(f, a, b, N, result)


    IMPLICIT NONE
    INTERFACE
      FUNCTION f(x)  ! Interface for the function to integrate
        REAL(8) :: f
        REAL(8), INTENT(IN) :: x
      END FUNCTION f
    END INTERFACE

    REAL(8), INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: N
    REAL(8), INTENT(OUT) :: result

    REAL(8) :: h, x_mid
    INTEGER :: i

    ! Check if N is valid
    IF (N <= 0) THEN
       PRINT *, "Error: N must be positive."
       STOP
    END IF

    ! Step size
    h = (b - a) / N

    ! Initialize the result
    result = zero

    ! Sum up the midpoint values
    DO i = 1, N
       x_mid = a + (i - half) * h
       result = result + f(x_mid)
    END DO

    ! Multiply by the step size to get the final result
    result = result * h

  END SUBROUTINE mid_point_integration

END MODULE mod_numerics  

!--------------------------------------------------------------------------------------------------
!     MODULE function_wrappers_1d
!>    @details wrapper for 1D functions
!!    @authors W. Villafana
!!    @date    Oct-28-2024
!--------------------------------------------------------------------------------------------------
MODULE mod_function_wrappers_1d

  USE CurrentProblemValues, ONLY: string_length, zero
  USE mod_print, ONLY: print_error
  IMPLICIT NONE
  CONTAINS

!--------------------------------------------------------------------------------------------------
!     FUNCTION evaluate_function
!>    @details This function selects and evaluates a predefined 1D function based on the specified 
!>             type. It acts as a central selector, allowing flexible choice among different 
!>             functions without modifying the calling code.
!!    @authors W. Villafana
!!    @date    Oct-28-2024
!     @input 
!       func_type - Character string specifying the function type ("bessel", "cylindrical", etc.).
!       x         - The input variable for the function, of type REAL(8).
!     @output
!       f_val     - The evaluated result of the function at the specified value x, of type REAL(8).
!--------------------------------------------------------------------------------------------------

  FUNCTION evaluate_function(func_type, x) RESULT(f_val)
    ! Selects a function to evaluate based on func_type
    CHARACTER(LEN=string_length), INTENT(IN) :: func_type
    REAL(8), INTENT(IN) :: x
    REAL(8) :: f_val
    CHARACTER(LEN=string_length) :: message, routine

    routine = 'evaluate_function'

    SELECT CASE (TRIM(func_type))
    CASE ("bessel_diffusion")
      f_val = bessel_diffusion(x)
    CASE ("cosh_sinh_diffusion")
      f_val = cosh_sinh_diffusion(x)
    CASE DEFAULT
      WRITE(message, '(A,A)') "Unknown function type. Only: 'bessel_diffusion', 'cosh_sinh_diffusion' are available. Received ",TRIM(func_type)
      CALL print_error(message,routine=routine)
    END SELECT
  END FUNCTION evaluate_function

  FUNCTION bessel_diffusion(x) RESULT(f_val)
    USE CurrentProblemValues, ONLY: Rw_bessel, bessel_first_zero_order_zero
    REAL(8), INTENT(IN) :: x
    REAL(8) :: f_val, alpha, cylindrical_volume_part
    alpha = bessel_first_zero_order_zero/Rw_bessel
    f_val = BESSEL_J0(x*alpha)
  END FUNCTION bessel_diffusion

  FUNCTION cosh_sinh_diffusion(y) RESULT(f_val)
    USE CurrentProblemValues, ONLY: Rw_bessel, bessel_first_zero_order_zero, Y_shift_exp_Y, Y0_exp
    REAL(8), INTENT(IN) :: y
    REAL(8) :: f_val, alpha, y_arg
    alpha = bessel_first_zero_order_zero/Rw_bessel
    y_arg = y-Y_shift_exp_Y
    f_val = COSH(alpha*y_arg)-COSH(alpha*Y0_exp)/SINH(alpha*Y0_exp)*SINH(alpha*y_arg)
  END FUNCTION cosh_sinh_diffusion

END MODULE mod_function_wrappers_1d
!--------------------------------------------------------------------------------------------------
!     MODULE array_utils
!>    @details Provides utility subroutines for array manipulation, including resizing arrays.
!!    @authors W. Villafana
!!    @date    Nov-19-2024
!--------------------------------------------------------------------------------------------------

MODULE array_utils

IMPLICIT NONE
CONTAINS

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE increase_array_size
!>    @details Dynamically increases the size of an allocatable array while preserving its content.
!>             If the array is not allocated, it initializes it to the specified size.
!!    @param[inout] array     The allocatable array to resize.
!!    @param[in]    new_size  The new size for the array.
!!    @authors      W. Villafana
!!    @date         Nov-19-2024
!--------------------------------------------------------------------------------------------------

  SUBROUTINE increase_array_size(array, new_size)
    INTEGER, ALLOCATABLE :: array(:)
    INTEGER :: new_size
    INTEGER, ALLOCATABLE :: temp(:)

    IF (.NOT. ALLOCATED(array)) THEN
       ALLOCATE(array(new_size))
    ELSE
       ALLOCATE(temp(SIZE(array)))
       temp = array  ! Store current values
       DEALLOCATE(array)
       ALLOCATE(array(new_size))
       array(:SIZE(temp)) = temp  ! Restore values
    ENDIF
  END SUBROUTINE increase_array_size

END MODULE array_utils
