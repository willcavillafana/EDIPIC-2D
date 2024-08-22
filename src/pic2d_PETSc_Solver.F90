!-----------------------------------
!
! This module initially was written by Janhunen Jani Salomon
! A number of changes listed below was made by DS (me)
!
module PETSc_Solver
  use ParallelOperationValues
  use CurrentProblemValues
  use BlockAndItsBoundaries
  use petsc
  implicit none
  private

#include <petsc/finclude/petsc.h>

  KSP ksp     ! solver object
  PC  pc      ! pre-conditioner object
  
  Mat Amat                ! Matrix type
  Vec bvec                ! Right hand side
  Vec xvec                ! Solution
  PetscInt ntoteq         ! total number of equations
  PetscInt nloceq         ! local number of equations
  
  integer :: parComm
  
  public SolverInitialization, Solve, InsertData, FetchData, SolverDestroy
  
contains
  
  subroutine Solve
    PetscErrorCode :: ierr

! sets KSP options from the options database
    call KSPSetFromOptions(ksp,ierr)
! solves linear system
    call KSPSolve(ksp, bvec, xvec, ierr)

    return
  end subroutine Solve
  
  subroutine SolverInitialization(comm)

    implicit none
    
    integer, intent(IN) ::  comm  !??? not used anywhere???

    PetscInt :: one, four, five, six

    character(20) :: petsc_config='petsc.rc'   ! PETSc database file
    PetscErrorCode :: ierr

    integer jbegin, jend, ibegin, iend
    PetscInt :: irow_global
    PetscInt :: jcolumn_global(1:6)
    PetscScalar :: value_at_jcol(1:6)
    integer i, j

    INTEGER nio, position_flag

    REAL(8), ALLOCATABLE :: eps_ishifted_j(:,:)
    REAL(8), ALLOCATABLE :: eps_i_jshifted(:,:)
    REAL(8) :: eps_shifted_quarter_2, eps_shifted_quarter, eps_shifted_half ! for points that are at a quarter of a segment
    INTEGER ALLOC_ERR
    REAL(8) :: factor_geom_cyl  ! additional factor for node on the right in the stencil
    REAL(8) :: factor_axis_geom_cyl  ! additional factor for node at the top and bottom in the stencil for axis.
    REAL(8) :: factor_geom_cyl_left, factor_geom_cyl_right
    REAL(8) :: dS1_dx, dS2_dx, dS3_dx, dS4_dx, rhs_coef, r_i ! geometrical coefs for Neuman boundary 
    LOGICAL :: neumann_flag ! to determine if current point is neumann or not

    ! By default we have no additional factors
    factor_geom_cyl = 1.0_8
    factor_axis_geom_cyl = 1.0_8
    factor_geom_cyl_left = 1.0_8
    factor_geom_cyl_right = 1.0_8
    rhs_coef = 1.0_8
   neumann_flag = .FALSE.
!    integer            :: m, n, nx, ny, i, j, k, ix, jy
!    PetscInt :: nrows, ncols, one=1, five=5, temp
!    PetscInt :: irow(5), jcol(5), tmpcol(5)
!    PetscScalar,parameter  :: wv(5)=(/ 0.25, 0.25, -1.0, 0.25, 0.25 /)
!    PetscScalar :: v(5)
!    PetscViewer :: viewer
!    PetscScalar :: tol
!    INTEGER local_x_max,local_x_min,local_y_min,local_y_max
!    INTEGER i_local,j_local

    one=1
    four=4
    five=5
    six = 6
    
! Initializes the petsc database and mpi
    call PetscInitialize(petsc_config, ierr)

! Remember the communicator for future reference
    parComm=PETSC_COMM_WORLD

! Creates a matrix where the type is determined from either a call to MatSetType() 
! or from the options database with a call to MatSetFromOptions()
! The default matrix type is AIJ
    call MatCreate(parComm, Amat, ierr)

    ntoteq=N_to_solve_total
    nloceq=block_N_of_nodes_to_solve   !######## ??????????????????????????? NEW ##########

! Set MPI distribution for A matrix [by S.J.]

! Sets the local and global sizes and checks to determine compatibility
! below we use PETSC_DECIDE instead of number of local rows/columns
!####??????    call MatSetSizes(Amat, PETSC_DECIDE, PETSC_DECIDE, ntoteq, ntoteq, ierr)
    call MatSetSizes(Amat, nloceq, nloceq, ntoteq, ntoteq, ierr)    !######## ??????????????????????????? NEW ##########

! Builds matrix object for a particular matrix type
! MATAIJ = "aij" - a matrix type to be used for sparce matrices
    call MatSetType(Amat, MATAIJ, ierr)

! Creates a matrix where the type is determined from the options database
    call MatSetFromOptions(Amat, ierr)
    five=5    !?????  why???

! For good matrix assembly performance the user should preallocate the matrix storage
! the second argument is the number of nonzeros per row (same for all rows)
    call MatSeqAIJSetPreallocation(Amat, six, PETSC_NULL_INTEGER, ierr)
    five=5    !????? why???

! Preallocates memory for a sparce parallel matrix in AIJ format (the default parallel petsc format)
! the second argument is the number of nonzeros per row in DIAGONAL portion of local submatrix
! the fourth argument is the number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix
    call MatMPIAIJSetPreallocation(Amat, six, PETSC_NULL_INTEGER, six, PETSC_NULL_INTEGER, ierr)
    one=1    !????  why???

! ny - the number of columns
! jcol - global indices of columns
! v - a logically two-dimensional array of values

!----------------------------------------
! Initialization of matrix coefficients written by DS
!

    jbegin = indx_y_min+1
    jend   = indx_y_max-1
    ibegin = indx_x_min+1
    iend   = indx_x_max-1

    IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
    IF (Rank_of_process_right.LT.0) iend   = indx_x_max
    IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
    IF (Rank_of_process_above.LT.0) jend   = indx_y_max

    ALLOCATE(eps_ishifted_j(ibegin:iend+1, jbegin:jend), STAT=ALLOC_ERR)   ! eps_ishifted_j(i,j) is between nodes {i-1,j} and {i,j}
    ALLOCATE(eps_i_jshifted(ibegin:iend, jbegin:jend+1), STAT=ALLOC_ERR)   ! eps_i_jshifted(i,j) is between nodes {i,j-1} and {i,j}

    eps_ishifted_j = 1.0_8
    eps_i_jshifted = 1.0_8

    DO j = jbegin, jend
       DO i = ibegin, iend+1
          CALL SET_EPS_ISHIFTED(i, j, eps_ishifted_j(i,j), 1 )
       END DO
    END DO

    DO j = jbegin, jend+1
       DO i = ibegin, iend
          CALL SET_EPS_JSHIFTED(i, j, eps_i_jshifted(i,j), 1 )
       END DO
    END DO

    irow_global = global_offset

!    j = indx_y_min !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
    IF (jbegin.EQ.indx_y_min) THEN
      j = indx_y_min
      ! boundary object along bottom border
      DO i = ibegin, iend
      ! print*,'value_bottom_is_1,i,j',i,j
         irow_global = irow_global + 1
         !          number_of_columns = 1
         jcolumn_global(1) = irow_global
         value_at_jcol(1) = 1.0_8
         call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
         ! If I have Neumann BCs, I need to compute coefficients 
         IF (block_has_neumann_bc_Y_bottom) THEN

         ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,j,neumann_flag)

            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN                    
               ! Left corner
               jcolumn_global(1) = irow_global + (iend-ibegin+1)    ! TOP
               ! this is at the left of the domain, along the BC
               IF ( i==indx_x_min .AND. ibegin==indx_x_min ) THEN

                  jcolumn_global(2) = irow_global                      ! CENTER         
                  jcolumn_global(3) = irow_global + 1                  ! RIGHT               
                  jcolumn_global(4) = irow_global + (iend-ibegin+1) + 1! TOP RIGHT    
                  
                  !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
                  ! Cartesian
                  dS1_dx = half
                  dS2_dx = half
                  rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS        
                  ! Cylindrical     
                  IF ( i_cylindrical==2 ) THEN
                     dS1_dx = pi*delta_x_m/4.0_8
                     dS2_dx = pi*delta_x_m/2.0_8
                     rhs_coef = pi*delta_x_m/4.0_8
                  END IF
                  CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) + 0.25_8, eps_shifted_quarter)   !right
                  CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j)+ 0.50_8, eps_shifted_quarter_2) !top                  

                  value_at_jcol(1) =    (eps_shifted_quarter_2*3.0_8/4.0_8*dS1_dx - eps_shifted_quarter*1.0_8/4.0_8*dS2_dx)/rhs_coef ! TOP
                  value_at_jcol(3) =  (- eps_shifted_quarter*1.0_8/4.0_8*dS2_dx + eps_shifted_quarter_2*3.0_8/4.0_8*dS1_dx)/rhs_coef ! RIGHT
                  value_at_jcol(4) =    (eps_shifted_quarter*1.0_8/4.0_8*dS2_dx + eps_shifted_quarter_2*1.0_8/4.0_8*dS1_dx)/rhs_coef ! TOP RIGHT
                  value_at_jcol(2) =  - (value_at_jcol(1)+value_at_jcol(3)+value_at_jcol(4))                                      ! CENTER

                  call MatSetValues(Amat, one, irow_global, 4, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 

               ! Left of the block. I might need to communicate with my left neighbor if it exists. 
               ELSE IF ( i==indx_x_min+1 ) THEN

                  ! No neighbor on the left, I can use my own node     ! LEFT
                  IF (indx_x_min==ibegin) THEN 
                     jcolumn_global(2) = irow_global - 1
                     jcolumn_global(6) = irow_global + (iend-ibegin+1) - 1! TOP LEFT
                  ELSE 
                     jcolumn_global(2) = process_left_bottom_right_inner_node + (j-indx_y_min-1) * process_left_solved_nodes_row_length  ! LEFT
                     jcolumn_global(6) = process_left_bottom_right_inner_node + (j-indx_y_min-1+1) * process_left_solved_nodes_row_length! TOP LEFT
                  END IF
                  jcolumn_global(3) = irow_global                      ! CENTER         
                  jcolumn_global(4) = irow_global + 1                  ! RIGHT
                  jcolumn_global(5) = irow_global + (iend-ibegin+1) + 1! TOP RIGHT
                                    

                  !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
                  ! Cartesian
                  dS1_dx = 1.0_8
                  dS2_dx = half
                  dS4_dx = half
                  rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS        
                  ! Cylindrical     
                  IF ( i_cylindrical==2 ) THEN
                     r_i = DBLE(i)*delta_x_m ! radius
                     dS1_dx = 2.0_8*pi*r_i*delta_x_m
                     dS2_dx = 2.0_8*pi*(r_i+delta_x_m)/2.0_8
                     dS4_dx = 2.0_8*pi*(r_i-delta_x_m)/2.0_8
                     rhs_coef = 2.0_8*pi*r_i
                  END IF
                  CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) + 0.25_8, eps_shifted_quarter)   !right
                  CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) + 0.25_8, eps_shifted_quarter_2) !left

                  value_at_jcol(1) =   eps_i_jshifted(i,j+1) - eps_shifted_quarter*1.0_8/4.0_8*dS2_dx/rhs_coef - eps_shifted_quarter_2*1.0_8/4.0_8*dS4_dx/ rhs_coef ! TOP
                  value_at_jcol(2) =   eps_shifted_quarter_2*3.0_8/4.0_8*dS4_dx/rhs_coef ! LEFT
                  value_at_jcol(4) =   eps_shifted_quarter*3.0_8/4.0_8*dS2_dx/rhs_coef ! RIGHT
                  value_at_jcol(5) =   eps_shifted_quarter*1.0_8/4.0_8*dS2_dx/rhs_coef ! TOP RIGHT
                  value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/4.0_8*dS4_dx/rhs_coef ! TOP LEFT
                  value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))                  

                  call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr)                   

               ! Right of the block. I might need to communicate with my right neighbor if it exists. 
               ELSE IF ( i==indx_x_max-1 ) THEN           
                  
                  jcolumn_global(2) = irow_global - 1                     ! LEFT
                  jcolumn_global(3) = irow_global                         ! CENTER         
                  ! No neighbor on the right, I can use my own node
                  IF ( indx_x_max==iend ) THEN
                     jcolumn_global(4) = irow_global + 1                  ! RIGHT  
                     jcolumn_global(5) = irow_global - (iend-ibegin+1) + 1! BOTTOM RIGHT
                  ! Neighbor on the right, I need to communicate
                  ELSE
                     jcolumn_global(4) = process_right_bottom_left_inner_node + (j-indx_y_min-1) * process_right_solved_nodes_row_length                 ! RIGHT
                     jcolumn_global(5) = process_right_bottom_left_inner_node + (j-indx_y_min-1+1) * process_right_solved_nodes_row_length! TOP RIGHT
                  END IF
                  jcolumn_global(6) = irow_global + (iend-ibegin+1) - 1! TOP LEFT                            

                  !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
                  ! Cartesian
                  dS1_dx = 1.0_8
                  dS2_dx = half
                  dS4_dx = half
                  rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS        
                  ! Cylindrical     
                  IF ( i_cylindrical==2 ) THEN
                     r_i = DBLE(i)*delta_x_m ! radius
                     dS1_dx = 2.0_8*pi*r_i
                     dS2_dx = 2.0_8*pi*(r_i+delta_x_m)/2.0_8
                     dS4_dx = 2.0_8*pi*(r_i-delta_x_m)/2.0_8
                     rhs_coef = 2.0_8*pi*r_i
                  END IF                  
                  CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) + 0.25_8, eps_shifted_quarter)   !right
                  CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) + 0.25_8, eps_shifted_quarter_2) !left
                  ! print*,'eps_shifted_quarter,eps_shifted_quarter_2',eps_shifted_quarter,eps_shifted_quarter_2

                  value_at_jcol(1) =   eps_i_jshifted(i,j+1) - eps_shifted_quarter*1.0_8/4.0_8*dS2_dx/rhs_coef - eps_shifted_quarter_2*1.0_8/4.0_8*dS4_dx/rhs_coef  ! TOP
                  value_at_jcol(2) =   eps_shifted_quarter_2*3.0_8/4.0_8*dS4_dx/rhs_coef ! LEFT
                  value_at_jcol(4) =   eps_shifted_quarter*3.0_8/4.0_8*dS2_dx/rhs_coef ! RIGHT
                  value_at_jcol(5) =   eps_shifted_quarter*1.0_8/4.0_8*dS2_dx/rhs_coef ! TOP RIGHT
                  value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/4.0_8*dS4_dx/rhs_coef ! TOP LEFT
                  value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))                  

                  call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr)                           

                  ! this is at the right of the domain, along the BC
               ELSE IF ( i==indx_x_max .AND. iend==indx_x_max ) THEN

                  jcolumn_global(2) = irow_global - 1               ! LEFT
                  jcolumn_global(3) = irow_global                   ! CENTER
                  jcolumn_global(4) = irow_global + (iend-ibegin+1) - 1! TOP LEFT                     

                  !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
                  ! Cartesian
                  dS1_dx = half
                  dS4_dx = half
                  rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS        
                  ! Cylindrical     
                  IF ( i_cylindrical==2 ) THEN
                     r_i = DBLE(i)*delta_x_m ! radius
                     dS1_dx = pi*(r_i-delta_x_m/4.0_8)
                     dS4_dx = pi*(r_i-delta_x_m/2.0_8)
                     rhs_coef = 2.0_8*pi*r_i
                  END IF             

                  ! Filling matrix
                  CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j)  + 0.25_8, eps_shifted_quarter)   !left
                  CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) + 0.5_8, eps_shifted_quarter_2)  !top       

                  value_at_jcol(1) =   (eps_shifted_quarter_2*3.0_8/4.0_8*dS1_dx - eps_shifted_quarter*1.0_8/4.0_8*dS4_dx)/rhs_coef ! TOP
                  value_at_jcol(2) =  (- eps_shifted_quarter_2*1.0_8/4.0_8*dS1_dx + eps_shifted_quarter*3.0_8/4.0_8*dS4_dx)/rhs_coef !LEFT
                  value_at_jcol(4) =   (eps_shifted_quarter*1.0_8/4.0_8*dS4_dx + eps_shifted_quarter_2*1.0_8/4.0_8*dS1_dx)/rhs_coef ! TOP LEFT
                  value_at_jcol(3) =   -(value_at_jcol(1)+value_at_jcol(2)+value_at_jcol(4)) ! CENTER


                  call MatSetValues(Amat, one, irow_global, 4, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr)                       

               ! This is a regular node at the top BC. Far from block boundaries and wall materials 
               ELSE

                  !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
                  ! Cartesian
                  dS1_dx = 1.0_8
                  dS2_dx = half
                  dS4_dx = half
                  rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS        
                  ! Cylindrical     
                  IF ( i_cylindrical==2 ) THEN
                     r_i = DBLE(i)*delta_x_m ! radius
                     dS1_dx = 2.0_8*pi*r_i
                     dS2_dx = 2.0_8*pi*(r_i+delta_x_m)/2.0_8
                     dS4_dx = 2.0_8*pi*(r_i-delta_x_m)/2.0_8
                     rhs_coef = 2.0_8*pi*r_i
                  END IF

                  jcolumn_global(2) = irow_global - 1                  ! LEFT
                  jcolumn_global(3) = irow_global                      ! CENTER         
                  jcolumn_global(4) = irow_global + 1                  ! RIGHT
                  jcolumn_global(5) = irow_global + (iend-ibegin+1) + 1! TOP RIGHT
                  jcolumn_global(6) = irow_global + (iend-ibegin+1) - 1! TOP LEFT
         
                  ! IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                  CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) + 0.25_8, eps_shifted_quarter)   !right
                  CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) + 0.25_8, eps_shifted_quarter_2) !left

                  value_at_jcol(1) =   eps_i_jshifted(i,j+1) - eps_shifted_quarter*1.0_8/4.0_8*dS2_dx/rhs_coef - eps_shifted_quarter_2*1.0_8/4.0_8*dS4_dx/ rhs_coef ! TOP
                  value_at_jcol(2) =   eps_shifted_quarter_2*3.0_8/4.0_8*dS4_dx/rhs_coef ! LEFT
                  value_at_jcol(4) =   eps_shifted_quarter*3.0_8/4.0_8*dS2_dx/rhs_coef ! RIGHT
                  value_at_jcol(5) =   eps_shifted_quarter*1.0_8/4.0_8*dS2_dx/rhs_coef ! TOP RIGHT
                  value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/4.0_8*dS4_dx/rhs_coef ! TOP LEFT

                  value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))   
                  call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr) 
               END IF ! IF loop over i index
            ENDIF ! If loop on neumann_flag
         END IF !IF (block_has_neumann_bc_Y_bottom) THEN     

         END DO
    END IF

!    j = indx_y_min+1 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    j = indx_y_min+1

    i = indx_x_min
    IF (ibegin.EQ.indx_x_min) THEN
      ! boundary object along left border
       irow_global = irow_global + 1
       IF ( .NOT.block_has_symmetry_plane_X_left .AND. .NOT.block_has_neumann_bc_X_left ) THEN ! this is a Cartesian or Cylindrical case with no symmetry (ie r_min>0 )
         ! Dirichlet (given potential) boundary
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       ELSE !This is either a Cartesian pr Cylindrical r-z case with symmetry axis
            ! the left border is a symmetry plane or Neumann BC
          IF (jbegin.EQ.indx_y_min) THEN                       ! BELOW
            ! boundary object along the bottom border
             jcolumn_global(1) = irow_global - (iend-ibegin+1)                          ! use the own node
             jcolumn_global(6) = irow_global - (iend-ibegin+1) + 1    ! BOTTOM RIGHT
          ELSE
            ! use a node from neighbor below
             jcolumn_global(1) = process_below_left_top_inner_node-1                    ! use a node from the neighbor below
             jcolumn_global(6) = process_below_left_top_inner_node-1 + 1    ! BOTTOM RIGHT
          END IF
          jcolumn_global(2) = irow_global                      ! CENTER
          jcolumn_global(3) = irow_global+1                    ! RIGHT
          jcolumn_global(4) = irow_global + (iend-ibegin+1)    ! ABOVE

            ! check whether the point is inside or at the surface of any inner object
          CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min, indx_y_min+1, nio, position_flag)

         ! note that we are at the left edge of the domain and the symmetry is applied here,
         ! so the boundary object must be symmetric relative to x=0 as well
         ! therefore we are either at the bottom surface, or inside, or at the top surface of the inner object
         SELECT CASE (position_flag)
            CASE (9)
            ! metal
               jcolumn_global(1) = irow_global
               value_at_jcol(1) = 1.0_8
               call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
               !  print*,'je passe par BC axis, metal'
               !             CASE (1,2)
!! dielectric surface above
!                value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
!             CASE (6,7)
!! dielectric surface below
!                value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr)              
            CASE DEFAULT 
! inside dielectric or plasma
               IF (i_cylindrical==2) THEN
                  factor_geom_cyl = two 
                  factor_axis_geom_cyl = one
                  rhs_coef = two
               END IF 
               IF (.NOT. block_has_neumann_bc_X_left ) THEN
                  value_at_jcol(1) =   eps_i_jshifted(i,j)*factor_axis_geom_cyl
                  value_at_jcol(2) = -(eps_i_jshifted(i,j)*factor_axis_geom_cyl + eps_i_jshifted(i,j+1)*factor_axis_geom_cyl + eps_ishifted_j(i+1,j)*factor_geom_cyl*rhs_coef + eps_ishifted_j(i+1,j)*factor_geom_cyl*rhs_coef)
                  value_at_jcol(3) =   (eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j))*factor_geom_cyl*rhs_coef ! this is 4 
                  value_at_jcol(4) =   eps_i_jshifted(i,j+1)*factor_axis_geom_cyl

                  call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) ! Same for both cylindrical r-z and Cartesian. left point= right point
               ELSE ! We have a Neumann BC. Cannot be cylindrical here 
                  ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
                  CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,j,neumann_flag)

                  ! This point is Neumann, I shall proceed
                  IF ( neumann_flag ) THEN        

                     jcolumn_global(5) = irow_global + (iend-ibegin+1) + 1    ! TOP RIGHT
                     
                     !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
                     ! Cartesian
                     dS1_dx = half
                     dS2_dx = 1.0_8
                     dS3_dx = half
                     rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS                          

                     CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) + 0.5_8, eps_shifted_quarter)   !top
                     CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) - 0.5_8, eps_shifted_quarter_2) !below    

                     value_at_jcol(1) =   3.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef ! BELOW
                     value_at_jcol(3) =  - 1.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef - 1.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef + eps_i_jshifted(i+1,j)  ! RIGHT
                     value_at_jcol(4) =   3.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef ! ABOVE*
                     value_at_jcol(5) =   1.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef ! TOP RIGHT
                     value_at_jcol(6) =   1.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef ! BOTTOM RIGHT

                     value_at_jcol(2) = -(value_at_jcol(1) + value_at_jcol(3) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6) ) ! CENTER
                     CALL MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr) 
                  ENDIF
               END IF                     
               !  print*,'axis,i,j',i,j
               !  print*,'value_at_jcol(1:4)',value_at_jcol(1:4)
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.25_8
               
         END SELECT

       END IF   !### IF (.NOT.block_has_symmetry_plane_X_left) THEN
    END IF      !### IF (ibegin.EQ.indx_x_min) THEN

!    i = indx_x_min+1

    i = indx_x_min+1
    irow_global = irow_global + 1

    IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
    ELSE

       IF (jbegin.EQ.indx_y_min) THEN                       ! BELOW
! boundary object along the bottom border
          jcolumn_global(1) = irow_global - (iend-ibegin+1)                          ! use the own node
       ELSE
! use a node from neighbor below
          jcolumn_global(1) = process_below_left_top_inner_node                      ! use a node from the neighbor below
       END IF
       IF (ibegin.EQ.indx_x_min) THEN                       ! LEFT
! boundary object along the left border
          jcolumn_global(2) = irow_global-1                                          ! use the own node
       ELSE
          jcolumn_global(2) = process_left_bottom_right_inner_node                   ! use a node from the left neighbor
       END IF
       jcolumn_global(3) = irow_global                      ! CENTER
       jcolumn_global(4) = irow_global+1                    ! RIGHT
       jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
       CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min+1, indx_y_min+1, nio, position_flag)

       SELECT CASE (position_flag)
          CASE (9)
! metal
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!          CASE (8)
!! dielectric surface on the right
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (4)
!! dielectric surface on the left
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (2)
!! dielectric surface above
!             value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (6)
!! dielectric surface below
!             value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
             IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
             value_at_jcol(1) =   eps_i_jshifted(i,j)
             value_at_jcol(2) =   eps_ishifted_j(i,j)
             value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
             value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
             value_at_jcol(5) =   eps_i_jshifted(i,j+1)
            !  print*,'inside,i,j',i,j
            !  print*,'value_inside',value_at_jcol(1:5)             
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.25_8
             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END SELECT
    END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

    IF (jbegin.EQ.indx_y_min) THEN
! boundary object along the bottom border
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
             jcolumn_global(2) = irow_global-1                    ! LEFT
             jcolumn_global(3) = irow_global                      ! CENTER
             jcolumn_global(4) = irow_global+1                    ! RIGHT
             jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, indx_y_min+1, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                   IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF    !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO    !### DO i = indx_x_min+2, indx_x_max-2

    ELSE   !### IF (jbegin.EQ.indx_y_min) THEN

! use a node from neighbor below
      DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = process_below_left_top_inner_node + (i-indx_x_min-1)  ! BELOW
             jcolumn_global(2) = irow_global-1                                         ! LEFT
             jcolumn_global(3) = irow_global                                           ! CENTER
             jcolumn_global(4) = irow_global+1                                         ! RIGHT
             jcolumn_global(5) = irow_global + (iend-ibegin+1)                         ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, indx_y_min+1, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                  IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO   !### DO i = indx_x_min+2, indx_x_max-2
    END IF   !### IF (jbegin.EQ.indx_y_min) THEN
 
    i = indx_x_max-1
    irow_global = irow_global + 1

    IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
    ELSE

       IF (jbegin.EQ.indx_y_min) THEN                       ! BELOW
! boundary object along the bottom border
          jcolumn_global(1) = irow_global - (iend-ibegin+1)                         ! use the own node
       ELSE
! use a node from neighbor below
          jcolumn_global(1) = process_below_left_top_inner_node + (i-indx_x_min-1)  ! use a node from the neighbor below
       END IF
       jcolumn_global(2) = irow_global-1                    ! LEFT
       jcolumn_global(3) = irow_global                      ! CENTER
       IF (iend.EQ.indx_x_max) THEN                         ! RIGHT
! boundary object along the right border
          jcolumn_global(4) = irow_global+1                                         ! use the own node
       ELSE
          jcolumn_global(4) = process_right_bottom_left_inner_node                  ! use a node from the right neighbor
       END IF
       jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
       CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_max-1, indx_y_min+1, nio, position_flag)

       SELECT CASE (position_flag)
          CASE (9)
! metal
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!          CASE (8)
!! dielectric surface on the right
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (4)
!! dielectric surface on the left
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (2)
!! dielectric surface above
!             value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (6)
!! dielectric surface below
!             value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
             IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
             value_at_jcol(1) =   eps_i_jshifted(i,j)
             value_at_jcol(2) =   eps_ishifted_j(i,j)
             value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
             value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
             value_at_jcol(5) =   eps_i_jshifted(i,j+1)
            !  print*,'value_at_jcol(1:4)_wil',value_at_jcol(1:4)
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.25_8
             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END SELECT
    END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

    i = indx_x_max

    IF (iend.EQ.indx_x_max) THEN
   ! boundary object along right border
       irow_global = irow_global + 1
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 

       ! If I have Neumann BCs, I need to compute coefficients. Also note that j = indx_y_min + 1 here ie, I might need to communicate with block below me if it exists
       IF (block_has_neumann_bc_X_right) THEN

         ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
         CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,j,neumann_flag)

         ! This point is Neumann, I shall proceed
         IF ( neumann_flag ) THEN         

            ! No neighbor below, I can use my own node
            IF ( indx_y_min==jbegin ) THEN
               jcolumn_global(1)   = irow_global - (iend-ibegin+1)    ! BELOW
               jcolumn_global(6)   = irow_global - (iend-ibegin+1) - 1   ! BOTTOM LEFT
            ! Ihave a neighbor, I must communicate
            ELSE            
               jcolumn_global(1)   = process_below_left_top_inner_node + (i-indx_x_min-1)    ! BELOW
               jcolumn_global(6)   = process_below_left_top_inner_node + (i-indx_x_min-1) - 1   ! BOTTOM LEFT
            END IF
            jcolumn_global(2) = irow_global - 1                  ! LEFT
            jcolumn_global(3) = irow_global                      ! CENTER         
            jcolumn_global(4) = irow_global + (iend-ibegin+1)    ! ABOVE
            jcolumn_global(5) = irow_global + (iend-ibegin+1) - 1    ! TOP LEFT
            

            !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
            ! Cartesian
            dS1_dx = half
            dS3_dx = half
            dS4_dx = 1.0_8
            rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS        
            ! Cylindrical     
            IF ( i_cylindrical==2 ) THEN
               r_i = DBLE(i)*delta_x_m ! radius
               dS1_dx = 2.0_8*pi*(r_i-delta_x_m/4.0_8)/2.0_8
               dS3_dx = 2.0_8*pi*(r_i-delta_x_m/4.0_8)/2.0_8
               dS4_dx = 2.0_8*pi*(r_i-delta_x_m/2.0_8)
               rhs_coef = 2.0_8*pi*(r_i-delta_x_m/4.0_8)
            END IF                  
            CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) + 0.5_8, eps_shifted_quarter)   !top
            CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) - 0.5_8, eps_shifted_quarter_2) !below    
            CALL SET_EPS_ISHIFTED(i, j, eps_shifted_half) !left

            value_at_jcol(1) =   3.0/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef  ! BELOW
            value_at_jcol(2) =   -eps_shifted_quarter*1.0_8/4.0_8*dS1_dx/rhs_coef - eps_shifted_quarter_2*1.0_8/4.0_8*dS3_dx/rhs_coef + eps_shifted_half*dS4_dx/rhs_coef ! LEFT
            value_at_jcol(4) =   eps_shifted_quarter*3.0_8/4.0_8*dS1_dx/rhs_coef ! ABOVE
            value_at_jcol(5) =   eps_shifted_quarter*1.0_8/4.0_8*dS1_dx/rhs_coef ! TOP LEFT
            value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/4.0_8*dS3_dx/rhs_coef ! BOTTOM LEFT
            value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))    ! CENTER 

            call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr) 
         END IF ! Neumann flag
       END IF         

    END IF

    DO j = indx_y_min+2, indx_y_max-2 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       i = indx_x_min

       IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
          irow_global = irow_global + 1
          IF (.NOT.block_has_symmetry_plane_X_left .AND. .NOT.block_has_neumann_bc_X_left) THEN
! Dirichlet (given potential) boundary
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
          ELSE
! the left border is a symmetry plane or Neumann BC
! use own nodes
             jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
             jcolumn_global(2) = irow_global                      ! CENTER
             jcolumn_global(3) = irow_global+1                    ! RIGHT
             jcolumn_global(4) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min, j, nio, position_flag)
       
! note that we are at the left edge of the domain and the symmetry is applied here,
! so the boundary object must be symmetric relative to x=0 as well
! therefore we are either at the bottom surface, or inside, or at the top surface of the inner object

            SELECT CASE (position_flag)
               CASE (9)
! metal
                  jcolumn_global(1) = irow_global
                  value_at_jcol(1) = 1.0_8
                  call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (1,2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = -1.0_8
!                   value_at_jcol(3) = 0.5_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
!                CASE (6,7)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = -1.0_8
!                   value_at_jcol(3) = 0.5_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr)              
               CASE DEFAULT 
! inside dielectric or plasma
                  IF (i_cylindrical==2) THEN
                     factor_geom_cyl = two 
                     factor_axis_geom_cyl = one
                     rhs_coef = two
                  END IF 
                  IF (.NOT. block_has_neumann_bc_X_left ) THEN
                     value_at_jcol(1) =   eps_i_jshifted(i,j)*factor_axis_geom_cyl
                     value_at_jcol(2) = -(eps_i_jshifted(i,j)*factor_axis_geom_cyl + eps_i_jshifted(i,j+1)*factor_axis_geom_cyl + eps_ishifted_j(i+1,j)*factor_geom_cyl*rhs_coef + eps_ishifted_j(i+1,j)*factor_geom_cyl*rhs_coef)
                     value_at_jcol(3) =   (eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j))*factor_geom_cyl*rhs_coef ! this is 4
                     value_at_jcol(4) =   eps_i_jshifted(i,j+1)*factor_axis_geom_cyl   
                     
                     call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
                  ELSE ! We have a Neumann BC. Cannot be cylindrical here 
                     
                     ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
                     CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,j,neumann_flag)

                     ! This point is Neumann, I shall proceed
                     IF ( neumann_flag ) THEN                             
                        jcolumn_global(5) = irow_global + (iend-ibegin+1) + 1    ! TOP RIGHT
                        jcolumn_global(6) = irow_global - (iend-ibegin+1) + 1    ! BOTTOM RIGHT
                        !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
                        ! Cartesian
                        dS1_dx = half
                        dS2_dx = 1.0_8
                        dS3_dx = half
                        rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS                          
      
                        CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) + 0.5_8, eps_shifted_quarter)   !top
                        CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) - 0.5_8, eps_shifted_quarter_2) !below    
      
                        value_at_jcol(1) =   3.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef ! BELOW
                        value_at_jcol(3) =  - 1.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef - 1.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef + eps_i_jshifted(i+1,j)  ! RIGHT
                        value_at_jcol(4) =   3.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef ! ABOVE*
                        value_at_jcol(5) =   1.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef ! TOP RIGHT
                        value_at_jcol(6) =   1.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef ! BOTTOM RIGHT
      
                        value_at_jcol(2) = -(value_at_jcol(1) + value_at_jcol(3) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6) ) ! CENTER
                        CALL MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr) 
                     ENDIF
                  END IF                                  
                  !  print*,'value_at_jcol(1:4)_wil2',value_at_jcol(1:4)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = -1.0_8
!                   value_at_jcol(3) = 0.5_8
!                   value_at_jcol(4) = 0.25_8
                  
            END SELECT

          END IF   !### IF (.NOT.block_has_symmetry_plane_X_left) THEN
       END IF      !### IF (ibegin.EQ.indx_x_min) THEN

!       i = indx_x_min+1

       i = indx_x_min+1
       irow_global = irow_global + 1

       IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
       ELSE

          jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
          IF (ibegin.EQ.indx_x_min) THEN                       ! LEFT
! boundary object along the left border
             jcolumn_global(2) = irow_global-1                                                                                     ! use the own node
          ELSE
             jcolumn_global(2) = process_left_bottom_right_inner_node + (j-indx_y_min-1) * process_left_solved_nodes_row_length     ! use a node from the left neighbor
          END IF
          jcolumn_global(3) = irow_global                      ! CENTER
          jcolumn_global(4) = irow_global+1                    ! RIGHT
          jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
          CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min+1, j, nio, position_flag)

          SELECT CASE (position_flag)
             CASE (9)
! metal
                jcolumn_global(1) = irow_global
                value_at_jcol(1) = 1.0_8
                call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!             CASE (8)
!! dielectric surface on the right
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(5) = 0.25_8
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (4)
!! dielectric surface on the left
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(5) = 0.25_8
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (2)
!! dielectric surface above
!                value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (6)
!! dielectric surface below
!                value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
               IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                value_at_jcol(1) =   eps_i_jshifted(i,j)
                value_at_jcol(2) =   eps_ishifted_j(i,j)
                value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
                value_at_jcol(5) =   eps_i_jshifted(i,j+1)
               !  print*,'value_at_jcol(1:4)_3',value_at_jcol(1:4)
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.25_8
                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          END SELECT
       END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
             jcolumn_global(2) = irow_global-1                    ! LEFT
             jcolumn_global(3) = irow_global                      ! CENTER
             jcolumn_global(4) = irow_global+1                    ! RIGHT
             jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, j, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                  IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
                  !  print*,'value_at_jcol(1:4)_wil4',value_at_jcol(1:4)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO   !### DO i = indx_x_min+2, indx_x_max-2

       i = indx_x_max-1
       irow_global = irow_global + 1

       IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
       ELSE

          jcolumn_global(1) = irow_global - (iend-ibegin+1)       ! BELOW
          jcolumn_global(2) = irow_global-1                       ! LEFT
          jcolumn_global(3) = irow_global                         ! CENTER
          IF (iend.EQ.indx_x_max) THEN                            ! RIGHT
! boundary object along the right border
             jcolumn_global(4) = irow_global+1                                                                                   ! use the own node
          ELSE
             jcolumn_global(4) = process_right_bottom_left_inner_node + (j-indx_y_min-1) * process_right_solved_nodes_row_length  ! use a node from the right neighbor
          END IF
          jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
          CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_max-1, j, nio, position_flag)

          SELECT CASE (position_flag)
             CASE (9)
! metal
                jcolumn_global(1) = irow_global
                value_at_jcol(1) = 1.0_8
                call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!             CASE (8)
!! dielectric surface on the right
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(5) = 0.25_8
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (4)
!! dielectric surface on the left
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(5) = 0.25_8
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (2)
!! dielectric surface above
!                value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (6)
!! dielectric surface below
!                value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
               IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                value_at_jcol(1) =   eps_i_jshifted(i,j)
                value_at_jcol(2) =   eps_ishifted_j(i,j)
                value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
                value_at_jcol(5) =   eps_i_jshifted(i,j+1)
               !  print*,'value_at_jcol(1:4)_wil5',value_at_jcol(1:4)
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.25_8
                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          END SELECT
       END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

       i = indx_x_max

       IF (iend.EQ.indx_x_max) THEN
         ! boundary object along right border
         irow_global = irow_global + 1
         jcolumn_global(1) = irow_global
         value_at_jcol(1) = 1.0_8
         call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 

         ! If I have Neumann BCs, I need to compute coefficients. 
         IF (block_has_neumann_bc_X_right) THEN

            ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,j,neumann_flag)

            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN            

               jcolumn_global(1)   = irow_global - (iend-ibegin+1)    ! BELOW
               jcolumn_global(6)   = irow_global - (iend-ibegin+1) - 1   ! BOTTOM LEFT
               jcolumn_global(2) = irow_global - 1                  ! LEFT
               jcolumn_global(3) = irow_global                      ! CENTER         
               jcolumn_global(4) = irow_global + (iend-ibegin+1)    ! ABOVE
               jcolumn_global(5) = irow_global + (iend-ibegin+1) - 1    ! TOP LEFT
               
               !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
               ! Cartesian
               dS1_dx = half
               dS3_dx = half
               dS4_dx = 1.0_8
               rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS        
               ! Cylindrical     
               IF ( i_cylindrical==2 ) THEN
                  r_i = DBLE(i)*delta_x_m ! radius
                  dS1_dx = 2.0_8*pi*(r_i-delta_x_m/4.0_8)/2.0_8
                  dS3_dx = 2.0_8*pi*(r_i-delta_x_m/4.0_8)/2.0_8
                  dS4_dx = 2.0_8*pi*(r_i-delta_x_m/2.0_8)
                  rhs_coef = 2.0_8*pi*(r_i-delta_x_m/4.0_8)
               END IF                  
               CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) + 0.5_8, eps_shifted_quarter)   !top
               CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) - 0.5_8, eps_shifted_quarter_2) !below    
               CALL SET_EPS_ISHIFTED(i, j, eps_shifted_half) !left
      
               value_at_jcol(1) =   3.0/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef  ! BELOW
               value_at_jcol(2) =   -eps_shifted_quarter*1.0_8/4.0_8*dS1_dx/rhs_coef - eps_shifted_quarter_2*1.0_8/4.0_8*dS3_dx/rhs_coef + eps_shifted_half*dS4_dx/rhs_coef ! LEFT
               value_at_jcol(4) =   eps_shifted_quarter*3.0_8/4.0_8*dS1_dx/rhs_coef ! ABOVE
               value_at_jcol(5) =   eps_shifted_quarter*1.0_8/4.0_8*dS1_dx/rhs_coef ! TOP LEFT
               value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/4.0_8*dS3_dx/rhs_coef ! BOTTOM LEFT
               value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))    ! CENTER 
      
               call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr) 
            END IF ! Neumann flag
         END IF         
   
       END IF

    END DO !### DO j = indx_y_min+2, indx_y_max-2

    j = indx_y_max-1 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    i = indx_x_min

    IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
       irow_global = irow_global + 1
       IF (.NOT.block_has_symmetry_plane_X_left .AND. .NOT.block_has_neumann_bc_X_left ) THEN
! Dirichlet (given potential) boundary
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       ELSE
! the left border is a symmetry plane or Neumann BC
          jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
          jcolumn_global(2) = irow_global                      ! CENTER
          jcolumn_global(3) = irow_global+1                    ! RIGHT
          IF (jend.EQ.indx_y_max) THEN                         ! ABOVE
! boundary object along the top border
             jcolumn_global(4) = irow_global + (iend-ibegin+1)                         ! use the own node
             jcolumn_global(6) = irow_global + (iend-ibegin+1) + 1    ! TOP RIGHT
          ELSE
             jcolumn_global(4) = process_above_left_bottom_inner_node-1                ! use a node from the neighbor above
             jcolumn_global(6) = process_above_left_bottom_inner_node-1 + 1    ! TOP RIGHT
          END IF

! check whether the point is inside or at the surface of any inner object
          CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min, indx_y_max-1, nio, position_flag)

! note that we are at the left edge of the domain and the symmetry is applied here,
! so the boundary object must be symmetric relative to x=0 as well
! therefore we are either at the bottom surface, or inside, or at the top surface of the inner object

         SELECT CASE (position_flag)
            CASE (9)
! metal
               jcolumn_global(1) = irow_global
               value_at_jcol(1) = 1.0_8
               call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!             CASE (1,2)
!! dielectric surface above
!                value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
!             CASE (6,7)
!! dielectric surface below
!                value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr)              
            CASE DEFAULT 
! inside dielectric or plasma
               IF (i_cylindrical==2) THEN
                  factor_geom_cyl = two 
                  factor_axis_geom_cyl = one
                  rhs_coef = two
               END IF                
               IF (.NOT. block_has_neumann_bc_X_left ) THEN
                  value_at_jcol(1) =   eps_i_jshifted(i,j)*factor_axis_geom_cyl
                  value_at_jcol(2) = -(eps_i_jshifted(i,j)*factor_axis_geom_cyl + eps_i_jshifted(i,j+1)*factor_axis_geom_cyl + eps_ishifted_j(i+1,j)*factor_geom_cyl*rhs_coef + eps_ishifted_j(i+1,j)*factor_geom_cyl*rhs_coef)
                  value_at_jcol(3) =   (eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j))*factor_geom_cyl*rhs_coef ! this is 4 
                  value_at_jcol(4) =   eps_i_jshifted(i,j+1)*factor_axis_geom_cyl

                  call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
               ELSE ! We have a Neumann BC. Cannot be cylindrical here 

                  ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
                  CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,j,neumann_flag)

                  ! This point is Neumann, I shall proceed
                  IF ( neumann_flag ) THEN                          
                     jcolumn_global(6) = irow_global - (iend-ibegin+1) + 1    ! BOTTOM RIGHT
                     
                     !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
                     ! Cartesian
                     dS1_dx = half
                     dS2_dx = 1.0_8
                     dS3_dx = half
                     rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS                          

                     CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) + 0.5_8, eps_shifted_quarter)   !top
                     CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) - 0.5_8, eps_shifted_quarter_2) !below    

                     value_at_jcol(1) =   3.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef ! BELOW
                     value_at_jcol(3) =  - 1.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef - 1.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef + eps_i_jshifted(i+1,j)  ! RIGHT
                     value_at_jcol(4) =   3.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef ! ABOVE*
                     value_at_jcol(5) =   1.0_8/4.0_8*eps_shifted_quarter*dS1_dx/rhs_coef ! TOP RIGHT
                     value_at_jcol(6) =   1.0_8/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef ! BOTTOM RIGHT

                     value_at_jcol(2) = -(value_at_jcol(1) + value_at_jcol(3) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6) ) ! CENTER
                     CALL MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr)                   
                  ENDIF
               END IF
               !  print*,'value_at_jcol(1:4)_wil6',value_at_jcol(1:4)
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.25_8
               
         END SELECT

       END IF
    END IF

!    i = indx_x_min+1

    i = indx_x_min+1
    irow_global = irow_global + 1

    IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
    ELSE

       jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
       IF (ibegin.EQ.indx_x_min) THEN                       ! LEFT
! boundary object along the left border
          jcolumn_global(2) = irow_global-1                                                                                      ! use the own node
       ELSE
          jcolumn_global(2) = process_left_bottom_right_inner_node + (j-indx_y_min-1) * process_left_solved_nodes_row_length     ! use a node from the left neighbor
       END IF
       jcolumn_global(3) = irow_global                      ! CENTER
       jcolumn_global(4) = irow_global+1                    ! RIGHT
       IF (jend.EQ.indx_y_max) THEN                         ! ABOVE
! boundary object along the top border
          jcolumn_global(5) = irow_global + (iend-ibegin+1)                         ! use the own node
       ELSE
          jcolumn_global(5) = process_above_left_bottom_inner_node                  ! use a node from the neighbor above
       END IF

! check whether the point is inside or at the surface of any inner object
       CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min+1, indx_y_max-1, nio, position_flag)

       SELECT CASE (position_flag)
          CASE (9)
! metal
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!          CASE (8)
!! dielectric surface on the right
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (4)
!! dielectric surface on the left
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (2)
!! dielectric surface above
!             value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (6)
!! dielectric surface below
!             value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
            IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
             value_at_jcol(1) =   eps_i_jshifted(i,j)
             value_at_jcol(2) =   eps_ishifted_j(i,j)
             value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
             value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
             value_at_jcol(5) =   eps_i_jshifted(i,j+1)
            !  print*,'value_at_jcol(1:4)_wil7',value_at_jcol(1:4)
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.25_8
             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END SELECT
    END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

    IF (jend.EQ.indx_y_max) THEN
! boundary object along the top border
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
             jcolumn_global(2) = irow_global-1                    ! LEFT
             jcolumn_global(3) = irow_global                      ! CENTER
             jcolumn_global(4) = irow_global+1                    ! RIGHT
             jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, indx_y_max-1, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                  IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
                  !  print*,'value_at_jcol(1:4)_wil8',value_at_jcol(1:4)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO   !### DO i = indx_x_min+2, indx_x_max-2
    ELSE   !### IF (jend.EQ.indx_y_max) THEN
! use a node from the neighbor above
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = irow_global - (iend-ibegin+1)                            ! BELOW
             jcolumn_global(2) = irow_global-1                                            ! LEFT
             jcolumn_global(3) = irow_global                                              ! CENTER
             jcolumn_global(4) = irow_global+1                                            ! RIGHT
             jcolumn_global(5) = process_above_left_bottom_inner_node + (i-indx_x_min-1)  ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, indx_y_max-1, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                  IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
                  !  print*,'value_at_jcol(1:4)_wil9',value_at_jcol(1:4)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO   !### DO i = indx_x_min+2, indx_x_max-2
    END IF   !### IF (jend.EQ.indx_y_max) THEN
 
    i = indx_x_max-1
    irow_global = irow_global + 1

    IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
    ELSE

       jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
       jcolumn_global(2) = irow_global-1                    ! LEFT
       jcolumn_global(3) = irow_global                      ! CENTER
       IF (iend.EQ.indx_x_max) THEN                         ! RIGHT
! boundary object along the right border
          jcolumn_global(4) = irow_global+1                                                                                    ! use the own node
       ELSE
          jcolumn_global(4) = process_right_bottom_left_inner_node + (j-indx_y_min-1) * process_right_solved_nodes_row_length  ! use a node from the right neighbor
       END IF
       IF (jend.EQ.indx_y_max) THEN                         ! ABOVE
! boundary object along the top border
          jcolumn_global(5) = irow_global + (iend-ibegin+1)                            ! use the own node
       ELSE
          jcolumn_global(5) = process_above_left_bottom_inner_node + (i-indx_x_min-1)  ! use a node from the neighbor above
       END IF

! check whether the point is inside or at the surface of any inner object
       CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_max-1, indx_y_max-1, nio, position_flag)

       SELECT CASE (position_flag)
          CASE (9)
! metal
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!          CASE (8)
!! dielectric surface on the right
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (4)
!! dielectric surface on the left
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (2)
!! dielectric surface above
!             value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (6)
!! dielectric surface below
!             value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
            IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
             value_at_jcol(1) =   eps_i_jshifted(i,j)
             value_at_jcol(2) =   eps_ishifted_j(i,j)
             value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
             value_at_jcol(4) =   eps_ishifted_j(i+1,j)*factor_geom_cyl
             value_at_jcol(5) =   eps_i_jshifted(i,j+1)
            !  print*,'value_at_jcol(1:4)_wil_10',value_at_jcol(1:4)
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.25_8
             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END SELECT
    END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

    i = indx_x_max

    IF (iend.EQ.indx_x_max) THEN
      ! boundary object along right border
      irow_global = irow_global + 1
      jcolumn_global(1) = irow_global
      value_at_jcol(1) = 1.0_8
      call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 

      ! If I have Neumann BCs, I need to compute coefficients. Also note that j = indx_y_max -1 1 here ie, I might need to communicate with block below me if it exists
      IF (block_has_neumann_bc_X_right) THEN

         ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
         CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,j,neumann_flag)

         ! This point is Neumann, I shall proceed
         IF ( neumann_flag ) THEN
            ! No neighbor below, I can use my own node
            jcolumn_global(1)   = irow_global - (iend-ibegin+1)    ! BELOW
            jcolumn_global(6)   = irow_global - (iend-ibegin+1) - 1   ! BOTTOM LEFT
            IF ( indx_y_max==jend ) THEN
               jcolumn_global(4) = irow_global + (iend-ibegin+1)    ! ABOVE
               jcolumn_global(5) = irow_global + (iend-ibegin+1) - 1    ! TOP LEFT            
            ! Ihave a neighbor, I must communicate
            ELSE            
               jcolumn_global(4) = process_above_left_bottom_inner_node + (i-indx_x_min-1)    ! ABOVE
               jcolumn_global(5) = process_above_left_bottom_inner_node + (i-indx_x_min-1) - 1    ! TOP LEFT
            END IF
            jcolumn_global(2) = irow_global - 1                  ! LEFT
            jcolumn_global(3) = irow_global                      ! CENTER         

            !!! I define geometrical coefs: 1= top, 2 = right, 3 = bottom, 4 = left
            ! Cartesian
            dS1_dx = half
            dS3_dx = half
            dS4_dx = 1.0_8
            rhs_coef = 1.0_8 ! we will have dx**2/2 for the volume RHS        
            ! Cylindrical     
            IF ( i_cylindrical==2 ) THEN
               r_i = DBLE(i)*delta_x_m ! radius
               dS1_dx = 2.0_8*pi*(r_i-delta_x_m/4.0_8)/2.0_8
               dS3_dx = 2.0_8*pi*(r_i-delta_x_m/4.0_8)/2.0_8
               dS4_dx = 2.0_8*pi*(r_i-delta_x_m/2.0_8)
               rhs_coef = 2.0_8*pi*(r_i-delta_x_m/4.0_8)
            END IF                  
            CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) + 0.5_8, eps_shifted_quarter)   !top
            CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) - 0.5_8, eps_shifted_quarter_2) !below    
            CALL SET_EPS_ISHIFTED(i, j, eps_shifted_half) !left

            value_at_jcol(1) =   3.0/4.0_8*eps_shifted_quarter_2*dS3_dx/rhs_coef  ! BELOW
            value_at_jcol(2) =   -eps_shifted_quarter*1.0_8/4.0_8*dS1_dx/rhs_coef - eps_shifted_quarter_2*1.0_8/4.0_8*dS3_dx/rhs_coef + eps_shifted_half*dS4_dx/rhs_coef ! LEFT
            value_at_jcol(4) =   eps_shifted_quarter*3.0_8/4.0_8*dS1_dx/rhs_coef ! ABOVE
            value_at_jcol(5) =   eps_shifted_quarter*1.0_8/4.0_8*dS1_dx/rhs_coef ! TOP LEFT
            value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/4.0_8*dS3_dx/rhs_coef ! BOTTOM LEFT
            value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))    ! CENTER 

            call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr) 
         END IF
      END IF           

    END IF

!    j = indx_y_max !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.

    j = indx_y_max

    IF (jend.EQ.indx_y_max) THEN
      ! boundary object along top border
       DO i = ibegin, iend
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 

       ! If I have Neumann BCs, I need to compute coefficients 
          IF (block_has_neumann_bc_Y_top) THEN

            ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,j,neumann_flag)

            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN
               ! Left corner
               jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
               ! this is at the left of the domain, along the BC
               IF ( i==indx_x_min .AND. ibegin==indx_x_min ) THEN

                  jcolumn_global(2) = irow_global                      ! CENTER         
                  jcolumn_global(3) = irow_global + 1                  ! RIGHT               
                  jcolumn_global(4) = irow_global - (iend-ibegin+1) + 1! BOTTOM RIGHT    
                  
                  ! I will fill matrix values now
                  IF ( i_cylindrical==2 ) factor_axis_geom_cyl = 2.0_8 ! I apply it twice because the 1/8 factor in cartesian comes from 1/4*1/2
                  CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) - 0.25_8, eps_shifted_quarter)   !right
                  CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j)- 0.50_8, eps_shifted_quarter_2) !below                  

                  value_at_jcol(1) =   (eps_shifted_quarter_2*3.0_8/8.0_8 - eps_shifted_quarter*1.0_8/8.0_8*factor_axis_geom_cyl)*factor_axis_geom_cyl
                  value_at_jcol(3) = (- eps_shifted_quarter*1.0_8/8.0_8*factor_axis_geom_cyl + eps_shifted_quarter_2*3.0_8/8.0_8)*factor_axis_geom_cyl!-(eps_i_jshifted(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                  value_at_jcol(4) =   (eps_shifted_quarter*1.0_8/8.0_8*factor_axis_geom_cyl + eps_shifted_quarter_2*1.0_8/8.0_8)*factor_axis_geom_cyl
                  value_at_jcol(2) =   -(value_at_jcol(1)+value_at_jcol(3)+value_at_jcol(4))

                  call MatSetValues(Amat, one, irow_global, 4, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 

               ! Left of the block. I might need to communicate with my left neighbor if it exists. 
               ELSE IF ( i==indx_x_min+1 ) THEN

                     ! No neighbor on the left, I can use my own node     ! LEFT
                     IF (indx_x_min==ibegin) THEN 
                        jcolumn_global(2) = irow_global - 1
                        jcolumn_global(6) = irow_global - (iend-ibegin+1) - 1! BOTTOM LEFT
                     ELSE 
                        jcolumn_global(2) = process_left_bottom_right_inner_node + (j-indx_y_min-1) * process_left_solved_nodes_row_length                  ! LEFT
                        jcolumn_global(6) = process_left_bottom_right_inner_node + (j-indx_y_min-1-1) * process_left_solved_nodes_row_length! BOTTOM LEFT
                     END IF
                     jcolumn_global(3) = irow_global                      ! CENTER         
                     jcolumn_global(4) = irow_global + 1                  ! RIGHT
                     jcolumn_global(5) = irow_global - (iend-ibegin+1) + 1! BOTTOM RIGHT
                                       

                     IF ( i_cylindrical==2 ) THEN
                        factor_geom_cyl = 2.0_8 ! I need this because 1/8 in cartesian comes from 1/4*1/2
                        factor_geom_cyl_right = 1.0_8 + 1.0_8/(2.0_8*DBLE(i))
                        factor_geom_cyl_left = 1.0_8 - 1.0_8/(2.0_8*DBLE(i))
                     END IF
                     CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) - 0.25_8, eps_shifted_quarter)   !right
                     CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) - 0.25_8, eps_shifted_quarter_2) !left
                     ! print*,'eps_shifted_quarter,eps_shifted_quarter_2',eps_shifted_quarter,eps_shifted_quarter_2
      
                     value_at_jcol(1) =   (eps_i_jshifted(i,j) - eps_shifted_quarter*1.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl - eps_shifted_quarter_2*1.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl )
                     ! value_at_jcol(1) =   (eps_i_jshifted(i,j) + eps_shifted_quarter*3.0_8/8.0_8*factor_geom_cyl_right + eps_shifted_quarter_2*3.0_8/8.0_8*factor_geom_cyl_left )
                     value_at_jcol(2) =   eps_shifted_quarter_2*3.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl
                     value_at_jcol(4) =   eps_shifted_quarter*3.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl
                     value_at_jcol(5) =   eps_shifted_quarter*1.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl
                     value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl
                     value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))                  

                     call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr)                   

               ! Right of the block. I might need to communicate with my right neighbor if it exists. 
               ELSE IF ( i==indx_x_max-1 ) THEN

                     IF ( i_cylindrical==2 ) THEN
                        factor_geom_cyl = 2.0_8 ! I need this because1/8 in cartesian comes from 1/4*1/2
                        factor_geom_cyl_right = 1.0_8 + 1.0_8/(2.0_8*DBLE(i))
                        factor_geom_cyl_left = 1.0_8 - 1.0_8/(2.0_8*DBLE(i))
                     END IF               
                  
                     jcolumn_global(2) = irow_global - 1                     ! LEFT
                     jcolumn_global(3) = irow_global                         ! CENTER         
                     ! No neighbor on the right, I can use my own node
                     IF ( indx_x_max==iend ) THEN
                        jcolumn_global(4) = irow_global + 1                  ! RIGHT  
                        jcolumn_global(5) = irow_global - (iend-ibegin+1) + 1! BOTTOM RIGHT
                     ! Neighbor on the right, I need to communicate
                     ELSE
                        jcolumn_global(4) = process_right_bottom_left_inner_node + (j-indx_y_min-1) * process_right_solved_nodes_row_length                 ! RIGHT
                        jcolumn_global(5) = process_right_bottom_left_inner_node + (j-indx_y_min-1-1) * process_right_solved_nodes_row_length! BOTTOM RIGHT
                     END IF
                     jcolumn_global(6) = irow_global - (iend-ibegin+1) - 1! BOTTOM LEFT                            

                     CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) - 0.25_8, eps_shifted_quarter)   !right
                     CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) - 0.25_8, eps_shifted_quarter_2) !left
                     ! print*,'eps_shifted_quarter,eps_shifted_quarter_2',eps_shifted_quarter,eps_shifted_quarter_2
      
                     value_at_jcol(1) =   (eps_i_jshifted(i,j) - eps_shifted_quarter*1.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl - eps_shifted_quarter_2*1.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl )
                     value_at_jcol(2) =   eps_shifted_quarter_2*3.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl
                     value_at_jcol(4) =   eps_shifted_quarter*3.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl
                     value_at_jcol(5) =   eps_shifted_quarter*1.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl
                     value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl
                     value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))                  

                     call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr)                           

                  ! this is at the right of the domain, along the BC
               ELSE IF ( i==indx_x_max .AND. iend==indx_x_max ) THEN

                     jcolumn_global(2) = irow_global - 1               ! LEFT
                     jcolumn_global(3) = irow_global                   ! CENTER
                     jcolumn_global(4) = irow_global - (iend-ibegin+1) - 1! BOTTOM LEFT                     

                     IF ( i_cylindrical==2 ) THEN
                        factor_geom_cyl = 2.0_8 ! I need this because1/8 in cartesian comes from 1/4*1/2
                        factor_geom_cyl_left = 1.0_8 - 1.0_8/(4.0_8*DBLE(i)-1.0_8)
                     END IF

                     ! Filling matrix
                     CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) - 0.25_8, eps_shifted_quarter)   !left
                     CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) - 0.5_8, eps_shifted_quarter_2) !below       

                     value_at_jcol(1) =   (eps_shifted_quarter_2*3.0_8/8.0_8 - eps_shifted_quarter*1.0_8/8.0_8*factor_geom_cyl_left)*factor_geom_cyl
                     value_at_jcol(2) = (- eps_shifted_quarter_2*1.0_8/8.0_8 + eps_shifted_quarter*3.0_8/8.0_8*factor_geom_cyl_left)*factor_geom_cyl!-(eps_i_jshifted(i,j) + eps_ishifted_j(i+1,j)*factor_geom_cyl)
                     value_at_jcol(4) =   (eps_shifted_quarter*1.0_8/8.0_8*factor_geom_cyl_left + eps_shifted_quarter_2*1.0_8/8.0_8)*factor_geom_cyl
                     value_at_jcol(3) =   -(value_at_jcol(1)+value_at_jcol(2)+value_at_jcol(4))


                     call MatSetValues(Amat, one, irow_global, 4, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr)                       

               ! This is a regular node at the top BC. Far from block boundaries and wall materials 
               ELSE

                  ! Remove corners as it was done previously
                  ! IF ( i==indx_x_min .OR. i==indx_x_max ) CYCLE
                  IF ( i_cylindrical==2 ) THEN
                     factor_geom_cyl = 2.0_8 ! I need this because1/8 in cartesian comes from 1/4*1/2
                     factor_geom_cyl_right = 1.0_8 + 1.0_8/(2.0_8*DBLE(i))
                     factor_geom_cyl_left = 1.0_8 - 1.0_8/(2.0_8*DBLE(i))
                  END IF  

                  ! jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
                  jcolumn_global(2) = irow_global - 1                  ! LEFT
                  jcolumn_global(3) = irow_global                      ! CENTER         
                  jcolumn_global(4) = irow_global + 1                  ! RIGHT
                  jcolumn_global(5) = irow_global - (iend-ibegin+1) + 1! BOTTOM RIGHT
                  jcolumn_global(6) = irow_global - (iend-ibegin+1) - 1! BOTTOM LEFT
         
                  ! IF ( i_cylindrical==2 ) factor_geom_cyl = DBLE(i+1)/DBLE(i)
                  CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) - 0.25_8, eps_shifted_quarter)   !right
                  CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) - 0.25_8, eps_shifted_quarter_2) !left
                  ! print*,'eps_shifted_quarter,eps_shifted_quarter_2',eps_shifted_quarter,eps_shifted_quarter_2

                  value_at_jcol(1) =   (eps_i_jshifted(i,j) - eps_shifted_quarter*1.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl - eps_shifted_quarter_2*1.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl )
                  value_at_jcol(2) =   eps_shifted_quarter_2*3.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl
                  value_at_jcol(4) =   eps_shifted_quarter*3.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl
                  value_at_jcol(5) =   eps_shifted_quarter*1.0_8/8.0_8*factor_geom_cyl_right*factor_geom_cyl
                  value_at_jcol(6) =   eps_shifted_quarter_2*1.0_8/8.0_8*factor_geom_cyl_left*factor_geom_cyl
                  value_at_jcol(3) = -(value_at_jcol(1) + value_at_jcol(2) + value_at_jcol(4) + value_at_jcol(5) + value_at_jcol(6))
         
                  call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr) 
               END IF ! IF loop over i index
            END IF ! IF ( neumann_flag )
          END IF !IF (block_has_neumann_bc_Y_top) THEN        
             
       END DO !DO i = ibegin, iend
    END IF

! end of initialization of matrix coefficients written by DS ---------------

    DEALLOCATE(eps_ishifted_j, STAT=ALLOC_ERR)   ! eps_ishifted_j(i,j) is between nodes {i-1,j} and {i,j}
    DEALLOCATE(eps_i_jshifted, STAT=ALLOC_ERR)   ! eps_i_jshifted(i,j) is between nodes {i,j-1} and {i,j}

! call MatSetValues(Amat, one, ix, ny, jcol, v, INSERT_VALUES,ierr) 
! inserts or adds a block of values into a matrix
! MatAssemblyBegin() and MatAssemblyEnd() must be called after all calls to matSetValues have been completed.
! Amat - the matrix
! one - the number of rows
! ix - global indices of rows
! ny - the number of columns
! jcol - global indices of columns
! v - a logically two-dimensional array of values
! INSERT_VALUES replaces existing entries with new values

! Commit changes to the matrix [S.J.]

! Begins assembling the matrix
    call MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY, ierr)
! Completes assembling the matrix
    call MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY, ierr)
    
! Save matrix for testing ....
!    if(.true.) then
!       call  PetscViewerBinaryOpen(parComm, 'Amat.bin', FILE_MODE_WRITE, viewer, ierr)
!       call  MatView(Amat, viewer, ierr)
!       call  PetscViewerDestroy(viewer,ierr)        
!    end if
    
! Create solver [S.J.]

! KSP is abstract petsc object that manages all Krylov methods. 
! It manages the linear solves (even when krylov accelerators are not used like in direct solvers

! Creates the default KSP context
! here ksp is location where to put KSP context
    call KSPCreate(parComm, ksp, ierr)

    call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)


! Sets the matrix associated with the linear system and a (possibly) different one associated with the preconditioner
! argument #3 is the matrix to be used in constructing the preconditioner, 
! usually the same as argument #2 (the matrix that defines the linear system)
    call KSPSetOperators(ksp, Amat, Amat, ierr)

! . . Set method (this set-up uses the LU-decomposition as the solver) [S.J.]

! returns a pointer ot the preconditioner context set with KSPSetPc
    call KSPGetPc(ksp, pc, ierr)

! sets the preconditioner to be used to calculate the application of the preconditioner on a vector
    call KSPSetUp(ksp, ierr)

! Create B and X vectors [S.J.]

! note that ### bvec ### is the vector of the right hand side and ### xvec ### is the solution vector

! get vector(s) compatible with the matrix, i.e. with the same parallel layout
    call MatCreateVecs(Amat, PETSC_NULL_VEC, bvec, ierr)

! configures the vector from the options database
    call VecSetFromOptions(bvec, ierr)

! creates a new vector of the same type as an existing vector
    call VecDuplicate(bvec, xvec, ierr)

   CALL VecSet(xvec,zero,ierr)
    
    ! Done with setting up vectors and matrices!
    return
  end subroutine SolverInitialization
  
!-------------------------------------
!
! this subroutine inserts dvec which in the actual call is the real(8) array rhs
! into vector bvec
  subroutine InsertData(nnodes,dvec)

    PetscInt, intent(IN) :: nnodes
    PetscScalar, intent(IN) :: dvec(:)

    PetscInt, save :: n, val
    PetscInt, allocatable, save :: jvec(:)

    integer ibegin, jbegin, iend, jend

    integer :: i, j !, k
!    PetscViewer :: viewer
!    character(20) :: filename='bvector.dat'
    PetscErrorCode :: ierr
    
    if(.not.allocated(jvec)) then

       allocate(jvec(nnodes))

! calculation of jvec modified by DS

       jbegin = indx_y_min+1
       jend   = indx_y_max-1
       ibegin = indx_x_min+1
       iend   = indx_x_max-1

       IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
       IF (Rank_of_process_right.LT.0) iend   = indx_x_max
       IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
       IF (Rank_of_process_above.LT.0) jend   = indx_y_max

       n=0
       do j = jbegin, jend  !indx_y_min+1, indx_y_max-1
          do i = ibegin, iend  !indx_x_min+1, indx_x_max-1
!          if(i<1.or.i>=global_maximal_i.or.j<1.or.j>=global_maximal_j) cycle
             n = n+1
             jvec(n)=global_offset + n  !(N_grid_block_x*N_grid_block_y)*Rank_of_process+N_grid_block_x*(j-indx_y_min-1)+i-indx_x_min-1 
          end do
       end do
      !  print*,'jvec,global_offset_wil',jvec,global_offset
      !  print*,'dvec_wil',dvec
    end if

    val = n 
    if(.not.(nnodes.eq.n)) print *,'Not conformant!'
    if(size(dvec)<n) print *,'dvec too small!', size(dvec), n

! inserts or adds values into certain locations of a vector
! bvec - vector to insert in
! val - number of elements to add
! jvec - indices where to add
! dvec - array of values
! INSERT_VALUES - replaces existing entries with new values
    call VecSetValues(bvec, val, jvec, dvec, INSERT_VALUES, ierr)
    
! begins assembling the vector, should be called after completing all calls to VecSetValues
    call VecAssemblyBegin(bvec, ierr)

! completes assembling the vector
    call VecAssemblyEnd(bvec, ierr)

!    call PetscViewerCreate(parComm,viewer,ierr)
!    call PetscViewerASCIIOpen(parComm,filename,viewer,ierr)
!!    call PetscViewerSetType(viewer,PETSCVIEWERMATLAB,ierr)
!    call VecView(bvec,viewer,ierr)
!    call PetscViewerDestroy(parComm,viewer,ierr)

    return
  end subroutine InsertData

  subroutine FetchData(nnodes, dvec)

    PetscInt, intent(IN) :: nnodes
    PetscScalar, intent(OUT) :: dvec(:)

    integer ibegin, jbegin, iend, jend

    integer :: n, i, j   !, k, val
    PetscInt, allocatable, save :: jvec(:)
!    PetscViewer :: viewer
!    character(20) :: filename='soln.mat'
    PetscErrorCode :: ierr
    
    if(.not.allocated(jvec)) then
       allocate(jvec(nnodes))

! calculation of jvec modified by DS

       jbegin = indx_y_min+1
       jend   = indx_y_max-1
       ibegin = indx_x_min+1
       iend   = indx_x_max-1

       IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
       IF (Rank_of_process_right.LT.0) iend   = indx_x_max
       IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
       IF (Rank_of_process_above.LT.0) jend   = indx_y_max

       n=0
       do j = jbegin, jend  !indx_y_min+1, indx_y_max-1
          do i = ibegin, iend  !indx_x_min+1, indx_x_max-1
!            if(i<1.or.i>=global_maximal_i.or.j<1.or.j>=global_maximal_j) cycle
             n = n+1
             jvec(n)=global_offset + n  !(N_grid_block_x*N_grid_block_y)*Rank_of_process+N_grid_block_x*(j-indx_y_min-1)+i-indx_x_min-1
          end do
       end do
    end if

! gets values from certain locations of a vector
    call VecGetValues(xvec, nnodes, jvec, dvec, ierr)
    
! begins assembling the vector
    call VecAssemblyBegin(xvec, ierr)
! completes assembling the vector
    call VecAssemblyEnd(xvec, ierr)

!  Get the local part of the solution
!    call PetscViewerCreate(parComm,viewer,ierr)
!    call PetscViewerASCIIOpen(parComm,filename,viewer,ierr)
!!    call PetscViewerSetType(viewer,PETSCVIEWERMATLAB,ierr)
!    call VecView(bvec,viewer,ierr)
!    call PetscViewerDestroy(parComm,viewer,ierr)
    return
  end subroutine FetchData
  
  subroutine SolverDestroy
    PetscErrorCode :: ierr
    call KSPDestroy(ksp,ierr)
    call MatDestroy(Amat,ierr)
    return
  end subroutine SolverDestroy
  
end module PETSc_Solver

!-------------------------------
! 7---6---5
! |       |
! 8   0   4
! |       |
! 1---2---3
!
SUBROUTINE FIND_INNER_OBJECT_CONTAINING_POINT(i, j, nio, position_flag)

  USE CurrentProblemValues

  IMPLICIT NONE

  INTEGER, INTENT(IN) ::  i, j  ! x,y indices of the point of interest
  INTEGER, INTENT(OUT) :: nio   ! number of inner object containing the point
                                ! it is meaningful only if position_flag>=0
  INTEGER, INTENT(OUT) :: position_flag   ! 9 for a metal object
                                          ! 0 inside a dielectric object
                                          ! 1 to 8 at the surface of the dielectric object
                                          ! -1 if there is no object containing the given point
  INTEGER n

  nio = -1
  position_flag = -1

  IF (N_of_inner_objects.EQ.0) RETURN

  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
       IF (i.LT.whole_object(n)%ileft) CYCLE
       IF (i.GT.whole_object(n)%iright) CYCLE
       IF (j.LT.whole_object(n)%jbottom) CYCLE
       IF (j.GT.whole_object(n)%jtop) CYCLE

! since we are here, the point is either at the surface or inside inner material object n

       nio = n  ! save the number of the inner object

       IF (whole_object(n)%object_type.EQ.METAL_WALL) THEN
! the material object is metal, stop scanning inner objects, dielectric properties are omitted (relevant if metal is inside dielectric)
          position_flag = 9
          RETURN
       END IF

       position_flag = 0
       IF (i.EQ.whole_object(n)%ileft) THEN
! point on the left side of a dielectric object
          IF (j.EQ.whole_object(n)%jbottom) THEN
! left bottom corner
             position_flag = 1
          ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! left top corner
             position_flag = 7
          ELSE
! left side, not a corner
             position_flag = 8
          END IF
       ELSE IF (i.EQ.whole_object(n)%iright) THEN
! point on the right side of a dielectric object
          IF (j.EQ.whole_object(n)%jbottom) THEN
! right bottom corner
             position_flag = 3
          ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! right top corner
             position_flag = 5
          ELSE
! right side, not a corner
             position_flag = 4
          END IF
       ELSE IF (j.EQ.whole_object(n)%jbottom) THEN
! point on the bottom side of a dielectric object
          position_flag = 2
       ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! point on the top side of a dielectric object
          position_flag = 6
       END IF
! note that we do not stop scanning inner material objects here since there may be metal objects inside the dielectric
    END DO

    RETURN

END SUBROUTINE FIND_INNER_OBJECT_CONTAINING_POINT

!-------------------------------
! 7---6---5
! |       |
! 8   0   4
! |       |
! 1---2---3
!
SUBROUTINE CHECK_IF_INNER_OBJECT_CONTAINS_POINT(myobject, i, j, position_flag)

  USE CurrentProblemValues

  IMPLICIT NONE

  TYPE(boundary_object), INTENT(IN) :: myobject
  INTEGER, INTENT(IN) :: i, j  ! x,y indices of the point of interest
  INTEGER, INTENT(OUT) :: position_flag   ! 9 for a metal object
                                          ! 0 inside a dielectric object
                                          ! 1 to 8 at the surface of the dielectric object
                                          ! -1 if there is no object containing the given point
  INTEGER n

  position_flag = -1

!  IF (N_of_inner_objects.EQ.0) RETURN

  IF (i.LT.myobject%ileft) RETURN
  IF (i.GT.myobject%iright) RETURN
  IF (j.LT.myobject%jbottom) RETURN
  IF (j.GT.myobject%jtop) RETURN

! since we are here, the point is either at the surface or inside inner material object n

  IF (myobject%object_type.EQ.METAL_WALL) THEN
! the material object is metal, stop scanning inner objects, dielectric properties are omitted (relevant if metal is inside dielectric)
     position_flag = 9
     RETURN
  END IF

  position_flag = 0
  IF (i.EQ.myobject%ileft) THEN
! point on the left side of a dielectric object
     IF (j.EQ.myobject%jbottom) THEN
! left bottom corner
        position_flag = 1
     ELSE IF (j.EQ.myobject%jtop) THEN
! left top corner
        position_flag = 7
     ELSE
! left side, not a corner
        position_flag = 8
     END IF
  ELSE IF (i.EQ.myobject%iright) THEN
! point on the right side of a dielectric object
     IF (j.EQ.myobject%jbottom) THEN
! right bottom corner
        position_flag = 3
     ELSE IF (j.EQ.myobject%jtop) THEN
! right top corner
        position_flag = 5
     ELSE
! right side, not a corner
        position_flag = 4
     END IF
  ELSE IF (j.EQ.myobject%jbottom) THEN
! point on the bottom side of a dielectric object
     position_flag = 2
  ELSE IF (j.EQ.myobject%jtop) THEN
! point on the top side of a dielectric object
     position_flag = 6
  END IF
! note that we do not stop scanning inner material objects here since there may be metal objects inside the dielectric

  RETURN

END SUBROUTINE CHECK_IF_INNER_OBJECT_CONTAINS_POINT

!--------------------------------------------------
!
REAL(8) FUNCTION Get_Surface_Charge_Inner_Object(i,j,position_flag, myobject)

  USE CurrentProblemValues
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  INTEGER i, j            ! x,y indices of the point of interest
  INTEGER position_flag   ! 9 for a metal object
                          ! 0 inside a dielectric object
                          ! 1 to 8 at the surface of the dielectric object
                          ! -1 if there is no object containing the given point
  TYPE(boundary_object) myobject

  INTEGER pos

  INTEGER i_left_top, i_right_top, i_right_bottom

  Get_Surface_Charge_Inner_Object = 0.0_8 

  IF (N_of_inner_objects.EQ.0) RETURN
!  IF (nio.LT.N_of_boundary_objects+1) RETURN
!  IF (nio.GT.N_of_boundary_and_inner_objects) RETURN
  IF (position_flag.LT.1) RETURN
  IF (position_flag.GT.8) RETURN

  pos=-1

! ilt   -  ---- irt
!  |             |
!  |             |
!  1 ilbb ---- irb
!
  i_left_top        =                  myobject%jtop   - myobject%jbottom + 1
  i_right_top       = i_left_top     + myobject%iright - myobject%ileft
  i_right_bottom    = i_right_top    + myobject%jtop   - myobject%jbottom
!  i_left_bottom_bis = i_right_bottom + myobject%iright - myobject%ileft - 1

  SELECT CASE (position_flag)
     CASE (1)
! left bottom corner
        pos=1
     CASE (8)
! left side
        pos = j - myobject%jbottom  + 1
     CASE (7)
! left top corner
        pos = i_left_top
     CASE (6)
! top side
        pos = i_left_top + i - myobject%ileft
     CASE (5) 
! right top corner
        pos = i_right_top
     CASE (4)
! right side
        pos = i_right_top + myobject%jtop - j
     CASE (3)
! right bottom corner
        pos = i_right_bottom
     CASE (2)
  ! bottom side
        pos = i_right_bottom + myobject%iright - i
  END SELECT

  IF ((pos.GE.1).AND.(pos.LE.myobject%N_boundary_nodes)) THEN
     Get_Surface_Charge_Inner_Object = myobject%surface_charge(pos)
     RETURN
  ELSE
     print '("Error in Get_Surface_Charge_Inner_Object")'
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END FUNCTION Get_Surface_Charge_Inner_Object

!--------------------------------------------------------------------------------------------------
!     DECIDE_NEUMANN_EXTERNAL_BOUNDARY
!>    @details Determine if point i,j at external boundary belongs to a Neumann BC or Metal. If it belongs to both, then metal will prevail 
!!    @authors W. Villafana
!!    @date    Jun-19-2023
!-------------------------------------------------------------------------------------------------- 
SUBROUTINE DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i, j, neumann_flag)  

   USE CurrentProblemValues, ONLY: string_length, whole_object, N_of_boundary_objects, METAL_WALL
   USE mod_print, ONLY: print_error
 
   IMPLICIT NONE
 
   INCLUDE 'mpif.h'
  
   !IN/OUT
   INTEGER, INTENT(IN) :: i, j  ! Current point
   LOGICAL, INTENT(OUT) :: neumann_flag ! 0=No Neumann, 1=Neumann
 
   !LOCAL 
   INTEGER                      :: n1 ! Boundary number 
   INTEGER                      :: last_metal_object ! last metal object 
   INTEGER                      :: i_found ! Determines if point belongs to an external boundary
   CHARACTER(LEN=string_length) :: message, routine

   routine = "DECIDE_NEUMANN_EXTERNAL_BOUNDARY"

   last_metal_object = 0 ! By default it is zero, ie, it does not belong to any metal object
   i_found = 0 ! By default point does not belong anywhere
   neumann_flag = .FALSE.

   ! Loop over external boundaries 
   DO n1 = 1,N_of_boundary_objects
      
      ! Determine if point i,j belongs to a boundary 
      IF ( i<whole_object(n1)%ileft   ) CYCLE
      IF ( i>whole_object(n1)%iright  ) CYCLE
      IF ( j<whole_object(n1)%jbottom ) CYCLE
      IF ( j>whole_object(n1)%jtop    ) CYCLE

      ! Point belongs to one boundary at least
      i_found = 1
      IF ( whole_object(n1)%object_type==METAL_WALL ) last_metal_object = n1
      
   ENDDO

   ! Check if point was found at any boundary (must be the case)
   ! IF ( i_found==0 ) THEN
   !    WRITE( message ,'(A,I4,A,I4,A)') "Grid point i=",i," and j=",j," does not belong to any external boundary. Double check configuration."
   !    CALL print_error ( message,routine )
   ! END IF

   ! Determine if point is metal or not. Last metal will prevail
   ! IF ( i_found==1 .AND. last_metal_object==0 ) neumann_flag = .TRUE.
   IF ( last_metal_object==0 ) neumann_flag = .TRUE.
 
 END SUBROUTINE DECIDE_NEUMANN_EXTERNAL_BOUNDARY

!-----------------------------------------------
!
! Precomputes all element of Poisson matrix. Returns weight for point on the left in the stencil
SUBROUTINE SET_EPS_ISHIFTED(i, j, eps, petsc_use)  ! here point {i,j} is between nodes {i-1,j} and {i,j}

  USE CurrentProblemValues

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  INTEGER, INTENT(IN) :: i, j
  REAL(8), INTENT(OUT) :: eps
  INTEGER, INTENT(IN), OPTIONAL :: petsc_use ! if we build the Petsc matrix, we put gemoetrical effect right now. For the external circuit we do it later

  INTEGER count  ! counts dielectric objects
  LOGICAL segment_is_inside_dielectric_object
  LOGICAL segment_is_on_surface_of_metal_object
  INTEGER n1
  REAL(8) :: factor_cyl ! Additional geometrical factor in matrix for cylindrical
  LOGICAL :: add_geom

  !!! Add a factor if we have cylindrical coordinates r-z. We compute the factor f_left=1-delta_r/(2rij) and f_right = 2-f_left for rij>0
  factor_cyl = one ! By default I do not have  a factor 
  add_geom = .FALSE.

   ! Add geometricla factor for Petsc use. Not pretty but will work. Use of petsc_use directly does not work. Seg fault I do not understand.
   IF (PRESENT ( petsc_use )) add_geom = .TRUE.

   ! Add geometrical coefficients only when building the Petsc matrix and only for cylindrical
   IF ( i_cylindrical==2 .AND. i>0 .AND. add_geom ) factor_cyl = one - one/(two*DBLE(i)) ! for i= 0 we are at the axis and this factor disepears in the scheme. All r-z simulations are assumed to include the axis

!   print*,'current i,j',i,j
! find all inner objects owning segment {i-1,j}-{i,j}
! assume that only two dielectric objects may own a common segment
! allow a metal object to be added on/axis top of that

  eps = 0.0_8
  count = 0
  segment_is_inside_dielectric_object = .FALSE.

! first, look through the dielectric objects only
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.EQ.METAL_WALL) CYCLE

! find an inner object containing point {i-1,j}

     IF (i-1.LT.whole_object(n1)%ileft) CYCLE
     IF (i-1.GT.whole_object(n1)%iright) CYCLE
     IF (j.LT.whole_object(n1)%jbottom) CYCLE
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (i.GT.whole_object(n1)%iright) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i-1,j}-{i,j} is owned by the object
     
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object

     IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is on the surface

        count = count + 1
        eps = eps + whole_object(n1)%eps_diel
      !   print*,'increment_eps,eps_object,i,j',eps,whole_object(n1)%eps_diel,i,j

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside
        segment_is_inside_dielectric_object = .TRUE.
        eps = whole_object(n1)%eps_diel*factor_cyl
      !   print*,'segm_inside,esp,i,j,factr',eps,i,j,factor_cyl
        count = 1
        EXIT

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! at this point
! if count==0 then the segment is either in/at metal object or in vacuum
! if count==1 then the segment either inside a dielectric object or on the surface of ONE (not two) dielectric object
! if count==2 then the segment is at the interface between two dielectrics
! if count==3 this is an error

! second, look through the metal objects only
! if segment belongs to a metal object 
! return zero if it is inside
! if it is on the surface
! return (#1) epsilon of adjacent dielectric object or (#2) 1 if the metal object is in vacuum or (#3) 0 if the metal object is attached to another metal object

  segment_is_on_surface_of_metal_object = .FALSE.

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.NE.METAL_WALL) CYCLE

! find an inner object containing point {i-1,j}

     IF (i-1.LT.whole_object(n1)%ileft) CYCLE
     IF (i-1.GT.whole_object(n1)%iright) CYCLE
     IF (j.LT.whole_object(n1)%jbottom) CYCLE
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (i.GT.whole_object(n1)%iright) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i-1,j}-{i,j} is owned by the object
 
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object
     
     IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is on the surface

        IF (segment_is_on_surface_of_metal_object) THEN
! #3
! another metal object already has this segment on its border, so the segment is between two metal objects
! return zero epsilon
           eps = 0.0_8
           RETURN
        END IF
! remember that a metal object was found
        segment_is_on_surface_of_metal_object = .TRUE.

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside

        eps = 0.0_8
        RETURN

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! note that cases when segment is inside a metal object or in between two attached metal objects are already processed above, followed by RETURN calls

  IF (segment_is_on_surface_of_metal_object) THEN
     IF (count.EQ.0) THEN
! segment on metal surface facing vacuum
        eps = 1.0_8*factor_cyl
     ELSE
! segment on metal surface facing dielectric
        eps = eps / count*factor_cyl
     END IF
     RETURN
  END IF

  IF (count.EQ.0) THEN
! segment is in vacuum
     eps = 1.0_8*factor_cyl
     RETURN
  ELSE IF (count.EQ.1) THEN
! segment is inside a dielectric object
     IF (segment_is_inside_dielectric_object) RETURN
! segment is on the surface of a dielectric object, facing vacuum
     eps = (1.0_8 + eps) / 2.0_8*factor_cyl
     RETURN
  ELSE IF (count.EQ.2) THEN
     eps = eps / DBLE(count)*factor_cyl
     RETURN
  ELSE
     PRINT '("Error-3 in SET_EPS_ISHIFTED for i/j ",2x,i4,2x,i4)', i, j
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE SET_EPS_ISHIFTED


!-----------------------------------------------
!
SUBROUTINE SET_EPS_JSHIFTED(i, j, eps, petsc_use)  ! here point {i,j} is between nodes {i,j-1} and {i,j}

  USE CurrentProblemValues
  USE BlockAndItsBoundaries, ONLY: block_has_symmetry_plane_X_left,indx_x_min
  USE ParallelOperationValues, ONLY: Rank_of_process

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  INTEGER, INTENT(IN) :: i, j
  REAL(8), INTENT(OUT) :: eps
  INTEGER, OPTIONAL, INTENT(IN) :: petsc_use ! if we build the Petsc matrix, we put gemoetrical effect right now. For the external circuit we do it later

  INTEGER count  ! counts dielectric objects
  LOGICAL segment_is_inside_dielectric_object
  LOGICAL segment_is_on_surface_of_metal_object
  INTEGER n1
  REAL(8) :: factor_cyl,factor_cyl_used ! Additional geometrical factor in matrix for cylindrical
  LOGICAL :: add_geom

  !!! Add a factor if we have cylindrical coordinates r-z. We compute the factor f_left=1-delta_r/(4rij) and f_right = 2-f_left for rij>0
  factor_cyl = one ! By default I do not have  a factor 
  factor_cyl_used = one ! By default I do not have  a factor 
  add_geom = .FALSE.

  ! Add geometricla factor for Petsc use. Not pretty but will work. Use of petsc_use directly does not work. Seg fault I do not understand.
  IF (PRESENT ( petsc_use )) add_geom = .TRUE.
  
   ! Add geometrical coefficients only when building the Petsc matrix and only for cylindrical
   IF ( i_cylindrical==2 .AND. i>0 .AND. add_geom ) factor_cyl = one - one/(four*DBLE(i)) ! for i= 0 we are at the axis and this factor disepears in the scheme. All r-z simulations are assumed to include the axis


! find all inner objects owning segment {i,j-1}-{i,j}
! assume that only two dielectric objects may own a common segment
! allow a metal object to be added on top of that

  eps = 0.0_8
  count = 0
  segment_is_inside_dielectric_object = .FALSE.

! first, look through the dielectric objects only
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.EQ.METAL_WALL) CYCLE

! find an inner object containing point {i,j-1}

     IF (i.LT.whole_object(n1)%ileft) CYCLE
     IF (i.GT.whole_object(n1)%iright) CYCLE
     IF (j-1.LT.whole_object(n1)%jbottom) CYCLE
     IF (j-1.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j-1} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i,j-1}-{i,j} is owned by the object
     
! find whether segment {i,j-1}-{i,j} is at the surface or inside
! it only can be on the surface if this is the left or the right edge of the object

     IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is on the surface

        count = count + 1
        
        ! For generic case (inside dielectric or Cartesian case)
        IF ( i_cylindrical==0 ) THEN
            eps = eps + whole_object(n1)%eps_diel
        ! Cylindrical case r-z
        ELSE IF ( i_cylindrical==2) THEN
            ! If dielectric is on the left side of my segment
            IF ( i==whole_object(n1)%iright ) THEN
               eps = eps + whole_object(n1)%eps_diel*factor_cyl
               factor_cyl_used = factor_cyl
            ! Segment is on the right side of my segment
            ELSE
               eps = eps + whole_object(n1)%eps_diel*(two-factor_cyl)
               factor_cyl_used = two - factor_cyl
            END IF
         ! Cylindrical z-theta or others (error for now)
         ELSE
            PRINT*,"Error, in SET_EPS_JSHIFTED. r-theta (i_cylindrical=1) or other geometries are not implemented"
            STOP
         END IF
         ! note that for cylindrical r-z with dielectric-vacuum vertical interface with epsilon=1, there is no change with respect to the Cartesian case

     ELSE   !### IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is inside
        segment_is_inside_dielectric_object = .TRUE.
        eps = whole_object(n1)%eps_diel
        count = 1
        EXIT

     END IF   !### ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! at this point
! if count==0 then the segment is either in/on metal object or in vacuum
! if count==1 then the segment is either inside a dielectric object or on the surface of ONE (not two) dielectric object
! if count==2 then the segment is at the interface between two dielectrics
! if count>=3 this is an error

! second, look through the metal objects only
! if segment belongs to a metal object 
! return zero if it is inside
! if it is on the surface
! return (#1) epsilon of adjacent dielectric object or (#2) 1 if the metal object is in vacuum or (#3) 0 if the metal object is attached to another metal object

  segment_is_on_surface_of_metal_object = .FALSE.

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.NE.METAL_WALL) CYCLE

! find an inner object containing point {i,j-1}

     IF (i.LT.whole_object(n1)%ileft) CYCLE
     IF (i.GT.whole_object(n1)%iright) CYCLE
     IF (j-1.LT.whole_object(n1)%jbottom) CYCLE
     IF (j-1.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i,j-1}-{i,j} is owned by the object
 
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object
     
     IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is on the surface

        IF (segment_is_on_surface_of_metal_object) THEN
! #3
! another metal object already has this segment on its border, so the segment is between two metal objects
! return zero epsilon
           eps = 0.0_8
           RETURN
        END IF
! remember that a metal object was found
        segment_is_on_surface_of_metal_object = .TRUE.

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside

        eps = 0.0_8
        RETURN

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! note that cases when segment is inside a metal object or in between two attached metal objects are already processed above, followed by RETURN calls

  IF (segment_is_on_surface_of_metal_object) THEN
     IF (count.EQ.0) THEN
! segment on metal surface facing vacuum
        eps = 1.0_8
     ELSE
! segment on metal surface facing dielectric
        eps = eps / count
     END IF
     RETURN
  END IF

  IF (count.EQ.0) THEN
! segment is in vacuum
     eps = 1.0_8
     RETURN
  ELSE IF (count.EQ.1) THEN
! segment is inside a dielectric object
     IF (segment_is_inside_dielectric_object) RETURN

     ! Segment is at surface but it is on symmetry axis: I can stop the computation
     IF ( block_has_symmetry_plane_X_left .AND. i==indx_x_min ) RETURN
     
     ! segment is on the surface of a dielectric object, facing vacuum
     ! Cartesian case 
     IF ( i_cylindrical==0 ) THEN
         eps = (1.0_8 + eps) / 2.0_8
     ELSE IF ( i_cylindrical==2 ) THEN
         eps = (one*factor_cyl_used + eps) / two ! BNote that if eps=1 then this is the same as Cartesian
     END IF
     RETURN
  ELSE IF (count.EQ.2) THEN
     eps = eps / DBLE(count)
     RETURN
  ELSE
     PRINT '("Error-3 in SET_EPS_JSHIFTED for i/j ",2x,i4,2x,i4)', i, j
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE SET_EPS_JSHIFTED

!-----------------------------------------------
! this is a simplified verions
! it is expected that point x is INSIDE some inner object or OUTSIDE, but not at the border of the object
!
SUBROUTINE GET_EPS_IN_POINT(x, y, eps)

  USE CurrentProblemValues

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  REAL(8), INTENT(IN) :: x, y
  REAL(8), INTENT(OUT) :: eps

  INTEGER count  ! counts dielectric objects
  LOGICAL point_is_inside_dielectric_object
  LOGICAL point_is_on_surface_of_metal_object
  INTEGER n1

! find all inner objects owning segment {i,j-1}-{i,j}
! assume that only two dielectric objects may own a common segment
! allow a metal object to be added on top of that

  eps = 0.0_8
  count = 0
  point_is_inside_dielectric_object = .FALSE.

! first, look through the dielectric objects only
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.EQ.METAL_WALL) CYCLE

! find an inner object containing point x, y

     IF (x.LT.whole_object(n1)%Xmin) CYCLE
     IF (x.GT.whole_object(n1)%Xmax) CYCLE
     IF (y.LT.whole_object(n1)%Ymin) CYCLE
     IF (y.GT.whole_object(n1)%Ymax) CYCLE

! since we are here, point x,y is either at the surface or inside inner object n1

     count = count + 1
     eps = eps + whole_object(n1)%eps_diel

     IF ( (x.GT.whole_object(n1)%Xmin).AND. &
        & (x.LT.whole_object(n1)%Xmax).AND. &
        & (y.GT.whole_object(n1)%Ymin).AND. &
        & (y.LT.whole_object(n1)%Ymax) ) THEN
! point inside
        point_is_inside_dielectric_object = .TRUE.
        EXIT
     END IF

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! at this point
! if count==0 then the point is either in/on metal object or in vacuum
! if count==1 then the point is either inside a dielectric object or on the surface of ONE (not two) dielectric object
! if count==2,3,4 then the point is at the interface between 2,3,4 dielectrics
! if count>=5 this is an error

! second, look through the metal objects only
! if point belongs to a metal object 
! return zero if it is inside
! if it is on the surface
! return (#1) epsilon of adjacent dielectric object or (#2) 1 if the metal object is in vacuum or (#3) 0 if the metal object is attached to another metal object

  point_is_on_surface_of_metal_object = .FALSE.

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.NE.METAL_WALL) CYCLE

! find an inner object containing point x, y

     IF (x.LT.whole_object(n1)%Xmin) CYCLE
     IF (x.GT.whole_object(n1)%Xmax) CYCLE
     IF (y.LT.whole_object(n1)%Ymin) CYCLE
     IF (y.GT.whole_object(n1)%Ymax) CYCLE

! since we are here, point x,y is either at the surface or inside inner object n1

     IF ( (x.EQ.whole_object(n1)%Xmin).OR. &
        & (x.EQ.whole_object(n1)%Xmax).OR. &
        & (y.EQ.whole_object(n1)%Ymin).OR. &
        & (y.EQ.whole_object(n1)%Ymax) ) THEN
! point is on the surface
        IF (point_is_on_surface_of_metal_object) THEN
! #3
! another metal object already has this segment on its border, so the segment is between two metal objects
! return zero epsilon
           eps = 0.0_8
           RETURN
        END IF
! remember that a metal object was found
        point_is_on_surface_of_metal_object = .TRUE.
     ELSE
! point is inside
        eps = 0.0_8
        RETURN
     END IF

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! note that cases when segment is inside a metal object or in between two attached metal objects are already processed above, followed by RETURN calls

  IF (point_is_on_surface_of_metal_object) THEN
     IF (count.EQ.0) THEN
! point on metal surface facing vacuum
        eps = 1.0_8
     ELSE
! point on metal surface facing dielectric
        eps = eps / count
     END IF
     RETURN
  END IF

  IF (count.EQ.0) THEN
! point is in vacuum
     eps = 1.0_8
     RETURN
  ELSE IF (count.EQ.1) THEN
! point is inside a dielectric object
     IF (point_is_inside_dielectric_object) RETURN
! segment is on the surface of a dielectric object, facing vacuum
     eps = (1.0_8 + eps) / 2.0_8
     RETURN
  ELSE IF (count.LT.4) THEN
     eps = eps / DBLE(count)
     RETURN
  ELSE
     PRINT '("Error-3 in GET_EPS_IN_POINT for x/y ",2x,f11.3,2x,f11.3)', x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE GET_EPS_IN_POINT

