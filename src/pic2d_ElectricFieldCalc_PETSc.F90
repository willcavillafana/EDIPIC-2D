
!--------------------------------------
!
! procedures InsertData, Solve, and FetchData used below were written by Janhunen Jani Salomon
!
!
SUBROUTINE SOLVE_POTENTIAL_WITH_PETSC

  USE PETSc_Solver
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries
  USE Diagnostics, ONLY : Save_probes_e_data_T_cntr, N_of_probes_block, Probe_params_block_list, probe_F_block

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  REAL(8), ALLOCATABLE :: queue(:)     ! array of phi-values aligned left-to-right--bottom-to-top
  REAL(8), ALLOCATABLE :: rhsvalue(:)  ! array of right-hand-side-values corresponding to the queue
  INTEGER ALLOC_ERR

  REAL(8) factor_rho

  INTEGER jbegin, jend, ibegin, iend
  INTEGER irow_global
  INTEGER nn

  INTEGER i, j, n, m, i_start, i_end, j_start, j_end, nobj

  INTEGER n1, n3
  REAL(8), ALLOCATABLE :: rbufer(:)

  INTEGER npb
  REAL(8) :: k_test_r, k_test_z,Lr,Lz,k_x_cart,k_y_cart, k_x_cart_neumann_right, k_y_cart_neumann_right, k_test_z_neumann_right, k_test_r_neumann_right, k_test_r_dirichlet, k_test_z_dirichlet ! Test for cylindrical r-z and neumann
  INTEGER :: index_maxi_r,index_maxi_z ! test r-z
  INTEGER :: ni0,position_flag ! test if I ma in an inner object
  REAL(8) :: alpha_r ! root of bessel function
  LOGICAL :: neumann_flag ! determines if point i,j is Neumann

  neumann_flag = .FALSE.

  ! I need collective operation because we do not know the exact size of the box... This is because of ghost points ...
!   CALL MPI_ALLREDUCE(indx_x_max, index_maxi_r,1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, stattus, ierr)
!   CALL MPI_ALLREDUCE(indx_y_max, index_maxi_z,1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, stattus, ierr)
  
!   Lr = (index_maxi_r)*delta_x_m !0.03009999915957451!0.030028263106942177!0.003025328740477562
!   Lz = (index_maxi_z)*delta_x_m!0.020099999383091927!0.020027175545692444!0.002025220077484846
!   print*,'Lr,Lz,index_maxi_r,index_maxi_z',Lr,Lz,index_maxi_r,index_maxi_z
!   k_test_r = 11./4*two*pi/Lr ! Basic test 
!   k_test_z = 5.0*two*pi/Lz ! Basic test
!   alpha_r = 18.071063967910924/Lr ! sixth root
   ! k_x_cart = (pi+5.0_8*pi)/Lr ! Test in Cartesian coordinates
   ! k_x_cart_neumann_right = (pi/2.0_8+7.0_8*pi)/Lr ! Test in Cartesian coordinates. Neumann right only
   ! k_y_cart_neumann_right = (pi+5.0_8*pi)/Lz ! Test in Cartesian coordinates. Neumann right only
   ! k_y_cart = (pi/2.0_8+5.0_8*pi)/Lz ! Test in Cartesian coordinates
   ! k_test_r = (pi/2.0_8+6.0_8*pi)/Lr ! Test in Cylindrical coordinates
   ! k_test_z = (pi/2.0_8+5.0_8*pi)/Lz ! Test in Cylindrical coordinates  
   ! k_test_r_neumann_right = (pi+6.0_8*pi)/Lr ! Test in Cylindrical coordinates  with neumann right and dirichlet otherzise  
   ! k_test_z_neumann_right = (pi+5.0_8*pi)/Lz ! Test in Cylindrical coordinates  with neumann right and dirichlet otherzise  
   ! k_test_r_dirichlet = (pi/2.0_8+6.0_8*pi)/Lr ! Test in Cylindrical coordinates  with dirichlet otherzise  
   ! k_test_z_dirichlet = (pi+5.0_8*pi)/Lz ! Test in Cylindrical coordinates  with dirichlet otherzise     

  ALLOCATE(   queue(1:block_N_of_nodes_to_solve), STAT = ALLOC_ERR)
  ALLOCATE(rhsvalue(1:block_N_of_nodes_to_solve), STAT = ALLOC_ERR)
 
  n=0
  rhsvalue = 0.0_8
  queue = 0.0_8
  phi = 0.0_8
 
  ! If poisson is not required, then let everything to be zero
  IF ( i_no_poisson==1 ) RETURN

  factor_rho = -one / N_of_particles_cell_dble

  jbegin = indx_y_min+1
  jend   = indx_y_max-1
  ibegin = indx_x_min+1
  iend   = indx_x_max-1

  IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
  IF (Rank_of_process_right.LT.0) iend   = indx_x_max
  IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
  IF (Rank_of_process_above.LT.0) jend   = indx_y_max

  irow_global = global_offset

  nn=0

!    j = indx_y_min !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
  IF (jbegin.EQ.indx_y_min) THEN
! boundary object along bottom border
     DO i = ibegin, iend
        nn = nn + 1
         IF ( .NOT.block_has_neumann_bc_Y_bottom ) THEN
            DO n = 1, N_of_local_object_parts_below
               m = index_of_local_object_part_below(n)
               i_start = local_object_part(m)%istart
               i_end   = local_object_part(m)%iend
               nobj   = local_object_part(m)%object_number
               IF ((i.GE.i_start).AND.(i.LE.i_end)) THEN
                  IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                     rhsvalue(nn) = whole_object(nobj)%phi
                  ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                     rhsvalue(nn) = whole_object(nobj)%phi_profile(i)
                  END IF
                  EXIT
               END IF
            END DO
         ELSE ! Neumann
            ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,indx_y_min,neumann_flag)

            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN                

               ! rhsvalue(nn) = -(k_x_cart**2+k_y_cart**2)*SIN(k_x_cart*DBLE(i)*delta_x_m)*(one)/ F_scale_V*delta_x_m**2/2.0 ! Cartesian. 
               ! IF (i==0) THEN
               !    rhsvalue(nn) = (-(k_test_r**2+k_test_z**2)*COS(k_test_r*delta_x_m*i)-k_test_r**2)/ F_scale_V*delta_x_m**2/two ! Cylindrical
               ! ELSE
               !    rhsvalue(nn) = (-(k_test_r**2+k_test_z**2)*COS(k_test_r*delta_x_m*i)-k_test_r/(DBLE(i)*delta_x_m)*SIN(k_test_r*delta_x_m*i))/ F_scale_V*delta_x_m**2/two ! Cylindrical
               ! END IF
               rhsvalue(nn) = factor_rho * (rho_i(i,indx_y_min) - rho_e(i,indx_y_min))/two ! Nodal volue at top is dx**2/2 in cartesian. In cylindrical, I keep the same convention, ie the RHS has a factoer dx**2/2, everythong else is put on the LHS              
            ELSE
               DO n = 1, N_of_local_object_parts_below
                  m = index_of_local_object_part_below(n)
                  i_start = local_object_part(m)%istart
                  i_end   = local_object_part(m)%iend
                  nobj   = local_object_part(m)%object_number
                  IF ((i.GE.i_start).AND.(i.LE.i_end)) THEN
                     IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                        rhsvalue(nn) = whole_object(nobj)%phi
                     ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                        rhsvalue(nn) = whole_object(nobj)%phi_profile(i)
                     END IF
                     EXIT
                  END IF
               END DO               
            ENDIF
         END IF
      END DO
  END IF

  DO j = indx_y_min+1, indx_y_max-1 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!     i = indx_x_min

     IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
        nn = nn + 1
         IF (.NOT.block_has_symmetry_plane_X_left) THEN
            IF (.NOT.block_has_neumann_bc_X_left) THEN
               DO n = 1, N_of_local_object_parts_left
                  m = index_of_local_object_part_left(n)
                  j_start = local_object_part(m)%jstart
                  j_end   = local_object_part(m)%jend
                  nobj   = local_object_part(m)%object_number
                  IF ((j.GE.j_start).AND.(j.LE.j_end)) THEN
                     IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                        rhsvalue(nn) = whole_object(nobj)%phi
                     ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                        rhsvalue(nn) = whole_object(nobj)%phi_profile(j)
                     END IF
                     EXIT
                  END IF
               END DO
            ELSE
               ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
               CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(indx_x_min,j,neumann_flag)

               ! This point is Neumann, I shall proceed
               IF ( neumann_flag ) THEN
                  rhsvalue(nn) = factor_rho * (rho_i(indx_x_min,j) - rho_e(indx_x_min,j))/two
                  ! rhsvalue(nn) = -(k_x_cart_neumann_right**2+k_y_cart_neumann_right**2)*(one)*SIN(k_y_cart_neumann_right*DBLE(j)*delta_x_m)/ F_scale_V*delta_x_m**2/two ! Test cartesian. Solution phi = cos(kx*x)*sin(ky*y). Dirichlet 0 at top/bottom/right and neumann at left
               ! Metal
               ELSE 
                  DO n = 1, N_of_local_object_parts_left
                     m = index_of_local_object_part_left(n)
                     j_start = local_object_part(m)%jstart
                     j_end   = local_object_part(m)%jend
                     nobj   = local_object_part(m)%object_number
                     IF ((j.GE.j_start).AND.(j.LE.j_end)) THEN
                        IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                           rhsvalue(nn) = whole_object(nobj)%phi
                        ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                           rhsvalue(nn) = whole_object(nobj)%phi_profile(j)
                        END IF
                        EXIT
                     END IF
                  END DO
               ENDIF                  
            ENDIF
        ELSE !IF (whole_object(nobj)%object_type.EQ.SYMMETRY_PLANE) THEN
           rhsvalue(nn) = factor_rho * (rho_i(indx_x_min,j) - rho_e(indx_x_min,j))
         !   rhsvalue(nn) = (-(k_test_r_neumann_right**2+k_test_z_neumann_right**2)*(one)*SIN(k_test_z_neumann_right*DBLE(j)*delta_x_m)-k_test_r_neumann_right**2*SIN(k_test_z_neumann_right*delta_x_m*j))/ F_scale_V*delta_x_m**2!/two ! Test cylidnrical. Solution phi = cos(kx*x)*sin(ky*y). Dirichlet 0 at bottom/top and neumann at right. symmetry axis also
         !   rhsvalue(nn) = (-(k_test_r**2+k_test_z**2)*COS(k_test_z*delta_x_m*j)-k_test_r**2*COS(k_test_z*delta_x_m*j))/ F_scale_V*delta_x_m**2  ! CYlindrical basic test. SOlution phi = cos(kr*r)*cos(kz*z). Left is Neumann, bottom is Neumann. Dirichlet top and right
         !   rhsvalue(nn) = -(k_test_r+(k_test_r**2+k_test_z**2))*SIN(k_test_z*delta_x_m*j)/ F_scale_V*delta_x_m**2!/two ! Ttest cylindrical with double Neumann
           !   rhsvalue(nn) = zero
         !   rhsvalue(nn) = -(k_test_r**2+(k_test_r**2+k_test_z**2))*SIN(k_test_z_dirichlet*delta_x_m*j)/ F_scale_V*delta_x_m**2 !phi = cos(kr*r)*sin(kz*z)
           ! Check if I am in inner dielecctric
         !   CALL FIND_INNER_OBJECT_CONTAINING_POINT(0,j,ni0,position_flag)
         !   IF ( position_flag<0 ) THEN ! in plasma at symmetry axis
         !       ! print*,'at symmetry axis in plasma',nn
            !   rhsvalue(nn) = -(k_test_r**2+(k_test_r**2+k_test_z**2))*SIN(k_test_z*delta_x_m*j)/ F_scale_V*delta_x_m**2/two
         !   ELSE IF ( position_flag==1 .OR. position_flag==7 ) THEN ! one of the corners at symmetry axis as well for dielectric object. Note i=indx_x_min
         !      print*,'at symmetry axis at surface inner object nn,i,j,indx_x_min',nn,i,j,indx_x_min
         !      !! If bottom dielectric
         !      IF ( j<index_maxi_z/2 ) rhsvalue(nn) = -(k_test_z*delta_x_m/two)/ F_scale_V ! Field will be zero inside dielectric
         !    !   IF ( j<index_maxi_z/2 ) rhsvalue(nn) = -((k_test_z+one*10000)*delta_x_m/two)/ F_scale_V ! Field will have Bessel*sinh variations inside dielectric
         !      !! If top dielectric  
         !      IF ( j>index_maxi_z/2 ) rhsvalue(nn) = (k_test_z*delta_x_m/two)/ F_scale_V ! Field will be zero inside dielectric
         !    !   IF ( j>index_maxi_z/2 ) rhsvalue(nn) = ((k_test_z+one*10000)*delta_x_m/two)/ F_scale_V ! Field will be zero inside dielectric
         !   ELSE ! I am inside inner object: no charge
         !    ! print*,'at symmetry axis in inner object nn,i,j',nn,i,j  
         !    rhsvalue(nn) = zero
         !   ENDIF
         END IF
     END IF
   !   print*,'valeur min index indx_x_min',indx_x_min
     DO i = indx_x_min+1, indx_x_max-1
        nn = nn + 1
        rhsvalue(nn) = factor_rho * (rho_i(i,j) - rho_e(i,j))
      !   rhsvalue(nn) = -(k_x_cart_neumann_right**2+k_y_cart_neumann_right**2)*COS(k_x_cart_neumann_right*DBLE(i)*delta_x_m)*SIN(k_y_cart_neumann_right*DBLE(j)*delta_x_m)/ F_scale_V*delta_x_m**2 ! Test cartesian. Solution phi = cos(kx*x)*sin(ky*y). Dirichlet 0 at top/bottom/right and neumann at left
      !   rhsvalue(nn) = -(k_test_r*SIN(delta_x_m*k_test_r_dirichlet*i)/(delta_x_m*i)+(k_test_r**2+k_test_z**2)*COS(delta_x_m*k_test_r_dirichlet*i))*SIN(k_test_z_dirichlet*delta_x_m*j)/ F_scale_V*delta_x_m**2 !Solution phi= cos(k_r*r)*sin(k_z*z). Dirichlet 
      !   rhsvalue(nn) = (-(k_test_r**2+k_test_z**2)*COS(k_test_r*delta_x_m*i)*COS(k_test_z*delta_x_m*j)-k_test_r/(DBLE(i)*delta_x_m)*SIN(k_test_r*delta_x_m*i)*COS(k_test_z*delta_x_m*j))/ F_scale_V*delta_x_m**2  ! CYlindrical. Basic test. SOlution phi = cos(kr*r)*cos(kz*z). Left is Neumann, bottom is Neumann. Dirichlet top and right
      !   rhsvalue(nn) = (-(k_test_r**2+k_test_z**2)*COS(k_test_r*delta_x_m*i)*SIN(k_test_z*delta_x_m*j)-k_test_r/(DBLE(i)*delta_x_m)*SIN(k_test_r*delta_x_m*i)*SIN(k_test_z*delta_x_m*j))/ F_scale_V*delta_x_m**2 
        !  rhsvalue(nn) = (-(k_test_r**2+k_test_z**2)*COS(k_test_r*delta_x_m*i)*SIN(k_test_z*delta_x_m*j)-k_test_r/(DBLE(i)*delta_x_m)*SIN(k_test_r*delta_x_m*i)*SIN(k_test_z*delta_x_m*j))/ F_scale_V*delta_x_m**2
      !   rhsvalue(nn) = -(k_x_cart_neumann_right**2+k_y_cart_neumann_right**2)*SIN(k_x_cart_neumann_right*DBLE(i)*delta_x_m)*SIN(k_y_cart_neumann_right*DBLE(j)*delta_x_m)/ F_scale_V*delta_x_m**2 ! Test cartesian. Solution phi = sin(kx*x)*sin(ky*y). Dirichlet 0 at left/bottom/top and neumann at right
      !   rhsvalue(nn) = (-(k_test_r_neumann_right**2+k_test_z_neumann_right**2)*COS(k_test_r_neumann_right*DBLE(i)*delta_x_m)*SIN(k_test_z_neumann_right*DBLE(j)*delta_x_m)-k_test_r_neumann_right/(DBLE(i)*delta_x_m)*SIN(k_test_r_neumann_right*delta_x_m*i)*SIN(k_test_z_neumann_right*delta_x_m*j))/ F_scale_V*delta_x_m**2 ! Test cylidnrical. Solution phi = cos(kx*x)*sin(ky*y). Dirichlet 0 at bottom/top and neumann at right. symmetry axis also
        !   rhsvalue(nn) = -(k_x_cart**2+k_y_cart**2)*SIN(k_x_cart*DBLE(i)*delta_x_m)*SIN(k_y_cart*DBLE(j)*delta_x_m)/ F_scale_V*delta_x_m**2 ! Test cartesian. Solution phi = sin(kx*x)*sin(ky*y). Dirichlet 0 at left/bottom/right and neumann at top
         ! rhsvalue(nn) = -(k_x_cart**2+k_y_cart**2)*SIN(k_x_cart*DBLE(i)*delta_x_m)*COS(k_y_cart*DBLE(j)*delta_x_m)/ F_scale_V*delta_x_m**2 ! Test cartesian. Solution phi = sin(kx*x)*cos(ky*y). Dirichlet 0 at left/top/right and neumann at bottom      
      !   rhsvalue(nn) = zero!-(k_test_r*SIN(delta_x_m*k_test_r*i)/(delta_x_m*i)+(k_test_r**2+k_test_z**2)*COS(delta_x_m*k_test_r*i))*SIN(k_test_z*delta_x_m*j)/ F_scale_V*delta_x_m**2 !Solution phi= cos(k_r*r)*sin(k_z*z). Dirichlet at BCs z=0, z=2cm,r=3cm. Divide by scaling factor
      !   CALL FIND_INNER_OBJECT_CONTAINING_POINT(i,j,ni0,position_flag)
      !   IF ( position_flag<0 ) THEN ! imposed right hand side to test solver in cyl coordinates
      !       rhsvalue(nn) = -(k_test_r*SIN(delta_x_m*k_test_r*i)/(delta_x_m*i)+(k_test_r**2+k_test_z**2)*COS(delta_x_m*k_test_r*i))*SIN(k_test_z*delta_x_m*j)/ F_scale_V*delta_x_m**2 !Solution phi= cos(k_r*r)*sin(k_z*z). Dirichlet at BCs z=0, z=2cm,r=3cm. Divide by scaling factor
      !   ELSE IF ( position_flag>0 .AND. position_flag<9 ) THEN ! Set appropriate surfce charge distribution when there is a dielectric object 
      !       ! print*,'on surface,nn,position_flag,i,j',nn,position_flag,i,j
      !       !! if object is on the right
      !       IF ( i>index_maxi_r/2 ) rhsvalue(nn) = (k_test_r*delta_x_m*SIN(k_test_z*delta_x_m*j))/ F_scale_V
      !       !! If bottom dielectric
      !       ! IF ( j<index_maxi_z/2 ) rhsvalue(nn) = -(k_test_z*delta_x_m*COS(k_test_r*delta_x_m*i))/ F_scale_V ! case eps=1 (should not mater) and field =0 inside dielectric 
      !       ! IF ( j<index_maxi_z/2 ) rhsvalue(nn) = -(k_test_z*delta_x_m*COS(k_test_r*delta_x_m*i) + 10000*BESSEL_J0(alpha_r*delta_x_m*i)*delta_x_m)/ F_scale_V ! eps = 4 (matters) and field is sinh(alpha*z)*J_0(alpha*r) in dielectric
      !       !! If top dielectric
      !       ! IF ( j>index_maxi_z/2 ) rhsvalue(nn) = (k_test_z*delta_x_m*COS(k_test_r*delta_x_m*i))/ F_scale_V ! case eps=1 (should not mater) and field =0 inside dielectric 
      !       ! IF ( j>index_maxi_z/2 ) rhsvalue(nn) = (k_test_z*delta_x_m*COS(k_test_r*delta_x_m*i)+ 10000*BESSEL_J0(alpha_r*delta_x_m*i)*delta_x_m)/ F_scale_V ! eps = 4 (matters) and field is sinh(alpha*z)*J_0(alpha*r) in dielectric
      !    ELSE
      !       ! print*,'inside',nn
      !       rhsvalue(nn) = zero
      !   ENDIF
        IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) rhsvalue(nn) = given_F_double_period_sys
         
     END DO
   !   print*,'maxi,mini',MAXVAL(rhsvalue),MINVAL(rhsvalue)

!     i = indx_x_max

     IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
        nn = nn + 1
        IF ( .NOT.block_has_neumann_bc_X_right ) THEN
            DO n = 1, N_of_local_object_parts_right
               m = index_of_local_object_part_right(n)
               j_start = local_object_part(m)%jstart
               j_end   = local_object_part(m)%jend
               nobj   = local_object_part(m)%object_number
               IF ((j.GE.j_start).AND.(j.LE.j_end)) THEN
                  IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                     rhsvalue(nn) = whole_object(nobj)%phi
                  ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                     rhsvalue(nn) = whole_object(nobj)%phi_profile(j)
                  END IF
                  EXIT
               END IF
            END DO
         ELSE ! Neuman right
            ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(indx_x_max,j,neumann_flag)

            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN            
               ! rhsvalue(nn) = -(k_x_cart_neumann_right**2+k_y_cart_neumann_right**2)*(-one)*SIN(k_y_cart_neumann_right*DBLE(j)*delta_x_m)/ F_scale_V*delta_x_m**2/two ! Test cartesian. Solution phi = sin(kx*x)*sin(ky*y). Dirichlet 0 at left/bottom/top and neumann at right            
               ! IF (i==0) THEN
               !    rhsvalue(nn) = (-(k_test_r_neumann_right**2+k_test_z_neumann_right**2)*(-one)*SIN(k_test_z_neumann_right*DBLE(j)*delta_x_m)-k_test_r_neumann_right**2*SIN(k_test_z_neumann_right*delta_x_m*j))/ F_scale_V*delta_x_m**2/two ! Test cylidnrical. Solution phi = cos(kx*x)*sin(ky*y). Dirichlet 0 at bottom/top and neumann at right. symmetry axis also
               ! ELSE
               !    rhsvalue(nn) = (-(k_test_r_neumann_right**2+k_test_z_neumann_right**2)*(-one)*SIN(k_test_z_neumann_right*DBLE(j)*delta_x_m)-k_test_r_neumann_right/(delta_x_m*DBLE(i))*SIN(k_test_r_neumann_right*delta_x_m*DBLE(i))*SIN(k_test_z_neumann_right*delta_x_m*j))/ F_scale_V*delta_x_m**2/two ! Test cylidnrical. Solution phi = cos(kx*x)*sin(ky*y). Dirichlet 0 at bottom/top and neumann at right. symmetry axis also
               ! ENDIF
               rhsvalue(nn) = factor_rho * (rho_i(indx_x_max,j) - rho_e(indx_x_max,j))/two ! Nodal volue at top is dx**2/2 in cartesian. In cylindrical, I keep the same convention, ie the RHS has a factoer dx**2/2, everythong else is put on the LHS  
               !METAL   
            ELSE
               DO n = 1, N_of_local_object_parts_right
                  m = index_of_local_object_part_right(n)
                  j_start = local_object_part(m)%jstart
                  j_end   = local_object_part(m)%jend
                  nobj   = local_object_part(m)%object_number
                  IF ((j.GE.j_start).AND.(j.LE.j_end)) THEN
                     IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                        rhsvalue(nn) = whole_object(nobj)%phi
                     ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                        rhsvalue(nn) = whole_object(nobj)%phi_profile(j)
                     END IF
                     EXIT
                  END IF
               END DO               
            ENDIF
         END IF
     END IF

  END DO   !### DO j = indx_y_min+1, indx_y_max-1

!    j = indx_y_max !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
  IF (jend.EQ.indx_y_max) THEN
! boundary object along top border
     DO i = ibegin, iend
        nn = nn + 1
        IF ( .NOT.block_has_neumann_bc_Y_top ) THEN
            DO n = 1, N_of_local_object_parts_above
               m = index_of_local_object_part_above(n)
               i_start = local_object_part(m)%istart
               i_end   = local_object_part(m)%iend
               nobj   = local_object_part(m)%object_number
               IF ((i.GE.i_start).AND.(i.LE.i_end)) THEN
                  IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                     rhsvalue(nn) = whole_object(nobj)%phi
                  ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                     rhsvalue(nn) = whole_object(nobj)%phi_profile(i)
                  END IF
                  EXIT
               END IF
            END DO
         ELSE
            ! Double check if current point is Neumann or not (a cluster could have both Neumann and metal)
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,indx_y_max,neumann_flag)

            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN                    
               ! Cylindrical
               ! IF (i==0) THEN
               !    rhsvalue(nn) =(-(k_test_r**2+k_test_z**2)*COS(k_test_r*delta_x_m*i)*SIN(k_test_z*delta_x_m*j)-k_test_r*SIN(k_test_z*delta_x_m*j))/ F_scale_V*delta_x_m**2/2.0!TEst cylindrical double neumann
               ! ELSE
               !    rhsvalue(nn) =(-(k_test_r**2+k_test_z**2)*COS(k_test_r*delta_x_m*i)*SIN(k_test_z*delta_x_m*j)-k_test_r/(DBLE(i)*delta_x_m)*SIN(k_test_r*delta_x_m*i)*SIN(k_test_z*delta_x_m*j))/ F_scale_V*delta_x_m**2/2.0!TEst cylindrical double neumann
               ! END IF
               !    rhsvalue(nn) = -(k_x_cart**2+k_y_cart**2)*SIN(k_x_cart*DBLE(i)*delta_x_m)*(-one)/ F_scale_V*delta_x_m**2/2.0 ! Cartesian
               rhsvalue(nn) = factor_rho * (rho_i(i,indx_y_max) - rho_e(i,indx_y_max))/two ! Nodal volue at top is dx**2/2 in cartesian. In cylindrical, I keep the same convention, ie the RHS has a factoer dx**2/2, everythong else is put on the LHS  
            ELSE
               DO n = 1, N_of_local_object_parts_above
                  m = index_of_local_object_part_above(n)
                  i_start = local_object_part(m)%istart
                  i_end   = local_object_part(m)%iend
                  nobj   = local_object_part(m)%object_number
                  IF ((i.GE.i_start).AND.(i.LE.i_end)) THEN
                     IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                        rhsvalue(nn) = whole_object(nobj)%phi
                     ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                        rhsvalue(nn) = whole_object(nobj)%phi_profile(i)
                     END IF
                     EXIT
                  END IF
               END DO               
            ENDIF
         END IF
     END DO
  END IF

! fool proof
  IF (nn.NE.block_N_of_nodes_to_solve) THEN
     PRINT '("proc ",i4," :: Error-1 in SOLVE_POTENTIAL_WITH_PETSC ",2(2x,i9))', Rank_of_process, nn, block_N_of_nodes_to_solve
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

! find metal inner objects 
  nn = 0
  DO j = jbegin, jend
     DO i = ibegin, iend
        nn = nn + 1
        DO nobj = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (whole_object(nobj)%object_type.NE.METAL_WALL) CYCLE
           IF (i.LT.whole_object(nobj)%ileft) CYCLE
           IF (i.GT.whole_object(nobj)%iright) CYCLE
           IF (j.LT.whole_object(nobj)%jbottom) CYCLE
           IF (j.GT.whole_object(nobj)%jtop) CYCLE
           rhsvalue(nn) = whole_object(nobj)%phi
!print *, i, j, rhsvalue(nn)
        END DO
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!print '("proc ",i4," before InsertData")', Rank_of_process

  call InsertData(nn, rhsvalue)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!print '("proc ",i4," before Solve")', Rank_of_process
!stop
!   print*,'test rhs',rhsvalue
  call Solve

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!print '("proc ",i4," before FetchData")', Rank_of_process
!stop

  call FetchData(nn, queue)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!print '("proc ",i4," after FetchData")', Rank_of_process
!stop

  nn = 0
  DO j = jbegin, jend
     DO i = ibegin, iend
        nn = nn+1
        phi(i,j) = queue(nn)
     END DO
  END DO

! exchange potential values in overlap nodes with neighbor processes (blocks)

  n1 = indx_y_max - indx_y_min + 1
  n3 = indx_x_max - indx_x_min + 1

  IF (WHITE) THEN  
! "white processes"

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

     IF (Rank_of_process_right.GE.0) THEN
! ## 1 ## send right potential in column left of the right edge
        rbufer(1:n1) = phi(indx_x_max-1, indx_y_min:indx_y_max)
        CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 2 ## send left potential in column right of the left edge
        rbufer(1:n1) = phi(indx_x_min+1, indx_y_min:indx_y_max)
        CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 3 ## receive from left potential along the left edge
        CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min, indx_y_min:indx_y_max) = rbufer(1:n1)   
     END IF

     IF (Rank_of_process_right.GE.0) THEN
! ## 4 ## receive from right potential along the right edge
        CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process_right, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_max, indx_y_min:indx_y_max) = rbufer(1:n1)
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

     IF (Rank_of_process_above.GE.0) THEN
! ## 5 ## send up potential in the row below the top edge
        rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_max-1)
        CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 6 ## send down potential in the row above the bottom edge
        rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_min+1)
        CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 7 ## receive from below potential in the bottom edge
        CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process_below, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min:indx_x_max, indx_y_min) = rbufer(1:n3)
     END IF

     IF (Rank_of_process_above.GE.0) THEN
! ## 8 ## receive from above potential in the top edge
        CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process_above, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min:indx_x_max, indx_y_max) = rbufer(1:n3)
     END IF

  ELSE
! "black" processes

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

     IF (Rank_of_process_left.GE.0) THEN
! ## 1 ## receive from left potential along the left edge
        CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min, indx_y_min:indx_y_max) = rbufer(1:n1)   
     END IF

     IF (Rank_of_process_right.GE.0) THEN
! ## 2 ## receive from right potential along the right edge
        CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process_right, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_max, indx_y_min:indx_y_max) = rbufer(1:n1)
     END IF

     IF (Rank_of_process_right.GE.0) THEN
! ## 3 ## send right potential in column left of the right edge
        rbufer(1:n1) = phi(indx_x_max-1, indx_y_min:indx_y_max)
        CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 4 ## send left potential in column right of the left edge
        rbufer(1:n1) = phi(indx_x_min+1, indx_y_min:indx_y_max)
        CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

     IF (Rank_of_process_below.GE.0) THEN
! ## 5 ## receive from below potential in the bottom edge
        CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process_below, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min:indx_x_max, indx_y_min) = rbufer(1:n3)
     END IF

     IF (Rank_of_process_above.GE.0) THEN
! ## 6 ## receive from above potential in the top edge
        CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process_above, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min:indx_x_max, indx_y_max) = rbufer(1:n3)
     END IF

     IF (Rank_of_process_above.GE.0) THEN
! ## 7 ## send up potential in the row below the top edge
        rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_max-1)
        CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 8 ## send down potential in the row above the bottom edge
        rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_min+1)
        CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

  END IF

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

! ################ diagnostics, electrostatic potential #################
  IF (T_cntr.EQ.Save_probes_e_data_T_cntr) THEN
     DO npb = 1, N_of_probes_block
        i = Probe_params_block_list(1,npb)
        j = Probe_params_block_list(2,npb)
        probe_F_block(npb) = phi(i,j)
     END DO
  END IF
!

  DEALLOCATE(rhsvalue, STAT=ALLOC_ERR)
  DEALLOCATE(queue,   stat=alloc_err)
 
END SUBROUTINE SOLVE_POTENTIAL_WITH_PETSC


!-------------------------------------------------------
!
SUBROUTINE CALCULATE_ELECTRIC_FIELD

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  REAL(8) factor2_E_from_F
  REAL(8) factor1_E_from_F

  INTEGER i, j, k, pos, n1, n3, pos1, pos2
  INTEGER shift, bufsize 

  REAL(8), ALLOCATABLE :: rbufer(:)
  REAL(8), ALLOCATABLE :: loc_EX(:,:)
  REAL(8), ALLOCATABLE :: loc_EY(:,:)
  INTEGER ALLOC_ERR
  LOGICAL :: neumann_flag

  neumann_flag = .FALSE.
  factor2_E_from_F = F_scale_V / (E_scale_Vm * 2.0_8 * delta_x_m)
  factor1_E_from_F = F_scale_V / (E_scale_Vm * delta_x_m)

  ! Deactivate electric field and poisson
   IF ( i_no_poisson==1 )   phi = zero

  IF (cluster_rank_key.EQ.0) THEN


! volume
     DO j = indx_y_min+1, indx_y_max-1
        DO i = indx_x_min+1, indx_x_max-1
           EX(i,j) = factor2_E_from_F * (phi(i-1,j)-phi(i+1,j))
           EY(i,j) = factor2_E_from_F * (phi(i,j-1)-phi(i,j+1))
        END DO
     END DO
! left edge if there is no neighbour
     IF (Rank_of_process_left.LT.0) THEN
        i = indx_x_min
        DO j = indx_y_min+1, indx_y_max-1
           EX(i,j) = factor1_E_from_F * (phi(i,j)-phi(i+1,j))
           EY(i,j) = factor2_E_from_F * (phi(i,j-1)-phi(i,j+1))           
        END DO
     END IF
! right edge if there is no neighbour
     IF (Rank_of_process_right.LT.0) THEN
        i = indx_x_max
        DO j = indx_y_min+1, indx_y_max-1
           EX(i,j) = factor1_E_from_F * (phi(i-1,j)-phi(i,j))
           EY(i,j) = factor2_E_from_F * (phi(i,j-1)-phi(i,j+1))           
        END DO
     END IF
! edge above if there is no neighbour
     IF (Rank_of_process_above.LT.0) THEN
        j = indx_y_max
        DO i = indx_x_min+1, indx_x_max-1
           EX(i,j) = factor2_E_from_F * (phi(i-1,j)-phi(i+1,j))
           EY(i,j) = factor1_E_from_F * (phi(i,j-1)-phi(i,j))           
        END DO
     END IF
! edge below if there is no neighbour
     IF (Rank_of_process_below.LT.0) THEN
        j = indx_y_min
        DO i = indx_x_min+1, indx_x_max-1
           EX(i,j) = factor2_E_from_F * (phi(i-1,j)-phi(i+1,j))
           EY(i,j) = factor1_E_from_F * (phi(i,j)-phi(i,j+1))           
        END DO
     END IF
! use ranks of neighbor processes (blocks) because only corners of clusters are identified 
! left bottom corner if it is surrounded by walls
     IF ((Rank_of_process_left.LT.0).AND.(Rank_of_process_below.LT.0)) THEN
        i = indx_x_min
        j = indx_y_min
        EX(i,j) = factor1_E_from_F * (phi(i,j)-phi(i+1,j))
        EY(i,j) = factor1_E_from_F * (phi(i,j)-phi(i,j+1))  
     END IF
! left top corner if it is surrounded by walls
     IF ((Rank_of_process_left.LT.0).AND.(Rank_of_process_above.LT.0)) THEN
        i = indx_x_min
        j = indx_y_max
        EX(i,j) = factor1_E_from_F * (phi(i,j)-phi(i+1,j))
        EY(i,j) = factor1_E_from_F * (phi(i,j-1)-phi(i,j))  
     END IF
! right bottom corner if it is surrounded by walls
     IF ((Rank_of_process_right.LT.0).AND.(Rank_of_process_below.LT.0)) THEN
        i = indx_x_max
        j = indx_y_min
        EX(i,j) = factor1_E_from_F * (phi(i-1,j)-phi(i,j))
        EY(i,j) = factor1_E_from_F * (phi(i,j)-phi(i,j+1))  
     END IF
! right top corner if it is surrounded by walls
     IF ((Rank_of_process_right.LT.0).AND.(Rank_of_process_above.LT.0)) THEN
        i = indx_x_max
        j = indx_y_max
        EX(i,j) = factor1_E_from_F * (phi(i-1,j)-phi(i,j))
        EY(i,j) = factor1_E_from_F * (phi(i,j-1)-phi(i,j))  
     END IF

! receive fields from field calculators
     DO k = 2, cluster_N_blocks

        shift = (field_calculator(k)%indx_x_max - field_calculator(k)%indx_x_min + 1) * &
              & (field_calculator(k)%indx_y_max - field_calculator(k)%indx_y_min + 1)
        bufsize = 2*shift
        ALLOCATE(rbufer(bufsize), STAT = ALLOC_ERR)
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)

        pos=1
        j = field_calculator(k)%indx_y_min
        IF (field_calculator(k)%indx_y_min.EQ.c_indx_y_min) THEN
           IF (field_calculator(k)%indx_x_min.EQ.c_indx_x_min) THEN
! left bottom corner of block is left bottom corner of cluster
              EX(c_indx_x_min,j) = rbufer(pos)
              EY(c_indx_x_min,j) = rbufer(pos+shift)
           END IF
           pos = pos + 1
! bottom edge of block is on bottom edge of cluster
           DO i = field_calculator(k)%indx_x_min+1, field_calculator(k)%indx_x_max-1
              EX(i,j) = rbufer(pos)
              EY(i,j) = rbufer(pos+shift)
              pos = pos + 1
           END DO
           IF (field_calculator(k)%indx_x_max.EQ.c_indx_x_max) THEN
! right bottom corner of block is right bottom corner of cluster
              EX(c_indx_x_max,j) = rbufer(pos)
              EY(c_indx_x_max,j) = rbufer(pos+shift)
           END IF
           pos = pos + 1
        END IF

        pos = field_calculator(k)%indx_x_max - field_calculator(k)%indx_x_min + 2
        DO j = field_calculator(k)%indx_y_min+1, field_calculator(k)%indx_y_max-1
           IF (field_calculator(k)%indx_x_min.EQ.c_indx_x_min) THEN
              EX(c_indx_x_min,j) = rbufer(pos)
              EY(c_indx_x_min,j) = rbufer(pos+shift)
           END IF
           pos = pos + 1
           DO i = field_calculator(k)%indx_x_min+1, field_calculator(k)%indx_x_max-1
              EX(i,j) = rbufer(pos)
              EY(i,j) = rbufer(pos+shift)
              pos = pos + 1
           END DO
           IF (field_calculator(k)%indx_x_max.EQ.c_indx_x_max) THEN
              EX(c_indx_x_max,j) = rbufer(pos)
              EY(c_indx_x_max,j) = rbufer(pos+shift)
           END IF
           pos = pos + 1      
        END DO

        j = field_calculator(k)%indx_y_max
        IF (field_calculator(k)%indx_y_max.EQ.c_indx_y_max) THEN
           IF (field_calculator(k)%indx_x_min.EQ.c_indx_x_min) THEN
! left top corner of block is left top corner of cluster
              EX(c_indx_x_min,j) = rbufer(pos)
              EY(c_indx_x_min,j) = rbufer(pos+shift)
           END IF
           pos = pos + 1
! top edge of block is on top edge of cluster
           DO i = field_calculator(k)%indx_x_min+1, field_calculator(k)%indx_x_max-1
              EX(i,j) = rbufer(pos)
              EY(i,j) = rbufer(pos+shift)
              pos = pos+1
           END DO
           IF (field_calculator(k)%indx_x_max.EQ.c_indx_x_max) THEN
! right top corner of block is right top corner of cluster
              EX(c_indx_x_max,j) = rbufer(pos)
              EY(c_indx_x_max,j) = rbufer(pos+shift)
           END IF
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)

     END DO
      ! enforce boundary condition EX=0 at the symmetry plane
      IF ( symmetry_plane_X_left )  EX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = 0.0_8
      IF ( neumann_Y_cluster_top ) THEN
         DO i=c_indx_x_min, c_indx_x_max
            
            ! Check if point i,j is Neumann 
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,c_indx_y_max,neumann_flag)
            
            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN            
               EY(i, c_indx_y_max) = zero
            ENDIF 
         ENDDO
      END IF
      IF ( neumann_Y_cluster_bottom ) THEN
         DO i=c_indx_x_min, c_indx_x_max
            
            ! Check if point i,j is Neumann 
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(i,c_indx_y_min,neumann_flag)
            
            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN            
               EY(i, c_indx_y_min) = zero
            ENDIF 
         ENDDO      
      ENDIF
      IF ( neumann_X_cluster_left ) THEN 
         DO j=c_indx_y_min, c_indx_y_max
            
            ! Check if point i,j is Neumann 
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(c_indx_x_min,j,neumann_flag)
            
            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN       
               EX(c_indx_x_min, j) = zero     
            ENDIF 
         ENDDO      
      END IF
         
      IF ( neumann_X_cluster_right ) THEN
         DO j=c_indx_y_min, c_indx_y_max
            
            ! Check if point i,j is Neumann 
            CALL DECIDE_NEUMANN_EXTERNAL_BOUNDARY(c_indx_x_max,j,neumann_flag)
            
            ! This point is Neumann, I shall proceed
            IF ( neumann_flag ) THEN       
               EX(c_indx_x_max, j) = zero     
            ENDIF 
         ENDDO      
      END IF              
      

! exchange boundary field values with neighbor masters

     n1 = c_indx_y_max - c_indx_y_min + 1
     n3 = c_indx_x_max - c_indx_x_min + 1

     IF (WHITE_CLUSTER) THEN  
! "white processes"

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:(n1+n1)), STAT=ALLOC_ERR)

        IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right electric fields in column left of the right edge
           rbufer(1:n1)           = EX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer((n1+1):(n1+n1)) = EY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left electric field in column right of the left edge
           rbufer(1:n1)           = EX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer((n1+1):(n1+n1)) = EY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left electric field along the left edge
           CALL MPI_RECV(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(1:n1)   
           EY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer((n1+1):(n1+n1))
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right electric field along the right edge
           CALL MPI_RECV(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(1:n1)
           EY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer((n1+1):(n1+n1))
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:(n3+n3)), STAT=ALLOC_ERR)

        IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up electric fields in the row below the top edge
           rbufer(1:n3)           = EX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer((n3+1):(n3+n3)) = EY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           CALL MPI_SEND(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down electric fields in the row above the bottom edge
           rbufer(1:n3)           = EX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer((n3+1):(n3+n3)) = EY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           CALL MPI_SEND(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below electric fields in the bottom edge
           CALL MPI_RECV(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(1:n3)
           EY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer((n3+1):(n3+n3))
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above electric fields in the top edge
           CALL MPI_RECV(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(1:n3)
           EY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer((n3+1):(n3+n3))
        END IF

     ELSE
! "black" processes

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:(n1+n1)), STAT=ALLOC_ERR)

        IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left electric field along the left edge
           CALL MPI_RECV(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(1:n1)   
           EY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer((n1+1):(n1+n1))
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right electric field along the right edge
           CALL MPI_RECV(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(1:n1)
           EY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer((n1+1):(n1+n1))
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right electric fields in column left of the right edge
           rbufer(1:n1)           = EX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer((n1+1):(n1+n1)) = EY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left electric field in column right of the left edge
           rbufer(1:n1)           = EX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer((n1+1):(n1+n1)) = EY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:(n3+n3)), STAT=ALLOC_ERR)

        IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below electric fields in the bottom edge
           CALL MPI_RECV(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(1:n3)
           EY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer((n3+1):(n3+n3))
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above electric fields in the top edge
           CALL MPI_RECV(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(1:n3)
           EY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer((n3+1):(n3+n3))
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up electric fields in the row below the top edge
           rbufer(1:n3)           = EX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer((n3+1):(n3+n3)) = EY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           CALL MPI_SEND(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down electric fields in the row above the bottom edge
           rbufer(1:n3)           = EX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer((n3+1):(n3+n3)) = EY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           CALL MPI_SEND(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

! ready to send complete field array to all members of the cluster

  ELSE

! we make the array of the size of the whole block, however, the edge values will be used in the field master
! only if the edge is where there is no neighbour
     ALLOCATE(loc_EX(indx_x_min:indx_x_max, indx_y_min:indx_y_max), STAT = ALLOC_ERR)
     ALLOCATE(loc_EY(indx_x_min:indx_x_max, indx_y_min:indx_y_max), STAT = ALLOC_ERR)

     DO j = indx_y_min+1, indx_y_max-1
        DO i = indx_x_min+1, indx_x_max-1
           loc_EX(i,j) = factor2_E_from_F * (phi(i-1,j)-phi(i+1,j))
           loc_EY(i,j) = factor2_E_from_F * (phi(i,j-1)-phi(i,j+1))
        END DO
     END DO
! left edge if there is no neighbour
     IF (Rank_of_process_left.LT.0) THEN
        i = indx_x_min
        DO j = indx_y_min+1, indx_y_max-1
           loc_EX(i,j) = factor1_E_from_F * (phi(i,j)-phi(i+1,j))
           loc_EY(i,j) = factor2_E_from_F * (phi(i,j-1)-phi(i,j+1))           
        END DO
     END IF
! right edge if there is no neighbour
     IF (Rank_of_process_right.LT.0) THEN
        i = indx_x_max
        DO j = indx_y_min+1, indx_y_max-1
           loc_EX(i,j) = factor1_E_from_F * (phi(i-1,j)-phi(i,j))
           loc_EY(i,j) = factor2_E_from_F * (phi(i,j-1)-phi(i,j+1))           
        END DO
     END IF
! edge above if there is no neighbour
     IF (Rank_of_process_above.LT.0) THEN
        j = indx_y_max
        DO i = indx_x_min+1, indx_x_max-1
           loc_EX(i,j) = factor2_E_from_F * (phi(i-1,j)-phi(i+1,j))
           loc_EY(i,j) = factor1_E_from_F * (phi(i,j-1)-phi(i,j))           
        END DO
     END IF
! edge below if there is no neighbour
     IF (Rank_of_process_below.LT.0) THEN
        j = indx_y_min
        DO i = indx_x_min+1, indx_x_max-1
           loc_EX(i,j) = factor2_E_from_F * (phi(i-1,j)-phi(i+1,j))
           loc_EY(i,j) = factor1_E_from_F * (phi(i,j)-phi(i,j+1))           
        END DO
     END IF
! use ranks of neighbor processes (blocks) because only corners of clusters are identified 
! besides, the process may belong to a different cluster (not the one where the master is its field master)
! left bottom corner if it is surrounded by walls
     IF ((Rank_of_process_left.LT.0).AND.(Rank_of_process_below.LT.0)) THEN
        i = indx_x_min
        j = indx_y_min
        loc_EX(i,j) = factor1_E_from_F * (phi(i,j)-phi(i+1,j))
        loc_EY(i,j) = factor1_E_from_F * (phi(i,j)-phi(i,j+1))  
     END IF
! left top corner if it is surrounded by walls
     IF ((Rank_of_process_left.LT.0).AND.(Rank_of_process_above.LT.0)) THEN
        i = indx_x_min
        j = indx_y_max
        loc_EX(i,j) = factor1_E_from_F * (phi(i,j)-phi(i+1,j))
        loc_EY(i,j) = factor1_E_from_F * (phi(i,j-1)-phi(i,j))  
     END IF
! right bottom corner if it is surrounded by walls
     IF ((Rank_of_process_right.LT.0).AND.(Rank_of_process_below.LT.0)) THEN
        i = indx_x_max
        j = indx_y_min
        loc_EX(i,j) = factor1_E_from_F * (phi(i-1,j)-phi(i,j))
        loc_EY(i,j) = factor1_E_from_F * (phi(i,j)-phi(i,j+1))  
     END IF
! right top corner if it is surrounded by walls
     IF ((Rank_of_process_right.LT.0).AND.(Rank_of_process_above.LT.0)) THEN
        i = indx_x_max
        j = indx_y_max
        loc_EX(i,j) = factor1_E_from_F * (phi(i-1,j)-phi(i,j))
        loc_EY(i,j) = factor1_E_from_F * (phi(i,j-1)-phi(i,j))  
     END IF

! send to its field master
     bufsize = 2 * (indx_x_max - indx_x_min + 1) * (indx_y_max - indx_y_min + 1)
     ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)
     pos1 = 1
     pos2 = indx_x_max - indx_x_min + 1
     DO j = indx_y_min, indx_y_max
        rbufer(pos1:pos2) = loc_EX(indx_x_min:indx_x_max,j)
        pos1 = pos2 + 1
        pos2 = pos2 + indx_x_max - indx_x_min + 1
     END DO
     DO j = indx_y_min, indx_y_max
        rbufer(pos1:pos2) = loc_EY(indx_x_min:indx_x_max,j)
        pos1 = pos2 + 1
        pos2 = pos2 + indx_x_max - indx_x_min + 1
     END DO
     CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_master, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
    
     DEALLOCATE(loc_EX, STAT=ALLOC_ERR)
     DEALLOCATE(loc_EY, STAT=ALLOC_ERR)
     DEALLOCATE(rbufer, STAT=ALLOC_ERR)

! ready to receive complete field array from the master of the cluster

  END IF

! master of the cluster distributes complete field arrays

  bufsize = (c_indx_x_max - c_indx_x_min + 1) * (c_indx_y_max - c_indx_y_min + 1)
  call mpi_barrier(mpi_comm_world, ierr)
  CALL MPI_BCAST(EX, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr) 
  call mpi_barrier(mpi_comm_world, ierr)
  CALL MPI_BCAST(EY, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)

! each member (master and non-master) in each cluster accumulates fields for ions in the whoole cluster domain 
! (in this case there is no need for communications at all)
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        acc_EX(i,j) = acc_EX(i,j) + EX(i,j)
     END DO
  END DO
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        acc_EY(i,j) = acc_EY(i,j) + EY(i,j)
     END DO
  END DO
  
END SUBROUTINE CALCULATE_ELECTRIC_FIELD

!------------------------------------------------------
!
SUBROUTINE CLEAR_ACCUMULATED_FIELDS

  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INTEGER i, j

  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        acc_EX(i,j) = 0.0_8
     END DO
  END DO
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        acc_EY(i,j) = 0.0_8
     END DO
  END DO
 
END SUBROUTINE CLEAR_ACCUMULATED_FIELDS
