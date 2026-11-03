! ***********************************************************************
!! Octree (Collisions, Barnes-Hut, Dehnen and optmized BH)
!
!> Objectives
!  This module generates an octree based in the positions of the particles. 
!  The octree can be used to detect collisions and to evaluate the 
!  gravitational interaction between the particles using one of three
!  methods:
!  (1) the classical Barnes-Hut approach (Barnes & Hut, 1986);
!  (2) the Dehnen very-fast algorithm (Dehnen, 2000) and the Marcello's
!      version, conserving total angular momentum (Marcello, 2017);
!  (3) an optimized BH method, using shared interaction lists (Barnes, 1990).
!
!> Initializing the octree
!  The octree needs to be pre initiated
!
!       CALL tree % pre_init(m, method, &
!                             collide (OPT), radii (OPT), &
!                             group_size (OPT))
!
!  The first is the mass, the second needs to be one of the available
!  methods (see below). The collide and radii parameters are optional and
!  are relative to the collision detection. The group_size is the parameter
!  relative to the method (3).
!
!  To init the tree, just run
!
!       CALL tree % init(x, y, z)
!
!> Using the methods
!  To use (1), pre init the octree using method = 10, 11 or 12, being,
!  respectively, the BH method with monopole, quadrupole and octupole.
!  To evaluate the forces, call the "forces" function for each particle:
!
!       DO p = 1, N
!           fs(:,p) = tree % forces(p, theta2, G, eps2)
!       END DO
!
!  This can be done in parallel (over this DO loop), and the output
!  of the "forces" function is an array of dimension 3.
!
!  To use (2), pre init the octree using method = 20, 21 or 22, being,
!  respectively, the Dehnen's method with monopole, with quadrupole
!  and using the Marcello (2017) correction to conserve angular
!  momentum. To evaluate the forces, call the "dehnen_eval" routine:
!
!       CALL tree % dehnen_eval(theta2, eps2, G, fs)
!
!  At the moment, this is not parallelizable, as long it does a dual
!  traversal and the parallelization of this is not imeddiate.
!
!  To use (3), pre-init the tree using method = 30, 31 or 32, being,
!  respectively, the BH method with monopole, quadrupole and octupole.
!  You also can specify the 'group_size' argument; the default is 64.
!  To evaluate the forces, call the "forces_group" routine:
!
!       CALL tree % forces_groups(theta2, G, eps2, fs, nt)
!
!  where 'nt' is the number of threads, so this method is invasive in
!  this sense.
!
!> References
!  * BARNES, Josh; HUT, Piet. A hierarchical O(N log N) force-calculation
!    algorithm. Nature, v. 324, n. 6096, p. 446–449, dez. 1986.
!  * DEHNEN, Walter. A Very Fast and Momentum-conserving Tree Code. The 
!    Astrophysical Journal, v. 536, n. 1, p. L39-L42, 10 jun. 2000.
!  * MARCELLO, Dominic C. A Very Fast and Angular Momentum Conserving 
!    Tree Code. The Astronomical Journal, v. 154, n. 3, p. 92, 11 ago. 2017.
!  * BARNES, Joshua E. A modified tree code: Don’t laugh; It runs. Journal 
!    of Computational Physics, v. 87, n. 1, p. 161–170, 1 mar. 1990. 
!
!> Modified
!  2026.11.03
!
!> Created
!  2026.06.15
!
!> Author
!  oap
!
MODULE octree_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC OctreeType

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

    TYPE :: OctreeType
        LOGICAL :: is_allocated = .FALSE.
        INTEGER :: method

        ! max depth of the tree
        INTEGER :: max_depth = 25
        ! amplificator of the size of the quadtree-root
        REAL(pf) :: side_amplificator = 1.2_pf

        ! timers
        REAL(pf) :: time_basic, time_multipole, time_group, time_construct

        ! multipole
        INTEGER :: multipole

        INTEGER :: most_depth
        INTEGER, ALLOCATABLE :: counter_for_each_level(:)

        ! to save_txt
        INTEGER :: save_txt

        ! masses, positions and number of bodies (N)
        REAL(pf), ALLOCATABLE :: m(:), x(:), y(:), z(:)
        INTEGER :: N

        ! nodes
        INTEGER :: number_of_nodes = 0
        INTEGER :: max_number_of_nodes
        REAL(pf), ALLOCATABLE :: ns_cx(:), ns_cy(:), ns_cz(:) ! centers
        REAL(pf), ALLOCATABLE :: ns_halfside(:), ns_L2(:)
        REAL(pf), ALLOCATABLE :: ns_mass(:), ns_qcm_x(:), ns_qcm_y(:), ns_qcm_z(:)
        INTEGER, ALLOCATABLE :: ns_particle(:)
        INTEGER, ALLOCATABLE :: ns_type(:)
        INTEGER, ALLOCATABLE :: ns_child(:,:)
        INTEGER, ALLOCATABLE :: ns_depth(:)

        REAL(pf), ALLOCATABLE :: ns_quad(:,:) ! quadrupole terms
        REAL(pf), ALLOCATABLE :: ns_oct(:,:)  ! octupole terms

        ! to detect collisions
        REAL(pf), ALLOCATABLE :: ns_max_radius(:)
        REAL(pf), ALLOCATABLE :: radii(:)

        ! for the Dehnen method
        REAL(pf), ALLOCATABLE :: ns_rmax(:), ns_force(:,:), ns_hess(:,:), ns_third(:,:)
        INTEGER :: ind_oct(3,3,3) ! indexes of the octupoles terms

        ! for the grouped Barnes-Hut method
        INTEGER :: group_size = 64
        INTEGER :: ngroups = 0
        INTEGER, ALLOCATABLE :: ns_count(:), ns_first(:), ns_last(:)
        INTEGER, ALLOCATABLE :: order(:), groups(:)
        REAL(pf), ALLOCATABLE :: gx(:), gy(:), gz(:), gm(:)
    CONTAINS
        PROCEDURE :: pre_init, init, clear
        PROCEDURE :: allocate_nodes, add_node, allocate_subnode, add_to_subnode, add
        PROCEDURE :: incorporate_particles

        ! for collisions
        PROCEDURE :: detect_collisions, sphere_intersects_node

        ! for tree methods
        PROCEDURE :: evaluate_multipole

        ! for the Barnes-Hut method
        PROCEDURE :: forces => bh_forces_over_p

        ! for the Dehnen method
        PROCEDURE :: dehnen_eval ! <- use this to eval the forces
        PROCEDURE :: dehnen_forces ! this and others are internal
        PROCEDURE :: evaluate_rmax
        PROCEDURE :: mutual_interaction_walk
        PROCEDURE :: evaluate_mutual_interaction_monopole
        PROCEDURE :: evaluate_mutual_interaction_quadrupole

        ! for the optmized (grouped) Barnes-Hut
        PROCEDURE :: prepare_groups, forces_groups
    END TYPE
CONTAINS

!*****************************************************************************************
! DEFAULT ROUTINES
! some routines that construct the octree, to use to detect collisions or to evalute the
! forces using the Barnes-Hut method or the Dehnen method.
!*****************************************************************************************
SUBROUTINE pre_init (self, m, method_par, collide_par, radii, group_size)
! this subroutine starts the octree and prepare it to use for some objective.
! the methods are:
! 10: barnes-hut (monopole)
! 11: barnes-hut (quadrupole)
! 12: barnes-hut (octupole)
! 20: dehnen (monopole)
! 21: dehnen (quadrupole)
! 30: grouped barnes-hut (monopole)
! 31: grouped barnes-hut (quadrupole)
! 32: grouped barnes-hut (octupole)
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: m(:) ! masses vector
    INTEGER,  INTENT(IN), OPTIONAL :: method_par
    LOGICAL,  INTENT(IN), OPTIONAL :: collide_par
    REAL(pf), INTENT(IN), OPTIONAL :: radii(:) ! radii vector
    INTEGER,  INTENT(IN), OPTIONAL :: group_size

    INTEGER :: method
    LOGICAL :: collide

    method = 0
    IF (PRESENT(method_par)) method = method_par
    self % method = method

    collide = .FALSE.
    IF (PRESENT(collide_par)) collide = collide_par

    IF (.NOT. self % is_allocated) THEN
        self % N = SIZE(m)
        ALLOCATE(self % m(self % N))

        ! spatial variables
        ALLOCATE(self % x(self % N))
        ALLOCATE(self % y(self % N))
        ALLOCATE(self % z(self % N))

        ! about the depth
        ALLOCATE(self % counter_for_each_level(self % max_depth+1))

        ! allocate the other vectors
        self % max_number_of_nodes = 3 * self % N
        self % number_of_nodes = 0
        CALL self % allocate_nodes()
        
        ! allocating specific vectors
        IF (method == 20 .OR. method == 21 .OR. method == 22) THEN
            ALLOCATE(self % ns_rmax(self % max_number_of_nodes))
            ALLOCATE(self % ns_force(3, self % max_number_of_nodes))

            ! indexes of the octupole moments, to use in the angular
            ! momentum conserving method (Marcello's method)
            self % ind_oct(1,1,:) = [1, 4, 7]
            self % ind_oct(1,2,:) = [4, 2, 10]
            self % ind_oct(1,3,:) = [7, 10, 3]
            self % ind_oct(2,1,:) = [4, 2, 10]
            self % ind_oct(2,2,:) = [2, 5, 8]
            self % ind_oct(2,3,:) = [10, 8, 6]
            self % ind_oct(3,1,:) = [7, 10, 3]
            self % ind_oct(3,2,:) = [10, 8, 6]
            self % ind_oct(3,3,:) = [3, 6, 9]
        ENDIF
        IF (method == 21 .OR. method == 22) THEN ! Dehnen with quadrupoles
            ALLOCATE(self % ns_hess(6, self % max_number_of_nodes))
            ALLOCATE(self % ns_third(10, self % max_number_of_nodes))
        ENDIF

        IF (collide) THEN
            ALLOCATE(self % ns_max_radius(self % max_number_of_nodes))
        ENDIF

        ! multipoles
        SELECT CASE (method)
            CASE (10, 20, 30) ! monopoles
                self % multipole = 1

            CASE (11, 21, 31) ! quadrupoles
                self % multipole = 4
                ALLOCATE(self % ns_quad(6, self % max_number_of_nodes))

            CASE (12, 22, 32) ! octupole
                self % multipole = 8
                ALLOCATE(self % ns_quad(6, self % max_number_of_nodes))
                ALLOCATE(self % ns_oct(14, self % max_number_of_nodes))
        END SELECT

        ! optimized Barnes-Hut
        IF (method == 30 .OR. method == 31 .OR. method == 32) THEN
            self % group_size = 64
            IF (PRESENT(group_size)) self % group_size = group_size
        ENDIF
    ENDIF

    CALL self % clear()

    ! saving the masses
    self % m = m
    IF (collide) self % radii = radii
    
    self % is_allocated = .TRUE.
END SUBROUTINE

SUBROUTINE clear (self)
    CLASS(OctreeType), INTENT(INOUT) :: self

    ! default vectors
    self % number_of_nodes = 0
    self % ns_cx = 0.0_pf
    self % ns_cy = 0.0_pf
    self % ns_cz = 0.0_pf
    self % ns_halfside = 0.0_pf
    self % ns_L2 = 0.0_pf
    self % ns_mass = 0.0_pf
    self % ns_qcm_x = 0.0_pf
    self % ns_qcm_y = 0.0_pf
    self % ns_qcm_z = 0.0_pf
    self % ns_particle = -1
    self % ns_type = 0
    self % ns_depth = 0
    self % most_depth = 0
    self % counter_for_each_level = 0
    
    ! collision
    IF (ALLOCATED(self % ns_max_radius)) self % ns_max_radius = 0.0_pf
    IF (self % multipole >= 4) self % ns_quad = 0.0_pf
    IF (self % multipole >= 8) self % ns_oct = 0.0_pf

    ! Dehnen method
    IF (self % method == 20 .OR. self % method == 21 .OR. self % method == 22) THEN
        self % ns_rmax = 0.0_pf
        self % ns_force = 0.0_pf
    ENDIF
    IF (self % method == 21 .OR. self % method == 22) THEN
        self % ns_hess = 0.0_pf
        self % ns_third = 0.0_pf
    ENDIF
END SUBROUTINE

SUBROUTINE init (self, x, y, z, save_txt)
! this subroutine inits the tree by allocating the global vectors and adding each particle
! in a node. if its the case it saves the root information too.
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: x(:), y(:), z(:)
    INTEGER,  INTENT(IN), OPTIONAL :: save_txt

    REAL(pf) :: infos_root(4)
    INTEGER  :: idx_root, p

    REAL(pf) :: t0, t1

    IF (.NOT. self % is_allocated) STOP "Not allocated! Run pre_init first!"

    CALL CPU_TIME(t0) ! timer
    ! initializing the general structure of the tree
        ! clear
        CALL self % clear()

        ! saving information
        self % x = x
        self % y = y
        self % z = z

        ! init the root
        infos_root = node_size_center(x, y, z)
        CALL self % add_node(infos_root(1), infos_root(2), infos_root(3), &
            self%side_amplificator*infos_root(4), 0, idx_root)

        ! if wants to save_txt
        self % save_txt = -1
        IF (PRESENT(save_txt)) THEN
            self % save_txt = save_txt
            WRITE (self % save_txt, *) self%ns_cx(1), self%ns_cy(1), self%ns_cz(1), self%ns_halfside(1)
        ENDIF
        
        ! add the children to root
        DO p = 1, self % N
            CALL self % add(idx_root, p)
        END DO

        ! add the physical information of the particles if use monopole.
        ! if use quadrupole or octupole, this is made in the evaluate_multipole
        IF (self % multipole == 1) CALL self % incorporate_particles()
    CALL CPU_TIME(t1)
    self % time_basic = t1 - t0

    ! if wants to use multipole
    self % time_multipole = 0.0_pf
    IF (self % multipole > 1) THEN
        CALL CPU_TIME(t0)
            CALL self % evaluate_multipole()
        CALL CPU_TIME(t1)
        self % time_multipole = t1 - t0
    ENDIF

    ! if using the optmized Barnes-Hut method, mount the groups
    self % time_group = 0.0_pf
    IF (self % method == 30 .OR. self % method == 31 .OR. self % method == 32) THEN
        CALL CPU_TIME(t0)
            CALL self % prepare_groups()
        CALL CPU_TIME(t1)
        self % time_group = t1 - t0
    ENDIF

    ! general timer
    self % time_construct = self % time_basic + self % time_multipole + self % time_group
END SUBROUTINE

SUBROUTINE allocate_nodes (self)
! this subroutine allocates the global vectors that stores informations about the nodes.
! if the tree already have nodes, so its a reallocation. in that case, the old information
! are stored in temp vectors and the vectors are reallocated with twice times the old size
! without lose the old information.
! @calledby init, add_node
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), ALLOCATABLE :: temp_real(:)
    INTEGER, ALLOCATABLE :: temp_int(:), temp_int_2(:,:)
    INTEGER :: old_size, new_size

    ! GENERAL VECTORS
    ! node positions
    ALLOCATE(self % ns_cx(self % max_number_of_nodes))
    ALLOCATE(self % ns_cy(self % max_number_of_nodes))
    ALLOCATE(self % ns_cz(self % max_number_of_nodes))
    ! node size
    ALLOCATE(self % ns_halfside(self % max_number_of_nodes))
    ALLOCATE(self % ns_L2(self % max_number_of_nodes))
    ! node mass and center of mass
    ALLOCATE(self % ns_mass(self % max_number_of_nodes))
    ALLOCATE(self % ns_qcm_x(self % max_number_of_nodes))
    ALLOCATE(self % ns_qcm_y(self % max_number_of_nodes))
    ALLOCATE(self % ns_qcm_z(self % max_number_of_nodes))
    ! node particle index
    ALLOCATE(self % ns_particle(self % max_number_of_nodes))
    ! node type, 1 - leaf, 2 - twig
    ALLOCATE(self % ns_type(self % max_number_of_nodes))
    ! node depth
    ALLOCATE(self % ns_depth(self % max_number_of_nodes))
    ! node children
    ALLOCATE(self % ns_child(8, self % max_number_of_nodes))

END SUBROUTINE

FUNCTION node_size_center (x, y, z) RESULT (infos)
! this subroutine gets the node size and center (in space xyz) based in the position
! vectors, giving the minimal square that contains every body.
! @calledby init
    REAL(pf), INTENT(IN) :: x(:), y(:), z(:)
    REAL(pf) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(pf) :: infos(4)

    xmin = MINVAL(x)
    xmax = MAXVAL(x)
    ymin = MINVAL(y)
    ymax = MAXVAL(y)
    zmin = MINVAL(z)
    zmax = MAXVAL(z)

    infos(1) = 0.5_pf * (xmin + xmax)
    infos(2) = 0.5_pf * (ymin + ymax)
    infos(3) = 0.5_pf * (zmin + zmax)
    infos(4) = MAXVAL((/ xmax - xmin, ymax - ymin, zmax - zmin /))
END FUNCTION

SUBROUTINE add_node (self, cx, cy, cz, side, depth, idx)
! this subroutine adds a node to the tree as a leaf. if its necessary, it reallocate the
! global node state vectors with twice the original size.
! @calledby init, allocate_subnode
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: cx, cy, cz, side
    INTEGER, INTENT(IN) :: depth
    INTEGER, INTENT(INOUT) :: idx

    self % number_of_nodes = self % number_of_nodes + 1
    idx = self % number_of_nodes
    
    IF (idx > self % max_number_of_nodes) THEN
        PRINT *, "INSUFFICIENT NUMBER OF NODES: ", self % max_number_of_nodes
        STOP 0
    END IF

    ! starts without children and being a leaf
    self % ns_child(:,idx) = -1
    self % ns_type(idx) = 1
    self % ns_depth(idx) = depth
    IF (depth > self % most_depth) self % most_depth = depth
    self % counter_for_each_level(depth+1) = self % counter_for_each_level(depth+1) + 1

    self % ns_cx(idx) = cx
    self % ns_cy(idx) = cy
    self % ns_cz(idx) = cz
    self % ns_halfside(idx) = side / 2.0_pf
    self % ns_L2(idx) = side * side
    self % ns_mass(idx) = 0.0_pf
    self % ns_qcm_x(idx) = 0.0_pf
    self % ns_qcm_y(idx) = 0.0_pf
    self % ns_qcm_z(idx) = 0.0_pf
    self % ns_particle(idx) = -1
END SUBROUTINE

SUBROUTINE allocate_subnode (self, node_idx, index)
! this subroutine allocates a subnode based in the index of the subnode wrt the node.
! the index is the default: 1-NNO, 2-NNE, 3-NSO, 4-NSE, 5-SNO, 6-SNE, 7-SSO, 8-SSE.
! @calledby add_to_subnode
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: index
    INTEGER :: subnode_idx, d
    REAL(pf) :: h, h_half, cx, cy, cz, cx_sub, cy_sub, cz_sub
    INTEGER :: sign_x, sign_y, sign_z
    
    h = self % ns_halfside(node_idx)
    h_half = h / 2.0_pf
    cx = self % ns_cx(node_idx)
    cy = self % ns_cy(node_idx)
    cz = self % ns_cz(node_idx)
    d = self % ns_depth(node_idx) + 1

    sign_x = MERGE(1,-1,BTEST(index-1,0))
    sign_y = MERGE(-1,1,BTEST(index-1,1))
    sign_z = MERGE(-1,1,BTEST(index-1,2))

    cx_sub = cx + sign_x*h_half
    cy_sub = cy + sign_y*h_half
    cz_sub = cz + sign_z*h_half

    ! create subnode
    CALL self % add_node(cx_sub, cy_sub, cz_sub, h, d, subnode_idx)
    ! allocate it as child of node
    self % ns_child(index, node_idx) = subnode_idx

    IF (self % save_txt .NE. -1) WRITE (self % save_txt, *) d, cx_sub, cy_sub, cz_sub
END SUBROUTINE

SUBROUTINE add_to_subnode (self, node_idx, p)
! given a node/twig and a particle, this adds the particle to the correct leaf based
! in its position.
! @calledby add
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: p
    INTEGER :: subnode
    INTEGER :: subnode_index

    INTEGER :: ix, iy, iz

    ix = MERGE(1,0,self % x(p) >= self % ns_cx(node_idx))
    iy = MERGE(1,0,self % y(p) < self % ns_cy(node_idx))
    iz = MERGE(1,0,self % z(p) < self % ns_cz(node_idx))

    subnode_index = 1 + ix + 2*iy + 4*iz
    ! CALL self % index_subnode(node_idx, p, subnode_index)

    subnode = self % ns_child(subnode_index, node_idx)

    ! if the subnode isnt associated
    IF (subnode == -1) THEN
        CALL self % allocate_subnode(node_idx, subnode_index)
    ENDIF

    subnode_index = self % ns_child(subnode_index, node_idx)

    ! now add
    CALL self % add(subnode_index, p)
END SUBROUTINE

SUBROUTINE add (self, node_idx, p)
! given a node and a particle, it adds the particle to the node as a particle if the
! node is empty and as a leaf if the node is already a leaf/particle.
! @calledby init
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: p ! particle index
    INTEGER :: old_p
    REAL(pf) :: pm, px, py, pz, old_mass
    
    ! an empty node become a particle
    IF (self % ns_particle(node_idx) == -1 .AND. self % ns_type(node_idx) == 1) THEN
        self % ns_particle(node_idx) = p
        RETURN
    ENDIF

    ! if isnt empty, it become a twig
    IF (self % ns_type(node_idx) == 1) THEN
        ! in this case, its now a twig
        self % ns_type(node_idx) = 2
        old_p = self % ns_particle(node_idx)
        self % ns_particle(node_idx) = -1

        ! if it is not at the deepest level, we can subdivide
        IF (self % ns_depth(node_idx) < self % max_depth) THEN
            CALL self % add_to_subnode(node_idx, old_p)
        ENDIF
    ENDIF

    ! now add the new particle if its not at the deepest level subdividing it
    IF (self % ns_depth(node_idx) < self % max_depth) THEN
        CALL self % add_to_subnode(node_idx, p)
    ENDIF
END SUBROUTINE

SUBROUTINE incorporate_particles (self)
! incorporate the physical information of the particles to the nodes. this is called 
! only if using monopole. in multipole>=4 cases, the evaluate_multipoles subroutine 
! do this work.
! @calledby init
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER :: node_idx, p, child_idx, child

    DO node_idx = self % number_of_nodes, 1, -1
        ! if its a leaf
        IF (self % ns_type(node_idx) == 1) THEN
            p = self % ns_particle(node_idx)
            IF (p == -1) CYCLE
            
            self % ns_mass(node_idx) = self % m(p)
            self % ns_qcm_x(node_idx) = self % x(p)
            self % ns_qcm_y(node_idx) = self % y(p)
            self % ns_qcm_z(node_idx) = self % z(p)
        
        ! if its a twig
        ELSE
            self % ns_mass(node_idx) = 0.0_pf
            self % ns_qcm_x(node_idx) = 0.0_pf
            self % ns_qcm_y(node_idx) = 0.0_pf
            self % ns_qcm_z(node_idx) = 0.0_pf

            DO child = 1, 8
                child_idx = self % ns_child(child, node_idx)
                IF (child_idx == -1) CYCLE

                self % ns_mass(node_idx) = self % ns_mass(node_idx) + &
                    self % ns_mass(child_idx)

                self % ns_qcm_x(node_idx) = self % ns_qcm_x(node_idx) + &
                    self % ns_mass(child_idx) * self % ns_qcm_x(child_idx)
                self % ns_qcm_y(node_idx) = self % ns_qcm_y(node_idx) + &
                    self % ns_mass(child_idx) * self % ns_qcm_y(child_idx)
                self % ns_qcm_z(node_idx) = self % ns_qcm_z(node_idx) + &
                    self % ns_mass(child_idx) * self % ns_qcm_z(child_idx)
            END DO

            self % ns_qcm_x(node_idx) = self % ns_qcm_x(node_idx) / self % ns_mass(node_idx)
            self % ns_qcm_y(node_idx) = self % ns_qcm_y(node_idx) / self % ns_mass(node_idx)
            self % ns_qcm_z(node_idx) = self % ns_qcm_z(node_idx) / self % ns_mass(node_idx)
        ENDIF
    END DO
END SUBROUTINE



!*****************************************************************************************
! COLLISIONS ROUTINES
! we can use the octree to detect collisions. for this the routine detect_collisions
! receives a particle index p and look for nodes or particles that intersects the sphere
! node. the list of nodes is traversed using the depth first traversal (DFS) algorithm.
!*****************************************************************************************
SUBROUTINE detect_collisions (self, p, collisions)
    CLASS(OctreeType), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: p
    INTEGER, INTENT(INOUT) :: collisions(:)
    
    REAL(pf) :: px, py, pz, dx, dy, dz, dist2, rsum
    REAL(pf) :: dpx, dpy, dpz
    INTEGER :: stack(self % number_of_nodes)
    INTEGER :: top, node_idx, child_idx, q, i, indice

    ! particle cache info
    px = self % x(p)
    py = self % y(p)
    pz = self % z(p)

    ! initialize
    ! to avoid recalculations, this was disabled
    ! colliders = 0
    ! collisions = 0
    top = 1
    stack(top) = 1

    ! dfs iterative
    DO WHILE (top > 0)

        node_idx = stack(top)
        top = top - 1

        ! empty node
        IF (self % ns_mass(node_idx) == 0.0_pf) CYCLE

        ! leaf
        IF (self % ns_type(node_idx) == 1) THEN
            q = self % ns_particle(node_idx)

            IF (q == -1) CYCLE ! empty
            IF (q == p)  CYCLE ! same particle

            IF (p > q) indice = (p-1)*(p-2)/2 + q
            IF (p < q) indice = (q-1)*(q-2)/2 + p
            IF (collisions(indice) .NE. 0) CYCLE
            collisions(indice) = -1

            dx = self % x(q) - px
            dy = self % y(q) - py
            dz = self % z(q) - pz

            dist2 = dx*dx + dy*dy + dz*dz
            rsum = self%radii(p) + self%radii(q)

            IF (dist2 <= rsum*rsum) THEN
                collisions(indice) = 1
            ENDIF

        ! if not intersects
        ELSE IF (.NOT. sphere_intersects_node(self, node_idx, px, py, pz, self%radii(p))) THEN
            CYCLE
        
        ! not leaf
        ELSE
            DO i = 1, 8
                child_idx = self % ns_child(i, node_idx)

                IF (child_idx .NE. -1) THEN
                    top = top + 1
                    stack(top) = child_idx
                ENDIF
            END DO
        ENDIF
    END DO
END SUBROUTINE

PURE FUNCTION sphere_intersects_node (self, node_idx, px, py, pz, r) RESULT(hit)
    CLASS(OctreeType), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: node_idx
    REAL(pf), INTENT(IN) :: px, py, pz, r
    LOGICAL :: hit

    REAL(pf) :: dx, dy, dz
    REAL(pf) :: cx, cy, cz, h
    REAL(pf) :: dist2
    REAL(pf) :: r2

    cx = self % ns_cx(node_idx)
    cy = self % ns_cy(node_idx)
    cz = self % ns_cz(node_idx)
    h  = self % ns_halfside(node_idx)

    dx = MAX(ABS(px - cx) - h, 0.0_pf)
    dy = MAX(ABS(py - cy) - h, 0.0_pf)
    dz = MAX(ABS(pz - cz) - h, 0.0_pf)

    dist2 = dx*dx + dy*dy + dz*dz

    r2 = r + self%ns_max_radius(node_idx)
    r2 = r2 * r2
    hit = (dist2 <= r2)
END FUNCTION



!*****************************************************************************************
! MULTIPOLE ROUTINE
! this subroutine evaluates the multipoles for the octree, and it can be up to quadrupoles
! or octupoles.
!*****************************************************************************************
SUBROUTINE evaluate_multipole (self)
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER :: i, child_idx, p, node_idx
    REAL(pf) :: pm, px, py, pz
    REAL(pf) :: dxi, dyi, dzi

    INTEGER  :: d, d_idx
    REAL(pf) :: x_sd, y_sd, z_sd

    ! multipole state vectors
    self % ns_quad = 0.0_pf
    IF (self % multipole > 4) self % ns_oct  = 0.0_pf

    DO node_idx = self % number_of_nodes, 1, -1
        ! if its a leaf, it doesnt have contributions
        IF (self % ns_type(node_idx) == 1) THEN
            p = self % ns_particle(node_idx)
            IF (p == -1) CYCLE
            
            ! add the physical information of the particle
            self % ns_mass(node_idx) = self % m(p)
            self % ns_qcm_x(node_idx) = self % x(p)
            self % ns_qcm_y(node_idx) = self % y(p)
            self % ns_qcm_z(node_idx) = self % z(p)
            CYCLE
        ENDIF

        ! if its a twig but it is in the deepest level, it doesnt have contributions also
        IF (self % ns_depth(node_idx) == self % max_depth) THEN
            CYCLE
        ENDIF

        ! adding the physical information (mass and com)
        self % ns_mass(node_idx) = 0.0_pf
        self % ns_qcm_x(node_idx) = 0.0_pf
        self % ns_qcm_y(node_idx) = 0.0_pf
        self % ns_qcm_z(node_idx) = 0.0_pf

        DO d = 1, 8
            d_idx = self % ns_child(d, node_idx)
            IF (d_idx == -1) CYCLE

            self % ns_mass(node_idx) = self % ns_mass(node_idx) + &
                self % ns_mass(d_idx)

            self % ns_qcm_x(node_idx) = self % ns_qcm_x(node_idx) + &
                self % ns_mass(d_idx) * self % ns_qcm_x(d_idx)
            self % ns_qcm_y(node_idx) = self % ns_qcm_y(node_idx) + &
                self % ns_mass(d_idx) * self % ns_qcm_y(d_idx)
            self % ns_qcm_z(node_idx) = self % ns_qcm_z(node_idx) + &
                self % ns_mass(d_idx) * self % ns_qcm_z(d_idx)
        END DO

        self % ns_qcm_x(node_idx) = self % ns_qcm_x(node_idx) / self % ns_mass(node_idx)
        self % ns_qcm_y(node_idx) = self % ns_qcm_y(node_idx) / self % ns_mass(node_idx)
        self % ns_qcm_z(node_idx) = self % ns_qcm_z(node_idx) / self % ns_mass(node_idx)

        ! if its a twig, we need to avaliate the daughters
        DO d = 1, 8
            d_idx = self % ns_child(d, node_idx)
            IF (d_idx == -1) CYCLE

            pm = self % ns_mass(d_idx)
            px = self % ns_qcm_x(d_idx)
            py = self % ns_qcm_y(d_idx)
            pz = self % ns_qcm_z(d_idx)

            dxi = px - self % ns_qcm_x(node_idx)
            dyi = py - self % ns_qcm_y(node_idx)
            dzi = pz - self % ns_qcm_z(node_idx)

            ! if its a twig, first add its contribution
            IF (self % ns_type(d_idx) == 2) THEN
                self % ns_quad(:, node_idx) = self % ns_quad(:,node_idx) + self % ns_quad(:,d_idx)
                IF (self % multipole > 4) THEN
                    self % ns_oct(:, node_idx) = self % ns_oct(:,node_idx) + self % ns_oct(:,d_idx)
                ENDIF
            ENDIF

            self % ns_quad(1,node_idx) = self % ns_quad(1,node_idx) + pm * dxi**2    ! mxi2
            self % ns_quad(2,node_idx) = self % ns_quad(2,node_idx) + pm * dyi**2    ! myi2
            self % ns_quad(3,node_idx) = self % ns_quad(3,node_idx) + pm * dzi**2    ! mzi2
            self % ns_quad(4,node_idx) = self % ns_quad(4,node_idx) + pm * dxi * dyi ! mxyi
            self % ns_quad(5,node_idx) = self % ns_quad(5,node_idx) + pm * dxi * dzi ! mxzi
            self % ns_quad(6,node_idx) = self % ns_quad(6,node_idx) + pm * dyi * dzi ! myzi

            IF (self % multipole > 4) THEN
                self % ns_oct(1,node_idx) = self % ns_oct(1,node_idx) + pm * dxi * dxi**2
                self % ns_oct(2,node_idx) = self % ns_oct(2,node_idx) + pm * dxi * dyi**2
                self % ns_oct(3,node_idx) = self % ns_oct(3,node_idx) + pm * dxi * dzi**2

                self % ns_oct(4,node_idx) = self % ns_oct(4,node_idx) + pm * dyi * dxi**2
                self % ns_oct(5,node_idx) = self % ns_oct(5,node_idx) + pm * dyi * dyi**2
                self % ns_oct(6,node_idx) = self % ns_oct(6,node_idx) + pm * dyi * dzi**2
                
                self % ns_oct(7,node_idx) = self % ns_oct(7,node_idx) + pm * dzi * dxi**2
                self % ns_oct(8,node_idx) = self % ns_oct(8,node_idx) + pm * dzi * dyi**2
                self % ns_oct(9,node_idx) = self % ns_oct(9,node_idx) + pm * dzi * dzi**2
                
                self % ns_oct(10,node_idx) = self % ns_oct(10,node_idx) + pm * dxi * dyi * dzi
            ENDIF

            IF (self % multipole > 4 .AND. self % ns_type(d_idx) == 2) THEN
                self % ns_oct(1, node_idx) = self % ns_oct(1, node_idx) + &
                    3.0_pf * dxi * self % ns_quad(1, d_idx)
                self % ns_oct(2, node_idx) = self % ns_oct(2, node_idx) + &
                    2.0_pf * dyi * self % ns_quad(4, d_idx) + dxi * self % ns_quad(2, d_idx)
                self % ns_oct(3, node_idx) = self % ns_oct(3, node_idx) + &
                    2.0_pf * dzi * self % ns_quad(5, d_idx) + dxi * self % ns_quad(3, d_idx)

                self % ns_oct(4, node_idx) = self % ns_oct(4, node_idx) + &
                    2.0_pf * dxi * self % ns_quad(4, d_idx) + dyi * self % ns_quad(1, d_idx)
                self % ns_oct(5, node_idx) = self % ns_oct(5, node_idx) + &
                    3.0_pf * dyi * self % ns_quad(2, d_idx)
                self % ns_oct(6, node_idx) = self % ns_oct(6, node_idx) + &
                    2.0_pf * dzi * self % ns_quad(6, d_idx) + dyi * self % ns_quad(3, d_idx)

                self % ns_oct(7, node_idx) = self % ns_oct(7, node_idx) + &
                    2.0_pf * dxi * self % ns_quad(5, d_idx) + dzi * self % ns_quad(1, d_idx)
                self % ns_oct(8, node_idx) = self % ns_oct(8, node_idx) + &
                    2.0_pf * dyi * self % ns_quad(6, d_idx) + dzi * self % ns_quad(2, d_idx)
                self % ns_oct(9, node_idx) = self % ns_oct(9, node_idx) + &
                    3.0_pf * dzi * self % ns_quad(3, d_idx)

                self % ns_oct(10, node_idx) = self % ns_oct(10, node_idx) + &
                    dxi * self % ns_quad(6, d_idx) + &
                    dyi * self % ns_quad(5, d_idx) + &
                    dzi * self % ns_quad(4, d_idx)
            ENDIF
        END DO

        IF (self % multipole > 4) THEN
            self % ns_oct(11, node_idx) = SUM(self % ns_oct(1:10, node_idx))
            self % ns_oct(12, node_idx) = SUM(self % ns_oct(1:3, node_idx))
            self % ns_oct(13, node_idx) = SUM(self % ns_oct(4:6, node_idx))
            self % ns_oct(14, node_idx) = SUM(self % ns_oct(7:9, node_idx))
        ENDIF
    END DO
END SUBROUTINE



!*****************************************************************************************
! BARNES-HUT ROUTINE
! this subroutine evaluates the forces over a particle of index p using the Barnes-Hut
! method, considering particle-cell interactions and therefore not preserving the Newton's
! 3rd law.
! given a particle and the parameters, this evaluates the forces over the particle
! using the Barnes-Hut criterion and optionally a multipole expansion (quadrupole or octupole).
! it uses the depth first traversal (DFS) algorithm to evaluate the forces in the tree.
!*****************************************************************************************
FUNCTION bh_forces_over_p (self, p, par_theta2, par_G, par_eps2) RESULT (forces)
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: p ! particle index
    REAL(pf), INTENT(IN), OPTIONAL :: par_theta2, par_G, par_eps2 ! parameters
    REAL(pf) :: eps2, theta2, G

    REAL(pf) :: pm, px, py, pz
    REAL(pf) :: forces(3)
    REAL(pf) :: dx, dy, dz, dist2, L2, f

    INTEGER :: stack(self % number_of_nodes)
    INTEGER :: top, node_idx, i, child_idx

    REAL(pf) :: rinv, invR, invR3, invR5
    REAL(pf) :: Mx, My, Mz
    REAL(pf) :: Qx, Qy, Qz
    REAL(pf) :: mxi2, myi2, mzi2, mxyi, mxzi, myzi

    REAL(pf) :: dx2, dy2, dz2, dxyz, invR2

    REAL(pf) :: invR7, invR9, Ox, Oy, Oz, Oc
    REAL(pf) :: Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz

    ! default values
    eps2 = 0.0_pf
    theta2 = 0.0_pf
    G = 1.0_pf

    ! replace if present
    IF (PRESENT(par_eps2))   eps2 = par_eps2
    IF (PRESENT(par_theta2)) theta2 = par_theta2
    IF (PRESENT(par_G))      G = par_G

    ! particle cache info
    pm = self % m(p)
    px = self % x(p)
    py = self % y(p)
    pz = self % z(p)

    ! initialize
    forces = 0.0_pf
    top = 1
    stack(top) = 1

    ! dfs iterative
    DO WHILE (top > 0)

        node_idx = stack(top)
        top = top - 1

        ! empty node
        IF (self % ns_mass(node_idx) == 0.0_pf) CYCLE

        ! self interaction
        IF (self % ns_type(node_idx) == 1) THEN
            IF (self % ns_particle(node_idx) == p) CYCLE
        ENDIF

        ! geometry
        dx = self % ns_qcm_x(node_idx) - px
        dy = self % ns_qcm_y(node_idx) - py
        dz = self % ns_qcm_z(node_idx) - pz
        dist2 = dx*dx + dy*dy + dz*dz
        L2 = self % ns_L2(node_idx)

        ! leaf (or twig at the deepest level) or bh criterion
        IF (self % ns_type(node_idx) == 1 & ! leaf
            .OR. L2 <= theta2 * dist2     & ! bh criterion
            .OR. self % ns_depth(node_idx) == self % max_depth & ! deepest level
        ) THEN
            IF (self % ns_type(node_idx) == 1 .OR. self % multipole == 1) THEN
                rinv = 1.0_pf / SQRT(dist2 + eps2)
                rinv = rinv * rinv * rinv
                f = G * pm * self % ns_mass(node_idx) * rinv
                forces(1) = forces(1) + f * dx
                forces(2) = forces(2) + f * dy
                forces(3) = forces(3) + f * dz
            ELSE
                invR = 1.0_pf / SQRT(dist2 + eps2)
                invR2 = 1.0_pf / (dist2 + eps2)
                invR3 = invR * invR2
                invR5 = invR3 * invR * invR
                invR7 = invR5 * invR * invR
                invR9 = invR7 * invR * invR

                Mx = dx * invR3 * self % ns_mass(node_idx)
                My = dy * invR3 * self % ns_mass(node_idx)
                Mz = dz * invR3 * self % ns_mass(node_idx)
                
                mxi2 = self % ns_quad(1, node_idx)
                myi2 = self % ns_quad(2, node_idx)
                mzi2 = self % ns_quad(3, node_idx)
                mxyi = self % ns_quad(4, node_idx)
                mxzi = self % ns_quad(5, node_idx)
                myzi = self % ns_quad(6, node_idx)

                dx2  = 15.0_pf * dx * dx * invR2 - 3.0_pf
                dy2  = 15.0_pf * dy * dy * invR2 - 3.0_pf
                dz2  = 15.0_pf * dz * dz * invR2 - 3.0_pf
                dxyz = 15.0_pf * dx * dy * dz * invR2
                
                Qx = (mxi2 * (dx2 - 6.0_pf) + myi2 * dy2 + mzi2 * dz2) * 0.5_pf * dx
                Qy = (myi2 * (dy2 - 6.0_pf) + mzi2 * dz2 + mxi2 * dx2) * 0.5_pf * dy
                Qz = (mzi2 * (dz2 - 6.0_pf) + mxi2 * dx2 + myi2 * dy2) * 0.5_pf * dz

                Qx = Qx + mxyi * dy * dx2 + mxzi * dz * dx2 + myzi * dxyz
                Qy = Qy + myzi * dz * dy2 + mxyi * dx * dy2 + mxzi * dxyz
                Qz = Qz + mxzi * dx * dz2 + myzi * dy * dz2 + mxyi * dxyz

                ! quadrupole
                IF (self % multipole == 4) THEN
                    forces(1) = forces(1) + G * pm * (Mx + invR5 * Qx)
                    forces(2) = forces(2) + G * pm * (My + invR5 * Qy)
                    forces(3) = forces(3) + G * pm * (Mz + invR5 * Qz)
                
                ! octupole
                ELSE
                    Ox = - 1.5_pf * invR5 * self%ns_oct(12, node_idx)
                    Oy = - 1.5_pf * invR5 * self%ns_oct(13, node_idx)
                    Oz = - 1.5_pf * invR5 * self%ns_oct(14, node_idx)

                    Oc = 7.5_pf * invR7 * (dx * self % ns_oct(12, node_idx) + &
                                        dy * self % ns_oct(13, node_idx) + &
                                        dz * self % ns_oct(14, node_idx))
                    Ox = Ox + dx * Oc
                    Oy = Oy + dy * Oc
                    Oz = Oz + dz * Oc

                    dx2 = dx*dx
                    dy2 = dy*dy
                    dz2 = dz*dz
                    Ox = Ox + 7.5_pf * invR7 * (dx2 * self % ns_oct(1, node_idx) + &
                                                dy2 * self % ns_oct(2, node_idx) + &
                                                dz2 * self % ns_oct(3, node_idx) + &
                                                2.0_pf * dx*dy * self % ns_oct(4, node_idx) + &
                                                2.0_pf * dx*dz * self % ns_oct(7, node_idx) + &
                                                2.0_pf * dy*dz * self % ns_oct(10, node_idx))
                    Oy = Oy + 7.5_pf * invR7 * (dx2 * self % ns_oct(4, node_idx) + &
                                                dy2 * self % ns_oct(5, node_idx) + &
                                                dz2 * self % ns_oct(6, node_idx) + &
                                                2.0_pf * dx*dy * self % ns_oct(2, node_idx) + &
                                                2.0_pf * dx*dz * self % ns_oct(10, node_idx) + &
                                                2.0_pf * dy*dz * self % ns_oct(8, node_idx))
                    Oz = Oz + 7.5_pf * invR7 * (dx2 * self % ns_oct(7, node_idx) + &
                                                dy2 * self % ns_oct(8, node_idx) + &
                                                dz2 * self % ns_oct(9, node_idx) + &
                                                2.0_pf * dx*dy * self % ns_oct(10, node_idx) + &
                                                2.0_pf * dx*dz * self % ns_oct(3, node_idx) + &
                                                2.0_pf * dy*dz * self % ns_oct(6, node_idx))

                    Oc = -17.5_pf * invR9 * ( &
                        dx2 * dx * self % ns_oct(1, node_idx) + &
                        dy2 * dy * self % ns_oct(5, node_idx) + &
                        dz2 * dz * self % ns_oct(9, node_idx) + &
                        3.0_pf * dx2 * dy * self % ns_oct(4, node_idx) + &
                        3.0_pf * dx2 * dz * self % ns_oct(7, node_idx) + &
                        3.0_pf * dy2 * dx * self % ns_oct(2, node_idx) + &
                        3.0_pf * dy2 * dz * self % ns_oct(8, node_idx) + &
                        3.0_pf * dz2 * dx * self % ns_oct(3, node_idx) + &
                        3.0_pf * dz2 * dy * self % ns_oct(6, node_idx) + &
                        6.0_pf * dx * dy * dz * self % ns_oct(10, node_idx))

                    Ox = Ox + dx * Oc
                    Oy = Oy + dy * Oc
                    Oz = Oz + dz * Oc

                    forces(1) = forces(1) + G * pm * (Mx + invR5 * Qx + Ox)
                    forces(2) = forces(2) + G * pm * (My + invR5 * Qy + Oy)
                    forces(3) = forces(3) + G * pm * (Mz + invR5 * Qz + Oz)
                ENDIF
            ENDIF
        ! if isnt leaf nor bh is valid, so go to children
        ELSE
            ! push children in the vec
            DO i = 1, 8
                child_idx = self % ns_child(I, node_idx)
                 IF (child_idx .NE. -1) THEN
                    top = top + 1
                    stack(top) = child_idx
                ENDIF
            END DO
        ENDIF
    END DO
END FUNCTION



!*****************************************************************************************
! DEHNEN ROUTINES
! the Dehnen algorithm consists in:
! (1) we evaluate the rmax for each (twig) node as being the maximum distance of a body
!     to the node com.
! (2) then we do a dual-traversal on the tree to evaluate the interactions between nodes.
!     this is done using a double DFS.
! (3) to evaluate the forces over the particles, we go over the tree starting by the root
!     and accumulating the forces on the childs. if the node is a leaf, so its a particle,
!     and as long the list of node indexes is ordered (the parent is always before the
!     child) all the forces on the particle was already evaluated, so we just apply it to
!     the particle.
!*****************************************************************************************
SUBROUTINE evaluate_rmax (self)
    CLASS(OctreeType), INTENT(INOUT) :: self
    
    REAL(pf) :: r_child_max
    INTEGER  :: node_idx, i, child_idx
    REAL(pf) :: dx, dy, dz, dist

    self % ns_rmax = 0.0_pf

    ! upward
    DO node_idx = self % number_of_nodes, 1, -1
        ! if its a leaf with (one) particle, its zero
        IF (self % ns_type(node_idx) == 1) THEN
            CYCLE

        ! if its an internal node, we combine the rmax of the childs
        ELSE
            r_child_max = 0.0_pf
            DO i = 1, 8
                child_idx = self % ns_child(i, node_idx)
                IF (child_idx == -1) CYCLE
                IF (self % ns_mass(child_idx) == 0.0_pf) CYCLE ! empty

                ! distance between the coms
                dx = self%ns_qcm_x(child_idx) - self%ns_qcm_x(node_idx)
                dy = self%ns_qcm_y(child_idx) - self%ns_qcm_y(node_idx)
                dz = self%ns_qcm_z(child_idx) - self%ns_qcm_z(node_idx)
                dist = SQRT(dx*dx + dy*dy + dz*dz)

                ! hmmmmm
                r_child_max = MAX(r_child_max, dist + self % ns_rmax(child_idx))
            END DO
            self % ns_rmax(node_idx) = r_child_max
        ENDIF
    END DO
END SUBROUTINE

SUBROUTINE mutual_interaction_walk (self, theta2, eps2, G)
! here we do a dual-traversal on the octree to determine which cells interact
! and how are the interactions
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: theta2, eps2, G

    INTEGER, ALLOCATABLE :: stack_a(:), stack_b(:)
    INTEGER :: top
    INTEGER :: a, b
    INTEGER :: i, j
    INTEGER :: child_a, child_b
    REAL(pf) :: dx, dy, dz, dist2, rsum2

    ! stacks of node pairs
    ALLOCATE(stack_a(36 * self % number_of_nodes))
    ALLOCATE(stack_b(36 * self % number_of_nodes))

    ! both starts at the root
    top = 1
    stack_a(top) = 1
    stack_b(top) = 1

    ! dual traversal
    DO WHILE (top > 0)

        a = stack_a(top)
        b = stack_b(top)
        top = top - 1

        ! empty nodes doesnt interact
        IF (self % ns_mass(a) == 0.0_pf) CYCLE
        IF (self % ns_mass(b) == 0.0_pf) CYCLE

        ! self-interaction
        IF (a == b) THEN
            ! a particle doesnt interact with itself
            IF (self % ns_type(a) == 1) CYCLE

            ! for a cell we interact the pairs
            DO i = 1, 8
                child_a = self % ns_child(i,a) 
                IF (child_a == -1) CYCLE
                DO j = i, 8
                    child_b = self % ns_child(j,a)
                    IF (child_b == -1) CYCLE
                    top = top + 1
                    stack_a(top) = child_a
                    stack_b(top) = child_b
                END DO
            END DO
            CYCLE
        ENDIF

        ! Dehnen MAC criterion
        dx = self % ns_qcm_x(a) - self % ns_qcm_x(b)
        dy = self % ns_qcm_y(a) - self % ns_qcm_y(b)
        dz = self % ns_qcm_z(a) - self % ns_qcm_z(b)
        dist2 = dx*dx + dy*dy + dz*dz
        rsum2 = (self % ns_rmax(a) + self % ns_rmax(b))**2

        ! well-separated
        IF (rsum2 <= theta2 * dist2) THEN
            ! we evaluate the coefficients of the expansion A <- B and B <- A
            IF (self % multipole == 1) THEN
                CALL evaluate_mutual_interaction_monopole(self, a, b, G, eps2)
            ELSE IF (self % multipole == 4 .OR. self % multipole == 8) THEN
                CALL evaluate_mutual_interaction_quadrupole(self, a, b, G, eps2)
            ELSE
                STOP "the multipole parameter needs to be 1, 4 or 8"
            ENDIF
            CYCLE
        ENDIF

        ! not well-separated, we subdivide the bigger node
        IF (self % ns_rmax(a) >= self % ns_rmax(b)) THEN
            ! opens A
            DO i = 1, 8
                child_a = self % ns_child(i,a)
                IF (child_a == -1) CYCLE
                top = top + 1
                stack_a(top) = child_a
                stack_b(top) = b
            END DO
        ELSE
            ! opens B
            DO i = 1, 8
                child_b = self % ns_child(i,b)
                IF (child_b == -1) CYCLE
                top = top + 1
                stack_a(top) = a
                stack_b(top) = child_b
            END DO
        ENDIF
    END DO
    DEALLOCATE(stack_a, stack_b)
END SUBROUTINE

SUBROUTINE evaluate_mutual_interaction_monopole (self, a, b, G, eps2)
! this evaluates the gravitational interaction between two nodes a and b using monopoles
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER,  INTENT(IN) :: a, b
    REAL(pf), INTENT(IN) :: G, eps2
    REAL(pf) :: dx, dy, dz, dist2, rinv3, factor

    ! vector point from A com to B com
    dx = self % ns_qcm_x(b) - self % ns_qcm_x(a)
    dy = self % ns_qcm_y(b) - self % ns_qcm_y(a)
    dz = self % ns_qcm_z(b) - self % ns_qcm_z(a)
    
    dist2 = dx*dx + dy*dy + dz*dz
    rinv3 = 1.0_pf / (dist2 + eps2)**1.5_pf
    factor = G * rinv3

    ! A <- B
    self % ns_force(1,a) = self % ns_force(1,a) + factor * self % ns_mass(b) * dx
    self % ns_force(2,a) = self % ns_force(2,a) + factor * self % ns_mass(b) * dy
    self % ns_force(3,a) = self % ns_force(3,a) + factor * self % ns_mass(b) * dz

    ! B <- A
    self % ns_force(1,b) = self % ns_force(1,b) - factor * self % ns_mass(a) * dx
    self % ns_force(2,b) = self % ns_force(2,b) - factor * self % ns_mass(a) * dy
    self % ns_force(3,b) = self % ns_force(3,b) - factor * self % ns_mass(a) * dz
END SUBROUTINE

SUBROUTINE evaluate_mutual_interaction_quadrupole (self, a, b, G, eps2)
! this evaluates the gravitational interaction between two nodes a and b using quadrupoles
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: a, b
    REAL(pf), INTENT(IN) :: G, eps2

    INTEGER :: i
    REAL(pf) :: Ma, Mb
    REAL(pf) :: dx, dy, dz
    REAL(pf) :: r2, rinv3, rinv5, rinv7
    REAL(pf) :: D1(3), D2(6), D3(10)
    REAL(pf) :: Qa(6), Qb(6)
    REAL(pf) :: CQa(3), CQb(3)

    REAL(pf) :: R(3), ci(3), CRR(3), Cjkl
    INTEGER  :: j, k, l
    REAL(pf) :: Fc(3)

    Ma = self%ns_mass(a)
    Mb = self%ns_mass(b)

    ! R = ZA - ZB
    dx = self%ns_qcm_x(a) - self%ns_qcm_x(b)
    dy = self%ns_qcm_y(a) - self%ns_qcm_y(b)
    dz = self%ns_qcm_z(a) - self%ns_qcm_z(b)

    r2 = dx*dx + dy*dy + dz*dz + eps2

    rinv3 = 1.0_pf / r2**1.5_pf
    rinv5 = 1.0_pf / r2**2.5_pf
    rinv7 = 1.0_pf / r2**3.5_pf

    ! considering g(r) = G / sqrt(r^2 + eps^2), Dehnen defines
    ! D^n = | (1/r d/dr)^n g(r) |_{r=|R|}

    ! D_i^(1) = R_i D^1 (eq. 7b)
    D1(1) = -G * dx * rinv3
    D1(2) = -G * dy * rinv3
    D1(3) = -G * dz * rinv3

    ! D_ij^(2) = delta_ij D^1 + R_i R_j D^2 (eq. 7c)
    D2(1) = G * (3.0_pf*dx*dx*rinv5 - rinv3) ! xx
    D2(2) = G * (3.0_pf*dy*dy*rinv5 - rinv3) ! yy
    D2(3) = G * (3.0_pf*dz*dz*rinv5 - rinv3) ! zz
    D2(4) = G * 3.0_pf * dx*dy*rinv5 ! xy
    D2(5) = G * 3.0_pf * dx*dz*rinv5 ! xz
    D2(6) = G * 3.0_pf * dy*dz*rinv5 ! yz

    ! D_ijk^(3) = (delta_ij R_k + delta_jk R_i + delta_ki R_j) D^2 + R_i R_j R_k D^3 (eq. 7d)
    D3(1)  = G * (-15.0_pf*dx*dx*dx*rinv7 + 9.0_pf*dx*rinv5) ! xxx
    D3(2)  = G * (-15.0_pf*dx*dx*dy*rinv7 + 3.0_pf*dy*rinv5) ! xxy
    D3(3)  = G * (-15.0_pf*dx*dx*dz*rinv7 + 3.0_pf*dz*rinv5) ! xxz
    D3(4)  = G * (-15.0_pf*dx*dy*dy*rinv7 + 3.0_pf*dx*rinv5) ! xyy
    D3(5)  = G * (-15.0_pf*dx*dy*dz*rinv7)                   ! xyz
    D3(6)  = G * (-15.0_pf*dx*dz*dz*rinv7 + 3.0_pf*dx*rinv5) ! xzz
    D3(7)  = G * (-15.0_pf*dy*dy*dy*rinv7 + 9.0_pf*dy*rinv5) ! yyy
    D3(8)  = G * (-15.0_pf*dy*dy*dz*rinv7 + 3.0_pf*dz*rinv5) ! yyz
    D3(9)  = G * (-15.0_pf*dy*dz*dz*rinv7 + 3.0_pf*dy*rinv5) ! yzz
    D3(10) = G * (-15.0_pf*dz*dz*dz*rinv7 + 9.0_pf*dz*rinv5) ! zzz

    ! specific quadrupoles of the nodes
    Qa(:) = self%ns_quad(:,a) / Ma
    Qb(:) = self%ns_quad(:,b) / Mb

    ! C_Q = Q_jk D_ijk^(3)
    ! C_Q,x
    CQa(1) = Qa(1)*D3(1) + 2.0_pf*Qa(4)*D3(2) + 2.0_pf*Qa(5)*D3(3) &
           + Qa(2)*D3(4) + 2.0_pf*Qa(6)*D3(5) + Qa(3)*D3(6)

    CQb(1) = Qb(1)*D3(1) + 2.0_pf*Qb(4)*D3(2) + 2.0_pf*Qb(5)*D3(3) &
           + Qb(2)*D3(4) + 2.0_pf*Qb(6)*D3(5) + Qb(3)*D3(6)

    ! C_Q,y
    CQa(2) = Qa(1)*D3(2) + 2.0_pf*Qa(4)*D3(4) + 2.0_pf*Qa(5)*D3(5) &
           + Qa(2)*D3(7) + 2.0_pf*Qa(6)*D3(8) + Qa(3)*D3(9)

    CQb(2) = Qb(1)*D3(2) + 2.0_pf*Qb(4)*D3(4) + 2.0_pf*Qb(5)*D3(5) &
           + Qb(2)*D3(7) + 2.0_pf*Qb(6)*D3(8) + Qb(3)*D3(9)

    ! C_Q,z
    CQa(3) = Qa(1)*D3(3) + 2.0_pf*Qa(4)*D3(5) + 2.0_pf*Qa(5)*D3(6) &
           + Qa(2)*D3(8) + 2.0_pf*Qa(6)*D3(9) + Qa(3)*D3(10)

    CQb(3) = Qb(1)*D3(3) + 2.0_pf*Qb(4)*D3(5) + 2.0_pf*Qb(5)*D3(6) &
           + Qb(2)*D3(8) + 2.0_pf*Qb(6)*D3(9) + Qb(3)*D3(10)


    ! coefficients of the local expansion in A
    ! dg_A (x) = Mb [D1 + 1/2 Qb D3 + D2 . x + 1/2 x . D3 x]
    self%ns_force(1,a) = self%ns_force(1,a) + Mb * (D1(1) + 0.5_pf*CQb(1))
    self%ns_force(2,a) = self%ns_force(2,a) + Mb * (D1(2) + 0.5_pf*CQb(2))
    self%ns_force(3,a) = self%ns_force(3,a) + Mb * (D1(3) + 0.5_pf*CQb(3))
    
    self%ns_hess(1:6,a) = self%ns_hess(1:6,a) + Mb * D2(:)
    self%ns_third(1:10,a) = self%ns_third(1:10,a) + Mb * D3(:)

    ! coefficients of the local expansion in B
    ! using that D1(-R) = -D1(R), D2(-R) =  D2(R), D3(-R) = -D3(R)
    self%ns_force(1,b) = self%ns_force(1,b) - Ma * (D1(1) + 0.5_pf*CQa(1))
    self%ns_force(2,b) = self%ns_force(2,b) - Ma * (D1(2) + 0.5_pf*CQa(2))
    self%ns_force(3,b) = self%ns_force(3,b) - Ma * (D1(3) + 0.5_pf*CQa(3))

    self%ns_hess(1:6,b) = self%ns_hess(1:6,b) + Ma * D2(:)
    self%ns_third(1:10,b) = self%ns_third(1:10,b) - Ma * D3(:)

    ! implementation of the Marcello's version, using the octupoles terms to
    ! vanish the torque (conserving angular momentum)
    IF (self % multipole >= 8) THEN
        R = [dx, dy, dz]
        ci = 0.0_pf
        CRR = 0.0_pf
        DO j = 1, 3
            DO k = 1, 3
                DO l = 1, 3
                    Cjkl = Mb * self % ns_oct(self % ind_oct(j,k,l), a) &
                        - Ma * self%ns_oct(self % ind_oct(j,k,l), b)
                    IF (k == l) ci(j) = ci(j) + Cjkl
                    CRR(j) = CRR(j) + Cjkl * R(k) * R(l)
                END DO
            END DO
        END DO

        ! F = 1/2 [ D2 c + D3 CRR ],  D2 = 3G/r^5,  D3 = -15G/r^7
        Fc = 0.5_pf * ( 3.0_pf*G*rinv5 * ci - 15.0_pf*G*rinv7 * CRR )

        self % ns_force(:,a) = self%ns_force(:,a) + Fc / Ma
        self % ns_force(:,b) = self%ns_force(:,b) - Fc / Mb
    ENDIF
END SUBROUTINE

SUBROUTINE dehnen_forces (self, forces)
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(OUT) :: forces(:,:) ! N x 3
    
    INTEGER :: node_idx, d, child_idx
    INTEGER :: p_num, p

    REAL(pf) :: f_child(3)
    REAL(pf) :: dx, dy, dz
    ! hessian terms
    REAL(pf) :: hxx, hyy, hzz, hxy, hxz, hyz
    ! third terms
    REAL(pf) :: txxx, txxy, txxz
    REAL(pf) :: txyy, txyz, txzz
    REAL(pf) :: tyyy, tyyz, tyzz
    REAL(pf) :: tzzz

    forces = 0.0_pf

    ! downward passage to propagate the forces from the root
    DO node_idx = 1, self % number_of_nodes
        ! a leaf dont propagate nothing
        IF (self % ns_type(node_idx) == 1) CYCLE
        
        IF (self % multipole >= 4) THEN
            ! hessian matrix
            hxx = self % ns_hess(1, node_idx)
            hyy = self % ns_hess(2, node_idx)
            hzz = self % ns_hess(3, node_idx)
            hxy = self % ns_hess(4, node_idx)
            hxz = self % ns_hess(5, node_idx)
            hyz = self % ns_hess(6, node_idx)
            ! third derivatives matrix
            txxx = self%ns_third(1, node_idx)
            txxy = self%ns_third(2, node_idx)
            txxz = self%ns_third(3, node_idx)
            txyy = self%ns_third(4, node_idx)
            txyz = self%ns_third(5, node_idx)
            txzz = self%ns_third(6, node_idx)
            tyyy = self%ns_third(7, node_idx)
            tyyz = self%ns_third(8, node_idx)
            tyzz = self%ns_third(9, node_idx)
            tzzz = self%ns_third(10, node_idx)
        ENDIF

        ! if its a twig, propagate the forces to its children
        DO d = 1, 8
            child_idx = self % ns_child(d, node_idx)
            IF (child_idx == -1) CYCLE
            IF (self % ns_mass(child_idx) == 0.0_pf) CYCLE

            ! forces over the child
            f_child = self % ns_force(:, child_idx)

            ! the parent g field is added to the child field
            ! monopole contribution
            f_child = f_child + self % ns_force(:, node_idx)

            IF (self % multipole >= 4) THEN
            ! distance between the child and the com of parent
                dx = self%ns_qcm_x(child_idx) - self%ns_qcm_x(node_idx)
                dy = self%ns_qcm_y(child_idx) - self%ns_qcm_y(node_idx)
                dz = self%ns_qcm_z(child_idx) - self%ns_qcm_z(node_idx)

            ! for the quadrupole, the contributions are given as
            ! a_child = monopole + H d + 0.5 d . T d + O(d^3)
                f_child(1) = f_child(1) &
                    + hxx * dx + hxy * dy + hxz * dz & ! hessian
                    + txxy*dx*dy + txxz*dx*dz + txyz*dy*dz &
                    + 0.5_pf * (txxx*dx*dx + txyy*dy*dy + txzz*dz*dz)
                
                f_child(2) = f_child(2) &
                    + hxy * dx + hyy * dy + hyz * dz & ! hessian
                    + txyy*dx*dy + txyz*dx*dz + tyyz*dy*dz &
                    + 0.5_pf * (txxy*dx*dx + tyyy*dy*dy + tyzz*dz*dz)

                f_child(3) = f_child(3) &
                    + hxz * dx + hyz * dy + hzz * dz & ! hessian
                    + txyz*dx*dy + txzz*dx*dz + tyzz*dy*dz &
                    + 0.5_pf * (txxz*dx*dx + tyyz*dy*dy + tzzz*dz*dz)
                    
            ! we also includes the contributions on the hessian and the third
            ! H_child = H_parent + Td + O(d^2)
                IF (self % ns_type(child_idx) /= 1) THEN
                    self % ns_hess(1, child_idx) = self % ns_hess(1, child_idx) &
                        + hxx + txxx*dx + txxy*dy + txxz*dz
                    self % ns_hess(2, child_idx) = self % ns_hess(2, child_idx) &
                        + hyy + txyy*dx + tyyy*dy + tyyz*dz
                    self % ns_hess(3, child_idx) = self % ns_hess(3, child_idx) &
                        + hzz + txzz*dx + tyzz*dy + tzzz*dz
                    self % ns_hess(4, child_idx) = self % ns_hess(4, child_idx) &
                        + hxy + txxy*dx + txyy*dy + txyz*dz
                    self % ns_hess(5, child_idx) = self % ns_hess(5, child_idx) &
                        + hxz + txxz*dx + txyz*dy + txzz*dz
                    self % ns_hess(6, child_idx) = self % ns_hess(6, child_idx) &
                        + hyz + txyz*dx + tyyz*dy + tyzz*dz
                ENDIF

            ! T_child = T_parent + O(d)
                self % ns_third(:, child_idx) = self % ns_third(:, child_idx) &
                    + self % ns_third(:, node_idx)
            ENDIF

            self % ns_force(:, child_idx) = f_child

            ! if its a particle, update the forces
            IF (self % ns_type(child_idx) == 1) THEN
                p = self % ns_particle(child_idx)
                IF (p == -1) CYCLE
                forces(:,p) = self % m(p) * self % ns_force(:, child_idx)
            ENDIF
        END DO
    END DO
END SUBROUTINE

SUBROUTINE dehnen_eval (self, theta2, potsoft2, G, forces)
! this subroutine do all the process to evaluate the Dehnen algorithm, assuming
! that the tree is already constructed
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN)  :: theta2, potsoft2, G
    REAL(pf), INTENT(OUT) :: forces(:,:)

    IF (self % multipole >= 4) THEN
        self % ns_hess = 0.0_pf
        self % ns_third = 0.0_pf
    ENDIF
    self % ns_force = 0.0_pf

    ! evaluating the max radius
    CALL self % evaluate_rmax()

    ! now we do the mutual interaction walk (dual traversal)
    CALL self % mutual_interaction_walk(theta2, potsoft2, G)

    ! and finally eval the forces
    CALL self % dehnen_forces(forces)
END SUBROUTINE



!*****************************************************************************************
! OPTMIZED BARNES-HUT ROUTINES
! The subroutines below are relative to the optimization suggested by Barnes (1990),
! using shared interaction lists. The results are expected to be a little different of
! the original implementation of Barnes and Hut (1986), maybe a little better.
!*****************************************************************************************
SUBROUTINE prepare_groups (self)
! this subroutine prepares the groups to evaluate the forces.
! 1. count particle per node
! 2. constructs an interval [first,last] of a vector 'order' for each node
! 3. copy the bodies in the tree order
! 4. choose the groups
    CLASS(OctreeType), INTENT(INOUT) :: self
    
    INTEGER :: node_idx, child_idx, i, k, p, top
    INTEGER :: stack(8 * (self % max_depth + 2))

    IF (.NOT. ALLOCATED(self % ns_count)) THEN
        ALLOCATE(self % ns_count(self % max_number_of_nodes))
        ALLOCATE(self % ns_first(self % max_number_of_nodes))
        ALLOCATE(self % ns_last(self % max_number_of_nodes))
        ALLOCATE(self % groups(self % max_number_of_nodes))
        ALLOCATE(self % order(self % N))
        ALLOCATE(self % gx(self % N), self % gy(self % N))
        ALLOCATE(self % gz(self % N), self % gm(self % N))
    ENDIF

    ! 1. count, from the children to the parent (childs have bigger index than parents)
    DO node_idx = self % number_of_nodes, 1, -1
        IF (self % ns_type(node_idx) == 1) THEN
            self % ns_count(node_idx) = MERGE(1, 0, self % ns_particle(node_idx) .NE. -1)
        ELSE
            k = 0
            DO i = 1, 8
                child_idx = self % ns_child(i, node_idx)
                IF (child_idx .NE. -1) k = k + self % ns_count(child_idx)
            END DO
            self % ns_count(node_idx) = k
        ENDIF
    END DO

    ! the root must have N bodies
    IF (self % ns_count(1) .NE. self % N) THEN
        STOP "[prepare_grous] body lost"
    ENDIF

    ! 2. intervals, from parent to children
    self % ns_first(1) = 1
    DO node_idx = 1, self % number_of_nodes
        k = self % ns_first(node_idx)
        IF (self % ns_type(node_idx) == 1) THEN
            IF (self % ns_count(node_idx) == 1) self % order(k) = self % ns_particle(node_idx)
        ELSE
            DO i = 1, 8
                child_idx = self % ns_child(i, node_idx)
                IF (child_idx == -1) CYCLE
                self % ns_first(child_idx) = k
                k = k + self % ns_count(child_idx)
            END DO
        ENDIF
        self % ns_last(node_idx) = self % ns_first(node_idx) + self % ns_count(node_idx) - 1
    END DO

    ! 3. bodies in the tree order (contiguous arrays)
    DO i = 1, self % N
        p = self % order(i)
        self % gm(i) = self % m(p)
        self % gx(i) = self % x(p)
        self % gy(i) = self % y(p)
        self % gz(i) = self % z(p)
    END DO

    ! 4. groups: DFS, stop in node with a number of particles <= group_size
    self % ngroups = 0
    top = 1
    stack(top) = 1
    DO WHILE (top > 0)
        node_idx = stack(top)
        top = top - 1
        IF (self % ns_count(node_idx) == 0) CYCLE

        IF (self % ns_count(node_idx) <= self % group_size .OR. self % ns_type(node_idx) == 1) THEN
            self % ngroups = self % ngroups + 1
            self % groups(self % ngroups) = node_idx
        ELSE
            DO i = 1, 8
                child_idx = self % ns_child(i, node_idx)
                IF (child_idx == -1) CYCLE
                top = top + 1
                stack(top) = child_idx
            END DO
        ENDIF
    END DO
END SUBROUTINE

SUBROUTINE forces_groups (self, theta2, G, eps2, forces)
! evaluate the forces with the BH method using interaction lists.
! its necessary to call self % prepare_groups() before
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN)  :: theta2, G, eps2
    REAL(pf), INTENT(OUT) :: forces(:,:) ! (3,N) original order

    IF (.NOT. ALLOCATED(self % ns_count) .OR. self % ngroups == 0) THEN
        CALL self % prepare_groups()
    ENDIF
    
    ! the lists of each threads stay in local variables of group_worker
    !$OMP PARALLEL DEFAULT(SHARED)
    CALL group_worker(self, theta2, G, eps2, forces)
    !$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE group_worker (self, theta2, G, eps2, forces)
! called inside a parallel region. the lists (locals) are privates of each
! thread and the !$OMP DO distribute the groups
    CLASS(OctreeType), INTENT(IN) :: self
    REAL(pf), INTENT(IN)  :: theta2, G, eps2
    REAL(pf), INTENT(INOUT) :: forces(:,:) ! (3,N) original order

    ! direct fonts
    REAL(pf), ALLOCATABLE :: sx(:), sy(:), sz(:), sm(:)
    ! cells (mass, com(3), quad(6), octu(14))
    REAL(pf), ALLOCATABLE :: cd(:,:)

    INTEGER  :: group_index, group_number, group_first, group_last
    REAL(pf) :: group_min(3), group_max(3)

    INTEGER :: top, stack(8 * (self % max_depth + 2))
    INTEGER :: node_idx, child_idx, i, j, k
    INTEGER :: a, b, ns, nc
    REAL(pf) :: cc(3), dist2, dd
    
    REAL(pf) :: xi, yi, zi, ax, ay, az, f
    REAL(pf) :: dx, dy, dz, invr, invr2, invr3, invr5, invr7, invr9
    REAL(pf) :: qx, qy, qz, dqd, trq, s
    REAL(pf) :: ox, oy, oz, oc

    ! allocating the direct fonts
    ALLOCATE(sx(2048), sy(2048), sz(2048), sm(2048))

    SELECT CASE (self % multipole)
        CASE (1)
            ! mass and com(3) = 4
            ALLOCATE(cd(4, 512))
        CASE (4)
            ! mass, com(3) and quad(6) = 10
            ALLOCATE(cd(10, 512))
        CASE (8)
            ! mass, com(3), quad(6) and octu(14) = 24
            ALLOCATE(cd(24, 512))
    END SELECT
    cd = 0.0_pf

    !$OMP DO SCHEDULE(DYNAMIC, 1)
    DO group_index = 1, self % ngroups
        group_number = self % groups(group_index)
        group_first  = self % ns_first(group_number)
        group_last   = self % ns_last(group_number)

        ! box that contains all the particles of the group
        group_min(1) = MINVAL(self % gx(group_first:group_last))
        group_min(2) = MINVAL(self % gy(group_first:group_last))
        group_min(3) = MINVAL(self % gz(group_first:group_last))
        group_max(1) = MAXVAL(self % gx(group_first:group_last))
        group_max(2) = MAXVAL(self % gy(group_first:group_last))
        group_max(3) = MAXVAL(self % gz(group_first:group_last))

        ! traversing the tree for all the group, making the lists
        ns = 0
        nc = 0
        top = 1
        stack(top) = 1
        DO WHILE (top > 0)
            node_idx = stack(top)
            top = top - 1

            IF (node_idx == group_number) CYCLE ! is the group itself
            IF (self % ns_mass(node_idx) == 0.0_pf) CYCLE ! empty node

            ! ancestor of the group: open ever (it should not be accepted by the MAC)
            IF (self % ns_first(node_idx) <= group_first .AND. self % ns_last(node_idx) >= group_last) THEN
                DO i = 1, 8
                    child_idx = self % ns_child(i, node_idx)
                    IF (child_idx == -1) CYCLE
                    top = top + 1
                    stack(top) = child_idx
                END DO
                CYCLE
            ENDIF

            ! leaf: direct sum (exact and less expensive than a cell)
            IF (self % ns_type(node_idx) == 1) THEN
                a = self % ns_first(node_idx)
                b = self % ns_last(node_idx)
                DO WHILE (ns + (b - a + 1) > SIZE(sx))
                    CALL grow_real_vec(sx); CALL grow_real_vec(sy)
                    CALL grow_real_vec(sz); CALL grow_real_vec(sm)
                END DO
                sm(ns+1:ns+b-a+1) = self % gm(a:b)
                sx(ns+1:ns+b-a+1) = self % gx(a:b)
                sy(ns+1:ns+b-a+1) = self % gy(a:b)
                sz(ns+1:ns+b-a+1) = self % gz(a:b)
                ns = ns + (b - a + 1)
                CYCLE
            ENDIF

            ! MAC: distance^2 from node com to the group of the group
            cc = [self % ns_qcm_x(node_idx), self % ns_qcm_y(node_idx), self % ns_qcm_z(node_idx)]
            dist2 = 0.0_pf
            DO i = 1, 3
                dd = MAX(group_min(i) - cc(i), 0.0_pf, cc(i) - group_max(i))
                dist2 = dist2 + dd * dd
            END DO

            IF (self % ns_L2(node_idx) <= theta2 * dist2 &
                .OR. self % ns_depth(node_idx) >= self % max_depth) THEN
                ! accepted cell
                nc = nc + 1
                IF (nc > SIZE(cd, 2)) CALL grow_real_mat(cd)
                ! saving the cell information
                cd(1, nc) = self % ns_mass(node_idx)
                cd(2:4, nc) = cc
                IF (self % multipole >= 4) cd(5:10, nc) = self % ns_quad(:, node_idx)
                IF (self % multipole >= 8) cd(11:24, nc) = self % ns_oct(:, node_idx)
            ELSE
                DO i = 1, 8
                    child_idx = self % ns_child(i, node_idx)
                    IF (child_idx == -1) CYCLE
                    top = top + 1
                    stack(top) = child_idx
                END DO
            ENDIF
        END DO

        ! each particle of the group go through the lists (contiguous loops)
        DO i = group_first, group_last
            xi = self % gx(i); yi = self % gy(i); zi = self % gz(i)
            ax = 0.0_pf; ay = 0.0_pf; az = 0.0_pf

            ! in the group itself (skip the particle itself)
            DO j = group_first, group_last
                IF (j == i) CYCLE
                dx = self % gx(j) - xi
                dy = self % gy(j) - yi
                dz = self % gz(j) - zi
                dist2 = dx*dx + dy*dy + dz*dz + eps2
                invr = 1.0_pf / SQRT(dist2)
                invr3 = invr * invr * invr
                f = self % gm(j) * invr3
                ax = ax + f * dx
                ay = ay + f * dy
                az = az + f * dz
            END DO

            ! leafs of the other regions
            !$OMP SIMD REDUCTION(+:ax,ay,az) PRIVATE(dx,dy,dz,dist2,invr,invr3,f)
            DO j = 1, ns
                dx = sx(j) - xi
                dy = sy(j) - yi
                dz = sz(j) - zi
                dist2 = dx*dx + dy*dy + dz*dz + eps2
                invr = 1.0_pf / SQRT(dist2)
                invr3 = invr * invr * invr
                f = sm(j) * invr3
                ax = ax + f * dx
                ay = ay + f * dy
                az = az + f * dz
            END DO

            ! accepts cells, with d = cm - pos
            IF (self % multipole == 1) THEN
                !$OMP SIMD REDUCTION(+:ax,ay,az) PRIVATE(dx,dy,dz,dist2,invr,invr3,f)
                DO j = 1, nc
                    dx = cd(2,j) - xi
                    dy = cd(3,j) - yi
                    dz = cd(4,j) - zi
                    dist2 = dx*dx + dy*dy + dz*dz + eps2
                    invr = 1.0_pf / SQRT(dist2)
                    invr3 = invr * invr * invr
                    f = cd(1,j) * invr3
                    ax = ax + f * dx
                    ay = ay + f * dy
                    az = az + f * dz
                END DO
            ELSE IF (self % multipole == 4) THEN
                !$OMP SIMD REDUCTION(+:ax,ay,az) &
                !$OMP& PRIVATE(dx,dy,dz,dist2,invr,invr3,invr5,invr7,qx,qy,qz,dqd,trq,s)
                DO j = 1, nc
                    dx = cd(2,j) - xi
                    dy = cd(3,j) - yi
                    dz = cd(4,j) - zi
                    dist2 = dx*dx + dy*dy + dz*dz + eps2
                    invr = 1.0_pf / SQRT(dist2)
                    invr3 = invr * invr * invr
                    invr5 = invr3 * invr * invr
                    invr7 = invr5 * invr * invr

                    qx = cd(5,j)*dx + cd(8,j)*dy + cd(9,j)*dz
                    qy = cd(8,j)*dx + cd(6,j)*dy + cd(10,j)*dz
                    qz = cd(9,j)*dx + cd(10,j)*dy + cd(7,j)*dz
                    dqd = dx*qx + dy*qy + dz*qz
                    trq = cd(5,j) + cd(6,j) + cd(7,j)

                    s = cd(1,j)*invr3 - 1.5_pf*trq*invr5 + 7.5_pf*dqd*invr7
                    ax = ax + s*dx - 3.0_pf*qx*invr5
                    ay = ay + s*dy - 3.0_pf*qy*invr5
                    az = az + s*dz - 3.0_pf*qz*invr5
                END DO  
            ELSE IF (self % multipole == 8) THEN
                !$OMP SIMD REDUCTION(+:ax,ay,az) &
                !$OMP& PRIVATE(dx,dy,dz,dist2,invr,invr3,invr5,invr7,invr9,qx,qy,qz,dqd,trq,ox,oy,oz,oc)
                DO j = 1, nc
                    dx = cd(2,j) - xi
                    dy = cd(3,j) - yi
                    dz = cd(4,j) - zi
                    dist2 = dx*dx + dy*dy + dz*dz + eps2
                    invr = 1.0_pf / SQRT(dist2)
                    invr3 = invr * invr * invr
                    invr5 = invr3 * invr * invr
                    invr7 = invr5 * invr * invr
                    invr9 = invr7 * invr * invr

                    qx = cd(5,j)*dx + cd(8,j)*dy + cd(9,j)*dz
                    qy = cd(8,j)*dx + cd(6,j)*dy + cd(10,j)*dz
                    qz = cd(9,j)*dx + cd(10,j)*dy + cd(7,j)*dz
                    dqd = dx*qx + dy*qy + dz*qz
                    trq = cd(5,j) + cd(6,j) + cd(7,j)
                    
                    ! monopole + quadrupole
                    s = cd(1,j)*invr3 - 1.5_pf*trq*invr5 + 7.5_pf*dqd*invr7
                    ax = ax + s*dx - 3.0_pf*qx*invr5
                    ay = ay + s*dy - 3.0_pf*qy*invr5
                    az = az + s*dz - 3.0_pf*qz*invr5

                    ! octupole
                    oc = 7.5_pf * invr7 * (dx*cd(22,j) + dy*cd(23,j) + dz*cd(24,j))
                    ox = -1.5_pf*invr5*cd(22,j) + dx*oc
                    oy = -1.5_pf*invr5*cd(23,j) + dy*oc
                    oz = -1.5_pf*invr5*cd(24,j) + dz*oc

                    ox = ox +  7.5_pf*invr7*(dx*dx*cd(11,j) + dy*dy*cd(12,j) + dz*dz*cd(13,j))
                    ox = ox + 15.0_pf*invr7*(dx*dy*cd(14,j) + dx*dz*cd(17,j) + dy*dz*cd(20,j))

                    oy = oy +  7.5_pf*invr7*(dx*dx*cd(14,j) + dy*dy*cd(15,j) + dz*dz*cd(16,j))
                    oy = oy + 15.0_pf*invr7*(dx*dy*cd(12,j) + dx*dz*cd(20,j) + dy*dz*cd(18,j))

                    oz = oz +  7.5_pf*invr7*(dx*dx*cd(17,j) + dy*dy*cd(18,j) + dz*dz*cd(19,j))
                    oz = oz + 15.0_pf*invr7*(dx*dy*cd(20,j) + dx*dz*cd(13,j) + dy*dz*cd(16,j))

                    oc = dx*dx*dx*cd(11,j) + dy*dy*dy*cd(15,j) + dz*dz*dz*cd(19,j)
                    oc = oc + 3.0_pf * (dx*dx*dy*cd(14,j) + dx*dx*dz*cd(17,j) + dx*dy*dy*cd(12,j))
                    oc = oc + 3.0_pf * (dy*dy*dz*cd(18,j) + dx*dz*dz*cd(13,j) + dy*dz*dz*cd(16,j))
                    oc = oc + 6.0_pf * dx * dy * dz * cd(20,j)
                    oc = -17.5_pf * invr9 * oc
                    
                    ox = ox + dx * oc
                    oy = oy + dy * oc
                    oz = oz + dz * oc

                    ax = ax + ox
                    ay = ay + oy
                    az = az + oz
                END DO
            ENDIF
            ! original order
            k = self % order(i)
            forces(1,k) = G * self % gm(i) * ax
            forces(2,k) = G * self % gm(i) * ay
            forces(3,k) = G * self % gm(i) * az
        END DO
    END DO
    !$OMP END DO

    DEALLOCATE(sx, sy, sz, sm, cd)
END SUBROUTINE

SUBROUTINE grow_real_vec (a)
! just grow a real vector
    REAL(pf), ALLOCATABLE, INTENT(INOUT) :: a(:)
    REAL(pf), ALLOCATABLE :: b(:)
    ALLOCATE(b(2 * SIZE(a)))
    b(1:SIZE(a)) = a
    CALL MOVE_ALLOC(b, a)
END SUBROUTINE

SUBROUTINE grow_real_mat (a)
! just grow a real matrix
    REAL(pf), ALLOCATABLE, INTENT(INOUT) :: a(:,:)
    REAL(pf), ALLOCATABLE :: b(:,:)
    ALLOCATE(b(SIZE(a,1), 2 * SIZE(a,2)))
    b(:, 1:SIZE(a,2)) = a
    CALL MOVE_ALLOC(b, a)
END SUBROUTINE

END MODULE 