! ************************************************************
!! Octree (Barnes-Hut)
!
!> Objectives
!  This module generates an octree based in the positions of the
!  particles and evaluates the forces using the Barnes-Hut 
!  approximation with a quadrupole expansion.
!
!> Modified
!  2026.09.20
!
!> Created
!  2026.06.15
!
!> Author
!  oap
!
MODULE octree_mod
USE tipos
IMPLICIT NONE
PRIVATE
PUBLIC OctreeType

TYPE :: OctreeType
    ! max depth of the tree
    INTEGER :: max_depth = 32
    ! amplificator of the size of the quadtree-root
    REAL(pf) :: side_amplificator = 1.2_pf

    ! multipole
    INTEGER :: multipole

    ! about the levels
    INTEGER :: most_depth
    INTEGER, ALLOCATABLE :: counter_for_each_level(:)

    ! to save_txt
    INTEGER :: save_txt

    ! masses, positions and number of bodies (N)
    REAL(pf), ALLOCATABLE :: real_m(:), m(:), x(:), y(:), z(:), radii(:)
    INTEGER :: N

    ! if its allocated
    LOGICAL :: is_allocated = .FALSE.

    ! if the simulation has collisions
    LOGICAL :: collide

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
    REAL(pf), ALLOCATABLE :: ns_max_radius(:)

    REAL(pf), ALLOCATABLE :: ns_quad(:,:) ! quadrupole terms
    REAL(pf), ALLOCATABLE :: ns_oct(:,:)  ! octupole terms
CONTAINS
    PROCEDURE :: pre_init, init, clear
    PROCEDURE :: allocate_nodes, add_node, allocate_subnode, add_to_subnode, add
    PROCEDURE :: forces => evaluate_forces_over_p
    PROCEDURE :: detect_collisions
    PROCEDURE :: evaluate_multipole
END TYPE

CONTAINS

SUBROUTINE pre_init (self, m, em, radii, multipole, collide_par)
! this subroutine pre initialize the tree by allocating all the self variables and clearing
! if its necessary
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: m(:) ! masses vector
    REAL(pf), INTENT(IN) :: radii(:) ! radii vector
    LOGICAL,  INTENT(IN) :: em   ! equal masses boolean
    INTEGER,  INTENT(IN) :: multipole
    LOGICAL,  INTENT(IN), OPTIONAL :: collide_par
    INTEGER :: N

    ! collide parameter
    self % collide = .FALSE.
    IF (PRESENT(collide_par)) self % collide = collide_par
    
    ! if it isnt allocated, so is necessary to allocate the basic self variables
    IF (.NOT. self % is_allocated) THEN
        self % N = SIZE(m)
        ALLOCATE(self % real_m(self % N))
        ALLOCATE(self % m(self % N))

        IF (self % collide) ALLOCATE(self % radii(self % N))

        ! spatial variables
        ALLOCATE(self % x(self % N))
        ALLOCATE(self % y(self % N))
        ALLOCATE(self % z(self % N))

        ! about the depth
        ALLOCATE(self % counter_for_each_level(self % max_depth+1))

        ! multipoles
        self % multipole = multipole

        ! allocate the other vectors
        self % max_number_of_nodes = 8 * self % N
        self % number_of_nodes = 0
        CALL self % allocate_nodes()
    ENDIF

    CALL self % clear()

    ! saving the masses
    self % real_m = m
    self % m = m
    IF (em) self % m = 1.0_pf
    IF (self % collide) self % radii = radii
END SUBROUTINE

SUBROUTINE clear (self)
    CLASS(OctreeType), INTENT(INOUT) :: self
    
    ! starting or restarting
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
    self % ns_max_radius = 0.0_pf
    self % ns_type = 0
    self % ns_depth = 0
    self % ns_child = -1
    IF (self % multipole > 1) THEN
        self % ns_quad = 0.0_pf
    ENDIF
    IF (self % multipole > 4) THEN
        self % ns_oct = 0.0_pf
    ENDIF
END SUBROUTINE

SUBROUTINE init (self, m, qs, save_txt)
! this subroutine inits the tree by allocating the global vectors and adding each particle
! in a node. if its the case it saves the root information too.
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: m(:), qs(:,:)
    INTEGER, OPTIONAL :: save_txt
    REAL(pf) :: infos_root(4)
    INTEGER :: p, idx_root

    ! clear
    CALL self % clear()

    ! saving particles information
    self % x = qs(:,1)
    self % y = qs(:,2)
    self % z = qs(:,3)

    ! about the depth
    self % most_depth = 0
    self % counter_for_each_level = 0

    ! init the root
    infos_root = node_size_center(self%x, self%y, self%z)
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

    ! if wants to use multipole
    IF (self % multipole > 1) THEN
        CALL self % evaluate_multipole()
    ENDIF
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

    ! if the tree already exists, so is the case of reallocation
    IF (self % number_of_nodes > 0) THEN
        old_size = self % max_number_of_nodes
        self % max_number_of_nodes = 2 * old_size
        new_size = self % max_number_of_nodes
        
        ! allocate the temp vectors
        ALLOCATE(temp_real(old_size))
        ALLOCATE(temp_int(old_size))
        ALLOCATE(temp_int_2(8,old_size))

        ! now deallocate and reallocate
        temp_real = self % ns_cx
        DEALLOCATE(self % ns_cx)
        ALLOCATE(self % ns_cx(new_size))
        self % ns_cx(1:old_size) = temp_real

        temp_real = self % ns_cy
        DEALLOCATE(self % ns_cy)
        ALLOCATE(self % ns_cy(new_size))
        self % ns_cy(1:old_size) = temp_real

        temp_real = self % ns_cz
        DEALLOCATE(self % ns_cz)
        ALLOCATE(self % ns_cz(new_size))
        self % ns_cz(1:old_size) = temp_real

        temp_real = self % ns_halfside
        DEALLOCATE(self % ns_halfside)
        ALLOCATE(self % ns_halfside(new_size))
        self % ns_halfside(1:old_size) = temp_real

        temp_real = self % ns_L2
        DEALLOCATE(self % ns_L2)
        ALLOCATE(self % ns_L2(new_size))
        self % ns_L2(1:old_size) = temp_real

        temp_real = self % ns_mass
        DEALLOCATE(self % ns_mass)
        ALLOCATE(self % ns_mass(new_size))
        self % ns_mass(1:old_size) = temp_real

        temp_real = self % ns_qcm_x
        DEALLOCATE(self % ns_qcm_x)
        ALLOCATE(self % ns_qcm_x(new_size))
        self % ns_qcm_x(1:old_size) = temp_real

        temp_real = self % ns_qcm_y
        DEALLOCATE(self % ns_qcm_y)
        ALLOCATE(self % ns_qcm_y(new_size))
        self % ns_qcm_y(1:old_size) = temp_real

        temp_real = self % ns_qcm_z
        DEALLOCATE(self % ns_qcm_z)
        ALLOCATE(self % ns_qcm_z(new_size))
        self % ns_qcm_z(1:old_size) = temp_real

        temp_int = self % ns_particle
        DEALLOCATE(self % ns_particle)
        ALLOCATE(self % ns_particle(new_size))
        self % ns_particle(1:old_size) = temp_int

        temp_int = self % ns_type
        DEALLOCATE(self % ns_type)
        ALLOCATE(self % ns_type(new_size))
        self % ns_type(1:old_size) = temp_int
        
        temp_int_2 = self % ns_child
        DEALLOCATE(self % ns_child)
        ALLOCATE(self % ns_child(8,new_size))
        self % ns_child(:,1:old_size) = temp_int_2

        temp_real = self % ns_max_radius
        DEALLOCATE(self % ns_max_radius)
        ALLOCATE(self % ns_max_radius(new_size))
        self % ns_max_radius(1:old_size) = temp_real

        DEALLOCATE(temp_real, temp_int, temp_int_2)
    ELSE
        ALLOCATE(self % ns_cx(self % max_number_of_nodes))
        ALLOCATE(self % ns_cy(self % max_number_of_nodes))
        ALLOCATE(self % ns_cz(self % max_number_of_nodes))
        ALLOCATE(self % ns_halfside(self % max_number_of_nodes))
        ALLOCATE(self % ns_L2(self % max_number_of_nodes))
        ALLOCATE(self % ns_mass(self % max_number_of_nodes))
        ALLOCATE(self % ns_qcm_x(self % max_number_of_nodes))
        ALLOCATE(self % ns_qcm_y(self % max_number_of_nodes))
        ALLOCATE(self % ns_qcm_z(self % max_number_of_nodes))
        ALLOCATE(self % ns_particle(self % max_number_of_nodes))
        ALLOCATE(self % ns_type(self % max_number_of_nodes))
        ALLOCATE(self % ns_depth(self % max_number_of_nodes))
        ALLOCATE(self % ns_child(8, self % max_number_of_nodes))
        ALLOCATE(self % ns_max_radius(self % max_number_of_nodes))

        ! quadrupole
        IF (self % multipole > 1) THEN
            ALLOCATE(self % ns_quad(self % max_number_of_nodes, 6))
            self % ns_quad = 0.0_pf
        ENDIF

        ! octupole
        IF (self % multipole > 4) THEN
            ALLOCATE(self % ns_oct(self % max_number_of_nodes, 14))
            self % ns_oct = 0.0_pf
        ENDIF
    ENDIF

    self % is_allocated = .TRUE.
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
        CALL self % allocate_nodes()
        PRINT *, '[DEBUG] REALLOCATING!'
    END IF

    ! starts without children and being a leaf
    self % ns_child(:, idx) = -1
    self % ns_type(idx) = 1 ! is a leaf
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
    self % ns_max_radius(idx) = 0.0_pf
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
    REAL(pf) :: pm, px, py, pz, pr, old_mass
    
    ! get particle information
    pm = self % real_m(p)
    px = self % x(p)
    py = self % y(p)
    pz = self % z(p)
    IF (self % collide) pr = self % radii(p)

    ! an empty node become a particle
    IF (self % ns_particle(node_idx) == -1 .AND. self % ns_type(node_idx) == 1) THEN
        self % ns_particle(node_idx) = p

        ! add directly the mass and center of mass
        self % ns_mass(node_idx) = pm
        self % ns_qcm_x(node_idx) = px
        self % ns_qcm_y(node_idx) = py
        self % ns_qcm_z(node_idx) = pz
        IF (self % collide) self % ns_max_radius(node_idx) = pr

        RETURN
    ENDIF

    ! update the mass and center of mass if the node isnt empty
    old_mass = self % ns_mass(node_idx)
    self % ns_mass(node_idx) = self % ns_mass(node_idx) + pm
    self % ns_qcm_x(node_idx) = (self % ns_qcm_x(node_idx) * old_mass + px * pm) / self % ns_mass(node_idx)
    self % ns_qcm_y(node_idx) = (self % ns_qcm_y(node_idx) * old_mass + py * pm) / self % ns_mass(node_idx)
    self % ns_qcm_z(node_idx) = (self % ns_qcm_z(node_idx) * old_mass + pz * pm) / self % ns_mass(node_idx)

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

        ! update the max radius
        IF (self % collide) THEN
            IF (self % ns_max_radius(node_idx) < pr) self % ns_max_radius(node_idx) = pr
        ENDIF
    ENDIF

    ! now add the new particle if its not at the deepest level subdividing it
    IF (self % ns_depth(node_idx) < self % max_depth) THEN
        CALL self % add_to_subnode(node_idx, p)
    ENDIF
END SUBROUTINE

SUBROUTINE evaluate_multipole (self)
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER :: i, child_idx, p, node_idx
    REAL(pf) :: pm, px, py, pz
    REAL(pf) :: dxi, dyi, dzi

    INTEGER :: level, counter, d, d_idx
    REAL(pf) :: x_sd, y_sd, z_sd

    INTEGER :: queue_current(self % number_of_nodes), queue_next(self % number_of_nodes)
    INTEGER :: kqc, keqc, keqn, remaining

    ! multipole state vectors
    self % ns_quad = 0.0_pf
    self % ns_oct  = 0.0_pf

    ! if the max depth was reached, we start at the deepest level
    IF (self % most_depth == self % max_depth) THEN
        level = self % most_depth

    ! if no, we can start at the almost depth because the lower only will have leafs
    ELSE
        level = self % most_depth - 1
    ENDIF
    counter = 0
    
    keqc = self % number_of_nodes ! key end queue current
    keqn = 0 ! key end queue next
    kqc = 0  ! key queue current

    ! start by the last node
    queue_current = [(self % number_of_nodes - i + 1, i=1, self%number_of_nodes)]

    DO WHILE (level >= 0)
        
        kqc = kqc + 1
        node_idx = queue_current(kqc)

        ! if isnt in the level, get the next
        IF (self % ns_depth(node_idx) < level) THEN
            keqn = keqn + 1
            queue_next(keqn) = node_idx
            CYCLE

        ! if its in the level, evaluate
        ELSE IF (self % ns_depth(node_idx) == level) THEN
            counter = counter + 1
            IF (counter == self % counter_for_each_level(level+1)) THEN
                level = level - 1
                counter = 0
                remaining = keqc - kqc

                IF (level >= 0) THEN
                    IF (kqc < keqc) THEN
                        queue_current(1:remaining) = queue_current(kqc+1:keqc)
                        queue_current(remaining+1:remaining+keqn) = queue_next(1:keqn)
                        keqc = remaining + keqn
                    ELSE
                        queue_current(1:keqn) = queue_next(1:keqn)
                        keqc = keqn
                    ENDIF
                ENDIF

                keqn = 0
                kqc = 0
            ENDIF

            ! if its a leaf, it doesnt have contributions
            IF (self % ns_type(node_idx) == 1) THEN
                CYCLE
            ENDIF

            ! if its a twig but it is in the deepest level, it doesnt have contributions also
            IF (self % ns_depth(node_idx) == self % max_depth) THEN
                CYCLE
            ENDIF

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
                    self % ns_quad(node_idx,:) = self % ns_quad(node_idx,:) + self % ns_quad(d_idx,:)
                    IF (self % multipole > 4) THEN
                        self % ns_oct(node_idx,:) = self % ns_oct(node_idx,:) + self % ns_oct(d_idx,:)
                    ENDIF
                ENDIF

                self % ns_quad(node_idx, 1) = self % ns_quad(node_idx, 1) + pm * dxi**2    ! mxi2
                self % ns_quad(node_idx, 2) = self % ns_quad(node_idx, 2) + pm * dyi**2    ! myi2
                self % ns_quad(node_idx, 3) = self % ns_quad(node_idx, 3) + pm * dzi**2    ! mzi2
                self % ns_quad(node_idx, 4) = self % ns_quad(node_idx, 4) + pm * dxi * dyi ! mxyi
                self % ns_quad(node_idx, 5) = self % ns_quad(node_idx, 5) + pm * dxi * dzi ! mxzi
                self % ns_quad(node_idx, 6) = self % ns_quad(node_idx, 6) + pm * dyi * dzi ! myzi

                IF (self % multipole > 4) THEN
                    self % ns_oct(node_idx, 1) = self % ns_oct(node_idx, 1) + pm * dxi * dxi**2
                    self % ns_oct(node_idx, 2) = self % ns_oct(node_idx, 2) + pm * dxi * dyi**2
                    self % ns_oct(node_idx, 3) = self % ns_oct(node_idx, 3) + pm * dxi * dzi**2

                    self % ns_oct(node_idx, 4) = self % ns_oct(node_idx, 4) + pm * dyi * dxi**2
                    self % ns_oct(node_idx, 5) = self % ns_oct(node_idx, 5) + pm * dyi * dyi**2
                    self % ns_oct(node_idx, 6) = self % ns_oct(node_idx, 6) + pm * dyi * dzi**2
                    
                    self % ns_oct(node_idx, 7) = self % ns_oct(node_idx, 7) + pm * dzi * dxi**2
                    self % ns_oct(node_idx, 8) = self % ns_oct(node_idx, 8) + pm * dzi * dyi**2
                    self % ns_oct(node_idx, 9) = self % ns_oct(node_idx, 9) + pm * dzi * dzi**2
                    
                    self % ns_oct(node_idx, 10) = self % ns_oct(node_idx, 10) + pm * dxi * dyi * dzi
                ENDIF

                IF (self % multipole > 4 .AND. self % ns_type(d_idx) == 2) THEN
                    self % ns_oct(node_idx, 1) = self % ns_oct(node_idx, 1) + &
                        3.0_pf * dxi * self % ns_quad(d_idx, 1)
                    self % ns_oct(node_idx, 2) = self % ns_oct(node_idx, 2) + &
                        2.0_pf * dyi * self % ns_quad(d_idx, 4) + dxi * self % ns_quad(d_idx, 2)
                    self % ns_oct(node_idx, 3) = self % ns_oct(node_idx, 3) + &
                        2.0_pf * dzi * self % ns_quad(d_idx, 5) + dxi * self % ns_quad(d_idx, 3)

                    self % ns_oct(node_idx, 4) = self % ns_oct(node_idx, 4) + &
                        2.0_pf * dxi * self % ns_quad(d_idx, 4) + dyi * self % ns_quad(d_idx, 1)
                    self % ns_oct(node_idx, 5) = self % ns_oct(node_idx, 5) + &
                        3.0_pf * dyi * self % ns_quad(d_idx, 2)
                    self % ns_oct(node_idx, 6) = self % ns_oct(node_idx, 6) + &
                        2.0_pf * dzi * self % ns_quad(d_idx, 6) + dyi * self % ns_quad(d_idx, 3)

                    self % ns_oct(node_idx, 7) = self % ns_oct(node_idx, 7) + &
                        2.0_pf * dxi * self % ns_quad(d_idx, 5) + dzi * self % ns_quad(d_idx, 1)
                    self % ns_oct(node_idx, 8) = self % ns_oct(node_idx, 8) + &
                        2.0_pf * dyi * self % ns_quad(d_idx, 6) + dzi * self % ns_quad(d_idx, 2)
                    self % ns_oct(node_idx, 9) = self % ns_oct(node_idx, 9) + &
                        3.0_pf * dzi * self % ns_quad(d_idx, 3)

                    self % ns_oct(node_idx, 10) = self % ns_oct(node_idx, 10) + &
                        dxi * self % ns_quad(d_idx, 6) + &
                        dyi * self % ns_quad(d_idx, 5) + &
                        dzi * self % ns_quad(d_idx, 4)
                ENDIF
            END DO

            IF (self % multipole > 4) THEN
                self % ns_oct(node_idx, 11) = SUM(self % ns_oct(node_idx, 1:10))
                self % ns_oct(node_idx, 12) = SUM(self % ns_oct(node_idx, 1:3))
                self % ns_oct(node_idx, 13) = SUM(self % ns_oct(node_idx, 4:6))
                self % ns_oct(node_idx, 14) = SUM(self % ns_oct(node_idx, 7:9))
            ENDIF
        ENDIF

    END DO
END SUBROUTINE

SUBROUTINE evaluate_forces_over_p (self, p, par_theta2, par_G, par_eps2, forces)
! given a particle and the parameters, this evaluates the forces over the particle
! using the Barnes-Hut criterion and optionally a multipole expansion (quadrupole or octupole).
! it uses the depth first traversal (DFS) algorithm to evaluate the forces in the tree.
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: p ! particle index
    REAL(pf), INTENT(IN), OPTIONAL :: par_theta2, par_G, par_eps2 ! parameters
    REAL(pf) :: eps2, theta2, G

    REAL(pf) :: pm, px, py, pz
    REAL(pf), INTENT(OUT) :: forces(3)
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
    IF (PRESENT(par_eps2))    eps2 = par_eps2
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
                
                mxi2 = self % ns_quad(node_idx, 1)
                myi2 = self % ns_quad(node_idx, 2)
                mzi2 = self % ns_quad(node_idx, 3)
                mxyi = self % ns_quad(node_idx, 4)
                mxzi = self % ns_quad(node_idx, 5)
                myzi = self % ns_quad(node_idx, 6)

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
                    Ox = - 1.5_pf * invR5 * self%ns_oct(node_idx, 12)
                    Oy = - 1.5_pf * invR5 * self%ns_oct(node_idx, 13)
                    Oz = - 1.5_pf * invR5 * self%ns_oct(node_idx, 14)

                    Oc = 7.5_pf * invR7 * (dx * self % ns_oct(node_idx, 12) + &
                                        dy * self % ns_oct(node_idx, 13) + &
                                        dz * self % ns_oct(node_idx, 14))
                    Ox = Ox + dx * Oc
                    Oy = Oy + dy * Oc
                    Oz = Oz + dz * Oc

                    dx2 = dx*dx
                    dy2 = dy*dy
                    dz2 = dz*dz
                    Ox = Ox + 7.5_pf * invR7 * (dx2 * self % ns_oct(node_idx, 1) + &
                                                dy2 * self % ns_oct(node_idx, 2) + &
                                                dz2 * self % ns_oct(node_idx, 3) + &
                                                2.0_pf * dx*dy * self % ns_oct(node_idx, 4) + &
                                                2.0_pf * dx*dz * self % ns_oct(node_idx, 7) + &
                                                2.0_pf * dy*dz * self % ns_oct(node_idx, 10))
                    Oy = Oy + 7.5_pf * invR7 * (dx2 * self % ns_oct(node_idx, 4) + &
                                                dy2 * self % ns_oct(node_idx, 5) + &
                                                dz2 * self % ns_oct(node_idx, 6) + &
                                                2.0_pf * dx*dy * self % ns_oct(node_idx, 2) + &
                                                2.0_pf * dx*dz * self % ns_oct(node_idx, 10) + &
                                                2.0_pf * dy*dz * self % ns_oct(node_idx, 8))
                    Oz = Oz + 7.5_pf * invR7 * (dx2 * self % ns_oct(node_idx, 7) + &
                                                dy2 * self % ns_oct(node_idx, 8) + &
                                                dz2 * self % ns_oct(node_idx, 9) + &
                                                2.0_pf * dx*dy * self % ns_oct(node_idx, 10) + &
                                                2.0_pf * dx*dz * self % ns_oct(node_idx, 3) + &
                                                2.0_pf * dy*dz * self % ns_oct(node_idx, 6))

                    Oc = -17.5_pf * invR9 * ( &
                        dx2 * dx * self % ns_oct(node_idx,1) + &
                        dy2 * dy * self % ns_oct(node_idx,5) + &
                        dz2 * dz * self % ns_oct(node_idx,9) + &
                        3.0_pf * dx2 * dy * self % ns_oct(node_idx,4) + &
                        3.0_pf * dx2 * dz * self % ns_oct(node_idx,7) + &
                        3.0_pf * dy2 * dx * self % ns_oct(node_idx,2) + &
                        3.0_pf * dy2 * dz * self % ns_oct(node_idx,8) + &
                        3.0_pf * dz2 * dx * self % ns_oct(node_idx,3) + &
                        3.0_pf * dz2 * dy * self % ns_oct(node_idx,6) + &
                        6.0_pf * dx * dy * dz * self % ns_oct(node_idx,10))

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
                child_idx = self % ns_child(i, node_idx)
                IF (child_idx .NE. -1) THEN
                    top = top + 1
                    stack(top) = child_idx
                ENDIF
            END DO
        ENDIF
    END DO
END SUBROUTINE

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
END MODULE