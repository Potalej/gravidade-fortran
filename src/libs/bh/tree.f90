MODULE octree_mod
USE tipos
IMPLICIT NONE
PRIVATE
PUBLIC OctreeType

TYPE :: OctreeType
    ! max depth of the tree
    INTEGER :: max_depth = 32
    ! amplificator of the size of the quadtree-root
    REAL(pf) :: side_amplificator = 1.2d0

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
CONTAINS
    PROCEDURE :: pre_init, init, clear
    PROCEDURE :: allocate_nodes, add_node, allocate_subnode, add_to_subnode, add
    PROCEDURE :: forces => evaluate_forces_over_p
    PROCEDURE :: detect_collisions
END TYPE

CONTAINS

SUBROUTINE pre_init (self, m, em, radii, collide_par)
! this subroutine pre initialize the tree by allocating all the self variables and clearing
! if its necessary
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: m(:) ! masses vector
    REAL(pf), INTENT(IN) :: radii(:) ! radii vector
    LOGICAL,  INTENT(IN) :: em   ! equal masses boolean
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

    infos(1) = 0.5d0 * (xmin + xmax)
    infos(2) = 0.5d0 * (ymin + ymax)
    infos(3) = 0.5d0 * (zmin + zmax)
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

    self % ns_cx(idx) = cx
    self % ns_cy(idx) = cy
    self % ns_cz(idx) = cz
    self % ns_halfside(idx) = side / 2.0d0
    self % ns_L2(idx) = side * side
    self % ns_mass(idx) = 0.0d0
    self % ns_qcm_x(idx) = 0.0d0
    self % ns_qcm_y(idx) = 0.0d0
    self % ns_qcm_z(idx) = 0.0d0
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
    
    h = self % ns_halfside(node_idx)
    h_half = h / 2.0d0
    cx = self % ns_cx(node_idx)
    cy = self % ns_cy(node_idx)
    cz = self % ns_cz(node_idx)
    d = self % ns_depth(node_idx) + 1

    ! NNO
    IF (index == 1) THEN
        cx_sub = cx - h_half
        cy_sub = cy + h_half
        cz_sub = cz + h_half
    ! NNE
    ELSEIF (index == 2) THEN
        cx_sub = cx + h_half
        cy_sub = cy + h_half
        cz_sub = cz + h_half
    ! NSO
    ELSEIF (index == 3) THEN
        cx_sub = cx - h_half
        cy_sub = cy - h_half
        cz_sub = cz + h_half
    ! NSE
    ELSEIF (index == 4) THEN
        cx_sub = cx + h_half
        cy_sub = cy - h_half
        cz_sub = cz + h_half
    ! SNO
    ELSEIF (index == 5) THEN
        cx_sub = cx - h_half
        cy_sub = cy + h_half
        cz_sub = cz - h_half
    ! SNE
    ELSEIF (index == 6) THEN
        cx_sub = cx + h_half
        cy_sub = cy + h_half
        cz_sub = cz - h_half
    ! SSO
    ELSEIF (index == 7) THEN
        cx_sub = cx - h_half
        cy_sub = cy - h_half
        cz_sub = cz - h_half
    ! SSE
    ELSE
        cx_sub = cx + h_half
        cy_sub = cy - h_half
        cz_sub = cz - h_half
    ENDIF

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
    LOGICAL :: top, right, ztop

    top   = (self % y(p) > self % ns_cy(node_idx))
    right = (self % x(p) > self % ns_cx(node_idx))
    ztop  = (self % z(p) > self % ns_cz(node_idx))

    IF (top) THEN
        IF (right) THEN
            subnode_index = 2
        ELSE
            subnode_index = 1
        ENDIF
    ELSE
        IF (right) THEN
            subnode_index = 4
        ELSE
            subnode_index = 3
        ENDIF
    ENDIF

    ! if its in the z-south
    IF (.NOT. ztop) subnode_index = subnode_index + 4

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
        ! we cannot go beyond the depth limit
        IF (self % ns_depth(node_idx) >= self % max_depth) THEN
            PRINT *, "BIG PROBLEM !!! MAX DEPTH !!!"
            STOP 0
        ENDIF

        ! in this case, its now a twig
        self % ns_type(node_idx) = 2

        ! add the old particle as a particle per si
        old_p = self % ns_particle(node_idx)
        self % ns_particle(node_idx) = -1
        CALL self % add_to_subnode(node_idx, old_p)

        ! update the max radius
        IF (self % collide) THEN
            IF (self % ns_max_radius(node_idx) < pr) self % ns_max_radius(node_idx) = pr
        ENDIF
    ENDIF

    ! now add the new particle
    CALL self % add_to_subnode(node_idx, p)
END SUBROUTINE

SUBROUTINE evaluate_forces_over_p (self, p, par_theta2, par_G, par_eps, forces)
! given a particle and the parameters, this evaluates the forces over the particle
! using the Barnes-Hut criterion. it uses the depth first traversal (DFS) algorithm 
! to evaluate the forces in the tree.
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: p ! particle index
    REAL(pf), INTENT(IN), OPTIONAL :: par_theta2, par_G, par_eps ! parameters
    REAL(pf) :: eps, theta2, G

    REAL(pf) :: pm, px, py, pz
    REAL(pf), INTENT(OUT) :: forces(3)
    REAL(pf) :: dx, dy, dz, dist2, L2, f, rinv

    INTEGER :: stack(self % number_of_nodes)
    INTEGER :: top, node_idx, i, child_idx

    ! default values
    eps = 0.0d0
    theta2 = 0.0d0
    G = 1.0d0

    ! replace if present
    IF (PRESENT(par_eps))    eps = par_eps
    IF (PRESENT(par_theta2)) theta2 = par_theta2
    IF (PRESENT(par_G))      G = par_G

    ! particle cache info
    pm = self % m(p)
    px = self % x(p)
    py = self % y(p)
    pz = self % z(p)

    ! initialize
    forces = 0.0d0
    top = 1
    stack(top) = 1

    ! dfs iterative
    DO WHILE (top > 0)

        node_idx = stack(top)
        top = top - 1

        ! empty node
        IF (self % ns_mass(node_idx) == 0.0d0) CYCLE

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

        ! leaf or bh criterion
        IF (self % ns_type(node_idx) == 1 .OR. L2 < theta2 * dist2) THEN
            rinv = 1.0_pf / SQRT(dist2 + eps*eps)
            rinv = rinv * rinv * rinv
            f = G * pm * self % ns_mass(node_idx) * rinv
            forces(1) = forces(1) + f * dx
            forces(2) = forces(2) + f * dy
            forces(3) = forces(3) + f * dz
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

SUBROUTINE detect_collisions (self, p, collisions, colliders)
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: p
    INTEGER, INTENT(INOUT) :: collisions(self % N - 1), colliders
    
    REAL(pf) :: pm, px, py, pz, dx, dy, dz, dist2, rsum
    INTEGER :: stack(self % number_of_nodes)
    INTEGER :: top, node_idx, child_idx, q, i

    ! particle cache info
    pm = self % real_m(p)
    px = self % x(p)
    py = self % y(p)
    pz = self % z(p)

    ! initialize
    colliders = 0
    collisions = 0
    top = 1
    stack(top) = 1

    ! dfs iterative
    DO WHILE (top > 0)

        node_idx = stack(top)
        top = top - 1

        ! empty node
        IF (self % ns_mass(node_idx) == 0.0_pf) CYCLE

        ! if not intersects
        IF (.NOT. sphere_intersects_node(self, node_idx, px, py, pz, self%radii(p))) CYCLE

        ! leaf
        IF (self % ns_type(node_idx) == 1) THEN
            q = self % ns_particle(node_idx)

            IF (q == -1) CYCLE ! empty
            IF (q == p)  CYCLE ! same particle

            dx = self % x(q) - px
            dy = self % y(q) - py
            dz = self % z(q) - pz

            dist2 = dx*dx + dy*dy + dz*dz
            rsum = self%radii(p) + self%radii(q)

            IF (dist2 <= rsum*rsum) THEN
                colliders = colliders + 1
                collisions(colliders) = q
            ENDIF
        
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

PURE FUNCTION sphere_intersects_node(self, node_idx, px, py, pz, r) RESULT(hit)
    CLASS(OctreeType), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: node_idx
    REAL(pf), INTENT(IN) :: px, py, pz, r
    LOGICAL :: hit

    REAL(pf) :: dx, dy, dz
    REAL(pf) :: cx, cy, cz, h
    REAL(pf) :: dist2

    cx = self % ns_cx(node_idx)
    cy = self % ns_cy(node_idx)
    cz = self % ns_cz(node_idx)
    h  = self % ns_halfside(node_idx)

    dx = MAX(ABS(px - cx) - h, 0.0_pf)
    dy = MAX(ABS(py - cy) - h, 0.0_pf)
    dz = MAX(ABS(pz - cz) - h, 0.0_pf)

    dist2 = dx*dx + dy*dy + dz*dz

    hit = (dist2 <= r*r)
END FUNCTION
END MODULE