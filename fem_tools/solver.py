def Solve(mesh_func, N, U_exact, F, boundary, size = False):


    # Set up computational domain
    E, Nodes, Mesh = mesh_func(N)

    # Identify boundary
    bn = boundary(Nodes)

    NN = N**2

    if size != False:
        NN = size


    M = Stiffness(E, Nodes)
    rhs = RHS(f = F, g = U_exact, elements = E, nodes = Nodes)

    # Extension of initial condition is just U itself.

    #Assemble full Stiffness matrix and then delete boundary data:
    A = np.zeros([NN,NN], dtype = "float")
    for e in E[:,0]:
        N_s = E[e,1::]
        M_loc = M.loc(e)
        for i in range(3):
            for j in range(3):
                A[N_s[i], N_s[j]] += M_loc[i,j]

    B = rhs.b(boundary = bn)

    # remove nodes corresponding to boundary data:

    # for the array it's a bit trickier
    # first define a list without the elements in bn
    nns = np.arange(0, len(Nodes[:,1]))

    not_bn = np.array([x for x in nns if x not in bn])

    #delete these rows and columns
    A_in = A[not_bn,:][:,not_bn]

    # Apply M to G
    if size == False:
        G_init = U_exact(Nodes[:,1], Nodes[:,2])
        u_D = G_init.copy()
    else:
        G_init = U_exact(Nodes[:,1], Nodes[:,2])
        u_D = G_init.copy()

    # Set G to zero except for the boundary.
    u_D[not_bn] = 0.
    AG = np.matmul(A[not_bn,:], u_D)

    B_in = B[not_bn] - AG


    # Finally invert
    #U_in = np.matmul(np.linalg.inv(A_in), B_in)
    U_in = np.linalg.solve(A_in, B_in)

    # Combine with boundary data:
    U  = np.zeros([NN,], dtype = "float")
    U[not_bn] = U_in
    U += u_D

    return U, G_init
