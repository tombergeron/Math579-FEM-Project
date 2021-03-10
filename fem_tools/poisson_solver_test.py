from square_triangulate import Rect_QT
from matrix_assembly import Poisson_Assemble

def U_exact_2(x,y):

    return  np.sin(3*x)*np.cos(y)


def F(x,y):

    return 10*np.sin(3*x)*np.cos(y)


def H1_norm(E, Nodes, U, which = None):

    #Assemble full Stiffness and Mass matrices
    N = len(Nodes[:,0])
    A = np.zeros([N,N], dtype = "float")
    M = Stiffness(E, Nodes)

    Mass = np.zeros([N,N], dtype = "float")
    Mass_loc = (1/12)*np.array([[2,1,1],[1,2,1],[1,1,2]])

    for e in E[:,0]:
        N_s = E[e,1::]
        [X,Y] = [Nodes[N_s, 1], Nodes[N_s,2]]
        J = abs((X[1] - X[0])*(Y[2]-Y[0]) - (X[2]-X[0])*(Y[1]-Y[0]))
        Area = 0.5*J
        mass_loc = Area*Mass_loc
        M_loc = M.loc(e)
        # Add entries to global stiffness and mass matrices
        for i in range(3):
            for j in range(3):
                A[N_s[i], N_s[j]] += M_loc[i,j]
                Mass[N_s[i], N_s[j]] += mass_loc[i,j]

    u = U.reshape([N,]).copy()

    if which == None:

        error = np.sqrt(np.matmul(u, np.matmul(Mass,u)) + np.matmul(u, np.matmul(A,u)))

        return error

    if which == "L2":
        error = np.sqrt(np.matmul(u, np.matmul(Mass, u)))

        return error


# Square Domain
H1_error = []
L2_error = []
L_inf = []


#initialize triangulation
S = Rect_QT()
lvls = [1,2,3,4,5]



for k in lvls:
    S.add_layer()
    P = Poisson_Assemble(S)
    A = P.ass_mat()





# for N in Ns[:-1]:
#     E, Nodes, Mesh = Square_Triangulate(N)
#     M = Stiffness(E, Nodes)
#     [U, G_init] = Solve(Square_Triangulate, N, U_exact_2, F, on_boundary)
#     X,Y = Nodes[E[0,1::],1], Nodes[E[0,1::],2]
#     J = abs((X[1] - X[0])*(Y[2]-Y[0]) - (X[2]-X[0])*(Y[1]-Y[0]))
#     Area.append(0.5*J)
#     L_dom = Linear_interp(elements = E, nodes = Nodes, U = U)
#     E_r, Nodes_r, Mesh_r = Square_Triangulate(Ns[i+1])
#
#     U_refine = L_dom.eval(Nodes_r[:,1], Nodes_r[:,2])
#
#     e_h1 = H1_norm(E_r, Nodes_r, U_refine-U_exact_2(Nodes_r[:,1], Nodes_r[:,2]))
#     e_l2 = H1_norm(E_r, Nodes_r, U_refine-U_exact_2(Nodes_r[:,1], Nodes_r[:,2]), which = "L2")
#
#     H1_error.append(e_h1)
#     L2_error.append(e_l2)
#     i+=1
#
#
# #Convergence plots
# fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (8,8))
# hs = np.array(Area)
#
#
# m = np.log(H1_error[-1]/H1_error[-2])/np.log((hs[-1]/hs[-2]))
# # print(m)
# two = 10*hs**2
# one = 10*hs
#
# ax[0].loglog(hs, H1_error, c = "k")
# ax[0].scatter(hs, H1_error, c = "k")
# ax[0].loglog(hs, one, c = 'k', linestyle = '--')
# ax[0].text(hs[2], one[3],r"$\mathcal{O}(h)$" %(m), fontsize = 18)
# ax[0].set_xlabel("Area of largest triangle", fontsize = 15)
# ax[0].set_ylabel(r"$||u_{exact} - u||_{H^1}$", fontsize = 15)
# ax[0].tick_params(axis='both', labelsize=15)
# ax[0].set_title(r"$H^1 \, error$", fontsize = 15)
# ax[0].grid()
# plt.savefig("Poisson_convergence.png")
#
#
# m = np.log(L2_error[-1]/L2_error[-2])/np.log((hs[-1]/hs[-2]))
# # print(m)
#
# ax[1].loglog(hs, L2_error, c = "k")
# ax[1].scatter(hs, L2_error, c = "k")
# ax[1].loglog(hs, two, c = 'k', linestyle = '--')
# ax[1].text(hs[2], two[3],r"$\mathcal{O}(h^2)$" %(m), fontsize = 18)
# ax[1].set_xlabel("Area of largest triangle", fontsize = 15)
# ax[1].set_ylabel(r"$||u_{exact} - u||_{L^2}$", fontsize = 15)
# ax[1].set_title(r"$L^2 \, error$", fontsize = 15)
# ax[1].tick_params(axis='both', labelsize=15)
# ax[1].grid()
# plt.subplots_adjust(wspace = 0.5)
# plt.savefig("Poisson_convergence_L2.png")
