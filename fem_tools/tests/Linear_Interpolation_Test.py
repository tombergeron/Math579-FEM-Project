import numpy as np
from Square_Triangulate import Square_Triangulate
from Interpolator import Linear_interp
import matplotlib.pyplot as plt

# Linear interpolation test
def uu_exact(x,y):

    return  np.sin(3*x)*np.cos(y)


err = []
Ns = np.array([4, 8, 16, 32, 64, 128])
two = 0.5/(Ns**2)
for N in Ns:

    E, Nodes, Mesh = Square_Triangulate(N)
    U = uu_exact(Nodes[:,1],Nodes[:,2])
    L_dom = Linear_interp(elements = E, nodes = Nodes, U = U)
    E_r, Nodes_r, Mesh_r = Square_Triangulate(2*N)
    U_refine = L_dom.eval(Nodes_r[:,1], Nodes_r[:,2])

    Error = abs(uu_exact(Nodes_r[:,1],Nodes_r[:,2]) - U_refine)
    err.append(np.max(Error))

fig = plt.figure(figsize = (7,7))
m = np.log(err[-1]/err[-2])/np.log((Ns[-1]/Ns[-2]))

plt.loglog(Ns, err, c = "k", linestyle = 'solid')
plt.scatter(Ns,err, c= 'k')
plt.loglog(Ns, two, c = "k", linestyle = "--")
plt.text(Ns[2], two[-3],r"$\mathcal{O}(N^{2})$" %(-m), fontsize = 15)
plt.xlabel("Number of points along each axis")
plt.ylabel(r"$||u_{exact} - u||_{L^{\infty}}$", fontsize = 15)


plt.savefig("interpolation_test.png")
