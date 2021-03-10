import numpy as np
import pdb
# Code to assemble elemental stiffness matrices:
class Poisson_Assemble(object):
    """
    Class storing a quadtree structure defined over a rectangle.

    Inputs:
        QT - Quadtree class

    ## this will likely contain other inputs in the future

    Methods:
    loc(Tri)- Function to assemble local stiffness matrix.
        Inputs: list indexing nodes

    """
    def __init__(self, QT):
        self.QT = QT
        return

    def loc(self, x_s, y_s):

        T = 0.5*abs(((x_s[1]- x_s[0])*(y_s[2]-y_s[0]) - \
                 (x_s[2]-x_s[0])*(y_s[1]-y_s[0])))

        Gg = np.linalg.inv(np.array([[1,1,1],x_s,y_s]))
        g = np.array([[0,0],[1,0],[0,1]])
        G = np.matmul(Gg, g)
        # M = T*np.matmul(G, G.T)
        M = (0.5)*abs(np.linalg.det(np.array([[1,1,1],x_s,y_s])))*np.matmul(G, G.T)

        return M

    def number_nodes(self):
        """
        N = starting number of nodes, i  = number of refinements
        returns an array of # num of nodes along each axis at each stage
        """
        ref = [9]
        for n in range(self.QT.layers):
            ref.append(2*ref[n]-1)

        return ref

    def ass_mat(self):
        """
        function to assemble entire mass matrix for the Poisson system.
        """
        #find the number of nodes
        NN = self.number_nodes()[self.QT.layers]

        A = np.zeros([NN,NN], dtype = "float")

        #obtain the connectivity of the last layer in the given quadtree
        T = self.QT.get_T(layer = self.QT.layers)

        #assemble full mass matrix
        XY = np.array(self.QT.xy.copy())

        for tri in T:
            pdb.set_trace()
            xy_s = XY[tri]
            x_s_temp = [x for x in xy_s[:,0]]
            y_s_temp = [x for x in xy_s[:,1]]

            M_loc = self.loc(x_s_temp, y_s_temp)

            for i in range(3):
                for j in range(3):
                    A[tri[i], tri[j]] += M_loc[i,j]


        return A

#------------------------------------------------------------------------------

# RHS assembly

class RHS(object):
    """
    Class to store the RHS of our problem,
    input is the forcing function (RHS) of the problem

    input: F[X,Y]: Omega \to R^N,
    N is the number of nodes, X,Y \in R^{N^2}.
    """

    def __init__(self, f, g, elements, nodes):
        # make sure these are callable
        self.f = f
        self.g = g
        self.E = elements
        self.Nodes = nodes
        return

    def b(self, boundary):
       # Calculates volume force at every node
        F_vol = np.zeros([len(self.Nodes[:,0]),])
        for e in range(len(self.E[:,0])):
            n_s = self.E[e, 1::]
            X,Y = self.Nodes[n_s, 1], self.Nodes[n_s,2]
            Area = 0.5*abs((X[1] - X[0])*(Y[2]-Y[0]) - (X[2]-X[0])*(Y[1]-Y[0]))
            i = 0
            for n in n_s:
                F_vol[n] += (1/3)*Area*self.f(X[i], Y[i])
                i +=1

        return F_vol



def on_boundary(Nodes):
    """
    Only for the square
    TODO: Generalize this function
    """

    bn = [] # Boundary nodes
    eps = 1e-11
    for i in range(len(Nodes[:,0])):
        n, x, y = Nodes[i,0], Nodes[i,1], Nodes[i,2]

        # This is just for the square

        if (x < eps) or x > (1 - eps) or y < eps or y > (1-eps):
            bn.append(int(n))

    return bn
