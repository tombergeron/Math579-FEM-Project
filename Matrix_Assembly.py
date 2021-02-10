# Code to assemble elemental stiffness matrices:
class Stiffness(object):

    def __init__(self, Elements, Nodes):
        self.E = Elements
        self.Nodes = Nodes
        return

    def loc(self, element = 0):
        # returns local stiffness matrix for a certain element
        N_s = self.E[element, :][1::]

        x_s = self.Nodes[N_s, 1]
        y_s = self.Nodes[N_s, 2]

        T = 0.5*abs(((x_s[1]- x_s[0])*(y_s[2]-y_s[0]) - \
                 (x_s[2]-x_s[0])*(y_s[1]-y_s[0])))

        Gg = np.linalg.inv(np.array([[1,1,1],x_s,y_s]))
        g = np.array([[0,0],[1,0],[0,1]])
        G = np.matmul(Gg, g)
        # M = T*np.matmul(G, G.T)
        M = (0.5)*abs(np.linalg.det(np.array([[1,1,1],x_s,y_s])))*np.matmul(G, G.T)

        return M

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

    bn = [] # Boundary nodes
    eps = 1e-11
    for i in range(len(Nodes[:,0])):
        n, x, y = Nodes[i,0], Nodes[i,1], Nodes[i,2]

        # This is just for the square

        if (x < eps) or x > (1 - eps) or y < eps or y > (1-eps):
            bn.append(int(n))

    return bn
