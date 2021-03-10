import scipy.spatial
import numpy as np
class Linear_interp(object):

    """
    Class to perform linear interpolation on a triangulated grid
    Nodes, Elements and Values at the vertices are used to initialize
    interpolator.
    """

    def __init__(self, elements, nodes, U):
        self.nodes = nodes.copy()
        self.E = elements.copy()
        self.u = U.copy()
        return

    def eval(self, x_q, y_q, point = False):
        """
        Inputs x_q, y_q are the query points in the form of 1D arrays.

        Output: array of points corresponding to a linear interpolation
        of the mesh function
        """

        xs = self.nodes[:,1]
        ys = self.nodes[:,2]
        XY_C = np.zeros([len(xs),2])
        XY_C[:,0] = xs
        XY_C[:,1] = ys

        mytree = scipy.spatial.cKDTree(XY_C)

        XY_q = np.zeros([len(x_q),2])
        XY_q[:,0] = x_q
        XY_q[:,1] = y_q


        # now query points in triangulation,
        # the k argument asks for 3 nearest neighbours
        dist, indices = mytree.query(x = XY_q, k = 3)

#-------Can now perform linear interpolation with these distances and indices------------

        # Solve linear system for the weights in barycentric coordinates
        # Loop over rows:
        i = 0
        u_ref = []

        for row in indices:
            u_p = np.array(self.u[row])
            A = np.array([np.array(xs[row]), np.array(ys[row]),[1,1,1]])

            if np.linalg.cond(A) < 1e11:
                B = np.linalg.inv(A)
                b = np.array([x_q[i], y_q[i],1])
                w = np.matmul(B,b)
                u_av = (1/np.sum(w))*(np.sum(u_p*w))
                u_ref.append(u_av)
            else:
                u_av = (1/np.sum(dist[i]))*np.sum((np.array(dist[i])*u_p))
                u_ref.append(u_av)

            i+=1

        return np.array(u_ref)
