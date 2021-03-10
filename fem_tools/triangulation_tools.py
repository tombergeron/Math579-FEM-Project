import numpy as np
import stripy
import pdb

class Node():
    """
    Base class for an element of a quadtree

    Attributes:
        pts (set): indices of vertices defining the triangle
        parent (Node object): "parent node". If root node then
        returns a string "root".
        children (list of node objects): list of four nodes forming
        the branch of current node
    """
    def __init__(self, pts, parent = "root", layer = 0):

        self.pts = pts
        self.parent = parent
        self.children = []
        self.divided = False
        self.layer = layer
        return



class plTriangulation():

    """
    Class based on the Triangulation() object included in the stripy
    package

    Attributes:
        x : array of floats, shape (n,)
            stored Cartesian x coordinates from input
        y : array of floats, shape (n,)
            stored Cartesian y coordinates from input

    """
    def __init__(self, x,y):

        self.T = stripy.Triangulation(x,y)

        return

    @property
    def x(self):
        return self.T.x

    @property
    def y(self):
        return self.T.y

    @property
    def simplices(self):
        return self.T.simplices



class Rectangle():
    """
    Rectangle

    """

    def __init__(self, bl = [0,0], tr = [1,1]):

        """
        Initialize base triangle formed of 8 uniform triangles
        -----
        |\  |
        | \ |
        |  \|
        -----
        """

        xs = [bl[0], (tr[0] + bl[0])/2, tr[0], bl[0], (tr[0] + bl[0])/2, tr[0],
              bl[0], (tr[0] + bl[0])/2, tr[0]]

        ys = [bl[1],bl[1], bl[1],
              (tr[1] + bl[1])/2,(tr[1] + bl[1])/2, (tr[1] + bl[1])/2,
            tr[1], tr[1], tr[1]]

        # preferable data structure subject to change though
        xy = [np.array([xs[i], ys[i]]) for i in range(9)]

        self.xy = xy.copy()
        #connectivity matrix
        # initial spatial grid

        self.T = [{0,1,3},{1,3,4},{1,2,4},{2,4,5},{3,4,6},{4,6,7},{4,5,7},{5,7,8}]
        # first layer of tree structure
        self.Tree = [Node(t) for t in self.T]
        self.nodes = len(xy)
        self.layers = 0

        return


    def branch(self, node):
        # single triangle input
        tri = node.pts
        # use pts to slice xy list

        self.xy.append(self.xy[tri[0]]/2 + self.xy[tri[1]]/2)
        self.xy.append(self.xy[tri[1]]/2 + self.xy[tri[2]]/2)
        self.xy.append(self.xy[tri[2]]/2 + self.xy[tri[0]]/2)

        # Create another 4 Node children
        i = self.nodes
        tt = [[tri[0],i,i+2], [i,tri[1],i+1], [i+1,tri[2],i+2], [i, i+1, i+2]]
        node.children = [Node(t, parent = node, layer = node.layer+1) for t in tt]
        self.nodes +=3

        return

    def gen_grid(self, layers):

        """
        function to output and ordered list of vertices and a connectivity list
        of triangles (lists of 3 indices).

        input:
            layers (int): how many layers of refinement from initialized grid.
        """
        T_con = [list(x) for x in self.T]

        t_out = []
        xys = self.xy.copy()

        for tri in T_con:
            xys.append(xys[tri[0]]/2 + xys[tri[1]]/2)
            xys.append(xys[tri[1]]/2 + xys[tri[2]]/2)
            xys.append(xys[tri[2]]/2 + xys[tri[0]]/2)



        return


    def refine(self, node):

        if node.layer < self.layers:
            for child in node.children:
                self.refine(child)
        else:
            self.branch(node)
            return
        #self.nodes = 2 + 10*(N+k)**2
        return


    def add_layer(self):

        for node in self.Tree:
            self.refine(node)

        self.layers +=1

        return

    def connectivity(self, N, k, T_con):
        """
        output list of connectivity elements for the base level of children
        """
        if N.layer < k:
            for child in N.children:
                self.connectivity(child, k, T_con)
        else:
            T_con.append(N.pts)
            return
        #self.nodes = 2 + 10*(N+k)**2
        return

    def get_T(self, layer):
        TT = []
        for nn in self.Tree:
            self.connectivity(nn, layer, TT)
        return TT


    def NN_search(self, Qs):
        """
        Function to perform nearest neighbour search on quadtree structure

        Input: A - Node object
               Q: [N,3] numpy array

        Output: Triangle containing point
        """
        self.dist_init(Qs)
        # first find root node triangle:

        #compute distance from nodes:

        return



import scipy.io
# Small test to produce data


# For Rectangle
# # S = Rect_QT()
#
# for k in [1,2,3,4]:
#     S.add_layer()

# For circle
# phis = np.linspace(0, 2*np.pi, 50)
# x,y = np.cos(phis), np.sin(phis)
# x = np.append(x, [0])
# y = np.append(y, [0])
#
# tri = plTriangulation(x,y)

# For circle
x = np.random.uniform(-1,1,10)
y = np.random.uniform(-1,1,10)
[X,Y] = np.meshgrid(x,y)


tri = plTriangulation(x,y)


# TT = []
# for nn in S.Tree:
#     S.connectivity(nn, 4, TT)

scipy.io.savemat("../../Plotting-Tools/Data_Plots/sample_triangulate.mat", {"pts_x": tri.x , "pts_y": tri.y, "T": tri.simplices})




    # def uniform_refine(self, node, k):
    #     """
    #     function to add k layers to tree stucture of T.
    #
    #
    #     G- np.array - 1D list of coordinates
    #     N defines how many refinements to perform on the grid G
    #     output: quadtree structure for refined mesh
    #
    #                /0\
    #               /___\    "ordering they are added"
    #              /\ 3/ \
    #             /_1\/ 2_\
    #
    #     0 is the first node in parent.
    #     k = depth of tree structure. 0 corresponds to first layer.
    #     """
    #     if node.layer < k:
    #         self.branch(node)
    #         for child in node.children:
    #             self.uniform_refine(child, k)
    #     else:
    #         return
    #     #self.nodes = 2 + 10*(N+k)**2
    #     return




# #
# def Square_Triangulate(N, x0 = [0,0], x1 = [1,1]):
#     """
#     Inputs:
#     N - Number of point along each axis
#     x0,x1 - x,y coordinates of bottom left and top right corners of square.
#     """
#
#     x_coords = np.linspace(x0[0], x1[0], N, endpoint = True)
#     y_coords = np.linspace(x0[1], x1[1], N, endpoint = True)
#
#     # Put this into a meshgrid for easier access to coordinate values.
#
#     [X,Y] = np.array(np.meshgrid(x_coords, y_coords))
#     #Put this in a more convenient form:
#     Nodes = np.zeros([N**2, 3])
#     Nodes[:,0], Nodes[:,1], Nodes[:,2] = np.arange(0,int(N**2), dtype = int), X.reshape([N**2,]), Y.reshape([N**2,])
#
#     # Label and create triangulated mesh data structure
#     E_u1 = []
#     E_u2 = []
#     E_u3 = []
#     for i in range(0, N-1):
#         b = i*N  #beginning in nodes
#         l_row = Nodes[b:b+N,0].copy()
#         u_row = Nodes[b+N:b+2*N,0].copy()
#         E_u1 += list(l_row[0:N-1])
#         E_u2 += list(l_row[1:N])
#         E_u3 += list(u_row[1:N])
#         E_u1 += list(l_row[0:N-1])
#         E_u2 += list(u_row[0:N-1])
#         E_u3 += list(u_row[1:N])
#
#
#     M = len(E_u1)
#     E = np.zeros([M, 4], dtype = int)
#     E[:,0] = np.arange(0, M, dtype = int) # "--------=----------"
#     E[0:len(E_u1),1] = E_u1
#     E[0:len(E_u2),2] = E_u2
#     E[0:len(E_u3),3] = E_u3
#
#     return E, Nodes, [X,Y]


# # Function to plot mesh
# def Mesh_Plot(E, Nodes, boundary = None):
#     """
#     Function to output plot of the mesh created with labelled elements
#     """
#     fig = plt.figure(figsize = (7,7))
#
#     #First create scatter plot with all the nodes
#
#     plt.scatter(Nodes[:,1], Nodes[:,2], c = "k", alpha = 0.5)
#
#     #Then plot lines and label elements
#     for row in E:
#         num = row[0]
#         xs = Nodes[row[1::],1]
#         ys = Nodes[row[1::], 2]
#         # calculate center of element
#         x_c, y_c = np.sum(xs)/3, np.sum(ys)/3
#         plt.text(x_c, y_c, s = "%d" %num, fontsize = 10)
#
#         # Add lines connecting the nodes for each element
#         plt.plot([xs[0], xs[1]], [ys[0], ys[1]], c= "k")
#         plt.plot([xs[1], xs[2]], [ys[1], ys[2]], c= "k")
#         plt.plot([xs[2], xs[0]], [ys[2], ys[0]], c= "k")
#
#     if boundary != None:
#         plt.scatter(Nodes[boundary,1], Nodes[boundary, 2], c = "b")
#
#     return
#

#
# def red_refinement(N,num):
#     """
#     N = starting number of nodes, i  = number of refinements
#     returns an array of # num of nodes along each axis at each stage
#     """
#     ref = [N]
#     for n in range(num-1):
#         ref.append(2*ref[n]-1)
#
#     return ref
#
# Ns = red_refinement(4, 2)

# E, Nodes, Mesh = Square_Triangulate(Ns[0])
# Mesh_Plot(E, Nodes)
#
# E, Nodes, Mesh = Square_Triangulate(Ns[1])
# Mesh_Plot(E, Nodes)
