import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from itertools import combinations
from sage.misc.prandom import randrange, uniform, choice
from copy import copy


class Blob:
    def __init__(self, t=15, n=30):
        self.t = t
        # number of triangles in neck must be even
        if mod(t,2) != 0:
            self.t = t + 1
        self.n = n

        # "height" of neck (should be less than 0.5)
        self.delz = 0.2
        self.neck = {}

        self.pts = np.array(self.make_neck_points() + self.make_random_points())
        self.delauney = Delaunay(self.pts)
        self.triangles = [(tuple(tr), np.array([list(self.pts[i]) for i in tr]))
                        for tr in self.delauney.convex_hull]

        for (tr, points) in self.triangles:
            if (tr[0] <= 2*self.t and tr[1] <= 2*self.t and tr[2] <= 2*self.t):
                self.neck[tr] = points

        fig = plt.figure()
        self.ax = fig.add_subplot(111, projection='3d')

    def is_in_neck(self, tr):
        return (tr in self.neck)

    def make_neck_points(self):
        r = sqrt(1-(self.delz**2)/4)

        pts = []
        pt1 = [0,0,self.delz/2]
        pt2 = [0,0,-self.delz/2]
        for i in range(self.t):
            theta = N((2*pi/self.t)*i)
            pt1[0] = N(r*cos(theta))
            pt1[1] = N(r*sin(theta))

            theta2 = N(theta + (pi/self.t))
            pt2[0] = N(r*cos(theta2))
            pt2[1] = N(r*sin(theta2))

            pts.append(copy(pt1))
            pts.append(copy(pt2))
        return pts

    def make_random_points(self):
        # make random points
        S = SphericalDistribution()
        rand_pts = [list(S.get_random_element()) for i in range(self.n)]
        rand_pts = [p for p in rand_pts if abs(p[2]) > self.delz*1.8]
        pts = [[0.0, 0.0, 0.0]] + rand_pts
        return pts

    def plot_triangle(self):
        for (indx, tpts) in self.triangles:
            viz_tri = Poly3DCollection([tpts])
            if self.is_in_neck(indx):
                viz_tri.set_color('r')
            viz_tri.set_edgecolor('k')
            self.ax.add_collection3d(viz_tri)

    def show(self):
         self.plot_triangle()
         self.ax.scatter(self.pts[:, 0], self.pts[:, 1], self.pts[:, 2], c='k')
         #plt.savefig("Delaunay.png", dpi=600)
         plt.show()

    # recursively try to nudge one vertex at a time until we get a good flip
    def try_perturbation(self):
        ch1 = sorted(self.triangles)
        pert_pts = move_vert(pts)
        tr2 = Delaunay(pert_pts)
        ch2 = np.sort(tr2.convex_hull, axis=0)
        if convex_hull_contains_flip(ch1, ch2):
            return pert_pts
        else:
            return pert_pts
            #return try_perturbation(pts)

    # choose random triangle, choose random neighbors, and move vertices until
    # we get a flip
    def make_flip(self):
        adj = make_adjacency_graph(self.triangles)
        tr1 = choice(adj.keys())
        tr2 = choice(adj[tr1])
        vs = []
        for ti in tr1:

        return tr1, tr2




# given a vertex, perturb it in a random direction (proportional to epsilon)
# the returned vertex will still be on unit sphere
def perturb_one(i, verts, epsilon = 0.05):
    rand_pt = [ uniform(-epsilon,epsilon),
                uniform(-epsilon,epsilon),
                uniform(-epsilon,epsilon)]
    perturb = [xi + eps for (xi, eps) in zip(verts[i], rand_pt)]
    norm = sum([i**2 for i in perturb])
    rand_pt = [i/norm for i in perturb]
    return rand_pt

# move one random vertex in a random direction by 1/10 the average inter-vertex
# radius
def move_vert(verts):
    n = len(verts)
    i = randrange(0, n)
    avg_rad = N(2/sqrt(n))
    eps = 0.1*avg_rad
    verts[i] = perturb_one(i, verts, eps)
    return verts

# scan two lists of triangles composing convex hull
# if perturbation has only created one flip, only two triangles will be
# different
def convex_hull_contains_flip(ch1, ch2):
    diffs = []
    # don't allow triangles to be created or destroyed
    if len(ch1) != len(ch2):
        return False
    difn = 0
    for (ti, tj) in zip(ch1, ch2):
        if not np.array_equal(ti,tj):
            if (difn == 2):
                return False
            diffs.append((ti,tj))
            difn += 1

    if len(diffs) == 2:
        quad1 = [ti for (ti,tj) in diffs]
        quad2 = [tj for (ti,tj) in diffs]
        if is_pairwise_flip(quad1, quad2) == True:
            return True
    else:
        return False

# given two pairs of triangles, return true if they are related by one flip
# if they are the same vertex set (four vertices) but different triangulation,
# they must be a flip
# https://en.wikipedia.org/wiki/Delaunay_triangulation#Visual_Delaunay_definition:_Flipping
def is_pairwise_flip(ts1, ts2):
    if  (Set([indx for ti in ts1 for indx in ti])
        ==  Set([indx for ti in ts2 for indx in ti])):
        return True
    else:
        return False

def are_neighbors(t1, t2):
    check = [0,0,0]
    for (n,i) in enumerate(t1):
        for j in t2:
            if i == j:
                check[n] = 1
    if sum(check) == 2:
        return True
    else:
        return False

# not optimized
def make_adjacency_graph(ts):
    tis = sorted([ti for (ti, pts) in ts])
    neighbors = {}
    for t in tis:
        if t not in neighbors:
            neighbors[t] = []
        for s in tis:
            if len(neighbors[t]) < 3 and are_neighbors(s, t):
                neighbors[t].append(s)
    return neighbors





# TODO: see if we can access neighbors of triangles
# def count_both_sides(conv_hull):
#     c_hull = sort(conv_hull)



def do_movement(n, pts):
    for i in range(n):
        pts = try_perturbation(pts)
    tri = Delaunay(pts)
    return tri, pts

b = Blob()
