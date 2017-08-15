import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from itertools import combinations
from sage.misc.prandom import randrange, uniform, choice
from copy import copy


class Blob:
    def __init__(self, n=45):
        self.t = int(n/3)
        # number of triangles in neck must be even
        if mod(self.t,2) != 0:
            self.t = self.t + 1
        self.n = n

        # "height" of neck (should be less than 0.5)
        self.delz = 0.2
        self.neck = {}

        self.pts = np.array(self.make_neck_points() + self.make_random_points())
        self.delauney = Delaunay(self.pts)
        self.triangles = {tuple(sorted(tr)):np.array([list(self.pts[i]) for i in tr])
                        for tr in self.delauney.convex_hull}
        self.adj = make_adjacency_graph(self.triangles)

        for tr in self.triangles.keys():
            if (tr[0] <= 2*self.t and tr[1] <= 2*self.t and tr[2] <= 2*self.t):
                self.neck[tr] = self.triangles[tr]

        fig = plt.figure()
        self.ax = fig.add_subplot(111, projection='3d')

    def update_triangulation(self, pts):
        self.pts = pts
        self.delauney = Delaunay(self.pts)
        self.triangles = {tuple(tr):np.array([list(self.pts[i]) for i in tr])
                        for tr in self.delauney.convex_hull}
        self.adj = make_adjacency_graph(self.triangles)

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
        for tr in self.triangles.keys():
            pts = self.triangles[tr]
            viz_tri = Poly3DCollection([pts])
            if self.is_in_neck(tr):
                viz_tri.set_color('r')
            viz_tri.set_edgecolor('k')
            self.ax.add_collection3d(viz_tri)

    def show(self):
         self.plot_triangle()
         self.ax.scatter(self.pts[:, 0], self.pts[:, 1], self.pts[:, 2], c='k')
         #plt.savefig("Delaunay.png", dpi=600)
         plt.show()

    # choose random triangle, choose random neighbors, and move vertices until
    # we get a flip
    def make_flip(self):
        tr1 = choice(self.adj.keys())
        tr2 = choice(self.adj[tr1])
        # case 1: reconfigure inside neck
        if tr1 in self.neck and tr2 in self.neck:
            v1, v2, shared = get_opposite_verts(tr1, tr2)
            if self.squish_verts(v1, v2):
                new_t1 = tuple(sorted([v1,v2,shared[0]]))
                new_t2 = tuple(sorted([v1,v2,shared[1]]))
                del self.neck[tr1]
                del self.neck[tr2]
                self.neck[new_t1] = self.triangles[new_t1]
                self.neck[new_t2] = self.triangles[new_t2]
                return True
        if tr1 not in self.neck and tr2 not in self.neck:
            v1, v2, shared = get_opposite_verts(tr1, tr2)
            return self.squish_verts(v1, v2)


    def squish_verts(self, v1, v2):
        p1, p2 = self.pts[v1], self.pts[v2]
        vect = p2-p1
        epsilon = 0.03
        p1_new = p1 + epsilon*vect
        p1_new = p1_new/(np.linalg.norm(p1_new))
        p2_new = p2 - epsilon*vect
        p2_new = p2_new/(np.linalg.norm(p2_new))
        pts_new = copy(self.pts)
        pts_new[v1], pts_new[v2] = p1_new, p2_new
        delau = Delaunay(pts_new)
        ch1 = set([tuple(sorted(t)) for t in self.delauney.convex_hull])
        ch2 = set([tuple(sorted(t)) for t in delau.convex_hull])

        if ch1 == ch2:
            self.update_triangulation(pts_new)
            self.squish_verts(v1, v2)
        else:
            if convex_hull_contains_flip(ch1, ch2):
                print "made a flip!"
                self.update_triangulation(pts_new)
                return True
            else:
                return False

# given two adjacent triangles, get non-connected vertices
def get_opposite_verts(tr1, tr2):
    vs = [0,0]
    shared = []
    for ti in tr1:
        if ti not in tr2:
            vs[0] = ti
        else:
            shared.append(ti)
    for ti in tr2:
        if ti not in tr1:
            vs[1] = ti
    return vs[0], vs[1], shared


# scan two lists of triangles composing convex hull
# if perturbation has only created one flip, only two triangles will be
# different
def convex_hull_contains_flip(ch1, ch2):
    diffs = []
    # don't allow triangles to be created or destroyed
    if len(ch1) != len(ch2):
        return False

    diff1 = [t for t in ch1 if t not in ch2]
    diff2 = [t for t in ch2 if t not in ch1]
    if len(diff1) == 2 and len(diff2) == 2:
        return is_pairwise_flip(diff1, diff2)
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
    tis = sorted(ts.keys())
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


b = Blob()

def do_movement(n):
    for i in range(n):
        res = b.make_flip()
        if res:
            b.show()
            plt.close("all")

      


    # recursively try to nudge one vertex at a time until we get a good flip
#    def try_perturbation(self):
#        ch1 = sorted(self.triangles)
#        pert_pts = move_vert(pts)
#        tr2 = Delaunay(pert_pts)
#        ch2 = np.sort(tr2.convex_hull, axis=0)
#        if convex_hull_contains_flip(ch1, ch2):
#            return pert_pts
#        else:
#            return pert_pts
#            #return try_perturbation(pts)

## given a vertex, perturb it in a random direction (proportional to epsilon)
## the returned vertex will still be on unit sphere
#def perturb_one(i, verts, epsilon = 0.05):
#    rand_pt = [ uniform(-epsilon,epsilon),
#                uniform(-epsilon,epsilon),
#                uniform(-epsilon,epsilon)]
#    perturb = [xi + eps for (xi, eps) in zip(verts[i], rand_pt)]
#    norm = sum([i**2 for i in perturb])
#    rand_pt = [i/norm for i in perturb]
#    return rand_pt
#
## move one random vertex in a random direction by 1/10 the average inter-vertex
## radius
#def move_vert(verts):
#    n = len(verts)
#    i = randrange(0, n)
#    avg_rad = N(2/sqrt(n))
#    eps = 0.1*avg_rad
#    verts[i] = perturb_one(i, verts, eps)
#    return verts
