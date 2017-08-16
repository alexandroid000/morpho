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

        self.k = len(self.neck)
        self.l = 0
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')

    def update_triangulation(self, pts):
        self.pts = pts
        self.delauney = Delaunay(self.pts)
        self.triangles = {tuple(sorted(tr)):np.array([list(self.pts[i]) for i in tr])
                        for tr in self.delauney.convex_hull}
        self.adj = make_adjacency_graph(self.triangles)

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
        rand_pts = [p for p in rand_pts if abs(p[2]) > self.delz*1.5]
        pts = [[0.0, 0.0, 0.0]] + rand_pts
        return pts

    def plot_triangle(self):
        for tr in self.triangles.keys():
            pts = self.triangles[tr]
            viz_tri = Poly3DCollection([pts])
            if tr in self.neck:
                viz_tri.set_color('r')
            viz_tri.set_edgecolor('k')
            self.ax.add_collection3d(viz_tri)

    def show(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.plot_triangle()
        self.ax.scatter(self.pts[:, 0], self.pts[:, 1], self.pts[:, 2], c='k')
        #plt.savefig("Delaunay.png", dpi=600)
        plt.show()
        plt.close()

    def flip_and_add_both(self,t1, t2):
        if self.flip_triangs(t1, t2):
            if t1 in self.neck:
                del self.neck[t1]
            if t2 in self.neck:
                del self.neck[t2]
            new_t1, new_t2 = expected_flip(t1, t2)
            self.neck[new_t1] = self.triangles[new_t1]
            self.neck[new_t2] = self.triangles[new_t2]
            return True
        else:
            return False

    def allow(self, case):
        return True

    def global_cache(self, case):
#        print self.l
        # adding a triangle to neck
        if case == 1 and self.l < 3:
            self.l = self.l + 1
            return True
        # reconfigure inside neck
        elif case == 2:
            return True
        # remove triangle from neck
        elif case == 3:
            self.l = self.l - 1
            return True
        else:
            return False


    # choose random triangle, choose random neighbors, and move vertices until
    # we get a flip
    def make_flip(self, rule):
        t1 = choice(self.adj.keys())
        t2 = choice(self.adj[t1])
        new_t1, new_t2 = expected_flip(t1, t2)
        # case zero: allow all flips that do not affect neck
        if (t1 not in self.neck) and (t2 not in self.neck):
            self.flip_triangs(t1, t2)
        # case one: one in, one out - allow all, both are now in neck
        elif ((t1 in self.neck) and (t2 not in self.neck)) and rule(1):
#            print "flip one into neck"
            #check that only one neighbor is in neck
            ns = [n for n in self.adj[t2] if n in self.neck and n != t1]
            if len(ns) == 0:
                if self.flip_and_add_both(t1, t2):
                    self.k = self.k + 1

        elif ((t2 in self.neck) and (t1 not in self.neck)) and rule(1):
            ns = [n for n in self.adj[t1] if n in self.neck and n != t2]
            if len(ns) == 0:
                if self.flip_and_add_both(t1, t2):
                    self.k = self.k + 1

        # case two: reconfigurations inside neck
        # two sub cases
        elif (t1 in self.neck) and (t2 in self.neck):
            ns1 = [n for n in self.adj[t1] if n in self.neck and n != t2]
            ns2 = [n for n in self.adj[t2] if n in self.neck and n != t1]
            if len(ns1) == 1 and len(ns2) == 1:
                # check if any vertices shared by neighbors
                common = list(set(ns1[0]).intersection(set(ns2[0])))
                # sub case one: flip keeps both in neck
                if len(common) == 0 and rule(2):
#                    print "flip inside neck"
                    self.flip_and_add_both(t1, t2)
                # sub case two: flip excludes triangle not containing common
                # vertex from the neck
                elif len(common) == 1 and rule(3):
                    if self.flip_triangs(t1, t2):
#                        print "flip one out of neck"
                        del self.neck[t1]
                        del self.neck[t2]
                        self.k = self.k - 1
                        if common[0] in new_t1:
                            self.neck[new_t1] = self.triangles[new_t1]
                        else:
                            self.neck[new_t2] = self.triangles[new_t2]
            else:
                print("more than two neighbors of a diamond are in the neck")

    def flip_triangs(self, t1, t2):
        v1, v2, shared = get_opposite_verts(t1, t2)
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
            return self.flip_triangs(t1, t2)
        else:
            # check if we flipped the triangles we expected to
            new_t1, new_t2 = expected_flip(t1, t2)
            success, tris = new_triangles(ch1, ch2)
            if success and set([new_t1, new_t2]) == set(tris):
#                print "made a flip!"
                self.update_triangulation(pts_new)
                return True
            else:
                return False

def expected_flip(t1, t2):
    v1, v2, shared = get_opposite_verts(t1, t2)
    new_t1 = tuple(sorted([v1,v2,shared[0]]))
    new_t2 = tuple(sorted([v1,v2,shared[1]]))
    return new_t1, new_t2

# given two adjacent triangles, get non-connected vertices
def get_opposite_verts(t1, t2):
    vs = [0,0]
    shared = []
    for ti in t1:
        if ti not in t2:
            vs[0] = ti
        else:
            shared.append(ti)
    for ti in t2:
        if ti not in t1:
            vs[1] = ti
    return vs[0], vs[1], shared


# scan two lists of triangles composing convex hull
# if perturbation has only created one flip, only two triangles will be
# different
def new_triangles(ch1, ch2):
    diffs = []
    # don't allow triangles to be created or destroyed
    if len(ch1) != len(ch2):
#        print "perturbation created or destroyed triangles, abandoning"
        return False, []

    diff1 = [t for t in ch1 if t not in ch2]
    diff2 = [t for t in ch2 if t not in ch1]
    if len(diff1) == 2 and len(diff2) == 2:
        if is_pairwise_flip(diff1, diff2):
            return True, diff2
        else:
            return False, []
    else:
#        print "flipped too many triangles"
        return False, []

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
    neighbors = {t:[] for t in tis}
    for t in tis:
        for s in tis:
            if are_neighbors(s, t):
                neighbors[t].append(s)
    return neighbors

# TODO: see if we can access neighbors of triangles
# def count_both_sides(conv_hull):
#     c_hull = sort(conv_hull)


b = Blob()
c = Blob(200)

def plot_count(xdat, ydat, count="", title="", label=""):
    plt.figure()
    plt.plot(xdat, ydat, '-o')
    plt.title(title)
    plt.xlabel(count)
    plt.savefig(label+".pdf", bbox_inches='tight')
    #plt.show()

def do_movement(blob, n, viz = False):
    neck_len = [0]*n
    for i in range(n):
        if i%100 == 0:
            print i
        res = blob.make_flip(blob.allow)
        neck_len[i] = blob.k

    plot_count(range(n), neck_len, "Neck length", label="neck_len_n"+str(n))
    if viz:
        blob.show()

      


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
