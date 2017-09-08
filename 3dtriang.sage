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
        # probability of allowing flips INTO top or out of bottom
        self.prob_top = 1.0
        # prob of allowing flips into bottom or out of top
        self.prob_bottom = 0.0
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
                self.neck[tr] = self.adj[tr]

        self.k = len(self.neck)
        self.l = 0
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.top, self.bottom = self.initialize_halves()

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
            if tr in self.top:
                viz_tri.set_color('g')
            if tr in self.bottom:
                viz_tri.set_color('b')
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

    def allow(self, case):
        return True

    def global_cache(self, case):
        # reconfigure inside neck
        if case == 1:
            return True
        # remove triangle from neck
        elif case == 2:
            self.l = self.l - 1
            return True
        # flip triangle into neck
        elif case == 3 and self.l < 3:
            self.l = self.l + 1
            return True
        else:
            return False


    # choose random triangle, choose random neighbors, and move vertices until
    # we get a flip
    def make_flip(self, rule):
        t1 = choice(self.adj.keys())
        t2 = choice(self.adj[t1])
        new_t1, new_t2 = expected_flip(t1, t2)

        # CASE 0: allow all flips that do not affect neck
        if (t1 not in self.neck) and (t2 not in self.neck):
            if self.flip_triangs(t1, t2):
                if t1 in self.top:
                    del self.top[t1]
                    del self.top[t2]
                    self.top[new_t1] = self.adj[new_t1]
                    self.top[new_t2] = self.adj[new_t2]
                else:
                    del self.bottom[t1]
                    del self.bottom[t2]
                    self.bottom[new_t1] = self.adj[new_t1]
                    self.bottom[new_t2] = self.adj[new_t2]

        # CASE 1: flip two triangles in neck
        # two sub cases
        elif (t1 in self.neck) and (t2 in self.neck):
            ns1 = [n for n in self.adj[t1] if n in self.neck and n != t2]
            ns2 = [n for n in self.adj[t2] if n in self.neck and n != t1]
            if len(ns1) == 1 and len(ns2) == 1:
                # check if any vertices shared by neighbors
                common = list(set(ns1[0]).intersection(set(ns2[0])))
                # SUBCASE 0: flip keeps both in neck
                if len(common) == 0 and rule(1):
#                    print "flip inside neck"
                    self.flip_and_add_both(t1, t2)
                # SUBCASE 1: flip one triangle out of neck
                elif len(common) == 1 and rule(2):
                    if ns1[0] in self.top and random() <= self.prob_top:
                        print "flip into top"
                        self.flip_out_of_neck(t1, t2, common)
                    elif   ((ns1[0] in self.bottom) and 
                            (random() <= self.prob_bottom)):
                        print "flip into bottom"
                        self.flip_out_of_neck(t1, t2, common)

            else:
                print("more than two neighbors of a diamond are in the neck")

        # CASE 2: only remaining case: one in out out, flip one into neck
        elif ((t1 in self.neck) and (t2 not in self.neck)) and rule(3):
            if t2 in self.top and random() <= self.prob_bottom:
                print "flip out of top"
                self.flip_into_neck(t1, t2)
            elif (t2 in self.bottom) and random() <= self.prob_top:
                print "flip out of bottom"
                self.flip_into_neck(t1, t2)
        elif ((t2 in self.neck) and (t1 not in self.neck)) and rule(3):
            if t1 in self.top and random() <= self.prob_bottom:
                self.flip_into_neck(t2, t1)
                print "flip out of top"
            elif (t1 in self.bottom) and random() <= self.prob_bottom:
                print "flip out of bottom"
                self.flip_into_neck(t2, t1)

    def flip_into_neck(self, nt, t):
        ns = [n for n in self.adj[t] if n in self.neck and n != nt]
        #check that only one neighbor is in neck
        if len(ns) == 0:
            if t in self.top:
                if self.flip_and_add_both(nt, t):
                    self.k = self.k + 1
                    del self.top[t]
            else:
                if self.flip_and_add_both(nt, t):
                    self.k = self.k + 1
                    del self.bottom[t]

    def flip_out_of_neck(self, t1, t2, common):
        if self.flip_triangs(t1, t2):
#            print "flip one out of neck"
            del self.neck[t1]
            del self.neck[t2]
            self.k = self.k - 1
            new_t1, new_t2 = expected_flip(t1,t2)
            if common[0] in new_t1:
                new_neck_triangle = new_t2
                new_body_triangle = new_t1
            else:
                new_neck_triangle = new_t1
                new_body_triangle = new_t2

            self.neck[new_neck_triangle] = self.adj[new_neck_triangle]
            ns = [n for n in self.adj[new_body_triangle]
                    if n not in self.neck]
            if ns[0] in self.top:
                self.top[new_body_triangle] = self.adj[new_body_triangle]
            else:
                self.bottom[new_body_triangle] = self.adj[new_body_triangle]


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

    def flip_and_add_both(self,t1, t2):
        if self.flip_triangs(t1, t2):
            if t1 in self.neck:
                del self.neck[t1]
            if t2 in self.neck:
                del self.neck[t2]
            new_t1, new_t2 = expected_flip(t1, t2)
            self.neck[new_t1] = self.adj[new_t1]
            self.neck[new_t2] = self.adj[new_t2]
            return True
        else:
            return False

    def bfs(self, queue, collection):
        if queue == []:
            return collection
        new_queue = []
        for tr in queue:
            ns = self.adj[tr]
            for n in ns:
                if  (n not in self.neck) and \
                    (n not in collection) and \
                    (n not in queue) and \
                    (n not in new_queue):
                    new_queue.append(n)
            collection[tr] = ns
        return self.bfs(new_queue, collection)


    def initialize_halves(self):
        n1 = next(iter(self.neck.keys()))
        t1 = [n for n in self.adj[n1] if n not in self.neck][0]
        n2 = [n for n in self.adj[n1] if n in self.neck][0]
        t2 = [n for n in self.adj[n2] if n not in self.neck][0]
        top_half = self.bfs([t1],{})
        bottom_half = self.bfs([t2],{})
        return top_half, bottom_half


#    def count_both_sides(self):
#        tris = [t for t in self.triangles.keys() if t not in self.neck]
#
#        return 10

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

def make_adjacency_graph(ts):
    tis = sorted(ts.keys())
    neighbors = {t:[] for t in tis}
    for t in tis:
        for s in tis:
            if are_neighbors(s, t):
                neighbors[t].append(s)
    return neighbors



b = Blob()
c = Blob(100)

def plot_count(xdat, ydats, count="", title="", label=""):
    plt.figure()
    colors = ['r','g','b']
    for ((y, labl), col) in zip(ydats, colors):
        plt.plot(xdat, y, '-o', label=labl, color=col)
    plt.title(title)
    plt.xlabel(count)
    legend = plt.legend(loc='upper left')
    plt.savefig(label+".pdf", bbox_inches='tight')
    #plt.show()

def do_movement(blob, n, viz = False):
    neck_len = [0]*n
    top_size = [0]*n
    bottom_size = [0]*n
    for i in range(n):
        if i%100 == 0:
            print i
        blob.make_flip(blob.allow)
        neck_len[i] = blob.k
        top_size[i] = len(blob.top)
        bottom_size[i] = len(blob.bottom)

    plot_count(range(n), [(neck_len, "neck"), (top_size, "top"), (bottom_size,
    "bottom")], "Triangles over time", label="neck_len_n"+str(n))
    if viz:
        blob.show()

   


