import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from itertools import combinations
from sage.misc.prandom import randrange, uniform
from copy import copy

# number of triangles in neck
# must be even!
t = 10

if mod(t,2) != 0:
    print("you done messed up")
    quit()

# number of other points to be distributed on rest of sphere
n = 20

# "height" of neck (should be less than 0.5)
delz = 0.2

# make triangles
r = sqrt(1-(delz**2)/4)

pts = []

pt1 = [0,0,delz/2]
pt2 = [0,0,-delz/2]
for i in range(t):
    theta = N((2*pi/t)*i)
    pt1[0] = N(r*cos(theta))
    pt1[1] = N(r*sin(theta))

    theta2 = N(theta + (pi/t))
    pt2[0] = N(r*cos(theta2))
    pt2[1] = N(r*sin(theta2))

    pts.append(copy(pt1))
    pts.append(copy(pt2))


# make random points
S = SphericalDistribution()
rand_pts = [list(S.get_random_element()) for i in range(n)]
rand_pts = [p for p in rand_pts if abs(p[2]) > delz/2 + delz]
pts = [[0.0, 0.0, 0.0]] + pts + rand_pts

pts = np.array(pts)
tri = Delaunay(pts)

# given a vertex, perturb it in a random direction (proportional to epsilon)
# the returned vertex will still be on unit sphere
def perturb_one(i, verts):
    epsilon = 0.05
    rand_pt = [ uniform(-epsilon,epsilon),
                uniform(-epsilon,epsilon),
                uniform(-epsilon,epsilon)]
    perturb = [xi + eps for (xi, eps) in zip(verts[i], rand_pt)]
    norm = sum([i**2 for i in perturb])
    rand_pt = [i/norm for i in perturb]
    return rand_pt

def move_vert(verts):
    n = len(verts)
    i = randrange(0, n)
    verts[i] = perturb_one(i, verts)
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

# TODO make this better
#def share_edge(t1, t2):
#    [i1, j1, k1] = t1
#    [i2, j2, k2] = t2
#    if (i1 == i2 or i1 == j2 or i1 == k2):
#        if 
#

# TODO
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

# recursively try until we get a good flip
def try_perturbation(pts):
    tr1 = Delaunay(pts)
    ch1 = np.sort(tr1.convex_hull, axis=0)
    pert_pts = move_vert(pts)
    tr2 = Delaunay(pert_pts)
    ch2 = np.sort(tr2.convex_hull, axis=0)
    if convex_hull_contains_flip(ch1, ch2):
        return pert_pts
    else:
        return try_perturbation(pts)

def is_in_neck(tr):
    if (tr[0] <= 2*t and tr[1] <= 2*t and tr[2] <= 2*t):
        return True
    else:
        return False

# TODO: see if we can access neighbors of triangles
# def count_both_sides(conv_hull):
#     c_hull = sort(conv_hull)


def plot_triangle(tri):
    pts = tri.points
    triangles = [(tr, np.array([list(pts[i]) for i in tr])) for tr in tri.convex_hull]
    for (indx, tpts) in triangles:
        viz_tri = Poly3DCollection([tpts])
        if is_in_neck(indx):
            viz_tri.set_color('r')
        viz_tri.set_edgecolor('k')
        ax.add_collection3d(viz_tri)
    return ax

# given a list of tetrahedra, plot collection in 3d
def plot_tetra(tetra, pts, color="green", lc="k", lw=1):
    combs = combinations(tetra, 3)
    for comb in combs:
        X = pts[comb, 0]
        Y = pts[comb, 1]
        Z = pts[comb, 2]
        verts = [zip(X, Y, Z)]
        triangle = Poly3DCollection(verts, facecolors=color, alpha=1.0)
        lines = Line3DCollection(verts, colors=lc, linewidths=lw)
        ax.add_collection3d(triangle)
        ax.add_collection3d(lines)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot_triangle(tri)
ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c='k')
#plt.savefig("Delaunay.png", dpi=600)
plt.show()

def do_movement(n, pts):
    objs = []
    tri = Delaunay(pts)
    objs.append(plot_triangle(tri))
    for i in range(n):
        pts = try_perturbation(pts)
        objs.append(plot_triangle(Delaunay(pts)))
    return objs


