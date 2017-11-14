import numpy as np
from scipy.spatial import Delaunay, ConvexHull
import matplotlib.pyplot as plt
from itertools import combinations
from functools import cmp_to_key
from copy import copy
import random

def random_points(n):
    random.seed()
    pts = [(random.random(), random.random()) for i in range(n)]
    return pts

def make_triangulation(n):
    pts = random_points(n)
    return Delaunay(pts)

def flatten(lofl):
    return [l for lelem in lofl for l in lelem]

# python library bookkeeping
# returns adjacency graph, but neighbors are unordered
def get_neighbor_vertices(d):
    (indices, indptr) = d.vertex_neighbor_vertices
    ns = {}
    for i in range(d.npoints):
        v = (i, tuple(d.points[i]))
        n_inds = indptr[indices[i]:indices[i+1]]
        n_pts = [(i, tuple(d.points[i])) for i in n_inds]
        ns[v] = n_pts
    return ns

def make_oriented_adj(d):
    pts = d.points
    ch = flatten(d.convex_hull)
    return get_neighbor_vertices(d)



# return true iff a is clockwise of b, with c as the center
# if points are collinear ones further from center will be "less"
# positive x axis is "zero"

# https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
def less(a, b, center):
    (ax, ay), (bx, by), (cx, cy) = a, b, center

    # check easy cases first
    if (ax-cx >= 0) and (bx-cx < 0):
        return True
    elif (ax-cx < 0) and (bx-cx >= 0):
        return False
    elif (ax-cx == 0.0) and (bx-cx == 0.0):
        if (ay-cy >= 0) or (by-cy >= 0):
            return ay > by
        return by > ay

    det = (ax-cx)*(by-cy)-(bx-cx)*(ay-cy)

    if det < 0.0:
        return True
    elif det > 0.0:
        return False

    d1 = (ax-cx)*(ax-cx) + (ay-cy)*(ay-cy)
    d2 = (bx-cx)*(bx-cx) + (by-cy)*(by-cy)
    return d1 > d2

# python sorting quirks - requires numerical value
def less_cmp(a, b, center):
    if less(a,b,center):
        return 1.0
    else:
        return -1.0


def sort_neighbors(v1, neighbors):
    v_cmp = lambda a,b: less_cmp(a,b,v1)
    ccw_neighbors = sorted(neighbors, key=cmp_to_key(v_cmp))
    return ccw_neighbors

def plot_tri(d):
    plt.triplot(d.points[:,0], d.points[:,1], d.simplices.copy())
    plt.plot(d.points[:,0], d.points[:,1], 'o')
    for j, p in enumerate(d.points):
        plt.text(p[0]-0.03, p[1]+0.03, j+1, ha='right')
    plt.show()

def write_adj(n):
    d = make_triangulation(n)
    adj = make_oriented_adj(d)

    print(adj)

#    alpha = 0
#    for i in range(n):
#        if i not in ch:
#            alpha = i
#            break
#
#    alpha += 1
#    gamma = adj[alpha-1][0] + 1
#
#    outfile = "NODECOUNT: "+str(n)+"\n"
#    outfile += "GEOMETRY: euclidean\n"
#    outfile += "ALPHA/BETA/GAMMA: "+str(alpha)+" "+str(gamma)+"\n"
#    outfile += "PACKNAME: test.p\n"
#    outfile += "FLOWERS:\n"
#    for i in range(n):
#        if i in ch:
#            outfile += str(i+1) + " " + str(len(adj[i])-1) + "\t"
#            outfile += "".join([str(j+1)+" " for j in adj[i]])
#        else:
#            outfile += str(i+1) + " " + str(len(adj[i])) + "\t"
#            outfile += "".join([str(j+1)+" " for j in adj[i]])+str(adj[i][0]+1)
#        outfile += "\n"
#    outfile += "\nEND"
#    print(outfile)
#    plot_tri(d)
