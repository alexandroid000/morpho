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

def rotate(l, i):
        return l[i:] + l[:i]


# return true iff a is clockwise of b, with c as the center
# if points are collinear ones further from center will be "less"

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
    (i, an), (j, bn), cn = a, b, center
    if less(an, bn, cn):
        return 1.0
    else:
        return -1.0

# sort_neighbors :: Point -> [(Int, Point)] -> [(Int, Point)]
def sort_neighbors(v1, neighbors):
    v_cmp = lambda a,b: less_cmp(a,b,v1)
    ccw_neighbors = sorted(neighbors, key=cmp_to_key(v_cmp))
    return ccw_neighbors

# python library bookkeeping
# returns adjacency graph, with neighbors sorted ccw
# Delaunay -> [Point, [(Int, Point)]
def make_oriented_adj(d):
    (indices, indptr) = d.vertex_neighbor_vertices
    ch = d.convex_hull
    ch_vis = np.unique(ch) # indices of ch verts
    ns = [0.0]*d.npoints

    for i in range(d.npoints):
        v = tuple(d.points[i])
        n_is = indptr[indices[i]:indices[i+1]] # indices of neighbors, unsorted
        n_pts = sort_neighbors(v, [(i, tuple(d.points[i])) for i in n_is])

        if (i in ch_vis) and (d.npoints != 3):
            print(i+1, " is in ch")
            ch_neighbors = [(j, tuple(d.points[j])) for j in n_is if j in ch_vis]
            ch_ns = [k for (k,pt) in sort_neighbors(d.points[i], ch_neighbors)]
            n1, n2 = ch_ns[0], ch_ns[-1]
            print(i+1, " neighbors are ", n1+1, n2+1)

            # rotate until "smaller" (more cw) point is first in list
            a, b = d.points[n1], d.points[n2]
            if less(a, b, d.points[i]):
                print("before rotation: ",[k for (k,pt) in n_pts])
                print("rotate by ", n1)
                n_pts = rotate(n_pts, n1)
                print("after rotation: ",[k for (k,pt) in n_pts])
            else:
                print("before rotation: ",[k for (k,pt) in n_pts])
                print("rotate by ", n2)
                n_pts = rotate(n_pts, n2)
                print("after rotation: ",[k for (k,pt) in n_pts])

        ns[i] = n_pts
    return ns


d = make_triangulation(10)
adj = make_oriented_adj(d)


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

    for i in range(n):
        ns = [ip+1 for (ip, pts) in adj[i]]
        print(i+1)
        print("\t",ns)

    plot_tri(d)
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
