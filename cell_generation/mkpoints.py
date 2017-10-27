import numpy as np
from scipy.spatial import Delaunay, ConvexHull
import matplotlib.pyplot as plt
from itertools import combinations
from copy import copy
import random

def random_points(n):
    random.seed()
    pts = [(random.random(), random.random()) for i in range(n)]
    return pts

def make_triangulation(n):
    pts = random_points(n)
    ch = ConvexHull(pts).vertices
    return Delaunay(pts), ch

def make_oriented_adj(n, d, ch):
#    (ind, indptr) = d.vertex_neighbor_vertices
#    adj = [indptr[ind[k]:ind[k+1]] for k in range(n)]
#    return adj
    adj = {}
    # initialize convex hull
    for (i,j,k) in list(zip(ch, ch[1:],ch[2:])) + [(ch[-2],ch[-1],ch[0]),
                                                   (ch[-1],ch[0],ch[1])]:
        adj[j] = [k,i]

    for [i,j,k] in d.simplices:
        if (i or j or k) not in ch:
            adj = check_triangle((i,j,k), adj)
            adj = check_triangle((j,k,i), adj)
            adj = check_triangle((k,i,j), adj)
    return adj

# TODO: there has to be some better algebraic structure here
# use CCW ordering of simplices to order neighbors
# since scipy unhelpfully discards this info
def check_triangle(tri, adj):
    (i,j,k) = tri
    if i not in adj:
        adj[i] = [j,k]
    else:
        neighbors = adj[i]
        if (j not in adj[i]) and (k not in adj[i]):
            adj[i] += [j,k]
        elif (j in adj[i]) and (k not in adj[i]):
            new_neighbors = []
            for n in adj[i]:
                if n == j: # match j -> insert k after
                    new_neighbors.append(j)
                    new_neighbors.append(k)
                else:
                    new_neighbors.append(n)
            adj[i] = new_neighbors
        elif (k in adj[i]) and (j not in adj[i]):
            new_neighbors = []
            for n in adj[i]:
                if n == k: # match k -> insert j before
                    new_neighbors.append(j)
                    new_neighbors.append(k)
                else:
                    new_neighbors.append(n)
            adj[i] = new_neighbors
    return adj


def plot_tri(d):
    plt.triplot(d.points[:,0], d.points[:,1], d.simplices.copy())
    plt.plot(d.points[:,0], d.points[:,1], 'o')
    for j, p in enumerate(d.points):
        plt.text(p[0]-0.03, p[1]+0.03, j+1, ha='right')

    plt.show()

def write_adj(n):
    d, ch = make_triangulation(n)
    adj = make_oriented_adj(n,d,ch)

    alpha = 0
    for i in range(n):
        if i not in ch:
            alpha = i
            break

    alpha += 1
    gamma = adj[alpha-1][0] + 1

    outfile = "NODECOUNT: "+str(n)+"\n"
    outfile += "GEOMETRY: euclidean\n"
    outfile += "ALPHA/BETA/GAMMA: "+str(alpha)+" "+str(gamma)+"\n"
    outfile += "PACKNAME: test.p\n"
    outfile += "FLOWERS:\n"
    for i in range(n):
        if i in ch:
            outfile += str(i+1) + " " + str(len(adj[i])-1) + "\t"
            outfile += "".join([str(j+1)+" " for j in adj[i]])
        else:
            outfile += str(i+1) + " " + str(len(adj[i])) + "\t"
            outfile += "".join([str(j+1)+" " for j in adj[i]])+str(adj[i][0]+1)
        outfile += "\n"
    outfile += "\nEND"
    print(outfile)
    plot_tri(d)
