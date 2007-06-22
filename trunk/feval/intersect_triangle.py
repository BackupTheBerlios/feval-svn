# fast intersection of a ray with a triangle
# code from http://jgt.akpeters.com/papers/MollerTrumbore97/code.html
# adapted to numpy
#
# @article{MollerTrumbore97,
#   author = "Tomas Möller and Ben Trumbore",
#   title = "Fast, Minimum Storage Ray-Triangle Intersection",
#   journal = "journal of graphics tools",
#   volume = "2",
#   number = "1",
#   pages = "21-28",
#   year = "1997",
# }


import numpy as N
EPSILON = 0.000001

def intersect_triangle(orig, dir, vert0, vert1, vert2):

    # find vectors for two edges sharing vert0
    edge1= vert1- vert0
    edge2= vert2- vert0

    # begin calculating determinant - also used to calculate U parameter
    pvec = N.cross(dir, edge2)

    # if determinant is near zero, ray lies in plane of triangle 
    det = N.dot(edge1, pvec)

    if (det > -EPSILON and det < EPSILON):
        return None
    inv_det = 1.0 / det

    # calculate distance from vert0 to ray origin 
    tvec = orig - vert0

    # calculate U parameter and test bounds
    u = N.dot(tvec, pvec) * inv_det
    if (u < 0.0 or u > 1.0):
        return None

    # prepare to test V parameter
    qvec = N.cross(tvec, edge1)

    # calculate V parameter and test bounds
    v = N.dot(dir, qvec) * inv_det
    if (v < 0.0 or u + v > 1.0):
        return None

    # calculate t, ray intersects triangle
    t = N.dot(edge2, qvec) * inv_det

    return t, u, v

if __name__ == "__main__":

    x1 = N.array((2, 0, 1))
    x2 = N.array((3, 5, 0))
    x3 = N.array((5, 3, 0))

    x0 = N.array((4, 3.5, 0))
    
    print intersect_triangle(N.zeros(3), N.array([3.5, 5, 0]), x1, x2, x3)
