# -*- coding: iso-8859-1 -*-

#============================================================
#
#           This file is part of FEval, a module for the
#           evaluation of Finite Element results
#
# Licencse: FEval is provided under the GNU General Public License (GPL)
#
# Authors:  Martin Lüthi, tinu@tnoo.net
#
# Homepage: http://feval.sourceforge.net
#
# History:  long, long, most of it in 2000
#           2001.09.21 (ml):  Code cleaned up for intial release
#
# Purpose:  Provide Finite Element shape functions
#
#============================================================


import numpy as N
try:
    import scipy.linalg as LA
except:
    print 'could not import scipy.linalg!'
    pass

## calculate coordinates and weights of Gauss points
## (cf. Numerical Recipies in C++, p.157)
##
## Results are the same as in Bathe 1982 for the first 7 digits 

#============================================================================
# try to use Psyco (psyco.sourceforge.net)
# if configured, this will speed up things considerably
# try:
#     import psyco
#     from psyco.classes import *
# except ImportError:
#     class _psyco:
#         def jit(self):      pass
#         def bind(self, f):  pass
#         def proxy(self, f): return f
#     psyco = _psyco()

class ShapeFunctionPrototype:
    """Defines the prototype of a interpolation function
    cornernodes defines the nodes at the geometrical corner
    We use MARC numbering, i.e. corner nodes first, anticlockwise
    """
    dim, nnodes   = 0, 0   # dimension and number of nodes
    nsides        = 0      # number of element sides
    cornernodes   = N.array([])     # list of corner nodes
    sidenodes     = N.array([])     # list of nodes on the side with index 
    nextnodes     = N.array([])     # list of nodes that are adjecent to the node with index 
    lcoordGauss   = None
    gaussShape    = None
    gaussShapeInv = None

    def __init__(self):
        self.f  = N.zeros(self.nnodes, N.float_)
        self.df = N.zeros( (self.dim, self.nnodes), N.float_)

    def __call__(self, args):
        return self.calcShape(args)

    # must be overrridden by the shape function
    def calcShape(self, lcoord):
        return None

    # must be overrridden by the shape function
    def calcShapeDeriv(self, lcoord):
        return None

    def nextPattern(self, lcoord):
        return [0]

    def calcGauss(self):
        """Calculate the inverse of the shape functions at the Gauss points"""
        a = []
        for lc in self.lcoordGauss:
            a.append(N.array(self.calcShape(lc)))
        self.gaussShape = N.take(N.array(a), self.cornernodes, 1)
        #self.gaussShapeInv = LA.inverse(self.gaussShape)
        self.gaussShapeInv = LA.pinv(self.gaussShape)


class ShapeFunction_Line2(ShapeFunctionPrototype):
    """Element function for linear element defined

    0-------1
    """

    name = 'Line2'
    dim, nnodes = 2, 2
    cornernodes = N.array([0,1])
    nsides      = 2
    sidetype    = 'Point1'
    sidenodes   = N.array(
                  [[0],
                   [1],
                   ])
    nextnodes   = N.array(
                  [[1],
                   [0],
                   ])
    gaussDist = 0.577350269189626  # 1./N.sqrt(3)
    lcoordGauss = N.array([ [-1.],
                            [ 1.],
                            ])*gaussDist
    
    def calcShape(self, lcoord):
        x = lcoord[0]
        return 0.5*N.array([
            1.0-x,
            1.0+x ])

    def calcShapeDeriv(self, lcoord):
        x = lcoord[0]
        return N.array([ -0.5, 0.5 ])

    def nextPattern(self, lcoord):
        x = lcoord * 1.01
        if   x >  1: return [1]
        elif x < -1: return [0]
        else:        return None



class ShapeFunction_Tri3(ShapeFunctionPrototype):
    """Element function for linear element defined in MARC Manual B2.4-1, 
    Taylor&Hughes (1981), p. 49

        2  
       / \  
      /   \ 
     /     \
    0-------1
    """

    name = 'Tri3'
    dim, nnodes = 2, 3
    cornernodes = N.array([0,1,2])
    nsides      = 3
    sidetype    = 'Line2'
    sidenodes   = N.array(
                  [[0,1],
                   [1,2],
                   [2,0],
                   ])
    nextnodes   = N.array(
                  [[1,2],
                   [0,2],
                   [0,1],
                   ])
    triangles   = N.array([[0,1,2]])
    #!!!! worng, this is still from quads
    gaussDist = 0.577350269189626  # 1./N.sqrt(3)
    lcoordGauss = N.array([ [-1., -1.],
                            [ 1., -1.],
                            [-1.,  1.],
                            [ 1.,  1.] ])*gaussDist
    
#     def calcShape(self, lcoord):
#         x, y = lcoord
#         return N.array([
#             1.-x-y,
#             x,
#             y ])
    def calcShape(self, lcoord):
        x, y = lcoord
        # 0.5*(x+1) [-1,1]  ->  x [0,1]
        x = 0.5*(x+1)
        y = 0.5*(y+1)
        return N.array([
            1.-x-y,
            x,
            y ])

    def calcShapeDeriv(self, lcoord):
        x, y = lcoord
        self.df[0,0] = -0.5
        self.df[0,1] =	0.5
        self.df[0,2] =	0.
        self.df[1,0] = -0.5
        self.df[1,1] =  0.
        self.df[1,2] =	0.5
        return self.df

    def nextPattern(self, lcoord):
        x,y = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x+y > 1: return [1,2]
        elif y < 0:   return [0,1]
        elif x < 0:   return [2,0]
        else:        return None


class ShapeFunction_Tri6(ShapeFunctionPrototype):
    """Element function for linear element defined in MARC Manual B2.4-1, 
    Taylor&Hughes (1981), p. 49

        2  
       / \  
      5   4 
     /     \
    0---3---1
    """

    name = 'Tri6'
    dim, nnodes = 2, 6
    cornernodes = N.array([0,1,2])
    nsides      = 3
    sidetype    = 'Line3'
    sidenodes   = N.array(
                  [[0,3,1],
                   [1,4,2],
                   [2,5,0],
                   ])
    nextnodes   = N.array(
                  [[1,2],
                   [0,2],
                   [0,1],
                   ])
    #triangles   = N.array([[0,1,2]])
    #!!!! worng, this is still from quads
    gaussDist = 0.577350269189626  # 1./N.sqrt(3)
    lcoordGauss = N.array([ [-1., -1.],
                            [ 1., -1.],
                            [-1.,  1.],
                            [ 1.,  1.] ])*gaussDist
    
    def calcShape(self, lcoord):
        xi1, xi2 = lcoord

        # 0.5*(x+1) [-1,1]  ->  x [0,1]
        y = 0.5*(xi1+1.)
        z = 0.5*(xi2+1.)
        x = 1. - y - z

        return N.array([
            2.*x*(x-0.5),
            2.*y*(y-0.5),
            2.*z*(z-0.5),
            4.*y*z,
            4.*z*x,
            4.*x*y,
            ])

    def calcShapeDeriv(self, lcoord):
        stop
        xi1, xi2 = lcoord

        # 0.5*(x+1) [-1,1]  ->  x [0,1]
        zeta1 = 0.5*(xi1+1.)
        zeta2 = 0.5*(xi2+1.)
        zeta0 = 1. - zeta1 - zeta2

        self.df[0,0] = 4.*zeta0-1.
        self.df[0,1] = 4.*zeta1-1.
        self.df[0,2] = 4.*zeta2-1.
        self.df[0,3] = 4.*zeta2-1.
        self.df[0,4] = 4.*zeta2-1.
        self.df[0,5] = 4.*zeta2-1.
        self.df[1,0] = -0.5
        self.df[1,1] =  0.
        self.df[1,2] =	0.5
        return self.df

    def nextPattern(self, lcoord):
        x,y = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x+y > 1: return [0,1]
        elif y < 0:   return [2,0]
        elif x < 0:   return [1,2]
        else:        return None


class ShapeFunction_Quad4(ShapeFunctionPrototype):
    """Element function for linear element defined in MARC Manual B2.4-1, 
    Taylor&Hughes (1981), p. 49

    3-------2
    |       |
    |       |
    |       |
    0-------1
    """

    name = 'Quad4'
    dim, nnodes = 2, 4
    cornernodes = N.array([0,1,2,3])
    nsides      = 4
    sidetype    = 'Line2'
    sidenodes   = N.array(
                  [[0,1],
                   [1,2],
                   [2,3],
                   [3,0],
                   ])
    nextnodes   = N.array(
                  [[1,3],
                   [0,2],
                   [1,3],
                   [0,2],
                   ])
    triangles   = N.array([[0,1,3],
                           [1,2,3]])
    gaussDist = 0.577350269189626  # 1./N.sqrt(3)
    lcoordGauss = N.array([ [-1., -1.],
                            [ 1., -1.],
                            [-1.,  1.],
                            [ 1.,  1.] ])*gaussDist
    
    def calcShape(self, lcoord):
        x, y = lcoord
        xy = x*y
        return 0.25*N.array([
            1.0-x-y+xy,
            1.0+x-y-xy,
            1.0+x+y+xy,
            1.0-x+y-xy ])

    def calcShapeDeriv(self, lcoord):
        x, y = lcoord
        return 0.25*N.array([
            [ -1.0+y,  1.0-y, 1.0+y, -1.0-y],
            [ -1.0+x, -1.0-x, 1.0+x, 1.0-x ]])

    def nextPattern(self, lcoord):
        x,y = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x >  1: return [1,2]
        elif x < -1: return [3,0]
        elif y >  1: return [2,3]
        elif y < -1: return [0,1]
        else:        return None


class ShapeFunction_Quad8(ShapeFunctionPrototype):
    """Element function for quadratic element defined in MARC Manual B2.7-1
    Taylor&Hughes (1981), p. 50
    Element nodes numbering is the same as for MARC

    3-----6-----2
    |(5) (6) (7)|
    |           |
    7(3)     (4)5
    |           |
    |(0) (1) (2)|
    0-----4-----1


    """

    name = 'Quad8'
    dim, nnodes = 2, 8
    cornernodes = [0,1,2,3]
    nsides      = 4
    sidetype    = 'Line3'
    sidenodes   = N.array(
                  [[0,4,1],
                   [1,5,2],
                   [2,6,3],
                   [3,7,0],
                   ])
    nextnodes   = N.array(
                  [[1,3],
                   [0,2],
                   [1,3],
                   [0,2],
                   ])
    triangles   = N.array([[7,0,4],
                           [4,1,5],
                           [5,2,6],
                           [6,3,7],
                           [7,4,5],
                           [5,6,7]])
    gaussDist = 0.774596669241483  # = N.sqrt(0.6)
    lcoordGauss = N.array([ [-1., -1.],
                            [ 0., -1.],
                            [ 1., -1.],
                            [-1.,  0.],
                            [ 1.,  0.],
                            [-1.,  1.],
                            [ 0.,  1.],
                            [ 1.,  1.] ])*gaussDist

    def calcShape(self, lcoord):
        x, y = lcoord
        xx, yy, xy = x*x, y*y, x*y
        xxy, xyy = xx*y, x*yy

        return 0.25*N.array([
        # the corner nodes
            (-1.0+xy+xx+yy-xxy-xyy),
            (-1.0-xy+xx+yy-xxy+xyy),
            (-1.0+xy+xx+yy+xxy+xyy),
            (-1.0-xy+xx+yy+xxy-xyy),
        # the mid-side nodes
            2.*(1.0-y-xx+xxy),
            2*(1.0+x-yy-xyy),
            2*(1.0+y-xx-xxy),
            2*(1.0-x-yy+xyy)])
		
    def calcShapeDeriv(self, lcoord):
        x, y = lcoord
        xx, yy, xy = x*x, y*y, x*y
        xxy, xyy, xy2 = xx*y, x*yy, xy*xy
        
        return 0.25*N.array([
            [
        # the corner nodes
            y+xx-xy2-yy,
            y+xx-xy2+yy,
            y+xx+xy2+yy,
            y+xx+xy2-yy,
        # the mid-side nodes
            (-x+xy)*4.,
            (1.0-yy)*2.,
            (-x-xy)*4.,
            (-1.0+yy)*2.,
            ],[
        # the corner nodes
            x+yy-xx-xy2,
            x+yy-xx+xy2,
            x+yy+xx+xy2,
            x+yy+xx-xy2,
        # the mid-side nodes
            (-1.0+xx)*2.,
            (-y-xy)*4.,
            (1.0-xx)*2.,
            (-y+xy)*4.]])
	
    def nextPattern(self, lcoord):
        x,y = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x >  1: return [1,2]
        elif x < -1: return [3,0]
        elif y >  1: return [2,3]
        elif y < -1: return [0,1]
        else:        return None


class ShapeFunction_Quad9(ShapeFunctionPrototype):
    """Element function for quadratic element defined in MARC Manual B2.7-1
    Taylor&Hughes (1981), p. 50
    Element nodes numbering is the same as for MARC

    3-----6-----2
    |(5) (6) (7)|
    |           |
    7(3)  8  (4)5
    |           |
    |(0) (1) (2)|
    0-----4-----1


    """

    name = 'Quad9'
    dim, nnodes = 2, 9
    cornernodes = [0,1,2,3]
    nsides      = 4
    sidetype    = 'Line3'
    sidenodes   = N.array(
                  [[0,4,1],
                   [1,5,2],
                   [2,6,3],
                   [3,7,0],
                   ])
    nextnodes   = N.array(
                  [[1,3],
                   [0,2],
                   [1,3],
                   [0,2],
                   ])
    triangles   = N.array([[7,0,4],
                           [4,1,5],
                           [5,2,6],
                           [6,3,7],
                           [7,4,5],
                           [5,6,7]])
    gaussDist = 0.774596669241483  # = N.sqrt(0.6)
    lcoordGauss = N.array([ [-1., -1.],
                            [ 0., -1.],
                            [ 1., -1.],
                            [-1.,  0.],
                            [ 1.,  0.],
                            [-1.,  1.],
                            [ 0.,  1.],
                            [ 1.,  1.] ])*gaussDist

    def calcShape(self, lcoord):
        print "not implemented correctly"
        stop

        
        x, y = lcoord
        xx, yy, xy = x*x, y*y, x*y
        xxy, xyy = xx*y, x*yy

        return 0.25*N.array([
        # the corner nodes
            (-1.0+xy+xx+yy-xxy-xyy),
            (-1.0-xy+xx+yy-xxy+xyy),
            (-1.0+xy+xx+yy+xxy+xyy),
            (-1.0-xy+xx+yy+xxy-xyy),
        # the mid-side nodes
            2.*(1.0-y-xx+xxy),
            2*(1.0+x-yy-xyy),
            2*(1.0+y-xx-xxy),
            2*(1.0-x-yy+xyy)])
		
    def calcShapeDeriv(self, lcoord):
        print "not implemented correctly"
        stop
        x, y = lcoord
        xx, yy, xy = x*x, y*y, x*y
        xxy, xyy, xy2 = xx*y, x*yy, xy*xy
        
        return 0.25*N.array([
            [
        # the corner nodes
            y+xx-xy2-yy,
            y+xx-xy2+yy,
            y+xx+xy2+yy,
            y+xx+xy2-yy,
        # the mid-side nodes
            (-x+xy)*4.,
            (1.0-yy)*2.,
            (-x-xy)*4.,
            (-1.0+yy)*2.,
            ],[
        # the corner nodes
            x+yy-xx-xy2,
            x+yy-xx+xy2,
            x+yy+xx+xy2,
            x+yy+xx-xy2,
        # the mid-side nodes
            (-1.0+xx)*2.,
            (-y-xy)*4.,
            (1.0-xx)*2.,
            (-y+xy)*4.]])
	
    def nextPattern(self, lcoord):
        x,y = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x >  1: return [1,2]
        elif x < -1: return [3,0]
        elif y >  1: return [2,3]
        elif y < -1: return [0,1]
        else:        return None



class ShapeFunction_Hex8(ShapeFunctionPrototype):
    """Element function for linear element defined in MARC Manual B2.4-1,
    Taylor&Hughes (1981), p. 49

    The integration points (in parentheses) are located at unexpected
    locations (for MARC)!

          7---------------6
         /|(6)        (7)/|
        / |             / |
       /  |            /  |
      /   |           /   |
     / (4)|       (5)/    |
    4---------------5     |
    |     |         |     |
    |     3---------|-----2
    |    / (2)      | (3)/
    |   /           |   /
    |  /            |  /
    | /             | /
    |/ (0)       (1)|/    
    0---------------1

       7-------6
      /|      /|
     / |     / |
    4-------5  | 
    |  3----|--2
    | /     | /
    |/      |/
    0-------1

    """

    name = 'Hex8'
    dim, nnodes = 3, 8
    cornernodes = [0,1,2,3,4,5,6,7]
    nsides      = 6
    sidetype    = 'Quad4'
    sidenodes   = N.array(
                  [[0,3,2,1],
                   [0,1,5,4],
                   [1,2,6,5],
                   [2,3,7,6],
                   [3,0,4,7],
                   [4,5,6,7],
                   ])
    nextnodes   = N.array(
                   [[1,3,4],
                   [0,2,5],
                   [1,3,6],
                   [0,2,7],
                   [0,5,7],
                   [1,4,6],
                   [2,5,7],
                   [3,4,6],
                   ])
    gaussDist = 0.577350269189626  # = 1./N.sqrt(3)
    lcoordGauss = N.array([ [-1., -1., -1.],
                            [ 1., -1., -1.],
                            [-1.,  1., -1.],
                            [ 1.,  1., -1.],
                            [-1., -1.,  1.],
                            [ 1., -1.,  1.],
                            [-1.,  1.,  1.],
                            [ 1.,  1.,  1.]])*gaussDist

    #     lcoordGauss = N.array([ [-1., -1., -1.],
    #                             [ 1., -1., -1.],
    #                             [ 1.,  1., -1.],
    #                             [-1.,  1., -1.],
    #                             [-1., -1.,  1.],
    #                             [ 1., -1.,  1.],
    #                             [ 1.,  1.,  1.],
    #                             [-1.,  1.,  1.]])*gaussDist

    def calcShape(self, lcoord):
        x, y, z = lcoord
        xy, xz, yz = x*y, x*z, y*z
        xyz = x*y*z
        
        return 0.125*N.array([
            1.0-x-y-z+xy+xz+yz-xyz,      # -1,-1,-1,
            1.0+x-y-z-xy-xz+yz+xyz,      #  1,-1,-1,
            1.0+x+y-z+xy-xz-yz-xyz,      #  1, 1,-1,
            1.0-x+y-z-xy+xz-yz+xyz,      # -1, 1,-1,
            1.0-x-y+z+xy-xz-yz+xyz,      # -1,-1, 1,
            1.0+x-y+z-xy+xz-yz-xyz,      #  1,-1, 1,
            1.0+x+y+z+xy+xz+yz+xyz,      #  1, 1, 1,
            1.0-x+y+z-xy-xz+yz-xyz])     # -1, 1, 1,

    def calcShapeDeriv(self, lcoord):
        x, y, z = lcoord
        xy, xz, yz = x*y, x*z, y*z

        self.df[0,0] = -1.0+y+z-yz
        self.df[1,0] = -1.0+x+z-xz
        self.df[2,0] = -1.0+x+y-xy
        self.df[0,1] =	1.0-y-z+yz
        self.df[1,1] = -1.0-x+z+xz
        self.df[2,1] = -1.0-x+y+xy
        self.df[0,2] =	1.0+y-z-yz
        self.df[1,2] =	1.0+x-z-xz
        self.df[2,2] = -1.0-x-y-xy
        self.df[0,3] = -1.0-y+z+yz
        self.df[1,3] =	1.0-x-z+xz
        self.df[2,3] = -1.0+x-y+xy
        self.df[0,4] = -1.0+y-z+yz
        self.df[1,4] = -1.0+x-z+xz
        self.df[2,4] =	1.0-x-y+xy
        self.df[0,5] =	1.0-y+z-yz
        self.df[1,5] = -1.0-x-z-xz
        self.df[2,5] =	1.0+x-y-xy
        self.df[0,6] =	1.0+y+z+yz
        self.df[1,6] =	1.0+x+z+xz
        self.df[2,6] =	1.0+x+y+xy
        self.df[0,7] = -1.0-y-z-yz
        self.df[1,7] =	1.0-x+z-xz
        self.df[2,7] =	1.0-x+y-xy
        self.df = self.df/8.0
        return self.df
    
    def nextPattern(self, lcoord):
        x,y,z = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x >  1: return self.sidenodes[2] #[1,2,6,5]
        elif x < -1: return self.sidenodes[4] #[0,4,7,3]
        elif y >  1: return self.sidenodes[3] #[2,3,7,6]
        elif y < -1: return self.sidenodes[1] #[0,1,5,4]
        elif z >  1: return self.sidenodes[5] #[4,5,6,7]
        elif z < -1: return self.sidenodes[0] #[0,3,2,1]
        else:        return None

class ShapeFunction_Hex20(ShapeFunctionPrototype):
    """Element function for linear element defined in MARC Manual B2.4-1,
    Taylor&Hughes (1981), p. 49

    Here we adopt the numbering from Libmesh, i.e. the second level
    of second order nodes comes befor the 3rd level

    The integration points (in parentheses) are located at unexpected
    locations (for MARC)!

#           7-------14------6
#          /|(6)        (7)/|
#         / |             / |
#        15 |            13 |
#       /  19           /   18
#      / (4)|       (5)/    |
#     4-------12------5     |
#     |     |         |     |
#     |     3------10-|-----2
#     |    / (2)      | (3)/
#    16   /           17  /
#     |  11           |  9
#     | /             | /
#     |/ (0)       (1)|/    
#     0-------8-------1

          7-------18------6
         /|(6)        (7)/|
        / |             / |
       19 |            17 |
      /  15           /   14
     / (4)|       (5)/    |
    4-------16------5     |
    |     |         |     |
    |     3------10-|-----2
    |    / (2)      | (3)/
   12   /           13  /
    |  11           |  9
    | /             | /
    |/ (0)       (1)|/    
    0-------8-------1


16 - 12
17 - 13
18 - 14
19 - 15
12 - 16
13 - 17
14 - 18
15 - 19

    """

    name = 'Hex20'
    dim, nnodes = 3, 20
    cornernodes = [0,1,2,3,4,5,6,7]
    nsides      = 6
    sidetype    = 'Quad8'
    sidenodes   = N.array(
                  [[0,3,2,1,11,10,9,8],       # side 0
                   [0,1,5,4,8, 13, 16, 12],   # side 1
                   [1,2,6,5,9, 14, 17, 13],   # side 2
                   [2,3,7,6,10, 15, 18, 14],  # side 3
                   [3,0,4,7,11, 12, 19, 15],  # side 4
                   [4,5,6,7,16, 17, 18, 19]   # side 5
                   ])
    nextnodes   = N.array(
                  [[1,3,4],
                   [0,2,5],
                   [1,3,6],
                   [0,2,7],
                   [0,5,7],
                   [1,4,6],
                   [2,5,7],
                   [3,4,6],
                   ])
    gaussDist = 0.774596669241483  # = N.sqrt(0.6)
    lcoordGauss = N.array([ [-1., -1., -1.],
                            [ 1., -1., -1.],
                            [-1.,  1., -1.],
                            [ 1.,  1., -1.],
                            [-1., -1.,  1.],
                            [ 1., -1.,  1.],
                            [-1.,  1.,  1.],
                            [ 1.,  1.,  1.]])*gaussDist

    def calcShape(self, lcoord):
        x, y, z = lcoord
        xy, xz, yz = x*y, x*z, y*z
        xx, yy, zz = x*x, y*y, z*z
        xyz, xxy, xxz, xyy = xy*z, xx*y, xx*z, x*yy
        yyz, xzz, yzz = yy*z, x*zz, y*zz
        xxyz, xyyz, xyzz = xxy*z, xyy*z, xyz*z

        self.f[0] = (x+y+z-xyz+xx+yy+zz-xxy-xyy-xxz-xzz-yyz-yzz+ \
                     xxyz+xyyz+xyzz-2.0)/8.0
        self.f[1] = (-x+y+z+xyz+xx+yy+zz-xxy+xyy-xxz+xzz-yyz-yzz+ \
                     xxyz-xyyz-xyzz-2.0)/8.0
        self.f[2] = (-x-y+z-xyz+xx+yy+zz+xxy+xyy-xxz+xzz-yyz+yzz- \
                     xxyz-xyyz+xyzz-2.0)/8.0
        self.f[3] = (x-y+z+xyz+xx+yy+zz+xxy-xyy-xxz-xzz-yyz+yzz- \
                     xxyz+xyyz-xyzz-2.0)/8.0
        self.f[4] = (x+y-z+xyz+xx+yy+zz-xxy-xyy+xxz-xzz+yyz-yzz- \
                     xxyz-xyyz+xyzz-2.0)/8.0
        self.f[5] = (-x+y-z-xyz+xx+yy+zz-xxy+xyy+xxz+xzz+yyz-yzz- \
                     xxyz+xyyz-xyzz-2.0)/8.0
        self.f[6] = (-x-y-z+xyz+xx+yy+zz+xxy+xyy+xxz+xzz+yyz+yzz+ \
                     xxyz+xyyz+xyzz-2.0)/8.0
        self.f[7] = (x-y-z-xyz+xx+yy+zz+xxy-xyy+xxz-xzz+yyz+yzz+ \
                     xxyz-xyyz-xyzz-2.0)/8.0

        self.f[8]  = (1.0-z-y+yz-xx+xxz+xxy-xxyz)/4.0
        self.f[9]  = (1.0-z+x-xz-yy+yyz-xyy+xyyz)/4.0
        self.f[10] = (1.0-z+y-yz-xx+xxz-xxy+xxyz)/4.0
        self.f[11] = (1.0-z-x+xz-yy+yyz+xyy-xyyz)/4.0
        self.f[16] = (1.0+z-y-yz-xx-xxz+xxy+xxyz)/4.0
        self.f[17] = (1.0+z+x+xz-yy-yyz-xyy-xyyz)/4.0
        self.f[18] = (1.0+z+y+yz-xx-xxz-xxy-xxyz)/4.0
        self.f[19] = (1.0+z-x-xz-yy-yyz+xyy+xyyz)/4.0
        self.f[12] = (1.0-y-x+xy-zz+yzz+xzz-xyzz)/4.0
        self.f[13] = (1.0-y+x-xy-zz+yzz-xzz+xyzz)/4.0
        self.f[14] = (1.0+y+x+xy-zz-yzz-xzz-xyzz)/4.0
        self.f[15] = (1.0+y-x-xy-zz-yzz+xzz+xyzz)/4.0
        return self.f

    def calcShapeDeriv(self, lcoord):
        x, y, z = lcoord
        xy, xz, yz = x*y, x*z, y*z
        xx, yy, zz = x*x, y*y, z*z
        xyz, xxy, xxz, xyy = xy*z, xx*y, xx*z, x*yy
        yyz, xzz, yzz = yy*z, x*zz, y*zz

        self.df[0, 0] =	 1.0-yz+2.0*x-2.0*xy-yy-2.0*xz-zz+2.0*xyz+yyz+yzz
        self.df[1, 0] =	 1.0-xz+2.0*y-xx-2.0*xy-2.0*yz-zz+xxz+2.0*xyz+xzz
        self.df[2, 0] =	 1.0-xy+2.0*z-xx-2.0*xz-yy-2.0*yz+xxy+xyy+2.0*xyz
        self.df[0, 1] = -1.0+yz+2.0*x-2.0*xy+yy-2.0*xz+zz+2.0*xyz-yyz-yzz
        self.df[1, 1] =	 1.0+xz+2.0*y-xx+2.0*xy-2.0*yz-zz+xxz-2.0*xyz-xzz
        self.df[2, 1] =	 1.0+xy+2.0*z-xx+2.0*xz-yy-2.0*yz+xxy-xyy-2.0*xyz
        self.df[0, 2] = -1.0-yz+2.0*x+2.0*xy+yy-2.0*xz+zz-2.0*xyz-yyz+yzz
        self.df[1, 2] = -1.0-xz+2.0*y+xx+2.0*xy-2.0*yz+zz-xxz-2.0*xyz+xzz
        self.df[2, 2] =	 1.0-xy+2.0*z-xx+2.0*xz-yy+2.0*yz-xxy-xyy+2.0*xyz
        self.df[0, 3] =	 1.0+yz+2.0*x+2.0*xy-yy-2.0*xz-zz-2.0*xyz+yyz-yzz
        self.df[1, 3] = -1.0+xz+2.0*y+xx-2.0*xy-2.0*yz+zz-xxz+2.0*xyz-xzz
        self.df[2, 3] =	 1.0+xy+2.0*z-xx-2.0*xz-yy+2.0*yz-xxy+xyy-2.0*xyz
        self.df[0, 4] =	 1.0+yz+2.0*x-2.0*xy-yy+2.0*xz-zz-2.0*xyz-yyz+yzz
        self.df[1, 4] =	 1.0+xz+2.0*y-xx-2.0*xy+2.0*yz-zz-xxz-2.0*xyz+xzz
        self.df[2, 4] = -1.0+xy+2.0*z+xx-2.0*xz+yy-2.0*yz-xxy-xyy+2.0*xyz
        self.df[0, 5] = -1.0-yz+2.0*x-2.0*xy+yy+2.0*xz+zz-2.0*xyz+yyz-yzz
        self.df[1, 5] =	 1.0-xz+2.0*y-xx+2.0*xy+2.0*yz-zz-xxz+2.0*xyz-xzz
        self.df[2, 5] = -1.0-xy+2.0*z+xx+2.0*xz+yy-2.0*yz-xxy+xyy-2.0*xyz
        self.df[0, 6] = -1.0+yz+2.0*x+2.0*xy+yy+2.0*xz+zz+2.0*xyz+yyz+yzz
        self.df[1, 6] = -1.0+xz+2.0*y+xx+2.0*xy+2.0*yz+zz+xxz+2.0*xyz+xzz
        self.df[2, 6] = -1.0+xy+2.0*z+xx+2.0*xz+yy+2.0*yz+xxy+xyy+2.0*xyz
        self.df[0, 7] =	 1.0-yz+2.0*x+2.0*xy-yy+2.0*xz-zz+2.0*xyz-yyz-yzz
        self.df[1, 7] = -1.0-xz+2.0*y+xx-2.0*xy+2.0*yz+zz+xxz-2.0*xyz-xzz
        self.df[2, 7] = -1.0-xy+2.0*z+xx-2.0*xz+yy+2.0*yz+xxy-xyy-2.0*xyz
        self.df[:, 0:8] = self.df[:, 0:8]/2.0
        
        self.df[0, 8] = -2.0*x+2.0*xz+2.0*xy-2.0*xyz
        self.df[1, 8] = -1.0+z+xx-xxz
        self.df[2, 8] = -1.0+y+xx-xxy
        self.df[0, 9] =	 1.0-z-yy+yyz
        self.df[1, 9] = -2.0*y+2.0*yz-2.0*xy+2.0*xyz
        self.df[2, 9] = -1.0-x+yy+xyy
        self.df[0, 10] = -2.0*x+2.0*xz-2.0*xy+2.0*xyz
        self.df[1, 10] =  1.0-z-xx+xxz
        self.df[2, 10] = -1.0-y+xx+xxy
        self.df[0, 11] = -1.0+z+yy-yyz
        self.df[1, 11] = -2.0*y+2.0*yz+2.0*xy-2.0*xyz
        self.df[2, 11] = -1.0+x+yy-xyy
        self.df[0, 16] = -2*x-2*xz+2*xy+2*xyz
        self.df[1, 16] = -1.0-z+xx+xxz
        self.df[2, 16] =  1.0-y-xx+xxy
        self.df[0, 17] =  1.0+z-yy-yyz
        self.df[1, 17] = -2*y-2*yz-2*xy-2*xyz
        self.df[2, 17] =  1.0+x-yy-xyy
        self.df[0, 18] = -2*x-2*xz-2*xy-2*xyz
        self.df[1, 18] =  1.0+z-xx-xxz
        self.df[2, 18] =  1.0+y-xx-xxy
        self.df[0, 19] = -1.0-z+yy+yyz
        self.df[1, 19] = -2*y-2*yz+2*xy+2*xyz
        self.df[2, 19] =  1.0-x-yy+xyy
        self.df[0, 12] = -1.0+y+zz-yzz
        self.df[1, 12] = -1.0+x+zz-xzz
        self.df[2, 12] = -2*z+2*yz+2*xz-2*xyz
        self.df[0, 13] =  1.0-y-zz+yzz
        self.df[1, 13] = -1.0-x+zz+xzz
        self.df[2, 13] = -2*z+2*yz-2*xz+2*xyz
        self.df[0, 14] =  1.0+y-zz-yzz
        self.df[1, 14] =  1.0+x-zz-xzz
        self.df[2, 14] = -2*z-2*yz-2*xz-2*xyz
        self.df[0, 15] = -1.0-y+zz+yzz
        self.df[1, 15] =  1.0-x-zz+xzz
        self.df[2, 15] =  -2*z-2*yz+2*xz+2*xyz
        self.df = self.df/4.0
        return self.df

    def nextPattern(self, lcoord):
        x,y,z = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x >  1: return [1,2,6,5]
        elif x < -1: return [0,3,7,4]
        elif y >  1: return [2,3,7,6]
        elif y < -1: return [0,1,5,4]
        elif z >  1: return [4,5,6,7]
        elif z < -1: return [0,1,2,3]
        else:        return None


        
class ShapeFunction_Hex27(ShapeFunctionPrototype):
    """Element function for linear element defined in MARC Manual B2.4-1,
    Taylor&Hughes (1981), p. 49

    Here we adopt the numbering from Libmesh, i.e. the second level
    of second order nodes comes before the 3rd level

    The integration points (in parentheses) are located at unexpected
    locations (for MARC)!

          7-------18------6
         /|(6)        (7)/|
        / |             / |
       19 |   [25]     17 |
      /  15    [23]   /   14     center node: 26
     / (4)|       (5)/    |
    4-------16------5     |
    | [24]|         | [22]|
    |     3------10-|-----2
    |    / (2)      | (3)/
   12   /  [21]     13  /
    |  11    [20]   |  9
    | /             | /
    |/ (0)       (1)|/    
    0-------8-------1

    """

    name = 'Hex27'
    dim, nnodes = 3, 27
    cornernodes = [0,1,2,3,4,5,6,7]
    nsides      = 6
    sidetype    = 'Quad9'
    sidenodes   = N.array([
        [0, 3, 2, 1, 11, 10,  9,  8, 20], # Side 0   (exodus: 5)   20 -> 22
        [0, 1, 5, 4,  8, 13, 16, 12, 21], # Side 1   (exodus: 1)   21 -> 26
        [1, 2, 6, 5,  9, 14, 17, 13, 22], # Side 2   (exodus: 2)   22 -> 25
        [2, 3, 7, 6, 10, 15, 18, 14, 23], # Side 3   (exodus: 3)   23 -> 27
        [3, 0, 4, 7, 11, 12, 19, 15, 24], # Side 4   (exodus: 4)   24 -> 24
        [4, 5, 6, 7, 16, 17, 18, 19, 25]  # Side 5   (exodus: 6)   25 -> 23
        ])
    nextnodes   = N.array(
                  [[1,3,4],
                   [0,2,5],
                   [1,3,6],
                   [0,2,7],
                   [0,5,7],
                   [1,4,6],
                   [2,5,7],
                   [3,4,6],
                   ])
    gaussDist = 0.774596669241483  # = N.sqrt(0.6)
    lcoordGauss = N.array([ [-1., -1., -1.],
                            [ 1., -1., -1.],
                            [-1.,  1., -1.],
                            [ 1.,  1., -1.],
                            [-1., -1.,  1.],
                            [ 1., -1.,  1.],
                            [-1.,  1.,  1.],
                            [ 1.,  1.,  1.]])*gaussDist

    def calcShape(self, lcoord):
        print 'not implemented'
        return None
#         x, y, z = lcoord
#         xy, xz, yz = x*y, x*z, y*z
#         xx, yy, zz = x*x, y*y, z*z
#         xyz, xxy, xxz, xyy = xy*z, xx*y, xx*z, x*yy
#         yyz, xzz, yzz = yy*z, x*zz, y*zz
#         xxyz, xyyz, xyzz = xxy*z, xyy*z, xyz*z
		
#         self.f[0] = (x+y+z-xyz+xx+yy+zz-xxy-xyy-xxz-xzz-yyz-yzz+ \
#                      xxyz+xyyz+xyzz-2.0)/8.0
#         self.f[1] = (-x+y+z+xyz+xx+yy+zz-xxy+xyy-xxz+xzz-yyz-yzz+ \
#                      xxyz-xyyz-xyzz-2.0)/8.0
#         self.f[2] = (-x-y+z-xyz+xx+yy+zz+xxy+xyy-xxz+xzz-yyz+yzz- \
#                      xxyz-xyyz+xyzz-2.0)/8.0
#         self.f[3] = (x-y+z+xyz+xx+yy+zz+xxy-xyy-xxz-xzz-yyz+yzz- \
#                      xxyz+xyyz-xyzz-2.0)/8.0
#         self.f[4] = (x+y-z+xyz+xx+yy+zz-xxy-xyy+xxz-xzz+yyz-yzz- \
#                      xxyz-xyyz+xyzz-2.0)/8.0
#         self.f[5] = (-x+y-z-xyz+xx+yy+zz-xxy+xyy+xxz+xzz+yyz-yzz- \
#                      xxyz+xyyz-xyzz-2.0)/8.0
#         self.f[6] = (-x-y-z+xyz+xx+yy+zz+xxy+xyy+xxz+xzz+yyz+yzz+ \
#                      xxyz+xyyz+xyzz-2.0)/8.0
#         self.f[7] = (x-y-z-xyz+xx+yy+zz+xxy-xyy+xxz-xzz+yyz+yzz+ \
#                      xxyz-xyyz-xyzz-2.0)/8.0

#         self.f[8]  = (1.0-z-y+yz-xx+xxz+xxy-xxyz)/4.0
#         self.f[9]  = (1.0-z+x-xz-yy+yyz-xyy+xyyz)/4.0
#         self.f[10] = (1.0-z+y-yz-xx+xxz-xxy+xxyz)/4.0
#         self.f[11] = (1.0-z-x+xz-yy+yyz+xyy-xyyz)/4.0
#         self.f[12] = (1.0+z-y-yz-xx-xxz+xxy+xxyz)/4.0
#         self.f[13] = (1.0+z+x+xz-yy-yyz-xyy-xyyz)/4.0
#         self.f[14] = (1.0+z+y+yz-xx-xxz-xxy-xxyz)/4.0
#         self.f[15] = (1.0+z-x-xz-yy-yyz+xyy+xyyz)/4.0
#         self.f[16] = (1.0-y-x+xy-zz+yzz+xzz-xyzz)/4.0
#         self.f[17] = (1.0-y+x-xy-zz+yzz-xzz+xyzz)/4.0
#         self.f[18] = (1.0+y+x+xy-zz-yzz-xzz-xyzz)/4.0
#         self.f[19] = (1.0+y-x-xy-zz-yzz+xzz+xyzz)/4.0
#         return self.f

    def calcShapeDeriv(self, lcoord):
        print 'not implemented'
        return None
#         x, y, z = lcoord
#         xy, xz, yz = x*y, x*z, y*z
#         xx, yy, zz = x*x, y*y, z*z
#         xyz, xxy, xxz, xyy = xy*z, xx*y, xx*z, x*yy
#         yyz, xzz, yzz = yy*z, x*zz, y*zz

#         self.df[0, 0] =	 1.0-yz+2.0*x-2.0*xy-yy-2.0*xz-zz+2.0*xyz+yyz+yzz
#         self.df[1, 0] =	 1.0-xz+2.0*y-xx-2.0*xy-2.0*yz-zz+xxz+2.0*xyz+xzz
#         self.df[2, 0] =	 1.0-xy+2.0*z-xx-2.0*xz-yy-2.0*yz+xxy+xyy+2.0*xyz
#         self.df[0, 1] = -1.0+yz+2.0*x-2.0*xy+yy-2.0*xz+zz+2.0*xyz-yyz-yzz
#         self.df[1, 1] =	 1.0+xz+2.0*y-xx+2.0*xy-2.0*yz-zz+xxz-2.0*xyz-xzz
#         self.df[2, 1] =	 1.0+xy+2.0*z-xx+2.0*xz-yy-2.0*yz+xxy-xyy-2.0*xyz
#         self.df[0, 2] = -1.0-yz+2.0*x+2.0*xy+yy-2.0*xz+zz-2.0*xyz-yyz+yzz
#         self.df[1, 2] = -1.0-xz+2.0*y+xx+2.0*xy-2.0*yz+zz-xxz-2.0*xyz+xzz
#         self.df[2, 2] =	 1.0-xy+2.0*z-xx+2.0*xz-yy+2.0*yz-xxy-xyy+2.0*xyz
#         self.df[0, 3] =	 1.0+yz+2.0*x+2.0*xy-yy-2.0*xz-zz-2.0*xyz+yyz-yzz
#         self.df[1, 3] = -1.0+xz+2.0*y+xx-2.0*xy-2.0*yz+zz-xxz+2.0*xyz-xzz
#         self.df[2, 3] =	 1.0+xy+2.0*z-xx-2.0*xz-yy+2.0*yz-xxy+xyy-2.0*xyz
#         self.df[0, 4] =	 1.0+yz+2.0*x-2.0*xy-yy+2.0*xz-zz-2.0*xyz-yyz+yzz
#         self.df[1, 4] =	 1.0+xz+2.0*y-xx-2.0*xy+2.0*yz-zz-xxz-2.0*xyz+xzz
#         self.df[2, 4] = -1.0+xy+2.0*z+xx-2.0*xz+yy-2.0*yz-xxy-xyy+2.0*xyz
#         self.df[0, 5] = -1.0-yz+2.0*x-2.0*xy+yy+2.0*xz+zz-2.0*xyz+yyz-yzz
#         self.df[1, 5] =	 1.0-xz+2.0*y-xx+2.0*xy+2.0*yz-zz-xxz+2.0*xyz-xzz
#         self.df[2, 5] = -1.0-xy+2.0*z+xx+2.0*xz+yy-2.0*yz-xxy+xyy-2.0*xyz
#         self.df[0, 6] = -1.0+yz+2.0*x+2.0*xy+yy+2.0*xz+zz+2.0*xyz+yyz+yzz
#         self.df[1, 6] = -1.0+xz+2.0*y+xx+2.0*xy+2.0*yz+zz+xxz+2.0*xyz+xzz
#         self.df[2, 6] = -1.0+xy+2.0*z+xx+2.0*xz+yy+2.0*yz+xxy+xyy+2.0*xyz
#         self.df[0, 7] =	 1.0-yz+2.0*x+2.0*xy-yy+2.0*xz-zz+2.0*xyz-yyz-yzz
#         self.df[1, 7] = -1.0-xz+2.0*y+xx-2.0*xy+2.0*yz+zz+xxz-2.0*xyz-xzz
#         self.df[2, 7] = -1.0-xy+2.0*z+xx-2.0*xz+yy+2.0*yz+xxy-xyy-2.0*xyz
#         self.df[:, 0:8] = self.df[:, 0:8]/2.0
        
#         self.df[0, 8] = -2.0*x+2.0*xz+2.0*xy-2.0*xyz
#         self.df[1, 8] = -1.0+z+xx-xxz
#         self.df[2, 8] = -1.0+y+xx-xxy
#         self.df[0, 9] =	 1.0-z-yy+yyz
#         self.df[1, 9] = -2.0*y+2.0*yz-2.0*xy+2.0*xyz
#         self.df[2, 9] = -1.0-x+yy+xyy
#         self.df[0, 10] = -2.0*x+2.0*xz-2.0*xy+2.0*xyz
#         self.df[1, 10] =  1.0-z-xx+xxz
#         self.df[2, 10] = -1.0-y+xx+xxy
#         self.df[0, 11] = -1.0+z+yy-yyz
#         self.df[1, 11] = -2.0*y+2.0*yz+2.0*xy-2.0*xyz
#         self.df[2, 11] = -1.0+x+yy-xyy
#         self.df[0, 12] = -2*x-2*xz+2*xy+2*xyz
#         self.df[1, 12] = -1.0-z+xx+xxz
#         self.df[2, 12] =  1.0-y-xx+xxy
#         self.df[0, 13] =  1.0+z-yy-yyz
#         self.df[1, 13] = -2*y-2*yz-2*xy-2*xyz
#         self.df[2, 13] =  1.0+x-yy-xyy
#         self.df[0, 14] = -2*x-2*xz-2*xy-2*xyz
#         self.df[1, 14] =  1.0+z-xx-xxz
#         self.df[2, 14] =  1.0+y-xx-xxy
#         self.df[0, 15] = -1.0-z+yy+yyz
#         self.df[1, 15] = -2*y-2*yz+2*xy+2*xyz
#         self.df[2, 15] =  1.0-x-yy+xyy
#         self.df[0, 16] = -1.0+y+zz-yzz
#         self.df[1, 16] = -1.0+x+zz-xzz
#         self.df[2, 16] = -2*z+2*yz+2*xz-2*xyz
#         self.df[0, 17] =  1.0-y-zz+yzz
#         self.df[1, 17] = -1.0-x+zz+xzz
#         self.df[2, 17] = -2*z+2*yz-2*xz+2*xyz
#         self.df[0, 18] =  1.0+y-zz-yzz
#         self.df[1, 18] =  1.0+x-zz-xzz
#         self.df[2, 18] = -2*z-2*yz-2*xz-2*xyz
#         self.df[0, 19] = -1.0-y+zz+yzz
#         self.df[1, 19] =  1.0-x-zz+xzz
#         self.df[2, 19] =  -2*z-2*yz+2*xz+2*xyz
#         self.df = self.df/4.0
#         return self.df

    def nextPattern(self, lcoord):
        x,y,z = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x >  1: return [1,2,6,5]
        elif x < -1: return [0,3,7,4]
        elif y >  1: return [2,3,7,6]
        elif y < -1: return [0,1,5,4]
        elif z >  1: return [4,5,6,7]
        elif z < -1: return [0,1,2,3]
        else:        return None
# all shape functions are registered here

shapeFunctions = {
    'Line2': ShapeFunction_Line2,
    'Tri3':  ShapeFunction_Tri3,
    'Tri6':  ShapeFunction_Tri6,
    'Quad4': ShapeFunction_Quad4,
    'Quad8': ShapeFunction_Quad8,
    'Quad9': ShapeFunction_Quad9,
    'Hex8' : ShapeFunction_Hex8, 
    'Hex20': ShapeFunction_Hex20,
    'Hex27': ShapeFunction_Hex27
   }

if __name__ == '__main__':

    sh6 = ShapeFunction_Tri6()
    sh3 = ShapeFunction_Tri3()

    def shape(zeta):
        zeta1, zeta2 = zeta
        zeta0 = 1. - zeta1 - zeta2
        return [2.*zeta0*(zeta0-0.5),
                2.*zeta1*(zeta1-0.5),
                2.*zeta2*(zeta2-0.5),
                4.*zeta1*zeta2,
                4.*zeta2*zeta0,
                4.*zeta0*zeta1,
                ]

    print shape([0.,0.])
    print sh6([-1.,-1.])
    print sh3([-1.,-1.])
    print '----------------------'
    print shape([1.,0.])
    print sh6([1.,-1.])
    print sh3([1.,-1.])
    print '----------------------'
    print shape([0.,1.])
    print sh6([-1.,1.])
    print sh3([-1.,1.])
    print '----------------------'
    print shape([0.5,0.5])
    print sh6([0.,0.])
    print sh3([0.,0.])
    print '----------------------'
    print shape([0.,0.5])
    print sh6([-1.,0.])
    print sh3([-1.,0.])
    print '----------------------'
    print shape([0.5,0.])
    print sh6([0.,-1.])
    print sh3([0.,-1.])
    print '----------------------'
    print shape([0.3,0.4])
    print sh6([-0.4,-0.2])
    print sh3([-0.4,-0.2])


#     for n, sf in shapeFunctions.items():
#         print '===== %s =====' % n
#         s = sf()
#         s.calcGauss()
#         print s.gaussShapeInv
    

