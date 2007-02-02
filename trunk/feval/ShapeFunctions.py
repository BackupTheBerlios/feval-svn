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

    def calcShape(self, lcoord):
        return None

    def calcShapeDeriv(self, lcoord):
        return None

    def evalShape(self, lcoord):
        self.calcShape(lcoord)
        return self.f

    def evalShapeDeriv(self, lcoord):
        self.calcShapeDeriv(lcoord)
        return self.df

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
    #!!!! worng
    gaussDist = 0.577350269189626  # 1./N.sqrt(3)
    lcoordGauss = N.array([ [-1., -1.],
                            [ 1., -1.],
                            [-1.,  1.],
                            [ 1.,  1.] ])*gaussDist
    
    def calcShape(self, lcoord):
        x, y = lcoord
        xy = x*y
        self.f[0] = 1.0-x-y+xy
        self.f[1] = 1.0+x-y-xy
        self.f[2] = 1.0+x+y+xy
        self.f[3] = 1.0-x+y-xy
        self.f = self.f*0.25
        return self.f

    def calcShapeDeriv(self, lcoord):
        x, y = lcoord
        self.df[0,0] = -1.0+y
        self.df[0,1] =	1.0-y
        self.df[0,2] =	1.0+y
        self.df[0,3] = -1.0-y
        self.df[1,0] = -1.0+x
        self.df[1,1] = -1.0-x
        self.df[1,2] =	1.0+x
        self.df[1,3] =	1.0-x
        self.df = self.df*0.25
        return self.df

    def nextPattern(self, lcoord):
        x,y = lcoord / max(N.absolute(lcoord)) * 1.01
        if   x >  1: return [1,2]
        elif x < -1: return [3,0]
        elif y >  1: return [2,3]
        elif y < -1: return [0,1]
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
        self.f[0] = 1.0-x-y+xy
        self.f[1] = 1.0+x-y-xy
        self.f[2] = 1.0+x+y+xy
        self.f[3] = 1.0-x+y-xy
        self.f = self.f*0.25
        return self.f

    def calcShapeDeriv(self, lcoord):
        x, y = lcoord
        self.df[0,0] = -1.0+y
        self.df[0,1] =	1.0-y
        self.df[0,2] =	1.0+y
        self.df[0,3] = -1.0-y
        self.df[1,0] = -1.0+x
        self.df[1,1] = -1.0-x
        self.df[1,2] =	1.0+x
        self.df[1,3] =	1.0-x
        self.df = self.df*0.25
        return self.df

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

        # the corner nodes
        self.f[0] = 0.25*(-1.0+xy+xx+yy-xxy-xyy)
        self.f[1] = 0.25*(-1.0-xy+xx+yy-xxy+xyy)
        self.f[2] = 0.25*(-1.0+xy+xx+yy+xxy+xyy)
        self.f[3] = 0.25*(-1.0-xy+xx+yy+xxy-xyy)
        # the mid-side nodes
        self.f[4] = 0.5*(1.0-y-xx+xxy)
        self.f[5] = 0.5*(1.0+x-yy-xyy)
        self.f[6] = 0.5*(1.0+y-xx-xxy)
        self.f[7] = 0.5*(1.0-x-yy+xyy)
        return self.f
		
    def calcShapeDeriv(self, lcoord):
        x, y = lcoord
        xx, yy, xy = x*x, y*y, x*y
        xxy, xyy, xy2 = xx*y, x*yy, xy*xy
        
        # the corner nodes
        self.df[0,0] = y+xx-xy2-yy
        self.df[0,1] =-y+xx-xy2+yy
        self.df[0,2] = y+xx+xy2+yy
        self.df[0,3] =-y+xx+xy2-yy
        self.df[1,0] = x+yy-xx-xy2
        self.df[1,1] =-x+yy-xx+xy2
        self.df[1,2] = x+yy+xx+xy2
        self.df[1,3] =-x+yy+xx-xy2
        self.df[:,0:4] = self.df[:,0:4]*0.25
        
        # the mid-side nodes
        self.df[0,4] = -x+xy
        self.df[0,5] = (1.0-yy)*0.5
        self.df[0,6] = -x-xy
        self.df[0,7] = (-1.0+yy)*0.5
        self.df[1,4] = (-1.0+xx)*0.5
        self.df[1,5] = -y-xy
        self.df[1,6] = (1.0-xx)*0.5
        self.df[1,7] = -y+xy
        return self.df
	
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
        
        self.f[0] = 1.0-x-y-z+xy+xz+yz-xyz	 # -1,-1,-1
        self.f[1] = 1.0+x-y-z-xy-xz+yz+xyz	 #  1,-1,-1
        self.f[2] = 1.0+x+y-z+xy-xz-yz-xyz	 #  1, 1,-1
        self.f[3] = 1.0-x+y-z-xy+xz-yz+xyz	 # -1, 1,-1
        self.f[4] = 1.0-x-y+z+xy-xz-yz+xyz	 # -1,-1, 1
        self.f[5] = 1.0+x-y+z-xy+xz-yz-xyz	 #  1,-1, 1
        self.f[6] = 1.0+x+y+z+xy+xz+yz+xyz	 #  1, 1, 1
        self.f[7] = 1.0-x+y+z-xy-xz+yz-xyz	 # -1, 1, 1
        self.f = self.f/8.0
        return self.f

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
        if   x >  1: return [1,2,6,5]
        elif x < -1: return [0,4,7,3]
        elif y >  1: return [2,3,7,6]
        elif y < -1: return [0,1,5,4]
        elif z >  1: return [4,5,6,7]
        elif z < -1: return [0,3,2,1]
        else:        return None


class ShapeFunction_Hex20(ShapeFunctionPrototype):
    """Element function for linear element defined in MARC Manual B2.4-1,
    Taylor&Hughes (1981), p. 49

    The integration points (in parentheses) are located at unexpected
    locations (for MARC)!

          7-------14------6
         /|(6)        (7)/|
        / |             / |
       15 |            13 |
      /  19           /   18
     / (4)|       (5)/    |
    4-------12------5     |
    |     |         |     |
    |     3------10-|-----2
    |    / (2)      | (3)/
   16   /           17  /
    |  11           |  9
    | /             | /
    |/ (0)       (1)|/    
    0-------8-------1

    """

    name = 'Hex20'
    dim, nnodes = 3, 20
    cornernodes = [0,1,2,3,4,5,6,7]
    nsides      = 6
    sidetype    = 'Quad8'
    sidenodes   = N.array(
                  [[0,11,3,10,2,9,1,8],
                   [0,8,1,17,5,12,4,16],
                   [1,9,2,18,6,13,5,17],
                   [2,10,3,19,7,14,6,18],
                   [3,11,0,16,4,15,7,19],
                   [4,12,5,13,6,14,7,15]
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
        self.f[12] = (1.0+z-y-yz-xx-xxz+xxy+xxyz)/4.0
        self.f[13] = (1.0+z+x+xz-yy-yyz-xyy-xyyz)/4.0
        self.f[14] = (1.0+z+y+yz-xx-xxz-xxy-xxyz)/4.0
        self.f[15] = (1.0+z-x-xz-yy-yyz+xyy+xyyz)/4.0
        self.f[16] = (1.0-y-x+xy-zz+yzz+xzz-xyzz)/4.0
        self.f[17] = (1.0-y+x-xy-zz+yzz-xzz+xyzz)/4.0
        self.f[18] = (1.0+y+x+xy-zz-yzz-xzz-xyzz)/4.0
        self.f[19] = (1.0+y-x-xy-zz-yzz+xzz+xyzz)/4.0
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
        self.df[0, 12] = -2*x-2*xz+2*xy+2*xyz
        self.df[1, 12] = -1.0-z+xx+xxz
        self.df[2, 12] =  1.0-y-xx+xxy
        self.df[0, 13] =  1.0+z-yy-yyz
        self.df[1, 13] = -2*y-2*yz-2*xy-2*xyz
        self.df[2, 13] =  1.0+x-yy-xyy
        self.df[0, 14] = -2*x-2*xz-2*xy-2*xyz
        self.df[1, 14] =  1.0+z-xx-xxz
        self.df[2, 14] =  1.0+y-xx-xxy
        self.df[0, 15] = -1.0-z+yy+yyz
        self.df[1, 15] = -2*y-2*yz+2*xy+2*xyz
        self.df[2, 15] =  1.0-x-yy+xyy
        self.df[0, 16] = -1.0+y+zz-yzz
        self.df[1, 16] = -1.0+x+zz-xzz
        self.df[2, 16] = -2*z+2*yz+2*xz-2*xyz
        self.df[0, 17] =  1.0-y-zz+yzz
        self.df[1, 17] = -1.0-x+zz+xzz
        self.df[2, 17] = -2*z+2*yz-2*xz+2*xyz
        self.df[0, 18] =  1.0+y-zz-yzz
        self.df[1, 18] =  1.0+x-zz-xzz
        self.df[2, 18] = -2*z-2*yz-2*xz-2*xyz
        self.df[0, 19] = -1.0-y+zz+yzz
        self.df[1, 19] =  1.0-x-zz+xzz
        self.df[2, 19] =  -2*z-2*yz+2*xz+2*xyz
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
        
# all shape functions are registered here

shapeFunctions = {
    'Tri3':  ShapeFunction_Tri3,
    'Quad4': ShapeFunction_Quad4,
    'Quad8': ShapeFunction_Quad8,
    'Hex8' : ShapeFunction_Hex8, 
    'Hex20': ShapeFunction_Hex20
   }

if __name__ == '__main__':
    for n, sf in shapeFunctions.items():
        print '===== %s =====' % n
        s = sf()
        s.calcGauss()
        print s.gaussShapeInv
    

