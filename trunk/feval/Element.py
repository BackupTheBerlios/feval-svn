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
# Purpose:  Container for Finite Element Model data
#
#============================================================

import numpy as N
import numpy.linalg as LA
import scipy.optimize

#============================================================================
# try to import Scientific Python (scipy.sourceforge.net)

# try:
#     import scipy.optimize._minpack
#     print 'using Scipy directly'

#     def minimize(element, globalcoord, accuracy=1.49012e-8):
#         """Find the local coordinates in the system of |element| if
#         the global coordinates |globalcoord| are given.
#         This is achieved with a Newton
#         procedure. The initial local coordinate is [0, 0].
#         The algorithm is evaluated until the accuracy of the result is
#         better than |accuracy| (in local coordinates).

#         We do this call the function directly, without the error checking in
#         scipy.optimize.fsolve
#         """
#         retval = scipy.optimize._minpack._hybrd(element.calcDev, \
#                                                 (0,)*element.dim,\
#                                                 (globalcoord,), 0, \
#                                                 accuracy, 200, -10, 10,\
#                                                 0.0, 100, None)
#         return retval[0]

#     print 'using Scipy normally'
#         # somewhat slower version with error checking
#     def minimize(element, globalcoord, accuracy=1.49012e-8):
#         retval = scipy.optimize.fsolve(element.calcDev, element.lcoord,
#                                        full_output=1, \
#                                        #fprime = element.calcDDev, \
#                                        args=(globalcoord,), xtol = accuracy)

#         return retval[0]


# except ImportError:
#     print 'Scipy not found'
#     raise ImportError

# except:
#     print 'problems defining scipy shortcut'
#     raise ImportError

def minimize(element, globalcoord, accuracy=1.49012e-8):
    if len(globalcoord) == element.shape.dim:
        retval = scipy.optimize.fsolve(element.calcDev, element.lcoord,
                                           full_output=1, \
                                           #fprime = element.calcDDev, \
                                           args=(globalcoord,), xtol = accuracy)
    # point embedded in a bigger space (2D element in 3D space)
    else:
        n1, n2, n3 = element.shape.cornernodes[:3]
        v1 = element.nodcoord[n2]-element.nodcoord[n1]
        v2 = element.nodcoord[n3]-element.nodcoord[n2]
        retval = scipy.optimize.fsolve(element.calcDevEmbedded, element.lcoord,
                                           full_output=1, \
                                           #fprime = element.calcDDev, \
                                           args=(globalcoord, v1, v2), xtol = accuracy)
    return retval[0]


class Element:
    """Holds the element properties:
    o shape is a shape function instance
    o dim and nnodes: number of (space) dimensions and number of nodes 
    o nodes is a list of node names (not instances, just the keys for the
	  model nodes dictionary)
    o name is the name of the element (the key in the model element dictionary)
    o nodcoord: array to keep the nodal coordinates
    o lcoord is to keep the current local coordinates
    """

    def __init__(self, nodes, shape, nodcoord, name = None):
        self.name     = name
        self.shape    = shape
        self.dim, self.nnodes = shape.dim, shape.nnodes
        self.lcoord   = N.zeros(self.dim, N.float_)
        # cut the list of the nodes to the actual nodes
        # (e.g. MARC sometimes has internal nodes)
        self.nodes    = nodes[:self.nnodes]
        #self.nodcoord = nodcoord[:self.nnodes,:self.shape.dim]
        self.nodcoord = nodcoord[:self.nnodes,:]
        self.N_ik = None

    def __repr__(self):
        return "Element %s" % str(self.name)

    def __str__(self):
        return "Element %s" % str(self.name)

    def mapNodalVar(self, nodalvar, lcoord=None, deriv=False):
        """Interpolate the nodal variables to the point given in
        |lcoord| using the element shape function
        if |deriv| is True calculate the spatial derivative of nodalvar"""
        if lcoord != None:
            self.lcoord = lcoord
        if deriv:
            sh = self.shape.calcShapeDeriv(self.lcoord) # the shape function derivatives
            inverse_jacobi = LA.inv(N.dot(self.nodcoord.transpose(), sh.transpose()))
            return N.dot(N.dot(nodalvar, sh.transpose()),inverse_jacobi)
        else:
            sh = self.shape.calcShape(self.lcoord) # the shape functions
            return N.dot(nodalvar, sh.transpose())

    def mapIntPointVar(self, intpointvar, lcoord=None, deriv=False):
        """Interpolate the interpolation point variables to the point given in
        |lcoord| using the element shape function"""
        if lcoord != None:
            self.lcoord = lcoord
        if self.shape.gaussShapeInv == None:
            self.shape.calcGauss()
        if deriv:
            sh = self.shape.calcShapeDeriv(self.lcoord) # the shape function derivatives
        else:
            sh = self.shape.calcShape(self.lcoord) # the shape functions
        nodvars = N.dot(self.shape.gaussShapeInv, intpointvar)
        return N.dot( sh, nodvars )

    def findGlobalCoord(self, lcoord=None):
        """Find the global coordinates if the local coordinates
        |lcoord| are given. This is straight-forward (and exactly
        the same as mapNodalVar(self.nodcoord, lcoord) ).
        """
        if lcoord != None:
            self.lcoord = lcoord
        return N.dot(self.shape.calcShape(self.lcoord), self.nodcoord )

    def calcDev(self, lc, globalcoord):
        """Return the difference of the global coordinate caluculated
        with the local coordinate |lc| and the shape functions to the
        global coordinate given in |globalcoord|""" 
        return N.dot( self.shape.calcShape(lc), self.nodcoord)-globalcoord

    def calcDevEmbedded(self, lc, globalcoord, v1, v2):
        """Return the difference of the global coordinate caluculated
        with the local coordinate |lc| and the shape functions to the
        global coordinate given in |globalcoord|""" 
        x = N.dot( self.shape.calcShape(lc), self.nodcoord)-globalcoord
        return N.dot(x, v1), N.dot(x, v2)

    def findLocalCoord(self, globalcoord, accuracy=1.49012e-8):
        """Find the local coordinates if the global coordinates
        |globalcoord| are given. The accuracy of the result is 
        better than |accuracy| (in local coordinates)."""

        # We use the pluggable function "minimize".
        #self.lcoord = minimize( self, globalcoord[:self.shape.dim], accuracy )
        self.lcoord = minimize( self, globalcoord, accuracy )
        return self.lcoord

    def containsPoint(self, globalcoord, accuracy = 1.e-5):
        """Test whether the element contains a point given in global
        coordinates
        """
        #lcoord = minimize( self, globalcoord[:self.shape.dim], accuracy )
        lcoord = minimize( self, globalcoord, accuracy )
        return (N.absolute(lcoord) <= 1.+accuracy).all()

    def setLocalCoord(self, lcoord):
        """set the local coorinate to |lcoord|"""
        self.lcoord = lcoord

if __name__ == '__main__':
    
    import ShapeFunctions
    sh = ShapeFunctions.shapeFunctions['Quad4']()
    nodes = [1, 2, 3, 4]
    #nodcoord = N.array(([0,0], [2,0], [3,2], [0,2]),dtype=N.float_)
    nodcoord = N.array(([0.2,0.], [10.,0.], [1,2], [0,3.]),dtype=N.float_)
    e = Element( nodes, sh, nodcoord, 1 )

    print e.findLocalCoord([1., 1.])
    print e.findGlobalCoord(e.findLocalCoord([10., 21.]))
    print e.findLocalCoord(e.findGlobalCoord([0., -1.]))
    print e.containsPoint([0., 2.], accuracy=0.01)

    print e.findGlobalCoord([-0., 0.])

    nodvar = e.nodcoord.transpose()
    
    print '---------'
    print e.mapNodalVar(nodvar,deriv=True)
