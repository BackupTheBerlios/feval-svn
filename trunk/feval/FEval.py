# -*- coding: iso-8859-1 -*-

#============================================================================
#
#           This file is part of FEval, a module for the
#           evaluation of Finite Element results
#
# Licencse: FEval is provided under the GNU General Public License (GPL)
#
# Author:   Martin Lüthi, tnoo@tnoo.net
#
# Homepage: http://www.sourceforge.net/projects/feval
#
# History:  2001.07.02 (ml):  Improved findNextElement. Performance 
#                             improved by 15% (2d) to 80% (3d) 
#           2001.09.21 (ml):  Code cleaned up for intial public release
#           2002.05.17 (ml):  Use Scipy instead of Multipack
#           2002.07.31 (ml):  Improved version of findModelBoundary
#           2002.08.08 (ml):  getIntPointVar works
#           2002.08.10 (ml):  all methods now allow to give coordinates as
#                             tuples or lists
#  
# Todo:     o findModelBoundary (make it work for outside points)
#           o getIntPointVar (there might still be some problems)
#           o handle case when scipy is not available (Newton algorithm)
#============================================================================

import numpy as N
import ShapeFunctions
from ModelData import *
import Element

#============================================================================
# try to use Psyco (psyco.sourceforge.net)
# if configured, this will speed up things considerably
try:
    import psyco
    from psyco.classes import *
except ImportError:
    class _psyco:
        def jit(self):      pass
        def bind(self, f):  pass
        def proxy(self, f): return f
    psyco = _psyco()


#  		# take the last local coordinates
#  		lc  = N.zeros(element.dim)
#  		sh  = element.shape.evalShape(lc)
#  		shd = element.shape.evalShapeDeriv(lc)
#  		print shd
#  		newgc = N.add.reduce( element.nodcoord * sh[:,N.newaxis] )
#  		print newgc
#  		dcoord = globalcoord - newgc
#  		it = 0; maxit = 100
#  		# FIXME: recalculate accuracy
#  		while (N.add.reduce(dcoord*dcoord) >= accuracy*accuracy) \
#  		      and (it < maxit):
#  			# the Jacobi-Matrix
#  			jacob = N.matrixmultiply( shd, element.nodcoord ) 
#  			ijacob = LA.inverse(jacob)
#  			lc = lc + N.add.reduce(ijacob * dcoord)

#  			# If point is outside of element:  return
#  			if it > 1 and not isInRange( lc, 1.0+accuracy ):
#  				return lc

#  			sh = element.shape.evalShape(lc)
#  			shd = element.shape.evalShapeDeriv(lc)
#  			# new global coordinates
#  			newgc = N.add.reduce( element.nodcoord * sh[:,N.newaxis] )
#  			dcoord = globalcoord - newgc
#  			it = it+1
#  		element.lcoord = lc
#  		return lc

# a convenience function, should be faster! 
def isInRange( val, bound ):
    """Returns true if all components of the array |val| are less than |bound|"""
    return N.alltrue( N.less_equal( N.absolute( val ), bound ) )

class FEModel(ModelData):

    def __init__(self, name='FE-Model', verbose=0):
        """
        |accuracy| is the accuracy of local coordinates
        |elem_tol| is the tolerance of the local coordinates and lets
        the user adjust the sloppyness of evaluation
        o if this value is exceeded, a new element is sought
        o if the element is at the boundary (no adjacent element
            exists), the element is returned, if the local coordinates
            within the elem_btol
        |elem_btol| is the tolerance of the model boundary,
        if the local coordinates + elem_btol is exceeded, the point
        is considered outside of the model
        """
        ModelData.__init__(self, name)
        self.elementMidPoints = None
        # reset all cached results
        self.nodeIndex = {}
        self.boundaryElems = {}
        self.boundaryNodes = {}

        self.accuracy = 1.49012e-8
        self.elem_tol = self.elem_btol = 0.0001
        self.makeModelCache()
        self.verbose = verbose

    def update(self):
        self.makeModelCache()

    def makeModelCache(self):
        """Create caches for model quantities to allow fast lookups 
        o an array with the mid-point coordinates of all elements
        o a dictionary (nodeIndex) with the nodes as key and the elements as value
          (this is the inverse of the connectivity self.Conn), but also contains
          loose nodes (not in any connectivity)

        Reset the element variables (this is useful when a new
        timestep is read)""" 
        # reset all cached results
        self.elementCache  ={None: None}  # dummy entry
        self.nodeIndex     = {}
        self.boundaryElems = {}
        self.boundaryNodes = {}
        midPoints = []
        self.removeUnusedNodes()
        try:
            mincoord = self.Coord.values()[0]
        except:
            mincoord = N.zeros(3)
        maxcoord = mincoord.copy()
        for elem in self.Conn.keys():	 # loop over elements
            coord = []
            for n in self.Conn[elem][1]:    # loop over the nodes
                coord.append( self.Coord[n] )
                mincoord = N.minimum(mincoord, self.Coord[n])
                maxcoord = N.maximum(maxcoord, self.Coord[n])
                self.nodeIndex.setdefault(n, []).append(elem)
            midPoints.append( N.asarray(coord).mean(axis=0) )

        self.elementMidPoints = N.asarray( midPoints, dtype=N.float_ )
        self.boundingbox      = N.asarray([mincoord, maxcoord])
        self.dirty['Conn'] = False
        self.dirty['Coord'] = False
        self.dirty['Var'] = False
        
    def removeUnusedNodes(self):
        """remove unused Nodes
        we should do this on creating the model cache
        maybe this should be called on reading
        """
        usednodes = set()
        for conn in self.Conn.values():
            usednodes.update(set((conn[1])))
        for node in self.Coord.keys():
            if not node in usednodes:
                self.Coord.pop(node)

    def findClosestNode(self, coord, recalc=True):
        """Find the name of the node closest to |coord|.
        the key of the node-dictionary is returned
        """
        allcoord = self.getAllCoordinatesAsArray()
        idx = N.sum((allcoord - coord)**2, axis=1).argmin()
        return self.Coord.keys()[idx]

    def findClosestElement(self, coord):
        """Find the name of the element closest to |coord|.
        the key of the element-dictionary is returned
        """
        dc = self.elementMidPoints - coord
        #eidx = N.argmin( N.sum( dc*dc, 1) )
        eidx = (dc*dc).sum(axis=1).argmin()
        ekey = self.Conn.keys()[eidx]
        return self.getElement(ekey)

    def findElement(self, coord, accuracy=1.49012e-8):
        """Find the element containing the global coordinate |coord|
        """
        n = self.findClosestNode(coord)
        for elem in self.nodeIndex[n]:
            e = self.getElement(elem)
            if e.containsPoint(coord, accuracy):
                return e
        return None

#     def findElement(self, coord):
#         """Find the element containing the global coordinate |coord|
#         """
#         e = self.findClosestElement(coord)
#         e_old = None
#         lc = e.findLocalCoord(coord, self.accuracy)
#         while e and not isInRange( lc, 1.0+self.elem_tol):
#             ee = self.findNextElement( e, lcoord = lc )
#             # if new element is within the model
#             # test for element ping-pong loop
#             if ee != None and e_old != None:
#                 print ee.name, e_old.name
#             if ee and ee != e_old:
#                 lc = ee.findLocalCoord(coord, self.accuracy)
#                 e_old = e
#                 e = ee
#             # if no new element exists (= boundary)
#             else:
#                 lc = e.findLocalCoord(coord, self.accuracy)
#                 if not isInRange( lc, 1.0+self.elem_btol):
#                     e = None
#                     if self.verbose:
#                         print "The point is not within the model: ", coord
#                 else:
#                     return e
#         return e

    def findElementsFromNodes(self, nodelist):
        """Find the elements containing two of the nodes in the |nodelist|,
        returns a list of (element, face-number) pairs, or []"""
        elems = {}
        for n in nodelist:
            ee = self.nodeIndex[n]
            for e in ee:
                e.setdefault(e, []).append(n)
                #                 if elems.has_key(e):
                #                     elems[e].append(n)
                #                 else:
                #                     elems[e] = [n]

        elemsides = []
        for e, nodes  in elems.items():
            if len(nodes) > 1:
                elem = self.getElement( e )
                # take the corner nodes
                cn = list( N.array(elem.nodes)[elem.shape.cornernodes] )
                cn.append(cn[0])
                for side in range(len(cn)-1):
                    if cn[side] in nodes and cn[side+1] in nodes:
                        elemsides.append( (e, side) )
        elemsides.sort()
        return elemsides

    def findNextElement(self, actelem, sidenr=0, lcoord=None):
        """Find the next element in the model, given the acutal
        element |actelem| and either
        o the side number |sidenr|
        o the local coordinats |lcoord| (one component is > 0)
        returns the next element or None"""
        # get the nodes on the side, either in lcoord, or directly
        if lcoord != None:
            nextnodes = actelem.shape.nextPattern(lcoord)
        else:
            nextnodes = list(
                set(actelem.shape.sidenodes[sidenr]).intersection(actelem.shape.cornernodes))
        # loop over the nodes in the list of the next nodes
        elems = []
        # find all elements that share the nextnodes
        for n in nextnodes:
            elems.extend(self.nodeIndex[actelem.nodes[n]])
        # delete the present element from that list
        elems.remove(actelem.name)
        # find the element with the same number of nodes on that side
        # was soll das. falsch!!!
        for e in elems:
            if elems.count(e) == len(nextnodes):
                return self.getElement( e )
        return None


    def getElement(self, elename):
        """Create an element instance from the connectivity information
        (the Element class holds a list of nodes)
        |elem| is the name (i.e. the number) of the element
        (This used to be ElementFactory)
        """
        # test whether the element already exists, shortcut if exists
        try:
            return self.elementCache[elename]
        except:
            pass
        conn = self.Conn[elename]
        nodes = conn[1]
        sh = ShapeFunctions.shapeFunctions[conn[0]]()

        # collect the nodal coordinates
        nodcoord = N.array([self.Coord[n] for n in nodes])

        # invoke an instance of the element
        e = Element.Element( nodes, sh, nodcoord, elename )
        self.elementCache[elename] = e
        return e

    def getAllCoordinatesAsArray(self):
        """Collect all coordinates and return them as numpy array
        """
        return N.array(self.Coord.values())

    def getVariables(self, element, nodvars=None, intpointvars=None, deriv=False):
        """Get a variable from the model for the |element| given
        The local coordinate of the |element| must be set to the correct value
        (as is the case after a findElement(point) or after explicitly
        setting through element.setLocalCoord(lcoord)).

        |nodvars| is the list of variable names defined on the nodes
        |intpointvars| is the list of variable names defined on the integration points

        The names of the variables can be queried through
        getNodVarInfo and getIntPointVarInfo
        """

        resNodvars, resIntPointVars = [], []

        # find the appropriate nodal variables
        if nodvars:
            for postvar in nodvars:
                idx = self.NodVarInfo.index(postvar)
                nodvar = []
                for node in element.nodes:
                    nodvar.append( self.NodVar[node][idx] )
                resNodvars.append( nodvar )
            resNodvars = N.asarray(resNodvars, dtype=N.float_)
            resNodvars = element.mapNodalVar(resNodvars, deriv=deriv)

        # find the appropriate integration point variable
        if intpointvars:
            idx = []
            resIntPointVars = self.IntPointVar[element.name]
            for postvar in intpointvars:
                idx.append(self.IntPointVarInfo.index(postvar))
            resIntPointVars = N.take( resIntPointVars, idx, 1 )
            resIntPointVars = element.mapIntPointVar(resIntPointVars, deriv=deriv)

        return resNodvars, resIntPointVars

    def getNodVar(self, point, postvars, deriv=False, accuracy=1.49012e-8):
        """Find the values of the variables |postvars| at
        the global coordinate |point|
        use |accuracy| to locate the points
        """
        point = N.asarray(point, dtype=N.float_)
        element = self.findElement(point, accuracy)
        # check wether the point is within the model
        if not element:
            return None
        resNodvars, resIntPointVars = self.getVariables(element, nodvars=postvars, deriv=deriv)
        return resNodvars

    def getIntPointVar(self, point, postvars, deriv=False):
        """Find the values of the integration point variables |postvars| at 
        the global coordinate |point|
        """
        point = N.asarray(point, dtype=N.float_)
        element = self.findElement(point)
        # check wether the point is within the model
        if not element:
            return None
        resNodvars, resIntPointVars = self.getVariables(element, intpointvars=postvars, deriv=deriv)
        return resIntPointVars

    def getAdjacentNodes(self, node, onboundary=False):
        """Get the nodes that are connected to |node| on an edge.
        If |node| is not a element corner node return None
        (this should happen only for higher order elements).
        """
        nextnodes = set()
        for elem in self.nodeIndex[node]:
            e = self.getElement(elem)
            enodes = N.asarray(e.nodes)
            if not node in enodes[e.shape.cornernodes]:
                return None
            nextnodes.update(enodes[e.shape.nextnodes[enodes == node]][0].tolist())
        if onboundary:
            nextnodes = set(self.findBoundaryNodes().keys()).intersection(nextnodes)
        return list(nextnodes)

    def findBoundaryNodes(self, recalc=False):
        """Find all the nodes on the model boundary.

        Returns a dict with node names as keys.

        The values are dicts with the elements that contain the node,
        and a list with the boundary sides of this element.
        """
        if not self.boundaryNodes or recalc:
            belems = self.findBoundaryElements(recalc=recalc)
            bnodes = {}
            for elename, sides in belems.items():
                e = self.getElement(elename)
                for side in sides:
                    nodes = N.array(e.nodes)[e.shape.sidenodes[side]]
                    for n in nodes:
                        bnodes.setdefault(n, {}).setdefault(elename, list()).append(side)
            self.boundaryNodes = bnodes
        return self.boundaryNodes

    def findBoundarySides(self, recalc=False):
        """Find all the nodes on the model boundary.

        Returns a dict with side numbers as keys.

        The values are lists with the node names that live on the sides.
        """
        # find the nodes on the boundary sides
        belems = self.findBoundaryElements(recalc=recalc)
        bsides = {}
        for elename, sides in belems.items():
            e = self.getElement(elename)
            for side in sides:
                nodes = N.array(e.nodes)[e.shape.sidenodes[side]]
                bsides.setdefault(side, set()).update(set(nodes))
        return bsides

    def findBoundaryElements(self, side=None, recalc=False):
        """Find all the elements on the model boundary. If |side| is given,
        find only those parts of the boundary that consist of a given element
        side (in the element local system; note that this might only be useful
        for regular quad grids).

        Returns a dict, where the keys are the element names and the values a
        list of all boundary sides
        """
        # use a dictionary (or a Set) to keep the boundary elements unique
        if not self.boundaryElems or not recalc:
            belems = {}
            if side != None:
                for elename in self.getElementNames():
                    e = self.getElement(elename)
                    if self.findNextElement(e, side) == None:
                        belems[elename] = [side]
            else:
                for elename in self.getElementNames():
                    e = self.getElement(elename)
                    for s in range(e.shape.nsides):
                        # sides on the boundary
                        if self.findNextElement(e, s) == None:
                            belems.setdefault(elename, []).append(s)
            self.boundaryElems = belems
        return self.boundaryElems

    def findModelBoundary(self, point, direction, savebmodel=True, culling=True):
        """Find the boundary of the model in |direction|, starting at the
        given |point|.  Return None if no boundary is found.
        |direction| and point| should be arrays, lists or tuples of the same
        |dimension as the model.  """

        #direction = direction / N.sqrt(N.dot(direction,direction))

        try:
            bmodel = self.bmodel_tri
        except:
            bmodel = self.extractBoundaryModel(triangles=True)
            if savebmodel:
                self.bmodel_tri = bmodel
        
        import intersect_triangle
        for ele in bmodel.getElementNames():
            e = bmodel.getElement(ele)
            res = intersect_triangle.intersect_triangle(point, direction, e.nodcoord)
            if res != None:
                # res also contains distance and local variables
                return res[-1]
        return None
            
                
        

    def findModelBoundaryOld(self, point, direction, accuracy=0.0001):
        """Find the boundary of the model in |direction|, starting at the
        given |point|.  Return None if no boundary is found.
        |direction| and point| should be arrays, lists or tuples of the same
        |dimension as the model.  """

        point, direction = N.asarray(point, dtype=N.float_), N.asarray(direction, dtype=N.float_)

        direction = direction / N.sqrt(N.dot(direction,direction))
        
        # set the accuracy to the given value
        elem_tol, elem_btol = self.elem_tol, self.elem_btol
        self.elem_tol, self.elem_btol = accuracy, accuracy
        
        e = self.findElement(point)
        # if we start outside the model: return find the next Element inside
        # the model
        if not e:
            #             node = self.getCoordinate( self.findClosestNode(point) )
            #             dist = N.sqrt(N.sum(node-point)**2)
            #             pt = point + dist*dir
            return None

        acc = 10.0
        while acc > accuracy:
            e = self.findElement(point)
            if not e:
                return None

            pt = point
            # find the boundary
            while 1:
                # dcoord is the half-size of the element in global coordinates
                # i.e. 0.5*(dx, dy, ..)
                dcoord = N.sum(N.dot( e.shape.calcShapeDeriv(N.zeros(e.dim)), e.nodcoord))
                pt = point + acc*direction*abs(dcoord)
                if self.verbose:
                    print 'element: ', e.name , ';  point  : ', pt
                ee = self.findElement(pt)
                if ee:
                    point = pt
                    e = ee
                else:
                    break

            # now we should be in some boundary element (probably not yet the correct..)
            lcoord = e.lcoord
            pt = point + acc*direction*abs(dcoord)
            lc = e.findLocalCoord(pt)
            if isInRange(lc, 1.):
                e.lcoord = lc
                point = pt
            else:
                e.lcoord = lcoord
                fact = max(2.5, max(lc)/2.5)
                acc = acc/fact  # 2.5 seems to be quite efficient

        # reset the tolerances
        self.elem_tol, self.elem_btol = elem_tol, elem_btol

        # go the exact element boundary
        e.lcoord = N.clip(e.lcoord, -1., 1.)
        return e.mapNodalVar(N.transpose(e.nodcoord))

    def sweepNodes(self, accuracy=1.e-12):
        """Remove duplicate nodes, throw away unused nodes, update the
        connectivity and the set definition

        |accuracy| gives the accuracy of the node coordinates

        This is a quite slow, but correct implementation.
        """
        self.makeModelCache()
        # remove all nodes (Coordinates) not in the index,
        # i.e. that are not part of an element connectivity
        removenodes = []
        for n in self.Coord:
            if not n in self.nodeIndex:
                removenodes.append(n)
        for n in removenodes:
            del self.Coord[n]
                
        # make a list of duplicate nodes
        nodecoords = N.array(self.Coord.values())
        nodenames  = N.array(self.Coord.keys())
        for name, coord in self.Coord.items():
            samenodes = nodenames[((coord- nodecoords)**2).sum(1) < accuracy]
            # if there is more than one node, retain the first and remove the rest
            if len(samenodes) > 1:
                for node in samenodes[1:]:
                    # maybe the node has been removed before. In that case we just go on
                    try:
                        del self.Coord[node]
                        # update the connectivity array
                        for elem in self.nodeIndex[node]:
                            elenodes = self.Conn[elem][1]
                            idx = list(elenodes).index(node)
                            elenodes[idx] = samenodes[0]
                    except:
                        pass
                    # update the sets: TODO

        self.makeModelCache()

    def rotate(self, phi, theta, psi):
        """rotate the mesh by the Euler angles
        |phi|, |theta|, |psi| should be given in degrees
        after this the update method should be called
        """
        p = -phi/180.*N.pi
        t = -theta/180.*N.pi
        s = -psi/180.*N.pi
        sp, cp = N.sin(p), N.cos(p)
        st, ct = N.sin(t), N.cos(t)
        ss, cs = N.sin(s), N.cos(s)

        for node, coord in self.Coord.items():
            if len(coord) == 3:
                x,y,z = coord
            elif len(coord) == 2:
                x,y = coord
                z   = 0.
            self.Coord[node] = N.array((
                ( cp*cs-sp*ct*ss)*x + ( sp*cs+cp*ct*ss)*y + (st*ss)*z,
                (-cp*ss-sp*ct*cs)*x + (-sp*st+cp*ct*cs)*y + (st*cs)*z,
                ( sp*st)*x          + (-cp*st)*y          + (ct)*z   ))

    def translate(self, dx):
        """translate the mesh by the amount |dx|
        which should be given as array of the correct length
        after this the update method should be called
        """
        dd = N.asarray(dx)
        for node, coord in self.Coord.items():
            self.Coord[node] = self.Coord[node] + dd

    def scale(self, factor):
        """scale the mesh with |factor|
        this just multiplies all coordinate values by |factor|
        """
        for node, coord in self.Coord.items():
            self.Coord[node] = self.Coord[node]*factor

    def extractBoundaryModel(self, triangles=False, faces=None):
        """extract a model that only contains the boundary nodes and sides
        if |faces| is given, only extract the elements with the given faces
        """
        bmodel = FEModel()
        belems = self.findBoundaryElements()

        for elename, sides in belems.items():
            e = self.getElement(elename)
            faces = faces or range(e.shape.nsides)
            for s in sides:
                if s in faces:
                    nodes = N.asarray(e.nodes)
                    shape = e.shape
                    if triangles:
                        for i, tri in enumerate(ShapeFunctions.shapeFunctions[shape.sidetype].triangles):
                            bmodel.setElement('%s-%elename-%d' % (str(elename), s, i), 'Tri3',
                                              nodes[shape.sidenodes[s]][tri])
                    else:
                        bmodel.setElement('%s-%elename' % (str(elename), s), shape.sidetype,
                                          nodes[shape.sidenodes[s]])
                    for c in nodes[shape.sidenodes[s]]:
                        bmodel.setCoordinate(c, self.getCoordinate(c))
        bmodel.renumberElements()
        bmodel.update()
        return bmodel

def TestFind():
    """Test function, mainly for profiling purposes"""
    for i in range(50):
        for n in m.getCoordinateNames():
            ci = m.getCoordinate(n)
            point = ci+N.array([0.0123,-0.0234])*i
            m.getNodVar(point,['v_x','f_y', 'v_z'])


def TestFindModelBoundary():
    """Test function, mainly for profiling purposes"""
    for c in m.Coord.items():
        coo = m.findModelBoundary( N.array([c[1][0],0.01]),
                                   N.array([1,-1]), accuracy=0.00001 )


# Tests
if __name__ == "__main__":
    from fecodes.marc.MarcPost import *
    from fecodes.marc.MarcFile import *

#     m = FEModel(verbose=1)
#     mf = MarcFile(m)
#     mf.readFile('/soft/numeric/feval/data/marc/test1.dat')

#     bednodes = m.Sets['node']['surf']    
#     elems = m.findElementsFromNodes(bednodes)
#     stop

#     print 'loading the model'
#     ## 2D-model test case
#     m = FEModel(verbose=1)
#     mf = MarcPostFile(m, '../data/marc/e7x1b.t16')
#     mf.readInc(1)

#     for n, c in m.Coord.items():
#         m.setNodVar(n, [c[0]*0.5, c[1]*2.])
#     m.setNodVarInfo(['X','Y'])
    
#     point = N.array([7,10,0])
#     print m.findElement(point).name
#     print m.getNodVar(point, m.NodVarInfo)
#     x =  m.getNodVar(point+N.array([1.0,0.1,0]), m.NodVarInfo, deriv=True)
#     print m.getNodVar(point, m.NodVarInfo, deriv=True)

#     print m.getAdjacentNodes(2)
#     #a = m.findBoundaryElements()

#     #     bmodel = m.extractBoundaryModel()
#     #     stop

#     m.setCoordinate(777, m.Coord[7])
#     m.Conn[23][1][3] = 777
#     m.update()

#     m.makeModelCache()
#     m.sweepNodes()
#     stop

#     e = m.findClosestElement([20, 20])
#     ee = m.findClosestElement([0, 0])

#     print m.findNextElement(e, lcoord = [0, -1]).name
#     print m.findNextElement(e, 0).name

#     stop

#     a = m.findBoundaryElements()
#     stop

#     point = m.getCoordinate(m.getCoordNames()[2])+N.array([-0.13,0.021])
#     print m.getNodVar(point,['d_x','d_y'])
#     point = m.getCoordinate(m.getCoordNames()[2])+N.array([-0.13,10.021])
#     boundary = m.findModelBoundary(point, N.array([5,-5]), accuracy=0.0000001)
#     print 'the boundary is at ', boundary

# #     print m.getIntPointVar(point,['t'])
 

#     stop

    # 3D-model
    try:
        m
    except:
        m  = FEModel()
        mf = MarcPostFile(m, '/home/tinu/projects/colle/marc/cgx8_dc0.80/cgx8dec.t16')
        mf.readInc(1)

    print m

    point = m.getCoordinate(1001)+N.array([-10,10,-50])
    direction = N.array([0,0,1])
    print m.findModelBoundaryOld(point, direction, accuracy=0.001)
    print 'no culling'
    print m.findModelBoundary(point, direction, culling=False)
    print 'culling'
    print m.findModelBoundary(point, direction, culling=True)

    print m
    m.bmodel_tri.renumberNodes(1)
    m.bmodel_tri.renumberElements(1)
    print m
    m.bmodel_tri.update()
    print m

    import feval.fecodes.gmv.GMVFile as gmv
    gf = gmv.GMVFile(m.bmodel_tri)
    gf.setWrite('gmvinput')
    gf.setWrite('nodes')
    gf.setWrite('cells')
    gf.setWrite('endgmv')
    gf.writeFile('A.gmv')




### profiling

#     print 'testing'
#     import profile

#     m.verbose=0
#     m.makeModelCache()
#     # profile.run('TestFind()')

#     psyco.bind(Element)
#     psyco.bind(FEModel)
#     psyco.bind(ShapeFunctions.ShapeFunction_Quad8)

#     m.makeModelCache()
#    profile.run('TestFindModelBoundary()')
# def TestFindModelBoundary():
#     """Test function, mainly for profiling purposes"""
#     for c in m.Coord.items():
#         coo = m.findModelBoundary( N.array([c[1][0],0.01]),
#                                    N.array([1,-1]), accuracy=0.00001 )


# # Tests
# if __name__ == "__main__":
#     from feval.fecodes.marc.MarcT16File import *
#     from feval.fecodes.marc.MarcFile import *

#     m = FEModel(verbose=1)
#     mf = MarcFile(m)
#     #mf.readFile('data/marc/test1.dat')
#     mf.readFile('/soft/numeric/feval/data/marc/test1.dat')

#     bednodes = m.Sets['node']['surf']    
#     elems = m.findElementsFromNodes(bednodes)
#     stop

#     print 'loading the model'
#     ## 2D-model test case
#     m = FEModel(verbose=1)
#     mf = MarcT16File(m, '../data/marc/e7x1b.t16', verbose=1)
#     mf.readInc(1)

#     point = m.getCoordinate(m.getCoordNames()[2])+N.array([-0.13,0.021])
#     print m.getNodVar(point,['d_x','d_y'])
#     point = m.getCoordinate(m.getCoordNames()[2])+N.array([-0.13,10.021])
#     boundary = m.findModelBoundary(point, N.array([5,-5]), accuracy=0.0000001)
#     print 'the boundary is at ', boundary
# #     print m.getIntPointVar(point,['t'])



#     ## 3D-model
# #     try:
# #         mf
# #     except:
# #         m = FEModel()
# #         m.verbose=1
# #         m.accuracy = 1.e-20
# #         ##mf = MarcT16File(m, '/home/tinu/numeric/feval/test/mm_22_51_22_3d_mid_g0.0-n1-a5.3.t16')
# #         mf = MarcT16File(m, '/home/tinu/projects/jako/marc/jako3db_polythermal.t16')
# #         mf.readInc(1)
# #     point = m.getCoordinate(m.getCoordinateNames()[1500])+N.array([-0.5,0.05,0.02])

# #     print 'nodvar ', m.getNodVar(point,['d_x','d_y'])
# #     print 'point accuracy: ', m.getIntPointVar(point,['x','y','z'])-point
# #     print m.getIntPointVar(point,['sigma_6'])
# #     print m.findModelBoundary(point, N.array([0,-5,-6]), accuracy=0.01)

# ### profiling

#     print 'testing'
#     import profile

#     m.verbose=0
#     m.makeModelCache()
#     # profile.run('TestFind()')

# #     psyco.bind(Element)
# #     psyco.bind(FEModel)
# #     psyco.bind(ShapeFunctions.ShapeFunction_Quad8)

#     m.makeModelCache()
# #    profile.run('TestFindModelBoundary()')

#     import feval.fecodes.gmv.GMVFile as gmv
#     datapath = '/home/tinu/fismo/trunk/visco3d/model-runs'
#     filename = 'with-tongue,n=1,etab=500,etaff=0.1.gmv'
#     m = FEModel()
#     m.verbose=0
#     m.accuracy = 1.e-20
#     mf = gmv.GMVFile(m)
#     mf.readFile(datapath+'/'+filename)
#     m.removeUnusedNodes()
    
#     for n, c in m.Coord.items():
#         m.setNodVar(n, c)
#     m.setNodVarInfo(['X', 'Y', 'Z'])
    
#     point = N.array([5000,100,0])
#     print m.findElement(point).name
#     print m.getNodVar(point, m.NodVarInfo, deriv=True)

