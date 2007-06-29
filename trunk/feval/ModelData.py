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

class ModelData:
    """Container for a reduced FE Model consisting of:
    - a connecitvity dictionary, each entry is of the form (elemtype, conn)
    - a dictionary of nodal coordinates
    - a list of nodal variables
    - a list of integration point variables
    - a dictionary of setstypes, each containing a dictionary of sets
	(lists) 

    The use of dictionaries for connectivity and nodes allows a flexible
    handling of non-continously numbered models.
    """

    DimNames = ('_x','_y','_z')

    def __init__(self, name='FE-Model'):
        self.name = name
        self.Conn = {}
        self.Coord = {}
        self.Intp = None
        self.Sets = {}
        self.initNodVar()
        self.initIntpVar()
        self.dirty = {'Coord': True, 'Conn': True, 'Var': False}

    def __repr__(self):
        """print a nice model summary"""
        statusFlag = {True: '*', False: ''}
        return "%s \n %6d Elements %s\n %6d Nodes %s\n %6d Nodal values %s"  \
               % (self.name,
                  len(self.getElementNames()), statusFlag[self.dirty['Conn']],
                  len(self.getNodeNames()), statusFlag[self.dirty['Coord']],
                  len(self.NodVarInfo), statusFlag[self.dirty['Var']],
                  )

    def update(self):
        pass

    def initNodVar(self):
        self.NodVarInfo = []
        self.NodVar = {}

    def initIntpVar(self):
        self.IntPointVar = {}
        self.IntPointVarInfo = []

    def makeModelCache(self):
        """dummy"""
        pass

    def setName(self, name):
        self.name = name

    def setConn(self, id, elemtype, conn):
        """deprecated, use setElement instead"""
        self.setElement(id, elemtype, conn)
        self.dirty['Conn'] = True

    def setElement(self, id, elemtype, conn):
        """set element with identifier |id|, type |elemtype|
        and the connectivity list |conn|"""
        # setting the connectivity as numpy array may change the type
        # this is a problem for integer and string names
        self.Conn[id] = (elemtype, conn)
        self.dirty['Conn'] = True

    def removeElement(self, id):
        """remove the element with name |id|"""
        del self.Conn[id]
        self.dirty['Conn'] = True

    def getElementNames(self):
        return self.Conn.keys()

    def getElementConn(self, id):
        return self.Conn[id]

    def setCoord(self, id, coord):
        """deprecated, use setCoordinate instead"""
        self.setCoordinate(id, coord)
        self.dirty['Coord'] = True

    def setCoordinate(self, id, coord):
        self.Coord[id] = N.asarray(coord, dtype=N.float_)
        self.dirty['Coord'] = True

    def removeCoordinate(self, id):
        """remove the element with name |id|"""
        del self.Coord[id]
        self.dirty['Coord'] = True

    def getCoordinateNames(self):
        return self.Coord.keys()

    def getIntPointVarInfo(self):
        """Get the integration point variable information
        """
        return self.IntPointVarInfo

    def getNodeNames(self):
        return self.Coord.keys()

    def getNodVarInfo(self):
        """Get the nodal variable information
        """
        return self.NodVarInfo

    def getCoordinate(self, id):
        return self.Coord[id]

    def getCoordNames(self):
        return self.Coord.keys()

    def setIntPointVar(self, id, intpointvar):
        self.IntPointVar[id] = intpointvar
        self.dirty['Var'] = True

    def setIntPointVarInfo(self, intpointvarinfo):
        """Wrapper to set the integration point variable information
        This is accessed from the file reader class (eg. MarcT16File)
        """
        self.IntPointVarInfo = intpointvarinfo
        self.dirty['Var'] = True

    def setNodVar(self, id, nodvar):
        """Set the nodal variables
        id is the node id
        """
        self.NodVar[id] = nodvar
        self.dirty['Var'] = True

        ## this is not needed, I guessk
        #     def getNodVar(self, id):
        #         """Get the nodal variables
        #         id is the node id
        #         """
        #         return self.NodVar[id]

    def getNodVarsAsArray(self):
        """Collect all nodal variables in an array.
        The array is sorted by the keys of the NodVar
        """
        import numpy as N
        keys = self.NodVar.keys()
        keys.sort()
        vars = []
        for k in keys:
            vars.append(self.NodVar[k])
        return N.asarray(vars, dtype=N.float_)

    def setNodVarInfo(self, nodvarinfo):
        """Wrapper to set the nodal variable information
        This is accessed from the file reader class (eg. MarcT16File)
        """
        self.NodVarInfo = nodvarinfo
        self.dirty['Var'] = True

    def setSet(self, type, name, set):
        """A dictionary of set-types, each containing a
        dictionary of sets (lists)
        Standard names for the set |type| are 'elem' and 'node'
        these get updated when the mesh is renumbered"""
        if not self.Sets.has_key(type):
            self.Sets[type] = {}
        self.Sets[type][name] = set

    def getSet(self, type, name):
        """get the set of type |type|, specified by |name|
        """
        try:
            return self.Sets[type][name]
        except:
            return None
            
    def getSetTypes(self):
        """get all types of sets
        """
        return self.Sets.keys()

    def getSetNames(self, type):
        """get all types of sets
        """
        try:
            return self.Sets[type].keys()
        except:
            return None

    def renumberNodes(self, base=1):
        """renumber all nodes starting with |base|
        also change the number in the connectivty and the set definition
        """
        newCoord = {}
        trans = {}
        newnr = base
        for nr, item in self.Coord.items():
            newCoord[newnr] = item
            trans[nr] = newnr
            newnr += 1
        self.Coord = newCoord
        for nr, item in self.Conn.items():
            newnodes = []
            self.setElement(nr, item[0],
                            [trans[n] for n in item[1]])
        if self.Sets.has_key('node'):
            for key, set in self.Sets['node'].items():
                self.Sets['node'][key] = [trans[n] for n in set]
        if len(self.NodVar) > 0:
            newNodVar = {}
            for n in self.NodVar.keys():
                newNodVar[trans[n]] = self.NodVar[n]
            self.NodVar = newNodVar
        self.makeModelCache()
        self.dirty['Coord'] = True
        

    def renumberElements(self, base=1):
        """renumber all elements starting with |base|
        also change the number in the set definition
        """
        newConn = {}
        trans = {}
        newnr = base
        for nr, item in self.Conn.items():
            newConn[newnr] = item
            trans[nr] = newnr
            newnr += 1
        self.Conn = newConn
        if self.Sets.has_key('elem'):
            for key, set in self.Sets['elem'].items():
                newelem = []
                for n in set:
                    newelem.append(trans[n])
                self.Sets['elem'][key] = newelem
        self.makeModelCache()
        self.dirty['Conn'] = True


if __name__ == '__main__':

    from feval.FEval import *

    m = FEModel()
    from feval.fecodes.marc.MarcFile import *  
    mf = MarcFile(m)
    mf.readFile('../data/marc/test1.dat')
    #mf.readInc(2)

    vars = m.getNodVarsAsArray()
    m.update()
    print m
    
    m.renumberElements(base=0)
    m.renumberNodes(base=0)
    
    print m
