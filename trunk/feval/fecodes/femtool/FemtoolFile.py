import re
import Numeric; N = Numeric
from feval.FEval import *
from feval.FETextFile import *

class FemtoolFile(FETextFile):
    """Parse a Femtool input-file
    Take care of the node pattern which is different from the MARC node
    pattern which  is applied throughout these routines.
    (The MARC node pattern is: first the corner nodes, then the side
    nodes, the Femtool pattern is corner-side-corner-side...)

    The file contains
    o an |Info| dictionary, which is updated by updateInfo()
    o dictionaries with |Dirichlet| and |Neumann| boundary conditions

    The node pattern is saved in the list |Variables|
    """
    type = 'femtool'

    # Dictionary of the Element types
    shapeFunctionDict = { 4: 'Quad4', \
                          8: 'Quad8' }
    dimensionDict     = { 'Quad4': 2, \
                          'Quad8': 2}

    # the node pattern is from Femtool to Marc
    nodePattern = {}
    nodePattern['Quad4'] = [0,1,2,3]
    nodePattern['Quad8'] = [0,2,4,6,1,3,5,7]

    # the inverse node pattern is from Marc to Femtool
    nodePatternInv = {}
    nodePatternInv['Quad4'] = [0,1,2,3]
    nodePatternInv['Quad8'] = [0,4,1,5,2,6,3,7]

    def __init__(self, model):
        FETextFile.__init__(self, model)
        self.Info = {}
        self.Info['ndim']   = 0
        self.Info['nshape'] = 0
        self.Info['nvar']   = 0
        self.Info['nele']   = 0
        self.Info['npts']   = 0
        self.Info['ndiri']  = 0
        self.Info['nneum']  = 0
        self.Info['ntrans'] = 0
        self.Info['nwhat']  = 0
        self.Dirichlet = {}
        self.Neumann   = {}
        self.Variables = {}

        # Pairs of characters used to mark comments
        self.Comment = [( '#', '\n' )] 

    def updateInfo(self):
        nele = len( self.model.Conn )
        self.Info['nele'] = nele
        self.Info['npts'] = len( self.model.Coord )
        if nele > 0:
            c = self.model.Conn[ self.model.Conn.keys()[0] ]
            self.Info['ndim'] = self.dimensionDict[c[0]]
        # count number of Dirichlet and Neumann boundary conditions
        diricount = [0]
        for actdiri in self.Dirichlet.values():
            diricount.append(len(actdiri))
        self.Info['ndiri'] = max(diricount)
        neumcount = [0]
        for actneum in self.Neumann.values():
            neumcount.append(len(actneum))
        self.Info['nneum'] = max(neumcount)


    def extract_initialize(self, linelist):
        words = string.split(linelist[1])
        ndim, nshape, nvar, nele, npts, ndiri, nneum, \
              ntrans, nwhat = map( string.atoi, words[:9] )
        self.Info['ndim']   = ndim
        self.Info['nshape'] = nshape
        self.Info['nvar']   = nvar
        self.Info['nele']   = nele
        self.Info['npts']   = npts
        self.Info['ndiri']  = ndiri
        self.Info['nneum']  = nneum
        self.Info['ntrans'] = ntrans
        self.Info['nwhat']  = nwhat

    def compose_initialize(self):
        lines = []
        lines.append('initialize\n')
        self.updateInfo()
        ndim = self.Info['ndim']   
        nshape = self.Info['nshape'] 
        nvar = self.Info['nvar']   
        nele = self.Info['nele']   
        npts = self.Info['npts']   
        ndiri = self.Info['ndiri']  
        nneum = self.Info['nneum']  
        ntrans = self.Info['ntrans']
        nwhat = self.Info['nwhat']
        lines.append( '%7i'*9 % (ndim, nshape, nvar, nele, npts, \
                                 ndiri, nneum, ntrans, nwhat) + '\n')
        return lines


    def extract_timeintegration(self, linelist):
        words = string.split(linelist[1])
        ntshape, ntnum, ntimp, ntover = map( string.atoi, words[:4] )
        self.Info['ntshape'] = ntshape
        self.Info['ntnum']   = ntnum
        self.Info['ntimp']   = ntimp
        self.Info['ntover']  = ntover
        self.times = N.array( map( string.atof, linelist[2:] ) )

    def compose_timeintegration(self):
        lines = []
        ntshape = self.Info['ntshape'] 
        ntnum = self.Info['ntnum']   
        ntimp = self.Info['ntimp']   
        ntover = self.Info['ntover']  
        line = '%8i'*8 % (ntshape, ntnum, ntimp, ntover)
        lines.append(line + '\n')
        for t in self.times:
            lines.append( str(t) + '\n' )
        return lines

    def extract_userdata(self, linelist):
        words = string.split(linelist[1])
        nupts, nuelem, nudata = map( string.atoi, words[:3] )
        self.Info['nupts']   = nupts
        self.Info['nuelem']  = nuelem
        self.Info['nudata']  = nudata
        ptr = 2
        if nupts > 0:
            ud = []
            dptr = self.Info['npts']
            ulines = linelist[ ptr:ptr+dptr ]
            for line in ulines:
                ud.append( map( string.atof, string.split(line[:-1]) ) )
            self.upts = N.array( ud )

        if nuelem > 0:
            ptr = ptr + dptr
            ud = []
            dptr = self.Info['nelem']
            ulines = linelist[ ptr:ptr+dptr ]
            for line in ulines:
                ud.append( map( string.atof, string.split(line[:-1]) ) )
            self.uelem = N.array( ud )

        if nudata > 0:
            ptr = ptr + dptr
            ud = []
            ulines = linelist[ptr:]
            for line in ulines:
                l = map( string.atof, string.split(line[:-1]) )
                for ll in l:
                    ud.append( ll )
            self.udata = N.array( ud )

    def compose_userdata(self):
        lines = []
        lines.append('userdata\n')
        nupts = self.Info['nupts']
        nuelem = self.Info['nuelem']  
        nudata = self.Info['nudata']
        lines.append( '%9i'*3 % (nupts, nuelem, nudata) + '\n')
        return lines


    def extract_points(self, linelist):
        linelist = linelist[1:]
        Coord = {}
        for line in linelist:
            words = string.split(line)
            id = string.atoi(words[0])
            c = N.array( map(string.atof, words[1:]))
            self.model.setCoord( id, c )

    def compose_points(self):
        lines = []
        lines.append('points\n')
        ckeys = self.model.Coord.keys()
        ckeys.sort()
        ndim = self.Info['ndim']
        for id in ckeys:
            coord = self.model.Coord[id]
            line = '%5i' % (id) + ' %12.4f'*ndim \
                   % tuple(coord[0:ndim].tolist()) + '\n'
            lines.append(line)
        return lines

    def extract_elements(self, linelist):
        linelist = linelist[1:]
        for line in linelist:
            words = map( string.atoi, string.split(line) )
            id = words[0]
            c = N.array(words[1:])
            try: 
                elemtype = self.shapeFunctionDict[c[0]]
            except:
                elemtype = ''
                print '**** Element type %i not defined!' % (c[0])
            nodes = N.take( c[1:], self.nodePattern[ elemtype ] )
            self.model.setElement( id, elemtype, nodes )

    def compose_elements(self):
        lines = []
        lines.append('elements\n')
        ckeys = self.model.Conn.keys()
        ckeys.sort()
        for id in ckeys:
            conn = self.model.Conn[id]
            nodes = conn[1]
            nodes = N.take(nodes, self.nodePatternInv[ conn[0] ] )
            nnod = len(nodes)
            line = '%5i%5i' % (id, nnod)
            ll = '%5i'*nnod % tuple(nodes)
            lines.append(line+ll+'\n')
        return lines

    def setDirichlet(self, key, value):
        self.Dirichlet[key] = value

    def getDirichlet(self, key):
        return self.Dirichlet[key]

    def extract_dirichlet(self, linelist):
        words = string.split(linelist[1])
        diricount = map( string.atoi, words )
        linelist = linelist[2:]
        for idiri in range(len(diricount)):
            ndiri = diricount[idiri]
            actdiri = {}
            for i in range(ndiri):
                line = linelist[i]
                words = string.split(line)
                node  = string.atoi(words[0])
                value = string.atof(words[1])
                actdiri[node] =  value
            self.setDirichlet(idiri+1, actdiri)
            linelist = linelist[ndiri:]
        self.updateInfo()

    def compose_dirichlet(self):
        lines = []
        lines.append('dirichlet\n')
        diricount = []
        dirikeys = self.Dirichlet.keys()
        dirikeys.sort()
        for key in dirikeys:
            diricount.append(len(self.getDirichlet(key)))
        lines.append( '%9i'*len(diricount) % tuple(diricount) + '\n')
        for key in dirikeys:
            actdiri = self.getDirichlet(key)
            keys = actdiri.keys()
            keys.sort()
            for k in keys:
                lines.append( '%9i %10f' % (k, actdiri[k]) + '\n')
        return lines


    def setNeumann(self, key, value):
        self.Neumann[key] = value

    def getNeumann(self, key):
        return self.Neumann[key]

    def extract_neumann(self, linelist):
        words = string.split(linelist[1])
        neumcount = map( string.atoi, words )
        linelist = linelist[2:]
        for ineum in range(len(neumcount)):
            nneum = neumcount[ineum]
            actneum = {}
            for i in range(nneum):
                line = linelist[i]
                words = string.split(line)
                node  = string.atoi(words[0])
                value = string.atof(words[1])
                actneum[node] =  value
            self.setNeumann(ineum+1, actneum)
            linelist = linelist[nneum:]
        self.updateInfo()

    def compose_neumann(self):
        lines = []
        lines.append('neumann\n')
        neumcount = []
        neumkeys = self.Neumann.keys()
        neumkeys.sort()
        for key in neumkeys:
            neumcount.append(len(self.getNeumann(key)))
        lines.append( '%9i'*len(neumcount) % tuple(neumcount) + '\n')
        for key in neumkeys:
            actneum = self.getNeumann(key)
            keys = actneum.keys()
            keys.sort()
            for k in keys:
                lines.append( '%9i %10f' % (k, actneum[k]) + '\n')
        return lines

    def setVariables(self, key, value):
        self.Variables[key] = value

    def getVariables(self, key):
        return self.Variables[key]

    def extract_variables(self, linelist):
        linelist = linelist[1:]
        ivar = 0
        for line in linelist:
            ivar = ivar+1
            self.setVariables(ivar, line[:-1])

    def compose_variables(self):
        lines = []
        lines.append('variables\n')
        varkeys = self.Variables.keys()
        varkeys.sort()
        for key in varkeys:
            lines.append(self.getVariables(key)+'\n')
        return lines

    def compose_calculate(self):
        return 'calculate\n'

### Test

if __name__ == '__main__':

    m = FEModel()
    ff = FemtoolFile(m)
    infilename  = os.path.join(feval.__path__[0],'data','femtool','test1.dat')
    outfilename = os.path.join(feval.__path__[0],'data','femtool','test1_out.dat')
    ff.readFile(infilename)

    ff.updateInfo()
    ff.setWrite('initialize')
    ff.setWrite('elements')
    ff.setWrite('points')

    ff.writeFile(outfilename)


