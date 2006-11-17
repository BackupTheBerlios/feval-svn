import re
import Numeric; N = Numeric

import feval.FEval
from feval.FETextFile import *

def marcatof( s ):
    """Convert a MARC Fortran String of the form
    '1.3345+02' into a floating number
    """
    ss = s[2:]
    ss=re.sub(r'(\+)','e+', ss)
    ss=re.sub(r'(-)','e-', ss)
    return string.atof(s[:2]+ss)

def marcftoa( x ):
    """Convert a MARC floating number into the form ' 1.3345+02' 
    """
    s = '%10.4e' %x
    s = '  ' + s[:-4]+s[-3:]
    return s[-10:]


class MarcFile(FETextFile):
    """Parse a MARC input-file
    """
    type = 'marc'

    # Dictionary of the Element types
    import feval.fecodes.marc.MarcShapeFunctions as MarcShapeFunctions
    shapeFunctionDict = MarcShapeFunctions.MarcShapeFunctionDict

    def __init__(self, model):
        FETextFile.__init__(self, model)
        self.elemtypes = {}
        # Pairs of characters used to mark comments
        self.Comment = [( '$', '\n' )]
        self.extended = False
        # Handling normal/extended precision
        self.connL, self.coordNameL, self.coordValueL = 5, 5, 10


    def extract_title(self, linelist):
        self.model.setName(linelist[0][6:])

    def extract_extended(self, linelist):
        self.extended = True
        self.connL, self.coordNameL, self.coordValueL = 10, 10, 20

    def compose_title(self):
        return ['title '+self.model.name]

    def extract_connectivity(self, linelist):
        """Empty lines are excluded in FETextFile, so we start with
        the next line """
        linelist = linelist[1:]
        # cope with fixed format
        cl = self.connL
        for line in linelist:
            id = int(line[:cl])
            et = int(line[cl:2*cl])
            eletype = self.shapeFunctionDict[et][0]
            nnodes  = self.shapeFunctionDict[et][2]
            self.elemtypes[eletype] = et

            # we do the nnodes trick to cut off pressure nodes (Hermann elements)
            # this is more elegant 
            c = []
            for i in range(nnodes):
                c.append( int(line[(2+i)*cl:(3+i)*cl]) )
            c = N.asarray(c)
            self.model.setElementConn(id, eletype, c)

    def compose_connectivity(self):
        lines = []
        lines.append('connectivity\n')
        lines.append('\n')
        ckeys = self.model.getElementNames()
        ckeys.sort()
        for id in ckeys:
            conn = self.model.Conn[id]
            line = '%5i%5i' % (id, self.elemtypes[conn[0]])
            nnod = len(conn[1])
            ll = '%5i'*nnod % tuple(conn[1])
            lines.append(line+ll+'\n')
        return lines

    def extract_coordinates(self, linelist):
        """fixed format"""
        linelist = linelist[2:]
        cn, cl = self.coordNameL, self.coordValueL
        for line in linelist:
            c = [line[0:cn], line[cn:cn+cl], line[cn+cl:cn+2*cl], line[cn+2*cl:cn+3*cl]]
            id = string.atoi(c[0])
            coord = N.array( map( marcatof, c[1:] ) )
            self.model.setCoord( id, coord )

    def compose_coordinates(self):
        """fixed format"""
        lines = []
        lines.append('coordinates\n')
        ckeys = self.model.getCoordinateNames()
        ckeys.sort()
        coord = self.model.Coord[ckeys[0]]
        ndir = len(coord)
        line = '%5i%5i' % (ndir, len(ckeys))
        lines.append(line + '\n')
        for id in ckeys:
            coord = self.model.Coord[id]
            ll = map( marcftoa, coord )
            line = '%5i' % (id) + string.join(ll, '') +'\n'
            lines.append(line)
        return lines

    def extract_define(self, linelist):
        """extract sets, simple and not yet complete"""
        words = string.split(linelist[0])
        type, name = words[1], words[3]
        if type == 'ndsq':
            type = 'node'
        linelist = linelist[1:]
        set = []
        for line in linelist:
            line = string.replace(line.lower(), 'c', '')
            words = map(string.atoi, (string.split(line)))
            set.extend(words)
        ## print type, name, len(linelist), len(set)
        self.model.setSet( type, name, set )


### Test

if __name__ == '__main__':

    infilename  = os.path.join( feval.__path__[0], 'data', 'marc', 'test1.dat' )
    outfilename = os.path.join( feval.__path__[0], 'data', 'marc', 'test1_out.dat' )

    m = feval.FEval.FEModel()
    mf = MarcFile(m)
    mf.readFile('/home/tinu/projects/jako/marc/jako3d_isothermal.dat')
    #mf.readFile(infilename)

    mf.setWrite('title')
    mf.setWrite('coordinates')
    mf.setWrite('connectivity')

    mf.writeFile(outfilename)

#      a= m.getSet('node','bed')
#      if a:
#          a.sort()
#      print a
