import re
import Numeric; N = Numeric
from feval.FEval import *
from feval.FETextFile import *

class TochnogFile(FETextFile):
    """Parse a Tochnog input-File
    """
    type = 'tochnog'

    # Dictionary of the Element types
    shapeFunctionDict = { '-hex8' : 'Hex8' , \
                          '-quad4': 'Quad4', \
                          '-quad8': 'Quad8', \
                          '-bar2' : 'Bar2', \
                          '-tet4' : 'Tet4', \
                          '-tria3': 'Tria3' }

    # inverse dictionary of the Element types
    invShapeFunctionDict = {}
    for k,v in shapeFunctionDict.items():
        invShapeFunctionDict[v] = k

    def __init__(self, model, mode='r'):
        FETextFile.__init__(self, model)
        self.nodes_already_written = None
        self.elements_already_written = None
        # Pairs of characters used to mark comments
        Comment = [( '(', ')' )] 

    def extract_echo(self, linelist):
        words = string.split(linelist[0])
        self.echo = string.strip(words[1])

    def compose_echo(self):
        lines = ['echo %s' % self.echo ]
        return lines

    def extract_number_of_space_dimensions(self, linelist):
        words = string.split(linelist[0])
        self.dim = string.atoi(words[1])

    def compose_number_of_space_dimensions(self):
        lines = ['number_of_space_dimensions %1i' % self.dim + '\n']
        self.nodes_already_written = None
        self.elements_already_written = None
        return lines

    def extract_node(self, linelist):
        """Reads one single node"""
        words = string.split(linelist[0])
        id = string.atoi(words[1])
        c = N.array( map(string.atof, words[2:]))
        self.model.setCoord( id, c )

    def compose_node(self):
        """writes all nodes at once"""
        self.nodes_already_written = 1
        lines = []
        ckeys = self.model.Coord.keys()
        ckeys.sort()
        # get the number of space dimensions
        ndim = len( self.model.Coord[ckeys[0]] )
        for id in ckeys:
            coord = self.model.Coord[id]
            line = 'node %5i' % (id) + ' %f'*ndim \
                   % tuple(coord.tolist()) + '\n'
            lines.append(line)
        lines.append('\n')
        return lines


    def extract_element(self, linelist):
        words = string.split(linelist[0])
        id = string.atoi(words[1])
        c = N.array( map( string.atoi, words[3:]) )
        try: 
            elemtype = self.shapeFunctionDict[words[2]]
        except:
            elemtype = ''
            print '**** Element type %i not defined!' % (c[0])
        self.model.setConn( id, elemtype, c)

    def compose_element(self):
        """writes all nodes at once"""
        if self.elements_already_written: return []
        self.elements_already_written = 1
        lines = []
        ckeys = self.model.Conn.keys()
        ckeys.sort()
        conn = self.model.Conn[ckeys[0]]
        nnode = len(conn[1])
        for id in ckeys:
            conn = self.model.Conn[id]
            line = 'element %5i  %s ' \
                   % (id, self.invShapeFunctionDict[conn[0]]) + \
                   ' %i'*nnode % tuple(conn[1].tolist()) + '\n'
            lines.append(line)
        lines.append('\n')
        return lines

### Test

if __name__ == '__main__':

    m = ModelData()
    tn = TochnogFile(m)

    infilename  = os.path.join( feval.__path__[0],'data','tochnog','axisym1.dat' )
    outfilename = os.path.join( feval.__path__[0],'data','tochnog','axisym1_out.dat' )

    tn.readFile(infilename)

    tn.setWrite('number_of_space_dimensions')
    tn.setWrite('echo')
    tn.setWrite('element')
    tn.setWrite('node')

    tn.writeFile(outfilename)

