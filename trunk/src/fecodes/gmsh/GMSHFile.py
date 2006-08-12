import re
import Numeric; N = Numeric
from feval.FEval import *
from feval.FETextFile import *

class GMSHFile(FETextFile):
    """Parse a GMSH input-file
    Take care of the node pattern which is different from the MARC node
    pattern which  is applied throughout these routines.
    (The MARC node pattern is: first the corner nodes, then the side
    nodes, the Femtool pattern is corner-side-corner-side...)

    The node pattern is saved in the list |Variables|
    """
    type = 'gmsh'

    # Dictionary of the Element types
    shapeFunctionDict = { 1: 'Line1', \
			  2: 'Tria3',
			  3: 'Quad4',
			  4: 'Tet4',
			  5: 'Hex8',
			  6: 'Prism6',
			  7: 'Pyram5' }

    # inverse dictionary of the Element types                               
    elementTypeDict = {}                                                    
    for k,v in shapeFunctionDict.items():                                   
        elementTypeDict[v] = k
        
    # the node pattern is from GMSH to Marc
    nodePattern = {}
    nodePattern['Quad4'] = [0,1,2,3]
    nodePattern['Tet4']  = [0,1,2,3]
    nodePattern['Quad8'] = [0,2,4,6,1,3,5,7]

    # the inverse node pattern is from Marc to GMSH
    nodePatternInv = {}
    nodePatternInv['Quad4'] = [0,1,2,3]
    nodePatternInv['Tet4']  = [0,1,2,3]
    nodePatternInv['Quad8'] = [0,4,1,5,2,6,3,7]

    def __init__(self, model):
	FETextFile.__init__(self, model)

	# Pairs of characters used to mark comments
	self.Comment = [( '#', '\n' )] 

    def isKeyword(self, line):
	"""check wether the line could contain a magic word
	return None if the line is not interesting
	In GMSH the magic words start with "$", so look for this
	"""
	if line[0] == '$':
	    return 1, string.lower(line[1:])
	else:
	    return 0, line

    def extract_nodes(self, linelist):
        print 'hallo extract'
        self.extract_nod(linelist)

    def extract_nod(self, linelist):
        """
        omit the first and the last line which contain $ELM and $ENDELM
        """
        linelist = linelist[2:]
        for line in linelist:
            words = string.split(line)
            id = string.atoi(words[0])
            c = N.array( map(string.atof, words[1:]))
            self.model.setCoord( id, c )

    def compose_nod(self):
        lines = []
        lines.append('$NOD\n')
        ckeys = self.model.Coord.keys()
        ckeys.sort()
        for id in ckeys:
            coord = self.model.Coord[id]
            ndim = len(coord)
            line = '%5i' % (id) + ' %f'*ndim \
                   % tuple(coord.tolist()) + '\n'
            lines.append(line)
        lines.append('$ENDNOD\n')
        return lines

    def extract_elements(self, linelist):
        self.extract_elm(linelist)


    def extract_elm(self, linelist):
	linelist = linelist[2:]
	for line in linelist:
	    words = map( string.atoi, string.split(line) )
	    id = words[0]
	    c = N.array(words[1:])
	    try: 
		elemtype = self.shapeFunctionDict[c[0]]
	    except:
		elemtype = ''
		print '**** Element type %i not defined!' % (c[0])
	    nodes = N.take( c[4:], self.nodePattern[ elemtype ] )
	    self.model.setElementConn( id, elemtype, nodes )

    def compose_elm(self):
	lines = []
	lines.append('$ELM\n')
	ckeys = self.model.Conn.keys()
	ckeys.sort()
	for id in ckeys:
	    conn = self.model.Conn[id]
	    nodes = conn[1]
  	    nodes = N.take(nodes, self.nodePatternInv[ conn[0] ] )
            type = self.elementTypeDict[conn[0]]
  	    nnod = len(nodes)
  	    line = '%5i%5i%5i%5i%5i' % (id, type, 1, 1, nnod)
	    ll = '%5i'*nnod % tuple(nodes)
	    lines.append(line+ll+'\n')
	lines.append('$ENDELM\n')
	return lines

### Test

if __name__ == '__main__':

    m = FEModel()
    ff = GMSHFile(m)
    infilename  = '/home/tinu/projects/wrangell/libmesh/sphere.msh'
    #infilename  = os.path.join(feval.__path__[0],'data','gmsh','simple.gmsh')
    outfilename = os.path.join(feval.__path__[0],'data','gmsh','simple_out.gmsh')
    ff.readFile(infilename)

    ff.setWrite('nod')
    ff.setWrite('elm')

    ff.writeFile(outfilename)


