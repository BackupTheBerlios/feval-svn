import re
import Numeric; N = Numeric
from feval.FEval import *

## note that this file is not based on FETextFile, but is a direct
## reader/writer which seems much simpler

class LibmeshFile(object):
    """Parse and write a XDR file
    This is XDA file type defined in libmesh
    """
    type = 'xdr'

    # Dictionary of the Element types (from enum_elem_type.h)
    shapeFunctionDict = {
                         '3':  'Tri3',
                         '6':  'Tri6',
                         '5':  'Quad4',
                         '6':  'Quad8', 
                         '7':  'Quad9', 
                         '10':  'Hex8', 
                         '11':  'Hex20', 
                         '12':  'Hex27', 
                         }

    # inverse dictionary of the Element types
    invShapeFunctionDict = {}
    for k,v in shapeFunctionDict.items():
        invShapeFunctionDict[v] = k

    def __init__(self, model):
        self.model = model
        self.elemNumber = 0
        self.nodeNumber = 0
        self.versionInfo = [0, 0]
        self.fileType    = ''
        # Pairs of characters used to mark comments
        self.Comment = [( '#', '\n' )]

    def readFile(self, filename):
        """Parse an input file
        """
        infile = file(filename)
        lines = infile.readlines()
        
        # first read the header
        self.fileType, self.versionInfo = lines[0].split()
        nelem      = int(lines[1].split()[0])      # number of elements
        nnodes     = int(lines[2].split()[0])      # number of nodes
        elemblocks = int(lines[6].split()[0])      # Num. Element Blocks.
        elemtypes  = lines[7].split()[:elemblocks] # Element types in each block.
        nlevels = len(lines[8].split('#')[0].split())/elemblocks
        neleblock  = map( int, lines[8].split()[:elemblocks]) # Num. of elements in each block.
        if not lines[10].startswith('Title String'):
            self.model.name = lines[10][:-1].strip()

        # and now we get the data
        # first comes the element data
        ll = 10
        for elem, nelem in zip(elemtypes, neleblock):
            if not self.shapeFunctionDict.has_key(elem):
                print 'No shape function found for element type', elem
            elemtype = self.shapeFunctionDict[elem]
            for l in range(nelem):
                ll += 1
                self.elemNumber += 1
                nodes = map(int, lines[ll].split())
                self.model.setElement(self.elemNumber, elemtype, nodes )

        # then the node data (coordinates)
        for node in range(nnodes):
            ll += 1
            self.nodeNumber += 1
            coord = map( float, lines[ll].split() )
            self.model.setCoordinate( self.nodeNumber, coord )

    def writeFile(self, filename, bcond = None):
        """Write a Libmesh xda file
        |bcond| should be a list of lists
        [[e1, s1, b1],
         [e2, s2, b2]]
        where e* are element names, s* side ids, and b* boundary condition ids
        """
        lines = []

        # go through all elements and see which ones have the same type
        # build a dict of elements
        weight = 0 # total weight
        elines = {}
        for ele in self.model.getElementNames():
            etype, conn = self.model.getElementConn(ele)
            nnod = len(conn)
            elines.setdefault(etype, []).append('%d '*nnod % tuple(conn) + '\n')
            weight += nnod

        # now we are ready to write the file
        lines.append('LIBM 0\n')
        lines.append('%d\t # Num. Elements\n' % len(self.model.Conn) )
        lines.append('%d\t # Num. Nodes\n'    % len(self.model.Coord) )
        lines.append('%d\t # Sum of Element Weights\n' % weight)
        if bcond:
            n_bound = len( bcond )
        else:
            n_bound = 0
        lines.append('%d\t # Num. Boundary Conds.\n' % n_bound)
        lines.append('65536\t # String Size (ignore)\n') 
        lines.append('%d\t # Num. Element Blocks.\n' % len(elines.keys()) )
        etypes, eblocks = '', ''
        for t in elines.keys():
            etypes  += '%s ' % self.invShapeFunctionDict[t]
            eblocks += '%s ' % len(elines[t])
        lines.append('%s\t # Element types in each block.\n' % etypes)
        lines.append('%s\t # Num. of elements in each block.\n' % eblocks)
        lines.append('Id String\n')
        title = 'Title String'
        if not self.model.name.startswith('FE-Model'):
            title = self.model.name.strip()
        lines.append(title+'\n')

        # add the element lines
        for t in elines.keys():
            lines.extend(elines[t])
        del elines

        # write all coordinates
        ckeys = sorted(self.model.getCoordinateNames())
        ndir = len(self.model.Coord[ckeys[0]])
        for id in ckeys:
            lines.append('%f '*ndir % tuple(self.model.getCoordinate(id)) + '\n')

        if bcond:
            for bc in bcond:
                lines.append('%d %d %d\n' % tuple(bc))

        file(filename,'w').writelines(lines)
        

### Test
if __name__ == '__main__':

    m = FEModel()
    mf = LibmeshFile(m)
    mf.readFile('/soft/numeric/libmesh/examples/ex10/mesh.xda')
    bc = [[1,5, 0],
          [2,6, 1000]]
    mf.writeFile('/soft/numeric/libmesh/examples/ex3/a2.xda', bcond = bc)
    
#     infilename  = os.path.join( feval.__path__[0], 'data', 'xdr', 'mesh.xda' )
#     outfilename = os.path.join( feval.__path__[0], 'data', 'xdr', 'ca_out.xda' )

#     from feval.fecodes.marc.MarcT16File import *  
    # inputfilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal_S10_E5.t16')
#     inputfilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal_S10_E5.t16')
#     mf = MarcT16File(m, inputfilename)
#     mf.readInc(3)

#    mf = MarcT16File(m, '/home/tinu/projects/colle/marc/cgx8_dc0.90/cgx8dec.t16')


#     from feval.fecodes.marc.MarcFile import *  
#     mf = MarcFile(m)
#     mf.readFile('../../../data/marc/test1.dat')
# #    mf.readFile('/home/tinu/projects/jako/marc/jako3d_isothermal.dat')

#     m.renumberNodes(base=0)
#     m.renumberElements(base=0)
#     bednodes = m.Sets['node']['bed']
#     elemsides = m.findElementsFromNodes(bednodes)
#     gf = XDRFile(m)

#     gf.setBoundaryElems(elemsides)
#     gf.setWrite('deal')
#     gf.setWrite('boundary')
#     gf.writeFile('jako.xdr')






