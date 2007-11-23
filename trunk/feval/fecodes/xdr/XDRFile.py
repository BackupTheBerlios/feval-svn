import re
import Numeric; N = Numeric
from feval.FEval import *
from feval.FETextFile import *

class XDRFile(FETextFile):
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
        FETextFile.__init__(self, model)
        self.elemNumber = 0
        self.nodeNumber = 0
        self.versionInfo = [0, 0]
        self.fileType    = ''
        # Pairs of characters used to mark comments
        self.Comment = [( '#', '\n' )]


    def readFile(self, filename):
        """Parse an input file, call the magic handlers
        """
        self.FileName = filename
        self.File = open(filename, 'r')
        linelist = self.File.readlines()
        if linelist[0].startswith('DEAL'):
            self.extract_deal(linelist)

    def extract_deal(self, linelist):
        """Extract all data from a file in DEAL fromat
        We do the whole parsing here, because there are no magic words.
        First we read the header
        Then come the data blocks
        """
        # here we are in the deal header
        self.versionInfo = map( int, linelist[0].split()[1].split(':') )
        self.fileType    = 'Deal'
        #nelem = int(linelist[1].split()[0])           # number of elements
        nnodes     = int(linelist[2].split()[0])           # number of nodes
        elemblocks = int(linelist[6].split()[0])      # Num. Element Blocks.
        elemtypes  = linelist[7].split()[:elemblocks] # Element types in each block.
        neleblock  = map( int, linelist[8].split()[:elemblocks]) # Num. of elements in each block.
        if not linelist[10].startswith('Title String'):
            self.model.name = linelist[10][:-1].strip()
        # and now we get the data
        #first comes the element data
        ll = 10
        for elem, nelem in zip(elemtypes, neleblock):
            if not self.shapeFunctionDict.has_key(elem):
                print 'No shape function found for element type', elem
            elemtype = self.shapeFunctionDict[elem]
            for l in range(nelem):
                ll += 1
                self.elemNumber += 1
                nodes = map(int, linelist[ll].split())
                self.model.setConn(self.elemNumber,
                                   elemtype,
                                   nodes )
        # then the node data (coordinates)
        for node in range(nnodes):
            ll += 1
            self.nodeNumber += 1
            coord = map( float, linelist[ll].split() )
            self.model.setCoord( self.nodeNumber, coord )

    def compose_deal(self):
        """Write a DEAL file"""
        lines = []

        # get all element keys and sort them
        ekeys = self.model.getElementNames()
        ##ekeys.sort()
        # go through all elements and see which ones have the same type
        weight = 0 # total weight
        elines = {}
        for id in ekeys:
            etype, conn = self.model.getElementConn(id)
            if not elines.has_key(etype):
                elines[etype] = []
            nnod = len(conn)
            elines[etype].append('%d '*nnod % tuple(conn) + '\n')
            weight += nnod

        # get all coordinate keys and sort them
        ckeys = self.model.getCoordinateNames()
        ##ckeys.sort()
        coord = self.model.Coord[ckeys[0]]
        ndir = len(coord)
        clines = []
        # write all coordinates
        for id in ckeys:
            clines.append('%f '*ndir % tuple(self.model.getCoordinate(id)) + '\n')

        # now we are ready to write the file
        lines.append('DEAL 003:003\n')
        lines.append('%d\t # Num. Elements\n' % len(ekeys) )
        lines.append('%d\t # Num. Nodes\n'    % len(ckeys) )
        lines.append('%d\t # Sum of Element Weights\n' % weight)
        n_bound = 0
        if self.magicID+'boundary' in self.verbatimData:
            n_bound = len( self.boundary_elems )
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
        if not self.model.name.startswith('FE-model'):
            title = self.model.name.strip()
        lines.append(title+'\n')

        # add the element lines
        for t in elines.keys():
            lines.extend(elines[t])
        lines.extend(clines)
        return lines

    def setBoundaryElems(self, boundary_elems):
        """Set the list of boundary elements with the boundary side
        [(elem, side)]"""
        self.boundary_elems = boundary_elems

    def compose_boundary(self):
        """Write the boundary"""
        lines = []
        lines.append('\n')    # empty line
        for be in self.boundary_elems:
            lines.append( '%d %d %d \n' % (be[0], be[1], 0) )
        return lines


### Test

if __name__ == '__main__':

    m = FEModel()
#     infilename  = os.path.join( feval.__path__[0], 'data', 'xdr', 'mesh.xda' )
#     outfilename = os.path.join( feval.__path__[0], 'data', 'xdr', 'ca_out.xda' )

#     from feval.fecodes.marc.MarcT16File import *  
    # inputfilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal_S10_E5.t16')
#     inputfilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal_S10_E5.t16')
#     mf = MarcT16File(m, inputfilename)
#     mf.readInc(3)

#    mf = MarcT16File(m, '/home/tinu/projects/colle/marc/cgx8_dc0.90/cgx8dec.t16')


    from feval.fecodes.marc.MarcFile import *  
    mf = MarcFile(m)
    mf.readFile('../../../data/marc/test1.dat')
#    mf.readFile('/home/tinu/projects/jako/marc/jako3d_isothermal.dat')

    m.renumberNodes(base=0)
    m.renumberElements(base=0)
    bednodes = m.Sets['node']['bed']
    elemsides = m.findElementsFromNodes(bednodes)
    gf = XDRFile(m)

    gf.setBoundaryElems(elemsides)
    gf.setWrite('deal')
    gf.setWrite('boundary')
    gf.writeFile('jako.xdr')






