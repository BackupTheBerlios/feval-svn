import re
import Numeric; N = Numeric
from feval.FEval import *
from feval.FETextFile import *

class UNVFile(FETextFile):
    """Parse a UNV file
    """
    type = 'unv'

    # Dictionary of the Element types
    # this is mainly based on the libmesh assignments
    shapeFunctionDict = {
                         41:  'Tri3',
                         91:  'Tri3',
                         42:  'Tri6',
                         92:  'Tri6',
                         44:  'Quad4',
                         94:  'Quad4',
                         45:  'Quad8',
                         95:  'Quad8',
                         56:  'Quad9',
                         111: 'Tet4',
                         112: 'Prism6',
                         115: 'Hex8', 
                         116: 'Hex20', 
                         118: 'Tet10', 
                         }

    # inverse dictionary of the Element types
    invShapeFunctionDict = {}
    for k,v in shapeFunctionDict.items():
        invShapeFunctionDict[v] = k

    # the node pattern is from UNV to Marc
    nodePattern = {}
    nodePattern['Hex20'] = [0,1,2,3, 4,5,6,7, 8,9,10,11, 16,17,18,19, 12,13,14,15]

    # the node pattern is from Marc to UNV
    nodePatternInv = {}
    nodePatternInv['Hex20'] = [0,1,2,3, 4,5,6,7, 8,9,10,11, 16,17,18,19, 12,13,14,15]

    def __init__(self, model):
        FETextFile.__init__(self, model)
        self.dataInfo = {}
        ## Block delimiter
        self.blockDelimiter = '-1'
        self.blockOpen = 0

    def isKeyword(self, line):
        """this overrides the standard keyword search
        if -1 then this is a new block
        We only have a new block after a blockDelimiter has already been found
        """
        if self.blockDelimiter == line.strip():
            if self.blockOpen:
                self.blockOpen = 0
                return 1, line
            else: 
                self.blockOpen = 1
                return 0, line
        else:
            return 0, line

    def findMagic(self, linelist):
        """this overwrites TextFile.findMagic 
        check whether a magic word occurs
        """
        if len(linelist) > 1:
            magicKey = linelist[2].strip()
            if magicKey in self.MagicWords:
                try:
                    fct = eval('self.extract_'+string.lower(magicKey))
                except:
                    fct = None
                if fct:
                    map( fct, [linelist] )

    def extract_2411(self, linelist):
        """read block 2411 (Nodes)"""
        for i in xrange((len(linelist)-3)/2):
            #  label,export,displ,color = linelist[3+2*i].split()
            label  = int( linelist[3+2*i].split()[0] )
            coord = map(float, linelist[4+2*i].replace('D','E').split())
            self.model.setCoord( label, N.asarray(coord) )

    def compose_2411(self):
        """compose block 2411 (Nodes)"""
        # get all coordinate keys and sort them
        ckeys = self.model.getCoordinateNames()
        if len(ckeys) == 0:
            return
        ckeys.sort()
        lines = []
        lines.append('    -1\n')
        lines.append('  2411\n')
        for id in ckeys:
            tmp = '%10d  '*4 + '\n'
            line = tmp % (id, 1, 1, 11)
            lines.append(line)
            line = '%25.16E%25.16E%25.16E\n' % tuple(self.model.Coord[id])
            lines.append(line)
        lines.append('    -1\n')
        return lines

    def extract_2412(self, linelist):
        """read block 2412 (Elements)
        reads only non-beam models (2 records per element)"""
        for i in xrange((len(linelist)-3)/2):
            label,type,physp,matp,color,nnodes = map(int, linelist[3+2*i].split() )
            nodes = map(int, linelist[4+2*i].split())
            type = self.shapeFunctionDict[type]
            self.model.setElement( label, type, N.asarray(nodes) )

    def compose_2412(self):
        """compose block 2412 (Elements)"""
        # get all coordinate keys and sort them
        ekeys = self.model.getElementNames()
        ekeys.sort()
        if len(ekeys) == 0:
            return
        lines = []
        lines.append('    -1\n')
        lines.append('  2412\n')
        for id in ekeys:
            elem = self.model.getElementConn(id)
            type = elem[0]
            ntype = self.invShapeFunctionDict[type]
            nnodes = len(elem[1])
            # find the node pattern for this element type (most are trivial)
            if self.nodePatternInv.has_key(type):
                nodepattern = N.asarray( self.nodePatternInv[type] )
            else:
                nodepattern = N.arange(nnodes)
            nodes = N.take(elem[1], nodepattern )
            # the last 
            tmp = '%8d  '*6 + '\n'
            line = tmp % (id, ntype, 0, 0, 1, nnodes)
            lines.append(line)
            tmp = '%8d  '*nnodes + '\n'
            line = tmp % tuple(nodes)
            lines.append(line)
        lines.append('    -1\n')
        return lines

    def extract_2414(self, linelist):
        """read block 2414 (Analysis data)
        """
        self.dataInfo['label'] = int(linelist[3])
        self.dataInfo['name']  = linelist[4]
        location = int(linelist[5])
        self.dataInfo['description'] = linelist[6:11]
        model,analy,datachar,restype,datatype,nvaldc = map(int, linelist[11].split())
        self.dataInfo['intdata1']  = map(int, linelist[12].split())
        self.dataInfo['intdata2']  = map(int, linelist[13].split())
        self.dataInfo['realdata1'] = map(float, linelist[14].split())
        self.dataInfo['realdata2'] = map(float, linelist[15].split())
        
        # conversion of datatype  (1: int, 2: real, 4: double)
        dataConv = (int,float,None,float)[datatype-1]
        # data at nodes
        if location == 1:    
            for i in xrange((len(linelist)-16)/2):
                id     = int(linelist[16+2*i])
                nodvar = N.asarray(map(dataConv, linelist[17+2*i].split()))
                self.model.setNodVar( id, nodvar )
        else:
            print 'location type %d not implemented\nOnly nodal datasets are implemented' % location
        self.model.setNodVarInfo( map(str, range(nvaldc)) )

            
    def compose_2414(self):
        """compose block 2414 (Analysis data)
        """
        ckeys = self.model.NodVar.keys()
        if len(ckeys) == 0:
            return
        ckeys.sort()
        lines = []
        lines.append('    -1\n')
        lines.append('  2414\n')
        lines.append(str( self.dataInfo.get('label', '     0') )+'\n')
        lines.append(self.dataInfo.get('name', 'FEval default title')+'\n')
        lines.append('     1\n')   # nodal variable (called location) is encoded by '1'
                                   # element values would be '2' (see documentation)
        lines.extend(self.dataInfo.get('description', ['FEval default description\n']*5))
        # datatype: '2': single precision, '4': double precision
        # nvaldc:   number of data values per node/element
        datatype = 2 
        nvaldc = len(self.model.getNodVarInfo())
        tmp = '%8d  '*6+'\n'
        line = tmp % (0, 0, 0, 0, datatype, nvaldc)
        lines.append(line)
        lines.append(self.dataInfo.get('intdata1', '     0'*8)+'\n')
        lines.append(self.dataInfo.get('intdata2', '     0'*8)+'\n')
        lines.append(self.dataInfo.get('realdata1', '  0.0000E+00'*6)+'\n')
        lines.append(self.dataInfo.get('realdata2', '  0.0000E+00'*6)+'\n')
        for id in ckeys:
            lines.append('    %d\n' % id)
            nodvar = self.model.NodVar[id]
            tmp = '%e  '*len(nodvar) + '\n'
            line = tmp % tuple(nodvar)
            lines.append(line)
        lines.append('    -1\n')
        return lines




### Test

if __name__ == '__main__':

    datafilename  = os.path.join( feval.__path__[0], 'data', 'unv', 'data_second_with_header_out.unv' )
    infilename  = os.path.join( feval.__path__[0], 'data', 'unv', 'pipe-mesh.unv' )
    outfilename = os.path.join( feval.__path__[0], 'data', 'unv', 'cgx8.unv' )

    m = FEModel()
    from feval.fecodes.marc.MarcT16File import *  
    ##poly_inputfilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal_S10_E5.t16')
    ##mf = MarcT16File(m, poly_inputfilename)
    mf = MarcT16File(m, '/home/tinu/projects/colle/marc/colle1_5o/colle1_5o.t16')
    ##mf = MarcT16File(m, '/home/tinu/numeric/marc/newsl_lafiua_tf4.0-tl4.5-tu2.5-g0.2-n4.-a5.3-v100-wl-0.06-0.0-550.0.t16')
    mf.readInc(-1)


#     from feval.fecodes.marc.MarcFile import *  
#     mf = MarcFile(m)
#     mf.readFile('/home/tinu/numeric/marc/ca_fm.dat')


    uf = UNVFile(m)
#     uf.readFile(infilename)
#     udf = UNVFile(m)
#     udf.readFile(datafilename)

    uf.setWrite('2414')
    uf.setWrite('2411')
    uf.setWrite('2412')
    #    uf.writeFile(outfilename)
    uf.writeFile('/home/tinu/projects/colle/marc/colle1_5o/colle1_5o.unv')

#     m.renumberElements()
#     m.renumberNodes()
#     from feval.fecodes.gmv.GMVFile import *  
#     gf = GMVFile(m)
#     gf.setWrite('gmvinput')
#     gf.setWrite('nodes')
#     gf.setWrite('cells')
#     gf.setWrite('variable')
#     gf.setWrite('endgmv')
    
#     gf.writeFile(outfilename)




