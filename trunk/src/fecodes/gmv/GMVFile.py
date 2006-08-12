import re
import Numeric; N = Numeric
from feval.FEval import *
from feval.FETextFile import *

class GMVFile(FETextFile):
    """Parse a GMV file
    """
    type = 'gmv'

    # Dictionary of the Element types
    shapeFunctionDict = {'hex 8':  'Hex8',
                         'hex8 8': 'Hex8',
                         'phex8 8': 'Hex8',
                         'tri 3':  'Tri3',
                         'quad 4': 'Quad4',
                         'pprism6': 'Prism6',
                         '8quad 8': 'Quad8',
                         'phex20 20': 'Hex20',
                         }

    # inverse dictionary of the Element types
    invShapeFunctionDict = {}
    for k,v in shapeFunctionDict.items():
        invShapeFunctionDict[v] = k

    def __init__(self, model):
        FETextFile.__init__(self, model)
        self.elemNumber = 0
        self.inVarBlock = 0
        self.nodVar     = []
        self.nodVarInfo = []
        # Pairs of characters used to mark comments
        self.Comment = [( '', '' )]


    def findMagic(self, linelist):
        """this overwrites TextFile.findMagic 
        check whether a magic word occurs
        """
        words = string.split(linelist[0])
        # try to match the magic word xxxx and
        # execute the corresponding handler (extract_xxxx)
        # if this is not possible: keep the input lines in a dictonary
        magicKey = string.lower(words[0])

        # test whether we are in the variable block,
        # i.e. between variable and endvars
        if magicKey == 'variable':
            self.inVarBlock = 1
        if magicKey == 'endvars':
            self.inVarBlock = 0
            self.finish_variable()
        
        if magicKey in self.MagicWords:
            try:
                fct = eval('self.extract_'+string.lower(magicKey))
            except:
                fct = None
            if fct:
                map( fct, [linelist] )
        elif self.inVarBlock:
            self.extract_variables(linelist)

    def extract_gmvinput(self, linelist):
        """read the header, do nothing"""
        pass 

    def compose_gmvinput(self):
        """print the header"""
        return ['gmvinput ascii\n\n']

    def compose_endgmv(self):
        """print the footer"""
        return ['\nendgmv\n']

    def extract_nodes(self, linelist):
        """Extract the nodes.
        The number of nodes is given in the first line of the input file.
        The nodes are numbered consecutively here."""
        nnodes = int(linelist[0].split()[1])
        coord = []
        linelist = linelist[1:]
        for line in linelist:
            coord.extend( map( float, line.split() ) )
        coord = N.asarray(coord)
        coord.shape = (len(coord)/nnodes, nnodes)
        for id in range(nnodes):
            self.model.setCoord( id+1, coord[:,id] )

    def compose_nodes(self):
        """Write nodes in increasing order, one line per coordinate"""
        lines = []
        # get all coordinate keys and sort them
        ckeys = self.model.getCoordinateNames()
        ckeys.sort()
        coord = self.model.Coord[ckeys[0]]
        ndir = len(coord)
        lines.append('nodes %i \n' % (len(ckeys)))
        # a loop per coordinate direction
        for d in range(ndir):
            for id in ckeys:
                lines.append('%f \n' % self.model.Coord[id][d])
        return lines

    def extract_element(self, linelist):
        """Extract the cells of any type
        The cells are numbered consecutively with self.elemNumber"""
        elemtype = linelist[0].strip()
        nodes = map(int, linelist[1].split())
        self.elemNumber += 1
        try: 
            elemtype = self.shapeFunctionDict[elemtype.lower()]
        except:
            print '**** Element type %i not defined!' % (elemtype)
            elemtype = ''
        self.model.setElementConn(self.elemNumber,
                                  elemtype,
                                  nodes )

    def extract_hex(self, linelist):
        """Extract the cells of type hex"""
        self.extract_element(linelist)

    def extract_hex8(self, linelist):
        """Extract the cells of type hex8"""
        self.extract_element(linelist)

    def extract_phex8(self, linelist):
        """Extract the cells of type hex8"""
        self.extract_element(linelist)

    def extract_phex20(self, linelist):
        """Extract the cells of type phex20"""
        self.extract_element(linelist)

    def extract_quad(self, linelist):
        """Extract the cells of type quad"""
        self.extract_element(linelist)

    def extract_8quad(self, linelist):
        """Extract the cells of type 8quad"""
        self.extract_element(linelist)

    def compose_cells(self):
        """Write elements in increasing order, one line per coordinate"""
        lines = []
        # get all element keys and sort them
        ckeys = self.model.getElementNames()
        ckeys.sort()
        lines.append('\ncells %i\n' % (len(ckeys)))
        # a loop per coordinate direction
        for id in ckeys:
            conn = self.model.Conn[id]
            nnod = len(conn[1])
            ll = '%i '*nnod % tuple(conn[1])
            line = '%s\n %s\n' % (self.invShapeFunctionDict[conn[0]], ll)
            lines.append(line)
        return lines
        
    def extract_velocity(self, linelist):
        """read the header, do nothing"""
        pass

    def extract_variables(self, linelist):
        """read the variable and collect the values"""
        varname, vartype = linelist[0].split()
        self.nodVarInfo.append(varname.strip())
        nodvar = ''.join(linelist[1:]).split()
        nodvar = map(float, nodvar)
        self.nodVar.append(nodvar)

    def finish_variable(self):
        if self.nodVar:
            nodvar = N.asarray(self.nodVar)
            for id in range(len(nodvar[0])):
                self.model.setNodVar( id+1, nodvar[:, id] )
            self.model.setNodVarInfo(self.nodVarInfo)

    def compose_variable(self, varNames=[], varNamesGMV = []):
        """Write variables (except velocity, which has its own block)"""
        lines = []
        # get all variables
        nodVarInfo = self.model.getNodVarInfo()
        if not varNames:
            varNames = nodVarInfo
        if varNames and len(varNames) > 0:
            # the nodal variables
            lines.append('\nvariable \n')
            allvars = self.model.getNodVarsAsArray()
            for var in varNames:
                idx = nodVarInfo.index(var)
                if var in nodVarInfo:
                    lines.append('%s 1\n' % var)
                    idx = nodVarInfo.index(var)
                    for v in allvars:
                        lines.append('%f ' % v[idx])
                    lines.append('\n')
            lines.append('endvars\n')

        return lines

### Test

if __name__ == '__main__':


#    infilename  = os.path.join( feval.__path__[0], 'data', 'gmv', 'test1.gmv' )
#    infilename  = os.path.join( '/soft/numeric/feval', 'data', 'gmv', 'test1.gmv' )
#     outfilename = os.path.join( feval.__path__[0], 'data', 'gmv', 'test1_out.gmv' )

    m = FEModel()
#     from feval.fecodes.marc.MarcT16File import *  
#     infilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal.t16')
#     #     outfilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal.gmv')
#     outfilename = os.path.expanduser('~/projects/colle/marc/cgx8_dc0.9.gmv')
#     #mf = MarcT16File(m, infilename)
#     # mf = MarcT16File(m, '/home/tinu/projects/colle/marc/cgx8_dc0.90/cgx8dec.t16')
#     mf = MarcT16File(m, '/home/tinu/numeric/marc/newsl_lafiua_tf4.0-tl4.5-tu2.5-g0.2-n4.-a5.3-v100-wl-0.06-0.0-550.0.t16')
#     mf.readInc(-1)

#     stop
#     from feval.fecodes.marc.MarcFile import *  
#     mf = MarcFile(m)
#     mf.readFile('/home/tinu/numeric/marc/ca_fm.dat')


    gf = GMVFile(m)
    #gf.readFile(infilename)
    gf.readFile('/home/tinu/projects/wrangell/libmesh/viscodens3d/out.gmv.001')
    #gf.readFile('/home/tinu/projects/lbie/src/A.gmv')

#     gf.setWrite('gmvinput')
#     gf.setWrite('nodes')
#     gf.setWrite('cells')
#     gf.setWrite('variable')
#     gf.setWrite('endgmv')

#     gf.writeFile(outfilename)




