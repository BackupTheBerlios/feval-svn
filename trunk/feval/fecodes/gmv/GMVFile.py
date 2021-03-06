import re
import numpy as N
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
                         '6tri 6':  'Tri6',
                         'quad 4': 'Quad4',
                         'pprism6 6': 'Prism6',
                         '8quad 8': 'Quad8',
                         'phex20 20': 'Hex20',
                         'phex27 27': 'Hex27',
                         }

    # inverse dictionary of the Element types
    invShapeFunctionDict = {}
    for k,v in shapeFunctionDict.items():
        invShapeFunctionDict[v] = k

    # the node pattern is from GMV to FEval
    nodePattern = {}
    nodePattern['Tri3']  = [0,1,2]
    nodePattern['Tri6']  = [0,1,2,3,4,5]
    nodePattern['Quad4'] = [0,1,2,3]
    nodePattern['Tet4']  = [0,1,2,3]
    nodePattern['Quad8'] = [0,1,2,3,4,5,6,7]
    nodePattern['Quad9'] = [0,1,2,3,4,5,6,7,8]
    nodePattern['Hex8']  = [0,1,2,3,4,5,6,7]
    nodePattern['Hex20'] = [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]
    nodePattern['Hex27'] = [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15, 24,20,21,22,23, 25,26]

    # the inverse node pattern is from FEval to GMV
    nodePatternInv = {}
    nodePatternInv['Tri3']  = [0,1,2]
    nodePatternInv['Tri6']  = [0,1,2,3,4,5]
    nodePatternInv['Quad4'] = [0,1,2,3]
    nodePatternInv['Tet4']  = [0,1,2,3]
    nodePatternInv['Quad8'] = [0,1,2,3,4,5,6,7]
    nodePatternInv['Quad9'] = [0,1,2,3,4,5,6,7,8]
    nodePatternInv['Hex8']  = [0,1,2,3,4,5,6,7]
    nodePatternInv['Hex20'] = [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]
    nodePatternInv['Hex27'] = [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15, 21,22,23,24,20, 25,26]

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

    def compose_material(self):
        """"""
        lines = []
        if not 'elem' in self.model.getSetTypes():
            return lines
        setnames = sorted(self.model.getSetNames('elem'))
        setnames.append('None')
        lines.append('material %i 0 \n' % (len(setnames)))
        for setname in setnames:
            lines.append(setname[:8]+'\n')
        setlines = ['%d\n' % (len(setnames))] * len(self.model.Conn)
        for j, setname in enumerate(setnames[:-1]):
            for ele in self.model.getSet('elem', setname):
                setlines[ele-1] = '%d\n' %(j+1)
        lines.extend(setlines)
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
        nnodes = N.take( nodes, self.nodePattern[ elemtype ] )
        self.model.setElement(self.elemNumber,
                              elemtype,
                              nnodes )

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

    def extract_tri(self, linelist):
        """Extract the cells of type quad"""
        self.extract_element(linelist)

    def extract_6tri(self, linelist):
        """Extract the cells of type quad"""
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
	    nodes = conn[1]
  	    nodes = N.take(nodes, self.nodePatternInv[ conn[0] ] )
            nnod = len(nodes)
            ll = '%i '*nnod % tuple(nodes)
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
                    lines.append('%s 1\n' % var.replace(' ', '_'))
                    idx = nodVarInfo.index(var)
                    for v in allvars:
                        lines.append('%f ' % v[idx])
                    lines.append('\n')
            lines.append('endvars\n')

        return lines

### Test

if __name__ == '__main__':

#     m = feval.FEval.FEModel()
#     import feval.fecodes.xdr.LibmeshFile as libm
#     lf = libm.LibmeshFile(m)
#     lf.readFile('/soft/numeric/libmesh/reference_elements/3D/one_hex27.xda')
#     m.update()

#     e = m.getElement(0)

#     nodes = N.take(e.nodes,e.shape.sidenodes[0])
#     print nodes

#     m.renumberElements(1)
#     m.renumberNodes(1)

#     for n in m.getCoordinateNames():
#         m.setNodVar(n, m.getCoordinate(n))
#     m.setNodVarInfo(['U','V','W'])

#     gf = GMVFile(m)

#     gf.setWrite('gmvinput')
#     gf.setWrite('nodes')
#     gf.setWrite('cells')
#     gf.setWrite('variable')
#     gf.setWrite('endgmv')
#     gf.writeFile('hex27.gmv')


# #    infilename  = os.path.join( feval.__path__[0], 'data', 'gmv', 'test1.gmv' )
# #    infilename  = os.path.join( '/soft/numeric/feval', 'data', 'gmv', 'test1.gmv' )
# #     outfilename = os.path.join( feval.__path__[0], 'data', 'gmv', 'test1_out.gmv' )

    m  = FEModel()
# #     from feval.fecodes.marc.MarcT16File import *  
# #     infilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal.t16')
# #     #     outfilename = os.path.expanduser('~/projects/jako/marc/jako3dd_polythermal.gmv')
# #     outfilename = os.path.expanduser('~/projects/colle/marc/cgx8_dc0.9.gmv')
# #     #mf = MarcT16File(m, infilename)
# #     # mf = MarcT16File(m, '/home/tinu/projects/colle/marc/cgx8_dc0.90/cgx8dec.t16')
# #     mf = MarcT16File(m, '/home/tinu/numeric/marc/newsl_lafiua_tf4.0-tl4.5-tu2.5-g0.2-n4.-a5.3-v100-wl-0.06-0.0-550.0.t16')
# #     mf.readInc(-1)

# #     stop
# #     from feval.fecodes.marc.MarcFile import *  
# #     mf = MarcFile(m)
# #     mf.readFile('/home/tinu/numeric/marc/ca_fm.dat')


    gf = GMVFile(m)
    gf.readFile('/scratch/tinu/species/final/hist_bg=0.006,T=0,ELA=2700,dELA=0/results/out.gmv.00000')




