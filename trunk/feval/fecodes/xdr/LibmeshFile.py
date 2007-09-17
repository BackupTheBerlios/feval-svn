import re
import numpy as N

class LibmeshFile(object):
    """Parse and write a XDR file
    This is XDA file type defined in libmesh
    Based on the documentation in /doc/latex/xda_format/xda_format.pdf
    """
    type = 'xdr/xda'

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
        self.model       = model
        self.versionInfo = [0, 0]
        self.fileType    = 'xda'
        self.bcond       = []

    def readFile(self, filename):
        """read either an ASCII or binary file 
        """
        if filename.endswith('.xda'):
            self.readFileASCII(filename)
        elif filename.endswith('.xdr'):
            self.readFileBinary(filename)
        else:
            print "don't know what to do"

    def writeFile(self, filename, bcond=None):
        """write either an ASCII or binary file 
        """
        if filename.endswith('.xda'):
            self.writeFileASCII(filename, bcond=bcond)
        elif filename.endswith('.xdr'):
            self.writeFileBinary(filename, bcond=bcond)
        else:
            print "don't know what to do"

    def readFileBinary(self, filename):
        """read the binary file
        """
        import struct
        f = file(filename,'rb')

        bytes = f.read(12)
        filetype, level = bytes[4:8], bytes[9]
        if not filetype == 'LIBM':
            print 'file type not (yet) supported', filetype
            stop

        nelems, nnodes, sumweight, nbcond, x, neleblocks = struct.unpack('!6i', f.read(24))
        eletypes = struct.unpack('!%si'%neleblocks, f.read(4*neleblocks))
        elenums  = struct.unpack('!%si'%neleblocks, f.read(4*neleblocks))
        # read the titles
        f.read(32)

        # read the element block
        conn = N.array(struct.unpack('!%si'%sumweight, f.read(4*sumweight)))
        import feval.ShapeFunctions
        cnt = 0
        for elem, nelem in zip(eletypes, elenums):
            if not self.shapeFunctionDict.has_key(str(elem)):
                print 'No shape function found for element type', elem
            elemtype = self.shapeFunctionDict[str(elem)]
            shape = feval.ShapeFunctions.shapeFunctions[elemtype]
            blocklen = nelem*(shape.nnodes+2)
            print blocklen, nelem
            for nodes in conn[cnt:cnt+blocklen].reshape(nelem, blocklen//nelem):
                self.model.setElement(nodes[-2], elemtype, nodes[:-2] )

            cnt += blocklen

        # read the node block
        coord = N.array(struct.unpack('!%sd'%(3*nnodes), f.read(8*3*nnodes)))
        coord.shape = ( len(coord)//3, 3)
        for n, c in enumerate(coord):
            self.model.setCoordinate(n, c)

        # read the boundary conditions
        if nbcond > 0:
            bcond = N.array(struct.unpack('!%si'%(3*nbcond), f.read(4*3*nnodes)))
            bcond.shape = (nbcond, 3)
            self.bcond = bcond
        f.close()

    def readFileBinary2(self, filename):
        """read the binary file
        """
        bytes = file(filename,'rb').read()

        filetype, level = bytes[4:8], bytes[9]
        if not filetype == 'LIBM':
            print 'file type not (yet) supported', filetype
            stop

        nelems, nnodes, sumweight, nbcond, x, neleblocks = N.fromstring(bytes[12:36], '>i')
        eletypes = N.fromstring(bytes[36:36+4*neleblocks], '>i')
        elenums  = N.fromstring(bytes[36+4*neleblocks:36+8*neleblocks], '>i')

        ptr = 36+8*neleblocks +32

        # read the element block
        conn = N.fromstring(bytes[ptr:ptr+4*sumweight], '>i')
        ptr += 4*sumweight
        import feval.ShapeFunctions
        cnt = 0
        for elem, nelem in zip(eletypes, elenums):
            if not self.shapeFunctionDict.has_key(str(elem)):
                print 'No shape function found for element type', elem
            elemtype = self.shapeFunctionDict[str(elem)]
            shape = feval.ShapeFunctions.shapeFunctions[elemtype]
            blocklen = nelem*(shape.nnodes+2)
            print blocklen, nelem
            for nodes in conn[cnt:cnt+blocklen].reshape(nelem, blocklen//nelem):
                self.model.setElement(nodes[-2], elemtype, nodes[:-2] )

            cnt += blocklen

        # read the node block
        coord = N.fromstring(bytes[ptr:ptr+8*3*nnodes], '>d')
        coord.shape = ( len(coord)//3, 3)
        for n, c in enumerate(coord):
            self.model.setCoordinate(n, c)
        ptr += 8*3*nnodes

        # read the boundary conditions
        if nbcond > 0:
            bcond = N.fromstring(bytes[ptr:ptr+4*3*nbcond], '>i')
            bcond.shape = (nbcond, 3)
            self.bcond = bcond


    def readFileASCII(self, filename):
        """Parse an ASCII input file
        """

        lines = file(filename).readlines()
        # remove empty lines
        lines = [line for line in lines if len(line) > 2]
        
        # read the header
        self.fileType, self.versionInfo = lines[0].split()
        nelem      = int(lines[1].split()[0])      # number of elements
        nnodes     = int(lines[2].split()[0])      # number of nodes
        nbcond     = int(lines[4].split()[0])      # number of boundary conditions
        elemblocks = int(lines[6].split()[0])      # Num. Element Blocks.
        elemtypes  = lines[7].split()[:elemblocks] # Element types in each block.
        nlevels    = len(lines[8].split('#')[0].split())/elemblocks
        neleblock  = map( int, lines[8].split()[:elemblocks]) # Num. of elements in each block.
        if not lines[10].startswith('Title String'):
            self.model.name = lines[10][:-1].strip()

        if self.fileType == 'LIBM':
            # the element data, connectivity followed by element id and parent id
            cnt = 10
            for elem, nelem in zip(elemtypes, neleblock):
                if not self.shapeFunctionDict.has_key(elem):
                    print 'No shape function found for element type', elem
                elemtype = self.shapeFunctionDict[elem]
                for l in range(nelem):
                    cnt += 1
                    nodes = map(int, lines[cnt].split())
                    self.model.setElement(nodes[-2], elemtype, nodes[:-2] )
        else: # fileType = 'DEAL' or 'MGF'
            # the element data, connectivity 
            cnt = 10
            elenum = 0
            for elem, nelem in zip(elemtypes, neleblock):
                if not self.shapeFunctionDict.has_key(elem):
                    print 'No shape function found for element type', elem
                elemtype = self.shapeFunctionDict[elem]
                for l in range(nelem):
                    cnt += 1
                    nodes = map(int, lines[cnt].split())
                    self.model.setElement(elenum, elemtype, nodes )
                    elenum +=1

        # then the node data (coordinates)
        nodeNumber = 0
        for node in range(nnodes):
            cnt += 1
            coord = map( float, lines[cnt].split() )
            self.model.setCoordinate( nodeNumber, coord )
            nodeNumber += 1

        # read the boundary conditions
        if nbcond > 0:
            cnt += 1
            bcond = []
            for n in range(nbcond):
                bcond.append(map(int, lines[cnt+n].split()))
            self.bcond = N.array(bcond)



    def writeFileASCII(self, filename, bcond = None):
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
            elines.setdefault(etype, []).append('%d '*nnod % tuple(conn) + '%d -1\n' % ele)
            weight += nnod
        weight += 2*len(self.model.Conn)
        
        # now we are ready to write the file

        # no refinement level: LIBM 0
        lines.append('LIBM 0\n')
        lines.append('%d\t # Num. Elements\n' % len(self.model.Conn) )
        lines.append('%d\t # Num. Nodes\n'    % len(self.model.Coord) )
        lines.append('%d\t # Sum of Element Weights\n' % weight )
        nbcond = 0
        if bcond != None:
            nbcond = len( bcond )
        lines.append('%d\t # Num. Boundary Conds.\n' % nbcond)
        lines.append('65536\t # String Size (ignore)\n') 
        lines.append('%d\t # Num. Element Blocks.\n' % len(elines.keys()) )
        etypes, eblocks = '', ''
        ekeys = sorted(elines.keys())
        for t in ekeys:
            etypes  += '%s ' % self.invShapeFunctionDict[t]
            eblocks += '%s ' % len(elines[t])
        lines.append('%s\t # Element types in each block.\n' % etypes)
        lines.append('%s\t # Num. of elements in each block at each level.\n' % eblocks)
        lines.append('Id String\n')
        title = 'Title String'
        if not self.model.name.startswith('FE-Model'):
            title = self.model.name.strip()
        lines.append(title+'\n')

        # add the element lines
        for t in ekeys:
            lines.extend(elines[t])
        del elines

        # write all coordinates
        ckeys = sorted(self.model.getCoordinateNames())
        for k in ckeys:
            lines.append('%f %f %f \n' % tuple(self.model.getCoordinate(k)))

        if bcond != None:
            #lines.append('\n')
            for bc in bcond:
                lines.append('%d %d %d\n' % tuple(bc))

        file(filename,'w').writelines(lines)
        

    def writeFileBinary(self, filename, bcond=None):
        """write the binary file
        """
        # go through all elements and see which ones have the same type
        # build a dict of elements
        weight = 0 # total weight
        etypes = {}
        for ele in self.model.getElementNames():
            etype, conn = self.model.getElementConn(ele)
            nnod = len(conn)
            etypes.setdefault(etype, []).append(ele)
            weight += nnod
        weight += 2*len(self.model.Conn)

        nbcond = 0
        if bcond != None:
            nbcond = len( bcond )

        # write the data
        import struct
        f = file(filename,'wb')

        f.write(struct.pack('!L', 6))
        f.write('LIBM 0\x00\x00')
        f.write(struct.pack('!6L',
                            len(self.model.Conn),
                            len(self.model.Coord),
                            weight,
                            nbcond,
                            65536,
                            len(etypes),
                            ))
        ekeys = sorted(etypes.keys())
        for k in ekeys:
            f.write(struct.pack('!L', int(self.invShapeFunctionDict[k])))
        for k in ekeys:
            f.write(struct.pack('!L', len(etypes[k])))
        s = 'Id String'
        f.write(struct.pack('!L',len(s)))
        f.write(struct.pack('!12s', s))
        s = 'Title String'
        f.write(struct.pack('!L',len(s)))
        f.write(struct.pack('!12s', s))

        # write the element block
        for k in ekeys:
            for ele in etypes[k]:
                etype, conn = self.model.getElementConn(ele)
                f.write(struct.pack('!%sL' %len(conn), * conn))
                f.write(struct.pack('!L', ele))
                f.write(struct.pack('!L', -1))
        
        # write the coordinate block
        ckeys = sorted(self.model.getCoordinateNames())
        for k in ckeys:
            f.write(struct.pack('!3d', * self.model.getCoordinate(k)))

        if bcond != None:
            for bc in bcond:
                f.write(struct.pack('!3L', * bc))

        f.close()

        
    def writeData(self, filename, ndata={}, edata={}):
        """write either an ASCII or binary file 
        """
        if filename.endswith('.xta'):
            self.writeDataASCII(filename, ndata=ndata, edata=edata)
        elif filename.endswith('.xtr'):
            self.writeDataBinary(filename, ndata=ndata, edata=edata)
        else:
            print "don't know what to do with file", filename


    def writeDataASCII(self, filename, ndata={}, edata={}):
        lines = []
        lines.append('\t# Data description\n')
        lines.append('REAL\t# type of values\n')
        lines.append('%d\t# No. of nodes for which data is stored\n' % len(ndata))
        lines.append('%d\t# No. of elements for which data is stored\n' % len(edata))
        for k in sorted(ndata):
            d = ndata[k]
            lines.append('%d\n' % k)
            lines.append('%d\n' % len(d))
            lines.append('%e  '*len(d) % tuple(d) + '\n')
        for k in sorted(edata):
            d = edata[k]
            lines.append('%d\n' % k)
            lines.append('%d\n' % len(d))
            lines.append('%e  '*len(d) % tuple(d) + '\n' )
        file(filename,'w').writelines(lines)

### Test
if __name__ == '__main__':

    import feval.FEval 

    m = feval.FEval.FEModel()
    mf = LibmeshFile(m)
#     mf.readFile('/soft/numeric/libmesh/examples/ex10/mesh.xda')
#     #mf.readFile("/soft/numeric/libmesh/reference_elements/3D/one_hex20.xda");
#     bc = [[0,1, 77],
#           [3,2, 1000]]

#     mf.writeFile('/home/tinu/projects/libmesh/ex2/a2.xda', bcond = bc)

    mf.readFile('/home/tinu/projects/libmesh/ex2/a2.xda')
    mf.writeFile('/home/tinu/projects/libmesh/ex2/a3.xda', bcond = mf.bcond)
    #mf.writeFile('/home/tinu/projects/libmesh/ex2/a3.xdr')
    mf.writeFile('/home/tinu/projects/libmesh/ex2/a3.xdr', bcond = mf.bcond)

    ndata = {1: [1.,99.],
             2: [2., 101.]}
    edata = {10: [10.],
             20: [20., 101.],
             3: [30., 101., 33333.],
             }
    mf.writeDataASCII('/home/tinu/projects/libmesh/ex2/a3.xta', ndata = ndata, edata=edata)
#     m2 = feval.FEval.FEModel()
#     mf2 = LibmeshFile(m2)
#     mf2.readFile('/home/tinu/projects/libmesh/ex2/a2.xdr')

