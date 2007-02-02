#
# 	$Id: MarcT16File.py 134 2006-04-04 15:29:57Z tnoo $	
#

import os, re, struct
from feval.FortranIO import FortranBinaryFile
#import Numeric; N = Numeric
import numpy as N

#============================================================================
# try to use Psyco (psyco.sourceforge.net)
# if configured, this will speed up things considerably
try:
    import psyco
    from psyco.classes import *
except ImportError:
    class _psyco:
        def jit(self):      pass
        def bind(self, f):  pass
        def proxy(self, f): return f
    psyco = _psyco()

# constants
MaxFileSizeInMemory = 10*1000*1000
IncPointerOffset    = -4  # this is the size of the Fortran record descriptor 

class MarcT16File(FortranBinaryFile):
    """MarcT16File: binary MARC output file
    Constructor: MarcT16File(filename, mode, endian)
    """

    # Dictionary of Analysis types: this gives the names of the nodal variables
    nodeVarNames     = {}
    nodeVarNames[1]  = ['d', 'xf']		      # displacements
    nodeVarNames[2]  = ['d', 'xf', 'r']	      # displacements w/r
    nodeVarNames[6]  = ['d', 'xf', 'v', 'r']  # rigid plastic w/r
    nodeVarNames[7]  = ['t']			      # heat transfer
    nodeVarNames[8]  = ['t', 'f']		      # heat transfer w/f
    nodeVarNames[42] = ['d', 'xf', 't', 'v', 'r','f']  # rigid plastic coupled w/r                

    intpVarNames      = {}
    intpVarNames[-3]  = 'z'
    intpVarNames[-2]  = 'y'
    intpVarNames[-1]  = 'x'
    intpVarNames[1]   = 'eps_1'
    intpVarNames[2]   = 'eps_2'
    intpVarNames[3]   = 'eps_3'
    intpVarNames[4]   = 'eps_4'
    intpVarNames[5]   = 'eps_5'
    intpVarNames[6]   = 'eps_6'
    intpVarNames[7]   = 'eps_II_pl'
    intpVarNames[8]   = 'eps_II_cr'
    intpVarNames[9]   = 't'
    intpVarNames[10]  = 'dt'
    intpVarNames[11]  = 'sigma_1'
    intpVarNames[12]  = 'sigma_2'
    intpVarNames[13]  = 'sigma_3'
    intpVarNames[14]  = 'sigma_4'
    intpVarNames[15]  = 'sigma_5'
    intpVarNames[16]  = 'sigma_6'
    intpVarNames[17]  = 'sigma_II'
    
    # all Elements with Hermann formulation
    # should be the list with all elements with extra nodes (not checked!!)
    ignorePressureNodeElements = [33,34,35,58,59,60,61,63,66,74,80,
                                  81,82,83,84,118,119,120,128,129,130,155,156,157]
    
    # Dictionary of the Element types
    import MarcShapeFunctions 
    shapeFunctionDict = MarcShapeFunctions.MarcShapeFunctionDict
    
    def __init__(self, model, filename, mode='r', endian = '@', verbose=0):
        self.model = model
        # big endian is only for files calculated on Suns or HP
        # otherwise use '<'
        FortranBinaryFile.__init__(self, filename, mode, endian, verbose=verbose)

        # dictionary of model properties, indexed with the
        # standard variable names of the MARC manual D
        self.Info = {}	  
        self.dataFiles = [self]
        self.findIncrements()

        # read title and non incremental data
        self._readBlock1()
        self._readBlocks2to13()	  

    def findIncrements(self):
        """Find the begin of Increments, marked with '****'"""
        # reset the file pointer, and look for '****'
        self.file.seek(0)		 
        rex = re.compile(r'[\*]{4}')	 
        self.incpointer = []
        oldpos = 0
        for i in range(self.filesize/MaxFileSizeInMemory + 1):
            data = self.file.read(MaxFileSizeInMemory) 
            m = rex.search(data)
            while m:
                pos = m.start() + IncPointerOffset
                self.incpointer.append(pos + i*MaxFileSizeInMemory)
                m = rex.search(data, pos + 10)
                oldpos = pos
        self.incpointer = N.array(self.incpointer)
        # reset the file pointer
        self.file.seek(0)		

    def IncrementInfo(self):
        """If increment info is not read: just do it"""
        try:
            return self.IncInfo
        except:					
            self.IncInfo = []
            f = self.dataFiles[-1]
            for ip in f.incpointer:
                f.file.seek(ip)
                f._readBlocks15to18()
                self.IncInfo.append( (f.Info['INC'], f.Info['TIME']) )
            self.IncInfo.sort()
            return self.IncInfo

    def readInc(self, inc):
        """Read the |inc|'th increment from the file. This is not the
        absolute increment number, but the increment at place inc in
        self.incpointer. The trick with self.dataFiles allows to read
        the results from several output files in one model. This is
        how decompsed data is read 
        """
        for f in self.dataFiles:
            if N.absolute(inc) in range(len(f.incpointer)):
                f.file.seek(f.incpointer[inc])
                f._readBlocks15to18()
                f._readBlocks19to27()
            else:
                pass
        # make a new model cache
        self.model.makeModelCache()

    def __readRecord(self, type=None):
        """Read a record from the binary file
        If the type is non-standard, a special interpretation of the
        record is invoked
        """
        record = FortranBinaryFile.readRecord(self, type)
        if record[0] <> None:
            return record
        else:
            # special record types
            record = record[1]
            if type == 'c4':			   
                lr = len(record)/4
                rr = ''
                for i in range(lr):
                    rr = rr+(struct.unpack(self.endian+'c', record[i*4]))[0]
                return rr
            elif type == 'postcode':
                id = (struct.unpack(self.endian+'i', record[0:4]))[0]
                record = record[4:]
                lr = len(record)/4
                rr = ''
                for i in range(lr):
                    rr = rr+(struct.unpack(self.endian+'c', record[i*4]))[0]
                return id, rr
            elif type == 'coord':
                id = (struct.unpack(self.endian+'i', record[0:4]))[0]
                record = record[4:]
                rr = self.fromstring(record, 'f')
                return id, rr

    # Block 1
    # ReadTitle
    def _readBlock1(self):
        title = self.__readRecord('c4')
        self.model.setName( title.strip() )
        stop

    def _readBlocks2to13(self):
        # Lists with the node and element names
        # These are used to correctly save the respective Quantities
        self.__NodeNames = []
        self.__ElementNames = []

        # Block 2
        # ReadHead, head1
        n = ('INUM', 'LNUM', 'MNUM', 'NDEG', 'NSTRES', 'INOD', 'IPSTCC', \
             'NADTIE', 'NCRD', 'NNODMX', 'IANTYP', 'ICOMPL', 'NBCTRA', \
             'POSTRV', 'NDISTL', 'NSET', 'NSPRNG', 'NDIE') 
        v = self.__readRecord('i')
        for i in range(len(n)):
            self.Info[n[i]] = v[i]

        # Block 3
        # ReadHead, head2
        n = ('NESETS', 'NNSETS', 'NISETS', 'NLSETS', 'NDSETS', 'NINSET', \
             'KELEM', 'KNODE', 'KINT', 'KLAYR', 'KDOF', 'KINC')
        v = self.__readRecord('i')
        for i in range(len(n)):
            self.Info[n[i]] = v[i]

        # Block 4
        # Dummy
        if self.Info['POSTRV'] >= 7:
            x = self.__readRecord()

        # Block 5
        # Domain decomposition ID
        (NPROCD, IDOMID) = (0, 0)
        if self.Info['POSTRV'] >= 7:
            (NPROCD, IDOMID) = self.__readRecord('i')
        self.Info['NPROCD'] = NPROCD
        self.Info['IDOMID'] = IDOMID

        # domain decomposition
        if NPROCD >=1 and IDOMID == 0:
            if len(self.dataFiles) == 1:	  # first time
                # loop over the domains
                for nproc in range(self.Info['NPROCD']):
                    path = os.path.dirname(self.filename)
                    file = os.path.basename(self.filename)
                    procfile = os.path.join(path, str(nproc+1)+file)
                    f = MarcT16File(self.model, procfile)
                    self.dataFiles.append( f )
	    return

        # Block 6
        # Element variables postcodes
        self.Post = []
        postcodes = []
        for i in range(self.Info['INUM']):
            p = self.__readRecord('postcode') 
            self.Post.append( p )
            if self.intpVarNames.has_key(p[0]):
                postcodes.append( self.intpVarNames[p[0]] )
            else:
                postcodes.append( p[0] )
        self.model.setIntPointVarInfo(postcodes)

        # Block 7
        # Element connectivities
        if (self.Info['MNUM'] >= 0 and self.Info['IPSTCC'] <> 0):
            for i in range(self.Info['MNUM']):
                c = self.__readRecord('i')
                id = c[0]
                # in decomposed models, some elements might be already defined
                if not self.model.Conn.has_key(id):
                    eletype = self.shapeFunctionDict[c[1]][0]
                    nnodes  = self.shapeFunctionDict[c[1]][2]
                    # we do the nnodes trick to cut off pressure nodes (Hermann elements)
                    # this is more elegant 
                    self.model.setConn(id, eletype, c[3:3+nnodes] )
                self.__ElementNames.append(id)


        # Block 8
        # Nodal coordinates
        if (self.Info['NCRD'] >= 0 and self.Info['IPSTCC'] <> 0):
            for i in range(self.Info['LNUM']):
                id, c = self.__readRecord('coord')
                self.model.setCoordinate(id, c)
                self.__NodeNames.append(id)

        # Block 9
        # Spring data
        if (self.Info['NSPRNG'] <> 0):
            self.__readRecord()

        # Block 10
        # Nodal transformations
        x = self.NodTrans = self.__readRecord('i')


        # Block 11
        # Ties due to adaptive meshing
        if (self.Info['NADTIE'] <> 0):
            self.TieData = self.__readRecord('i')

        # Block 12
        # Transformations
        if (self.Info['NBCTRA'] <> 0 and self.Info['IPSTCC'] <> 0):
            self.TransData = self.__readRecord('f')

        # Block 13
        # Set definition
        self.Sets = {}
        if (self.Info['NSET'] <> 0 and self.Info['IPSTCC'] <> 0):
            for i in range(self.Info['NSET']):
                setname = self.__readRecord('c4').strip()
                # read the number of set elements
                settype = self.__readRecord('i')[1]
                self.Sets[setname] = (settype, self.__readRecord('i'))

        # Block 14
        # Contact Geometric Data
        # not implemented

    def _readBlocks15to18(self):
        # Block 15
        # Begin Increment Indicator
        a = self.file.read(12)	 # record size in bytes

        # Block 16
        # Loadcase title
        if self.Info['POSTRV'] >= 7:
            self.loadcasetitle = self.__readRecord('c4')

        # Block 17
        # Integer Increment verification data
        n = ('NEWCC', 'INC', 'SUBINC', 'JANTYP',  'KNOD', 'IDMY1') 
        v = self.__readRecord('i')
        for i in range(len(n)):
            self.Info[n[i]] = v[i]
        # calculate JNODE (Manual D 9-23)
        self.Info['JNODE']=self.Info['KNOD']/self.Info['NDEG']

        # Block 18
        # Real Increment verification data
        n = ('TIME', 'FREQ', 'GMAS', 'DMY2', 'DMY3', 'DMY4')
        v = self.__readRecord('f')
        for i in range(len(n)):
            self.Info[n[i]] = v[i]

    def _readBlocks19to27(self):
        # Block 19
        # New non-incremental data
        if self.Info['NEWCC'] <> 0:
            self._readBlocks2to13()

        # Block 20
        # Magnitude of distributed loads
        if (self.Info['NDISTL'] <> 0 and self.Info['JANTYP'] <> 60 \
            and self.Info['JANTYP'] <> 61): 
            self.MagDist = self.__readRecord('f')

        # Block 21
        # Magnitude of spring forces
        if (self.Info['NSPRNG'] <> 0 or self.Info['JANTYP'] == 60 \
            or self.Info['JANTYP'] == 61):
            x =self.__readRecord('f')

        # Block 22
        # Magnitude of die forces
        if (self.Info['NDIE'] <> 0 and self.Info['JANTYP'] <> 60 \
            and self.Info['JANTYP'] <> 61):
            x = self.__readRecord('f')

        # Block 23
        # Values of element integration point variables
        if (self.Info['INUM'] <> 0 and self.Info['MNUM'] <> 0 \
            and self.Info['JANTYP'] >= 100):
            for i in xrange(self.Info['MNUM']):		 #	elements
                IntPointVar = []
                for j in range(self.Info['NSTRES']): #	integration points
                    IntPointVar.append( self.__readRecord('f') )
                IntPointVar = N.array(IntPointVar)
                id = self.__ElementNames[i]
                self.model.setIntPointVar(id, IntPointVar)

        # Block 24
        # Values of nodal variables
        if (self.Info['KNOD'] <> 0 and self.Info['JANTYP'] <> 60 \
            and self.Info['JANTYP'] <> 61):
            for i in xrange(self.Info['LNUM']):
                id = self.__NodeNames[i]
                self.model.setNodVar( id, self.__readRecord('f') )
            try:
                typelist = self.nodeVarNames[self.Info['IANTYP']]
            except:
                print '!!!!! ERROR: Analysis Type %i not implemented !!!!!' \
                      % (self.Info['IANTYP'])

            # set the names of the nodal variables
            dim = self.Info['NDEG']   # Maximum number of DoF's
            # (Manual D 9-5)
            types = []
            if dim > 1:
                for t in typelist:
                    for d in range(dim):
                        types.append(t+self.model.DimNames[d])
            else:
                types = typelist
            self.model.setNodVarInfo( types )

        # Block 25
        # Response gradients

        # Block 26
        # Response

        # Block 27
        # Design


class MarcT16FileEq(MarcT16File):
    """MarcT16FileEq: binary MARC output file with equidistant increments
    Is considerably faster on initialization, due to the fact that not the
    whole file has to be scanned for the position of the increments.
    Constructor: MarcT16FileEq(filename, mode, endian)
    """

    def findIncrements(self):
	"""Find the begin of Increments, marked with '****'"""
	self.file.seek(0)
	rex = re.compile(r'[\*]{4}')
	self.incpointer = []
	ninc = 0
	for i in range(self.filesize/MaxFileSizeInMemory + 1):
	    data = self.file.read(MaxFileSizeInMemory) # Anzahl Bytes,
	    m = rex.search(data)
	    while m:
		pos = m.start()	 + IncPointerOffset
		self.incpointer.append(pos + i*MaxFileSizeInMemory)
		m = rex.search(data, pos + 10)
		ninc = ninc+1
		if ninc == 2:
		    ip0 = self.incpointer[0]
		    dip = self.incpointer[-1] - self.incpointer[-2]
		    nincs = (self.filesize-ip0)/dip + 1
		    self.incpointer = N.arange(nincs)*dip + ip0
		    self.file.seek(0)
		    return None
	self.incpointer = N.array(self.incpointer)
	self.file.seek(0)		# reset the file position


class MarcT16FileEq2(MarcT16File):
    """MarcT16FileEq2: binary MARC output file with equidistant increments
    fuer Gwendo's files
    Constructor: MarcT16FileEq(filename, mode, endian)
    """

    def findIncrements(self):
        """Find the begin of Increments, marked with '****'"""
        self.file.seek(0)
        rex = re.compile(r'[\*]{4}')
        self.incpointer = []
        ninc = 0
        for i in range(self.filesize/MaxFileSizeInMemory + 1):
            data = self.file.read(MaxFileSizeInMemory) # Anzahl Bytes,
            m = rex.search(data)
            while m:
                pos = m.start()	 + IncPointerOffset
                self.incpointer.append(pos + i*MaxFileSizeInMemory)
                m = rex.search(data, pos + 10)
                if ninc == 5:
                    ip0 = self.incpointer[3]
                    dip = self.incpointer[5] - self.incpointer[3]
                    self.dip_short = self.incpointer[5] - self.incpointer[4]
                    nincs = (self.filesize-ip0)/dip
                    newincs = N.arange(nincs)*dip + ip0
                    oldincs = N.array(self.incpointer[0:2])
                    self.incpointer = N.concatenate( (oldincs, newincs) ) 
                    self.file.seek(0)
                    return None
                ninc = ninc+1
        self.incpointer = N.array(self.incpointer)
        self.file.seek(0)		# reset the file position

    def readInc(self, inc):
        if inc > 1:
            self.file.seek(self.incpointer[inc]-self.dip_short)
            self._readBlocks15to18()
            self._readBlocks19to27()

        self.file.seek(self.incpointer[inc])
        self._readBlocks15to18()
        self._readBlocks19to27()
        # make a new model cache
        self.model.makeModelCache()


### Test
if __name__ == '__main__':

    from feval.FEval import *


    m  = ModelData()
    mf = MarcT16File(m, '../../../data/marc/e7x1b.t16')
    mf.readInc(1)

    print 'test finished'
    



    
