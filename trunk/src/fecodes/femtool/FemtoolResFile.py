# Programm zum Einlesen von FEMTOOL (*.res) Dateien 
# Originalroutine ist printout.f (in FEMTOOL)                          *

import re, struct
from FortranIO import *
import Numeric; N = Numeric

class FemtoolResFile(FortranBinaryFile):
	"""Femtool binary output file
	Constructor: FemtoolResFile(filename, mode, endian)
	"""

	formDict = {  2: 'Quad', \
                      3: 'Hex'}
        
        analysisTypeDict  = {}
	analysisTypeDict[1] = ['p','v_x','v_y']

	def __init__(self, model, filename, mode='r', endian = '>'):
		self.model = model
		# big endian '>' is only for files calculated on Suns or HP
		# otherwise use '<'
		FortranBinaryFile.__init__(self, filename, mode, endian)
                self.readInfo()
                self.readNodes()
                self.readElements()
                self.readVariables()

        def readInfo(self):
            """read the header and define the Info quantities
               mrec        maximum record length in output file
               npara       number of parameters to be printed (is 10)
               nstart      (is 100)
               nicoord     file record number of icoord   = 2
               nmelpoint   file record number of melpoint = 3
               nmelem      file record number of melem    = 4
               ntcoord     ntcoord=irec (any value..)
               nscoord     nscoord=ntcoord+ntall (any value..)
               nrvec       nrvec=nscoord+nsall 
               ndum        dummy
               ndumx       dummy
               nvar        number of variables of the problem
               nrpts       total number of nodes of the problem
               ndim        number of dimensions of the problem
               nrvtot
               nelem       total number of elements of the problem
               nshape      
               nelpts      integration points per element"""

            self.Info = {}	  # dictionary of model properties
            self.file.seek(0)
            header = self.readBytes(25, 'l')
            self.recLength = header[0]
            n = ('mrec', 'npara', 'nstart', 'nicoord', 'nmelpoint', 'nmelem', \
                 'ntcoord', 'nscoord', 'nrvec', 'ndum', 'ndumx', 'nvar', \
                 'nrpts', 'ndim', 'nrvtot', 'nelem', 'nshape', 'nelpts', \
                 'timeshape', 'ntimeelems', 'ntimeimpl',  'ntimeprev', \
                 'ngtype1', 'ngtype2', 'itoffxx')
            for i in range(len(n)):
                self.Info[n[i]] = header[i]
            
        def readDataRecord(self, recnr, type, len=None):
            """Read a data record"""
            pos = (recnr-1)*self.recLength*4
            if not len:
                len = self.recLength
            return self.readBytes(len, type, offset=pos)
            

        def readNodes(self):
            """icoord is a pointer array for node number"""
            ndim, nrpts = self.Info['ndim'], self.Info['nrpts'] 
            self.__NodeNames = self.readDataRecord( self.Info['nicoord'], 'l')
            self.__NodeNames = self.__NodeNames[:nrpts]
            coord = self.readDataRecord( \
                self.Info['nscoord'], 'f', len = ndim*nrpts)
            coord.shape = (ndim, nrpts)
            for i in range(nrpts):
                self.model.setCoord( self.__NodeNames[i], coord[:,i])
        
        def readElements(self):
            """melpoint is a pointer array for the element"""
            nelem  = self.Info['nelem']
            nmelem = self.Info['nmelem']
            nelpts = self.Info['nelpts']
            dim    = self.Info['ndim']
            self.melpoint = self.readDataRecord( self.Info['nmelpoint'], 'l')
            self.melpoint = self.melpoint[:nelem]

            mellen=nelem*(nelpts+4)      # number of nodes + 4 info
            conn=self.readDataRecord( self.Info['nmelpoint']+1, 'l', \
                                      len=mellen)
            conn.shape = (nelem, nelpts+4)
            for c in conn[:,]:
                eletype = self.formDict[dim] + str(c[3])
                self.model.setConn(c[0], eletype, c[4:] )

        def readVariables(self):
            """  """
            nvar, nrpts = self.Info['nvar'], self.Info['nrpts'] 
            vars = self.readDataRecord( self.Info['nrvec'], 'f', \
                                        len = nvar*nrpts)
            self.vars = vars
            vars.shape = (nvar, nrpts)
            for i in range(nrpts):
                id = self.__NodeNames[i]
                self.model.setNodVar(id, vars[:,i])
            self.model.setNodVarInfo( self.analysisTypeDict[1] )

### Test

if __name__ == '__main__':
    
    from feval.FEval import *
    m = FEModel()
#      f = FemtoolResFile(m, '~/fem/compress/ca/ca_gm1a/ca_gm1_d1.0_f10.res' )
    f = FemtoolResFile(m, '~/fem/compress/ca/ca.res' )
