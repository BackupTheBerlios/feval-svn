# use py_post from MSC.MARC to access MARC Result files

# set the path and import the MARC library
import sys
#sys.path.append('../../../lib')
sys.path.append('/soft/numeric/feval/trunk/lib')
import py_post as M
import numpy as N
import MarcShapeFunctions

class MarcPostFile(object):
    """
    """

    def __init__(self, model, filename):
        self.model = model
        self.filename = filename
        try:
            self.file = M.post_open(filename)
        except:
            print "could not open postfile %s" %filename
            sys.exit(0)

    def readInc(self, inc):
        self.file.extrapolation('linear')  # linear extrapolation of integration point data
        try:
            self.file.moveto(inc+1)
        except:
            print "could not find increment %d" % inc

        # read the nodal data
        nodvarinfo = []
        nnodvar  = self.file.node_scalars()
        nnodvect = self.file.node_vectors()
        nodvarinfo = [self.file.node_scalar_label(i)
                      for i in range(nnodvar)]

        # read the integration point data
        # this data is extrapolated to the nodal data,
        intptvarinfo = []
        nintptvar = self.file.element_scalars()
        intptvarinfo = [self.file.element_scalar_label(i)
                        for i in range(nintptvar)]
        nodvarinfo.extend(intptvarinfo)
        self.model.setNodVarInfo( nodvarinfo )

        # read the nodes
        for n in range(self.file.nodes()):
            node = self.file.node(n)
            self.model.setCoordinate(node.id, N.array([node.x, node.y, node.z]))
            nodvar = [self.file.node_scalar(n, i) for i in range(nnodvar)]
            self.model.setNodVar( node.id, N.array(nodvar) )

        # read the elements
        # create a variable to store the integration point data
        # which gets extrapolated to the nodes
        newvals = {}
        for n in range(self.file.elements()):
            elem = self.file.element(n)
            sh = MarcShapeFunctions.MarcShapeFunctionDict[elem.type]
            # if not sh[2] == elem.len:
            #     print 'cutting element nodes from %d to %d ' %(elem.len, sh[2])
            self.model.setElement(elem.id, sh[0], elem.items[:sh[2]] )
            for i in range(nintptvar):
                sca = self.file.element_scalar(n, i)
                v = newvals.setdefault(i,{})
                for m in range(len(sca)):
                    a =  sca[m]
                    v.setdefault(a.id, []).append(a.value)
            
        # collect the integration point variable values and average them
        intptvars = {}
        for var in newvals:
            for node, vals in newvals[var].items():
                intptvars.setdefault(node,[]).append(N.array(vals).mean())
                
        for node, vals in intptvars.items():
            nodvars = self.model.NodVar[node]
            self.model.setNodVar( node, N.hstack((nodvars, vals)) )

        self.model.update()

    def IncrementInfo(self):
        """Plot information about increments"""
        try:
            return self.IncInfo
        except:					
            self.IncInfo = []
            for inc in range(self.file.increments()-1):
                self.file.moveto(inc+1)
                self.IncInfo.append((self.file.increment, self.file.time))
            return self.IncInfo

if __name__ == '__main__':
    from feval.FEval import *
    
    m  = FEModel()
    mf = MarcPostFile(m, '../../../data/marc/e7x1b.t16')
#    mf = MarcPostFile(m, '/home/tinu/projects/gorner/model/marc-riesen/gorner_nsl.t16')
#    mf = MarcPostFile(m, '/home/tinu/projects/colle/marc/ca8_t.t16')
    mf = MarcPostFile(m, '/home/soft/numeric/feval/trunk/data/marc/1_10_long.t16')
    mf.readInc(1)
    f = mf.file

    stop
    print 'read finished'

    m.renumberNodes()
    m.renumberElements()
    m.update()
    
    import feval.fecodes.gmv.GMVFile as gmv
    gf = gmv.GMVFile(m)

    gf.setWrite('gmvinput')
    gf.setWrite('nodes')
    gf.setWrite('cells')
    gf.setWrite('variable')
    gf.setWrite('endgmv')

    gf.writeFile('/home/tinu/projects/gorner/model/marc-riesen/gorner_nsl.gmv')
