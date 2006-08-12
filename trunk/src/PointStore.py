"""
Provide a container for element storage

This is especially useful if values have to be evaluated
over and over on the same model grid 
"""

import numpy as N

class PointStore:
    """
    A container to store elements and local coordinates for a set of points.
    Points are given as (a list of) global coordinates.
    """
    def __init__(self, model, name="", comment="", verbose=0):
        """provide a |model| to be evaluated and eventually a |name|"""
        if name:
            self.name = name
        elif model.name:
            self.name = model.name
        else:
            self.name = 'noname'
        self.comment = comment
        self.store = []
        self.verbose = verbose
        self.setModel( model )

    def __str__(self):
        return 'PointStore instance "%s"\n%d points stored' % (self.name, len(self.store))

    def __repr__(self):
        return self.__str__()

    def setModel(self, model, name=None):
        self.model = model

    def addPoint(self, point, verbose=0):
        """add a point to the store if it exists in the model"""
        point = N.asarray(point)
        e = self.model.findElement(point)
        if e:
            self.store.append((point, e.name, e.lcoord))
        else:
            self.store.append((point, None, None))
            if verbose or self.verbose:
                print "Pointstore: point is not within the model", point

    def addPoints(self, pointList, verbose=0):
        """add a list of Points"""
        for point in pointList:
            self.addPoint(point, verbose)

    def getVariables(self, nodvars=None, intpointvars=None, novalue=None ):
        """get the variables for all points in the store
        the variables are retured as the pair (nv, iv)
        where nv is a list of nodal variables and iv is a list of
        integration point variables
        """
        allnodvars, allintpointvars = [], []
        for pt, e, lcoord in self.store:
            if e:
                element = self.model.getElement(e)
                element.setLocalCoord(lcoord)
                nv, iv = self.model.getVariables( element,
                                                  nodvars=nodvars,
                                                  intpointvars=intpointvars )
                allnodvars.append(nv)
                allintpointvars.append(iv)
            else:
                allnodvars.append(novalue)
                allintpointvars.append(novalue)
        return allnodvars, allintpointvars

    def getPoints(self):
        allpoints = []
        for pt, e, lcoord in self.store:
            allpoints.append(pt)
        return allpoints

    def save(self, filename):
        """Save the point list to a file"""
        import cPickle
        file = open(filename,'w')
        try:
            file.write('!! created by PointStore 0.1, do not edit !!\n')
            file.write('Name:    '+self.name+'\n')
            file.write('Comment: '+self.comment+'\n')
            file.write('-------------------\n')
            cPickle.dump(self.name, file)
            cPickle.dump(self.comment, file)
            cPickle.dump(self.store, file)
        finally:
            file.close()

    def load(self, filename):
        """Load the point list from the file"""
        import cPickle
        file = open(filename)
        # read the header
        for i in range(4):
            file.readline()
        # read the data
        try:
            name = cPickle.load(file).strip()
            comment = cPickle.load(file).strip()
            self.store = cPickle.load(file)
            if name:
                self.name = name
            if comment:
                self.comment = comment
        finally:
            file.close()
            
if __name__ == '__main__':
    from feval.fecodes.marc.MarcT16File import *

    from feval.FEval import *

    print 'loading the model'
    ## 2D-model test case
    ## load the model as usual
    m = FEModel(verbose=0)
    mf = MarcT16File(m, 'data/marc/e7x1b.t16')
    mf.readInc(2)

    ## now create a point store
    store = PointStore(m, name='e7x1b', verbose=0)
    ## add some points to the point store
    store.addPoint([1,1], verbose=0)
    store.addPoints([[10.1,1],
                     [20,1],
                     [5.2,0.1]], verbose=0)    

    print '\n============\nOld Store \n============\n '
    ## get some model variables from the model at the points in the store
    nv, iv = store.getVariables(nodvars=['d_x','d_y','d_z'])
    for n in  nv:
        print n
    ## save the pointstore for later use
    store.name    = 'e7x1b test nodes'
    store.comment = 'A test-file for the PointStore'
    store.save('data/test/test.store')

    ## load a pointstore (probably in a later model run)
    newStore = PointStore(m, 'New Store')
    newStore.load('data/test/test.store')
    print '\n============\nNew Store \n============\n '
    print 'Name    : %s' % newStore.name
    print 'Comment : %s\n' % newStore.comment
    ## get some model variables from the model at the points in the store
    nv, iv = store.getVariables(nodvars=['d_x','d_y','d_z'])
    for n in  nv:
        print n
        
    
