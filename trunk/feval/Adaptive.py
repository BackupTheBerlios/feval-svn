# subdivide a Hex element
# ideas from LBIE

import numpy as N
import feval.FEval
import random

def subdivideHex(model, elem, lnodes):
    """subdivide a Hex8 element |elem|
    around the nodes |lnodes|
    remove |elem| from the model |model| and
    insert the subelements
    """
    assert(elem.shape.name == 'Hex8')

    # one edge needs to be refined
    if len(lnodes) == 1:
        corneridx = [[ 0,  1,  5,  4, 16, 17, 21, 20],
                     [ 1,  3, 15,  5, 17, 51, 63, 21],
                     [ 4,  5, 15, 12, 20, 21, 63, 60],
                     [ 16, 17, 21, 20, 48, 51, 63, 60]
                     ]
    elif len(lnodes) == 2:
        corneridx = [
                      [ 0,  1,  5,  4, 16, 17, 21, 20],
                      [ 1,  2,  6,  5, 17, 18, 22, 21],
                      [ 5,  6, 10,  9, 21, 22, 42, 41],
                      [ 2,  3,  7,  6, 18, 19, 23, 22],
                      [17, 18, 22, 21, 33, 34, 42, 41],
                      [ 4,  5,  9, 12, 20, 21, 41, 60],
                      [ 6,  7, 15, 10, 22, 23, 63, 42],
                      [ 9, 10, 15, 12, 41, 42, 63, 60],
                      [16, 17, 21, 20, 48, 33, 41, 60],
                      [18, 19, 23, 22, 34, 51, 63, 42],
                      [33, 34, 42, 41, 48, 51, 63, 60],
                     ]
    elif len(lnodes) == 4:
        corneridx = [
                      [ 0,  1,  5,  4, 16, 17, 21, 20],
                      [ 1,  2,  6,  5, 17, 18, 22, 21],
                      [ 2,  3,  7,  6, 18, 19, 23, 22],
                      [ 4,  5,  9,  8, 20, 21, 25, 24],
                      [ 5,  6, 10,  9, 21, 22, 26, 25],
                      [ 6,  7, 11, 10, 22, 23, 27, 26],
                      [ 8,  9, 13, 12, 24, 25, 29, 28],
                      [ 9, 10, 14, 13, 25, 26, 30, 29],
                      [10, 11, 15, 14, 26, 27, 31, 30],
                      [17, 18, 22, 21, 33, 34, 72, 71],
                      [22, 23, 27, 26, 72, 39, 43, 73],
                      [25, 26, 30, 29, 74, 73, 46, 45],
                      [20, 21, 25, 24, 36, 71, 74, 40],
                      [21, 22, 26, 25, 71, 72, 73, 74],
                      [16, 17, 21, 20, 48, 33, 71, 36],
                      [18, 19, 23, 22, 34, 51, 39, 72],
                      [26, 27, 31, 30, 73, 43, 63, 46],
                      [24, 25, 29, 28, 40, 74, 45, 60],
                      [71, 72, 73, 74, 33, 34, 46, 45],
                      [36, 71, 74, 40, 48, 33, 45, 60],
                      [72, 39, 43, 73, 34, 51, 63, 46],
                      [33, 34, 46, 45, 48, 51, 63, 60],
                    ]
        coord = elem.findGlobalCoord(N.array([1., 1., 1.5])/3.*2.-1.)
        m.setCoordinate(elem.name+'71', coord)
        coord = elem.findGlobalCoord(N.array([2., 1., 1.5])/3.*2.-1.)
        m.setCoordinate(elem.name+'72', coord)
        coord = elem.findGlobalCoord(N.array([2., 2., 1.5])/3.*2.-1.)
        m.setCoordinate(elem.name+'73', coord)
        coord = elem.findGlobalCoord(N.array([1., 2., 1.5])/3.*2.-1.)
        m.setCoordinate(elem.name+'74', coord)
    elif len(lnodes) == 8:
        corneridx = [
                      [ 0,  1,  5,  4, 16, 17, 21, 20],
                      [ 1,  2,  6,  5, 17, 18, 22, 21],
                      [ 2,  3,  7,  6, 18, 19, 23, 22],
                      [ 4,  5,  9,  8, 20, 21, 25, 24],
                      [ 5,  6, 10,  9, 21, 22, 26, 25],
                      [ 6,  7, 11, 10, 22, 23, 27, 26],
                      [ 8,  9, 13, 12, 24, 25, 29, 28],
                      [ 9, 10, 14, 13, 25, 26, 30, 29],
                      [10, 11, 15, 14, 26, 27, 31, 30],
                      [16, 17, 21, 20, 32, 33, 37, 36],
                      [17, 18, 22, 21, 33, 34, 38, 37],
                      [18, 19, 23, 22, 34, 35, 39, 38],
                      [20, 21, 25, 24, 36, 37, 41, 40],
                      [21, 22, 26, 25, 37, 38, 42, 41],
                      [22, 23, 27, 26, 38, 39, 43, 42],
                      [24, 25, 29, 28, 40, 41, 45, 44],
                      [25, 26, 30, 29, 41, 42, 46, 45],
                      [26, 27, 31, 30, 42, 43, 47, 46],
                      [32, 33, 37, 36, 48, 49, 53, 52],
                      [33, 34, 38, 37, 49, 50, 54, 53],
                      [34, 35, 39, 38, 50, 51, 55, 54],
                      [36, 37, 41, 40, 52, 53, 57, 56],
                      [37, 38, 42, 41, 53, 54, 58, 57],
                      [38, 39, 43, 42, 54, 55, 59, 58],
                      [40, 41, 45, 44, 56, 57, 61, 60],
                      [41, 42, 46, 45, 57, 58, 62, 61],
                      [42, 43, 47, 46, 58, 59, 63, 62],
                    ]
    else:
        print 'hallo'
        pass

    for i, idx in enumerate(corneridx):
        for n in idx:
            if n < 64 and not n in [0,3,12,15,48,51,60,63]:
                lcoord = N.array((n%4, n//4 %4, n//16))/3.*2. -1.
                coord  = elem.findGlobalCoord(lcoord)
                m.setCoordinate(elem.name+'%d' % n, coord)
        m.setElementConn(elem.name + '_%s' % i, 'Hex8', [elem.name+'%d' % n for n in idx])
    m.removeElement(e.name)



if __name__ == '__main__':
    m = feval.FEval.FEModel()
    m.setCoordinate(1, [0.,0.,0.])
    m.setCoordinate(2, [9.,0.,0.])
    m.setCoordinate(3, [9.,6.,0.])
    m.setCoordinate(4, [0.,6.,0.])
    m.setCoordinate(5, [0.,0.,3.])
    m.setCoordinate(6, [9.,0.,3.])
    m.setCoordinate(7, [9.,6.,3.])
    m.setCoordinate(8, [0.,6.,3.])
    m.setConn('x', 'Hex8', range(1,9))
    m.makeModelCache()
    m.renumberNodes(base = 0)
    
    e =  m.findClosestElement(N.asarray([0., 1., 0.5]))

    subdivideHex(m, e, [0,1,2,3, 4, 5, 6, 7])

    m.renumberNodes()
    m.renumberElements()

    import feval.fecodes.gmv.GMVFile as gmv
    mf = gmv.GMVFile(m)

    mf.setWrite('gmvinput')
    mf.setWrite('nodes')
    mf.setWrite('cells')
    mf.setWrite('variable')
    mf.setWrite('endgmv')

    mf.writeFile('adaptive.gmv.000')



