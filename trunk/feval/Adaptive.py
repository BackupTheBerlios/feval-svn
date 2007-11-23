# subdivide a Hex element

import numpy as N
import feval.FEval
import random

def subdivideHex(model, elem, globalrefinenodes, refine=3):
    """subdivide a Hex8 element |elem|
    around the nodes |refinenodes| (in global node numbering)
    remove |elem| from the model |model| and
    insert the subelements
    |refine| should be 3, except for refinement around one node
    when it can also be 2
    """

    assert(elem.shape.name == 'Hex8')
    elename = str(elem.name)

    # swap the element corners depending on the refine node / edge
    swapcorners = [
        [0, 1, 2, 3, 4, 5, 6, 7],   # u 0
        [1, 2, 3, 0, 5, 6, 7, 4],   # u 1
        [2, 3, 0, 1, 6, 7, 4, 5],   # u 2
        [3, 0, 1, 2, 7, 4, 5, 6],   # u 3
        [4, 7, 6, 5, 0, 3, 2, 1],   # o 0
        [5, 4, 7, 6, 1, 0, 3, 2],   # o 1
        [6, 5, 4, 7, 2, 1, 0, 3],   # o 2
        [7, 6, 5, 4, 3, 2, 1, 0],   # o 3
        [0, 4, 5, 1, 3, 7, 6, 2],   # s 0
        [1, 5, 6, 2, 0, 4, 7, 3],   # s 1
        [2, 6, 7, 3, 1, 5, 4, 0],   # s 2
        [7, 3, 0, 4, 6, 2, 1, 5],   # s 3
        ]
    # useful default value
    try:
        globalrefinenodes[0]
    except:
        globalrefinenodes = [globalrefinenodes]
    refinenodes = [elem.nodes.index(n) for n in globalrefinenodes]
    refinenode = 0
    n_refinenodes = len(refinenodes)


    # TODO: Maybe a better strategy to handle the cases of two and three nodes
    # if the number is not explicitly given, use full refinement
    if not n_refinenodes in [1,2,4,8]:
        n_refinenodes = 8

    # The new nodes for refinement are given as points on a grid of 64
    # point 0 is lower left, point 64 is upper right
    
    # only one edge needs to be refined
    if n_refinenodes == 1:
        corneridx = [[ 0,  1,  5,  4, 16, 17, 21, 20],
                     [ 1,  3, 15,  5, 17, 51, 63, 21],
                     [ 4,  5, 15, 12, 20, 21, 63, 60],
                     [16, 17, 21, 20, 48, 51, 63, 60]
                     ]
        refinenode = refinenodes[0]

    elif n_refinenodes == 2:
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
        edgeindex = {
                      (0, 1): 0,
                      (1, 2): 1,
                      (2, 3): 2,
                      (0, 3): 3,
                      (4, 5): 5,
                      (5, 6): 6,
                      (6, 7): 7,
                      (4, 7): 4,
                      (0, 4): 8,
                      (1, 5): 9,
                      (2, 6): 10,
                      (3, 7): 11
                     }
        refinenodes.sort()
        if edgeindex.has_key(tuple(refinenodes)):
            refinenode = edgeindex[tuple(refinenodes)]
        else:
            # if the specified edges dont exist refine the whole element
            # this means that two edges are refined, not a side
            n_refinenodes = 8

    elif n_refinenodes == 4:
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
        edgeindex = {
                      (0, 1, 2, 3): 0,    # unten
                      (4, 5, 6, 7): 5,    # oben
                      (0, 1, 4, 5): 8,    # !!!!!!!! falsche Transfromation
                      (1, 2, 6, 5): 9,    # rechts
                      (2, 3, 6, 7): 10,   # hinten
                      (0, 3, 4, 7): 11,   # links
                     }
        refinenodes.sort()
        if edgeindex.has_key(tuple(refinenodes)):
            refinenode = edgeindex[tuple(refinenodes)]
        else:
            # if the specified edges dont exist refine the whole element
            # this means that two edges are refined, not a side
            n_refinenodes = 8

    # refine around all nodes (i.e. refine the whole element)
    elif n_refinenodes == 8:
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

    # change the order of the corner nodes of the element according to the
    # refinement nodes 
    # this is no problem since the element gets removed at the end
    elem.nodes    = [elem.nodes[x] for x in swapcorners[refinenode]]
    elem.nodcoord = [elem.nodcoord[x] for x in swapcorners[refinenode]]

    # factor where to place the new nodes
    fact = 2./float(refine)

    # additional nodes on a 2.5th level, only for 4 node side
    if len(refinenodes) == 4:
        coord = elem.findGlobalCoord(N.array([1., 1., 1.5])*fact-1.)
        model.setCoordinate(elename+'_71', coord)
        coord = elem.findGlobalCoord(N.array([2., 1., 1.5])*fact-1.)
        model.setCoordinate(elename+'_72', coord)
        coord = elem.findGlobalCoord(N.array([2., 2., 1.5])*fact-1.)
        model.setCoordinate(elename+'_73', coord)
        coord = elem.findGlobalCoord(N.array([1., 2., 1.5])*fact-1.)
        model.setCoordinate(elename+'_74', coord)

    # the points at the location of the original corner nodes 
    cornerlist = [0,3,15,12,48,51,63,60]
    for i, corners in enumerate(corneridx):
        for n in corners:
            if n < 64 and not n in cornerlist:
                lcoord = N.array((n%4, n//4 %4, n//16))*fact -1.
                print lcoord
                coord  = elem.findGlobalCoord(lcoord)
                model.setCoordinate(elename+'_%d' % n, coord)
        nodenames = [elename+'_%d' % n for n in corners]
        # the order in the list is crucial!
        for n, nn in zip(elem.nodes, [elename+'_%d' % n for n in cornerlist]):
            try:
                nodenames[nodenames.index(nn)] = n
            except:
                pass

        # bring the nodes in the initial order
        # so that all faces are oriented as before refinement
        nodenames = [nodenames[j] for j in N.argsort(swapcorners[refinenode])]
        model.setElement(elename + '_%s' % i, 'Hex8', nodenames)
    model.removeElement(elem.name)

if __name__ == '__main__':
    m = feval.FEval.FEModel()
#     m.setCoordinate(1, [0.,0.,0.])
#     m.setCoordinate(2, [1.,0.,0.])
#     m.setCoordinate(3, [1.,1.,0.])
#     m.setCoordinate(4, [0.,1.,0.])
#     m.setCoordinate(5, [0.,0.,1.])
#     m.setCoordinate(6, [1.,0.,1.])
#     m.setCoordinate(7, [1.,1.,1.])
#     m.setCoordinate(8, [0.,1.,1.])
#     m.setConn('x', 'Hex8', range(1,9))
#     m.makeModelCache()
#     #m.renumberNodes(base = 1)
    
#     e =  m.findClosestElement(N.asarray([0., 1., 0.5]))

#     #subdivideHex(m, e, [0,1,2,3, 4, 5, 6, 7])
#     #subdivideHex(m, e, [0,4,1,5])
#     #subdivideHex(m, e, [7], fact=1.)
#     subdivideHex(m, e, 1, fact=1.)

#     m.renumberNodes()
#     m.renumberElements()

#     import feval.fecodes.gmsh.GMSHFile
#     gf = feval.fecodes.gmsh.GMSHFile.GMSHFile(m)
#     gf.setWrite('nod')
#     gf.setWrite('elm')
#     gf.writeFile('X.gmsh')


    import feval.Hexmesh
    print 'creating/reading file'

    # create Hexmesh
    feval.Hexmesh.makeCube(m,
                 0., 2.,
                 0., 3.,
                 0., 4.,
                 6, 6, 6)
    m.renumberNodes()
    m.renumberElements()
    
    refinenode = 284
    elems = m.nodeIndex[refinenode]
    for ele in elems:
        subdivideHex(m, m.getElement(ele), refinenode)

    m.sweepNodes()
    m.renumberNodes()
    m.renumberElements()
    
    import feval.fecodes.gmsh.GMSHFile
    gf = feval.fecodes.gmsh.GMSHFile.GMSHFile(m)
    gf.setWrite('nod')
    gf.setWrite('elm')
    gf.writeFile('X.gmsh')


