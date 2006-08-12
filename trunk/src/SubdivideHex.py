# subdivide a Hex element
# from LBIE

def subdivideHex(model, elem, lnodes):
    """subdivide a Hex8 element |elem|
    around the nodes |lnodes|
    remove |elem| from the model |model| and
    insert the subelements
    """
    assert(elem.ShapeFunction.name == 'Hex')



