#!/usr/bin/env python

####
#
#  Python script for the generation of a Femtool input file from a Marc file
#
#  Boundary conditions are calculated 
#
####

## import modules
import os
from feval.ModelData import *
from feval.fecodes.marc.MarcFile import *
from feval.fecodes.femtool.FemtoolFile import *

## file names
inputfilename = os.path.expanduser('~/projects/illimani/marc/ia.dat')
outfilename   = os.path.expanduser('~/projects/illimani/fem/ia.dat')

## initalize a FE-Model (data container with "intelligence")
m = ModelData()

## read the Marc input file
mf = MarcFile(m)
mf.readFile(inputfilename)

## initailise the Femtool File
ff = FemtoolFile(m)

## calculate and set boundary conditions in the Femtool file
## print available node sets in the Model with
##     print m.Sets['node'].keys()

print m.getSetTypes()

## set surface boundary conditions (p)
boundnodes = []
try:
    boundnodes.extend( m.getSet('node','surf_all') )
    boundnodes.extend( m.getSet('node','left_all') )
    boundnodes.extend( m.getSet('node','right_all') )
except:
    print "Warning: problems with node sets"

boundnodes.sort()
bounddiri = {}
for n in boundnodes:
    bounddiri[n] = 0.0
ff.setDirichlet(1, bounddiri)  # for the variable 1, i.e. p

## set bed boundary conditions (u,v)
boundnodes = m.Sets['node']['bed_all']
boundnodes.sort()
bounddiri = {}
for n in boundnodes:
    bounddiri[n] = 0.0
ff.setDirichlet(2, bounddiri)  # for the variable 2, i.e. u
ff.setDirichlet(3, bounddiri)  # for the variable 3, i.e. v

## set the variables and shape functions
ff.setVariables(1, '1 0 1 0 1 0 1 0')
ff.setVariables(2, '1 1 1 1 1 1 1 1')
ff.setVariables(3, '1 1 1 1 1 1 1 1')

ff.Info['nvar'] = 3
ff.Info['nshape'] = 2

## update the info block and
## set the blocks which are written to the file 
ff.updateInfo()
ff.setWrite('initialize')
ff.setWrite('elements')
ff.setWrite('points')
ff.setWrite('dirichlet')
ff.setWrite('variables')
ff.setWrite('calculate')

## write the file
ff.writeFile(outfilename)

